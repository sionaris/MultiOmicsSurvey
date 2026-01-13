#!/usr/bin/env Rscript

# =========================
# KLIC sample perturbations
# =========================

# ---- MOSEK env bootstrap (before loading Rmosek) ----
ensure_mosek_env <- function() {
  wd <- Sys.getenv("SLURM_SUBMIT_DIR", getwd())
  lic <- Sys.getenv("MOSEKLM_LICENSE_FILE", "")
  if (!nzchar(lic)) {
    cand <- file.path(wd, "mosek.lic")
    if (file.exists(cand)) {
      Sys.setenv(MOSEKLM_LICENSE_FILE = cand, MOSEK_LICENSE_FILE = cand)
    }
  } else {
    Sys.setenv(MOSEK_LICENSE_FILE = Sys.getenv("MOSEK_LICENSE_FILE", lic))
  }

  plat <- Sys.getenv("MOSEK_PLATFORM_DIR", "")
  bind <- Sys.getenv("MOSEK_BINDIR", "")

  if (!nzchar(plat) && !nzchar(bind)) {
    cand_plat <- file.path(Sys.getenv("HOME"), "mosek", "10.2", "tools", "platform", "linux64x86")
    if (dir.exists(cand_plat)) {
      plat <- cand_plat
      Sys.setenv(MOSEK_PLATFORM_DIR = plat)
    }
  }

  if (!nzchar(bind) && nzchar(plat)) {
    bind <- file.path(plat, "bin")
    Sys.setenv(MOSEK_BINDIR = bind)
  }

  if (nzchar(bind)) {
    Sys.setenv(PATH = paste(bind, Sys.getenv("PATH"), sep = .Platform$path.sep))
    Sys.setenv(LD_LIBRARY_PATH = paste(bind, Sys.getenv("LD_LIBRARY_PATH"), sep = .Platform$path.sep))
  }
}

ensure_mosek_env()

# ---- packages ----
packages <- c(
  "data.table", "dplyr", "openxlsx", "peakRAM",
  "cluster", "mclust", "Rfast",
  "SNFtool", "coca", "klic"
)

invisible(lapply(packages, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cran.r-project.org")
  }
}))

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(openxlsx)
  library(peakRAM)
  library(cluster)
  library(mclust)
  library(Rfast)
  library(SNFtool)
  library(coca)
  library(klic)
})

ensure_rmosek <- function() {
  # Install the MOSEK-provided R interface for your MOSEK major version (10.2 here)
  if (!requireNamespace("Rmosek", quietly = TRUE)) {
    install.packages("Rmosek", repos = "https://download.mosek.com/R/10.2")
  }

  ns <- asNamespace("Rmosek")

  # The actual solver entrypoint must exist; otherwise you've likely installed the CRAN helper/meta package.
  if (!exists("mosek", envir = ns, inherits = FALSE)) {
    stop(
      "Rmosek is installed but does not export `mosek()`. ",
      "This usually means you have the CRAN helper/meta package, not the MOSEK solver interface.\n",
      "Install the MOSEK R interface matching your MOSEK major version (e.g. from https://download.mosek.com/R/10.2), ",
      "or use MOSEK's builder script as described in their Rmosek installation docs."
    )
  }

  # Newer Rmosek uses mosek_version(); older scripts sometimes used mosek.version().
  ver <- NA_character_
  ver_err <- NA_character_

  if (exists("mosek_version", envir = ns, inherits = FALSE)) {
    ver <- tryCatch(Rmosek::mosek_version(), error = function(e) { ver_err <<- conditionMessage(e); NA_character_ })
  } else if (exists("mosek.version", envir = ns, inherits = FALSE)) {
    ver <- tryCatch(get("mosek.version", envir = ns)(), error = function(e) { ver_err <<- conditionMessage(e); NA_character_ })
  }

  if (!is.na(ver)) {
    message("Detected MOSEK library via Rmosek: ", ver)
  } else if (nzchar(ver_err)) {
    message("Rmosek loaded, but version query failed: ", ver_err)
    message("Continuing; solver availability will be exercised by the KLIC run.")
  } else {
    message("Rmosek loaded; no version function found (continuing).")
  }

  invisible(TRUE)
}

ensure_rmosek()

# ---- reproducibility ----
RNGversion("4.2.2")
set.seed(123)

algorithm <- "KLIC"

# ---- fixed KLIC settings ----
klic_final_kvals <- c(RNAseq = 4L, CNV = 2L, Methylation = 5L, miRNA = 5L, SNPs = 4L)
klic_globalK <- as.integer(Sys.getenv("KLIC_GLOBALK", "2"))
klic_iter <- as.integer(Sys.getenv("KLIC_ITER", "100"))

maxK <- as.integer(Sys.getenv("KLIC_MAXK", "5"))
B <- as.integer(Sys.getenv("KLIC_B", "250"))
pItem <- as.numeric(Sys.getenv("KLIC_PITEM", "0.8"))

input_path <- "Resources/mm_input.rds"
subsets_path <- "Resources/Sample_perturbations/sample_subsets.tsv.gz"
custom_fn_path <- "Resources/custom_performance_functions.R"

job_out <- Sys.getenv("BENCH_OUT_DIR", "")
if (nzchar(job_out)) {
  out_root <- job_out
} else {
  job_id <- Sys.getenv("SLURM_JOB_ID", "")
  out_root <- file.path("Results", "Performance", "Sample_perturbations", algorithm,
                        if (nzchar(job_id)) paste0("job_", job_id) else "local")
}
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(custom_fn_path)) stop("Missing: ", custom_fn_path)
source(custom_fn_path)  # must define coca_cc_mod()

# ---- inputs ----
input_full <- readRDS(input_path)

ground_truth_labels <- openxlsx::read.xlsx(
  file.path("Resources/algorithms", paste0(algorithm, "_ground_truth_labels.xlsx"))
)
if (!all(c("Sample.ID", "Cluster") %in% colnames(ground_truth_labels))) {
  stop("Ground truth file must contain columns: Sample.ID, Cluster")
}
ground_truth_labels$Sample.ID <- gsub("\\.", "-", ground_truth_labels$Sample.ID)
colnames(ground_truth_labels)[colnames(ground_truth_labels) == "Cluster"] <- "Cluster_gt"

sample_subsets <- data.table::fread(subsets_path)
need_cols <- c("fraction", "replicate", "sample_id")
if (!all(need_cols %in% colnames(sample_subsets))) {
  stop("sample_subsets.tsv.gz must contain columns: fraction, replicate, sample_id")
}

# ---- helpers ----
subset_input_cols <- function(input_obj, keep_ids) {
  out <- input_obj
  for (nm in names(out)) {
    if (!is.null(out[[nm]])) {
      kk <- keep_ids[keep_ids %in% colnames(out[[nm]])]
      out[[nm]] <- out[[nm]][, kk, drop = FALSE]
      if (ncol(out[[nm]]) == 0L || nrow(out[[nm]]) == 0L) out[[nm]] <- NULL
    }
  }
  out
}

align_by_common_samples <- function(X_list, ref_mod = NULL) {
  ids_list <- lapply(X_list, colnames)
  common <- Reduce(intersect, ids_list)
  if (length(common) == 0L) stop("No common sample IDs across modalities.")
  if (!is.null(ref_mod) && ref_mod %in% names(X_list)) {
    common <- common[common %in% colnames(X_list[[ref_mod]])]
  }
  lapply(X_list, function(m) m[, common, drop = FALSE])
}

ari_to_ground_truth <- function(gt_df, pred_df) {
  merged <- merge(gt_df, pred_df, by = "Sample.ID")
  if (nrow(merged) == 0) return(NA_real_)
  mclust::adjustedRandIndex(merged[["Cluster_gt"]], merged[["Cluster_pred"]])
}

# KLIC needs samples x features (transpose from features x samples)
to_sample_feature <- function(mat, make_binary = FALSE) {
  x <- t(as.matrix(mat))
  if (make_binary) {
    x <- apply(x, 2, as.numeric)
    x <- as.matrix(x)
    x[!is.na(x) & x != 0] <- 1
  }
  x
}

run_klic_fixed <- function(X_feat_by_samp) {
  continuous <- c("RNAseq", "CNV", "Methylation", "miRNA")
  categorical <- c("SNPs")
  mods_needed <- c(continuous, categorical)

  if (!all(mods_needed %in% names(X_feat_by_samp))) {
    missing <- setdiff(mods_needed, names(X_feat_by_samp))
    stop("Missing modalities for KLIC: ", paste(missing, collapse = ", "))
  }

  input_dists <- lapply(continuous, function(m) {
    x <- as.matrix(X_feat_by_samp[[m]])
    SNFtool::dist2(x, x)
  })
  names(input_dists) <- continuous
  input_dists[["SNPs"]] <- as.matrix(dist(as.matrix(X_feat_by_samp[["SNPs"]]), method = "binary"))

  nSamples <- nrow(X_feat_by_samp[[continuous[1]]])
  nDatasets <- length(mods_needed)

  allCM <- vector("list", length = maxK - 1L)
  for (k in 2:maxK) {
    CM <- array(NA_real_, dim = c(nSamples, nSamples, nDatasets))
    for (ii in seq_len(nDatasets)) {
      mod <- mods_needed[ii]
      temp_mat <- coca_cc_mod(
        dist         = input_dists[[mod]],
        clMethod     = "hclust",
        B            = B,
        K            = k,
        pItem        = pItem,
        hclustMethod = "average"
      )
      CM[, , ii] <- klic::spectrumShift(temp_mat, verbose = FALSE)
    }
    allCM[[k - 1L]] <- CM
  }

  kvals <- as.integer(klic_final_kvals[mods_needed])
  CM_final <- array(0, dim = c(nSamples, nSamples, nDatasets))
  for (ii in seq_len(nDatasets)) {
    idxInAllCM <- kvals[ii] - 2L + 1L
    CM_final[, , ii] <- allCM[[idxInAllCM]][, , ii]
  }

  local_params <- list(iteration_count = klic_iter, cluster_count = klic_globalK)
  res <- klic::lmkkmeans(CM_final, local_params)

  WKM <- matrix(0, nrow = nSamples, ncol = nSamples)
  for (j in seq_len(nDatasets)) {
    WKM <- WKM + (res$Theta[, j] %*% t(res$Theta[, j])) * CM_final[, , j]
  }

  list(res = res, WKM = WKM)
}

mods_required <- c("RNAseq", "CNV", "Methylation", "miRNA", "SNPs")
non_null_mods <- names(input_full)[!vapply(input_full, is.null, logical(1))]
if (!all(mods_required %in% non_null_mods)) {
  stop("mm_input.rds is missing required modalities for KLIC: ",
       paste(setdiff(mods_required, non_null_mods), collapse = ", "))
}

valid_ids <- colnames(input_full[["RNAseq"]])
sample_subsets <- sample_subsets[sample_id %in% valid_ids]
if (nrow(sample_subsets) == 0) stop("No sample IDs in sample_subsets match the input colnames.")

pairs <- unique(sample_subsets[, .(fraction, replicate)])
setorder(pairs, fraction, replicate)

p_by_mod_full <- sapply(mods_required, function(m) if (!is.null(input_full[[m]])) nrow(input_full[[m]]) else NA_integer_)
p_total_full <- sum(p_by_mod_full, na.rm = TRUE)

perf_rows <- vector("list", nrow(pairs))

for (i in seq_len(nrow(pairs))) {
  frac <- as.integer(pairs$fraction[i])
  repi <- as.integer(pairs$replicate[i])

  message(sprintf("[%s] %s sample subset %d%% rep %d", Sys.time(), algorithm, frac, repi))

  keep_ids <- unique(sample_subsets[fraction == frac & replicate == repi, sample_id])
  input_sub <- subset_input_cols(input_full, keep_ids)

  X_use <- input_sub[mods_required]
  X_use <- Filter(function(m) !is.null(m) && nrow(m) > 0L && ncol(m) > 0L, X_use)
  if (length(X_use) < 5L) {
    next
  }

  X_use <- align_by_common_samples(X_use, ref_mod = "RNAseq")
  sample_ids <- colnames(X_use[["RNAseq"]])
  n_samples <- length(sample_ids)

  X_sf <- list(
    RNAseq      = to_sample_feature(X_use[["RNAseq"]]),
    CNV         = to_sample_feature(X_use[["CNV"]]),
    Methylation = to_sample_feature(X_use[["Methylation"]]),
    miRNA       = to_sample_feature(X_use[["miRNA"]]),
    SNPs        = to_sample_feature(X_use[["SNPs"]], make_binary = TRUE)
  )

  klic_out <- NULL
  gc()
  pr <- peakRAM::peakRAM({
    klic_out <<- run_klic_fixed(X_sf)
  })
  elapsed_s <- as.numeric(pr$Elapsed_Time[1])
  peak_mib <- as.numeric(pr$Peak_RAM_Used_MiB[1])

  res <- klic_out$res
  WKM <- klic_out$WKM

  clusters_df <- data.frame(
    Sample.ID = gsub("\\.", "-", sample_ids),
    Cluster_pred = as.integer(res$clustering),
    stringsAsFactors = FALSE
  )

  ari_gt <- ari_to_ground_truth(ground_truth_labels, clusters_df)

  avg_width <- NA_real_
  if (length(unique(res$clustering)) > 1L) {
    sil <- cluster::silhouette(as.integer(res$clustering), as.dist(1 - WKM))
    avg_width <- summary(sil)$avg.width
  }

  clust_file <- file.path(out_root, sprintf("%s_clusters_samples_%dpct_rep%02d.tsv.gz", algorithm, frac, repi))
  data.table::fwrite(as.data.table(clusters_df), clust_file, sep = "\t", quote = FALSE, compress = "gzip")

  perf_rows[[i]] <- data.table(
    Algorithm = algorithm,
    Sample_Fraction = frac,
    Replicate = repi,
    n_samples = as.integer(n_samples),
    p_total = as.integer(p_total_full),
    p_RNAseq = as.integer(p_by_mod_full[["RNAseq"]]),
    p_CNV = as.integer(p_by_mod_full[["CNV"]]),
    p_Methylation = as.integer(p_by_mod_full[["Methylation"]]),
    p_miRNA = as.integer(p_by_mod_full[["miRNA"]]),
    p_SNPs = as.integer(p_by_mod_full[["SNPs"]]),
    n_clusters = as.integer(klic_globalK),
    KLIC_maxK = as.integer(maxK),
    KLIC_B = as.integer(B),
    KLIC_pItem = as.numeric(pItem),
    KLIC_iter = as.integer(klic_iter),
    KLIC_final_kvals = paste(klic_final_kvals, collapse = ","),
    ARI_to_Ground_Truth = ari_gt,
    Average_Silhouette_Width = avg_width,
    KLIC_Elapsed_Seconds = elapsed_s,
    KLIC_PeakRAM_MiB = peak_mib
  )

  rm(input_sub, X_use, X_sf, klic_out, res, WKM, clusters_df)
  gc()
}

perf_dt <- data.table::rbindlist(perf_rows, use.names = TRUE, fill = TRUE)
perf_file <- file.path(out_root, sprintf("%s_sample_perturbations_performance.tsv", algorithm))
data.table::fwrite(perf_dt, perf_file, sep = "\t", quote = FALSE, na = "NA")

meta_root <- Sys.getenv("BENCH_META_DIR", "")
si_path <- if (nzchar(meta_root)) file.path(meta_root, sprintf("%s_session_info.txt", algorithm))
else file.path(out_root, sprintf("%s_session_info.txt", algorithm))
writeLines(capture.output(sessionInfo()), si_path)

message("Wrote: ", perf_file)

