#!/usr/bin/env Rscript

# ==========================
# KLIC feature perturbations
# ==========================

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
  if (!requireNamespace("Rmosek", quietly = TRUE)) {
    install.packages("Rmosek", repos = "https://download.mosek.com/R/10.2")
  }
  ns <- asNamespace("Rmosek")
  if (!exists("mosek", envir = ns, inherits = FALSE)) {
    stop(
      "Rmosek is installed but does not export `mosek()`. ",
      "This usually means you have the CRAN helper/meta package, not the MOSEK solver interface.\n",
      "Install the MOSEK R interface matching your MOSEK major version (e.g. from https://download.mosek.com/R/10.2)."
    )
  }
  ver <- NA_character_
  ver_err <- NA_character_
  if (exists("mosek_version", envir = ns, inherits = FALSE)) {
    ver <- tryCatch(Rmosek::mosek_version(), error = function(e) { ver_err <<- conditionMessage(e); NA_character_ })
  } else if (exists("mosek.version", envir = ns, inherits = FALSE)) {
    ver <- tryCatch(get("mosek.version", envir = ns)(), error = function(e) { ver_err <<- conditionMessage(e); NA_character_ })
  }
  if (!is.na(ver)) message("Detected MOSEK library via Rmosek: ", ver)
  invisible(TRUE)
}

ensure_rmosek()

# ---- reproducibility ----
RNGversion("4.2.2")
set.seed(123)

algorithm <- "KLIC"
centiles <- c(0.10, 0.20, 0.50, 0.75, 0.90)

# ---- fixed KLIC settings ----
klic_final_kvals <- c(RNAseq = 4L, CNV = 2L, Methylation = 5L, miRNA = 5L, SNPs = 4L)
klic_globalK <- as.integer(Sys.getenv("KLIC_GLOBALK", "2"))
klic_iter <- as.integer(Sys.getenv("KLIC_ITER", "100"))

maxK <- as.integer(Sys.getenv("KLIC_MAXK", "5"))
B <- as.integer(Sys.getenv("KLIC_B", "250"))
pItem <- as.numeric(Sys.getenv("KLIC_PITEM", "0.8"))

input_path <- "Resources/mm_input.rds"
feat_rank_path <- "Resources/Feature_perturbations/feature_rankings.tsv.gz"
custom_fn_path <- "Resources/custom_performance_functions.R"

job_out <- Sys.getenv("BENCH_OUT_DIR", "")
if (nzchar(job_out)) {
  out_root <- job_out
} else {
  job_id <- Sys.getenv("SLURM_JOB_ID", "")
  out_root <- file.path("Results", "Performance", "Feature_perturbations", algorithm,
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

feat_rank <- data.table::fread(feat_rank_path)
if (!all(c("modality", "rank", "feature_id") %in% colnames(feat_rank))) {
  stop("feature_rankings.tsv.gz must contain columns: modality, rank, feature_id")
}

mods_required <- c("RNAseq", "CNV", "Methylation", "miRNA", "SNPs")
modalities_ranked <- intersect(unique(feat_rank$modality), names(input_full))
if (length(modalities_ranked) == 0) stop("No overlap between ranked modalities and input modalities.")

if (!all(mods_required %in% names(input_full))) {
  stop("mm_input.rds is missing required modalities for KLIC: ",
       paste(setdiff(mods_required, names(input_full)), collapse = ", "))
}

# ---- helpers ----
get_top_features <- function(feat_rank_dt, centile) {
  mods <- unique(feat_rank_dt$modality)
  out <- vector("list", length(mods)); names(out) <- mods
  for (mod in mods) {
    mod_dt <- feat_rank_dt[modality == mod][order(rank)]
    top_n <- max(1L, ceiling(centile * nrow(mod_dt)))
    out[[mod]] <- mod_dt$feature_id[seq_len(top_n)]
  }
  out
}

subset_input_features <- function(input_obj, top_features) {
  out <- input_obj
  for (mod in names(top_features)) {
    if (!is.null(out[[mod]])) {
      keep <- intersect(top_features[[mod]], rownames(out[[mod]]))
      if (length(keep) == 0L) {
        stop("After feature subsetting, 0 features remain for modality ", mod,
             ". Likely mismatch between feature_rankings feature_id and mm_input rownames.")
      }
      out[[mod]] <- out[[mod]][keep, , drop = FALSE]
      if (nrow(out[[mod]]) == 0L) out[[mod]] <- NULL
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

run_klic_fixed <- function(X_sxf) {
  continuous <- c("RNAseq", "CNV", "Methylation", "miRNA")
  categorical <- c("SNPs")
  mods_needed <- c(continuous, categorical)

  input_dists <- lapply(continuous, function(m) {
    x <- as.matrix(X_sxf[[m]])
    SNFtool::dist2(x, x)
  })
  names(input_dists) <- continuous
  input_dists[["SNPs"]] <- as.matrix(dist(as.matrix(X_sxf[["SNPs"]]), method = "binary"))

  nSamples <- nrow(X_sxf[[continuous[1]]])
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
    idxInAllCM <- kvals[ii] - 1L  # because allCM index is (k - 1)
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

# ---- main loop ----
perf_rows <- vector("list", length(centiles))

for (i in seq_along(centiles)) {
  centile <- centiles[i]
  pct <- as.integer(round(centile * 100))
  message(sprintf("[%s] %s feature subset %d%%", Sys.time(), algorithm, pct))

  top_feats <- get_top_features(feat_rank[modality %in% modalities_ranked], centile)
  input_sub <- subset_input_features(input_full, top_feats)

  X_use <- input_sub[mods_required]
  X_use <- Filter(function(m) !is.null(m) && nrow(m) > 0L && ncol(m) > 0L, X_use)
  if (length(X_use) < 5L) stop("After subsetting, KLIC has <5 modalities at ", pct, "%.")

  X_use <- align_by_common_samples(X_use, ref_mod = "RNAseq")
  sample_ids <- colnames(X_use[["RNAseq"]])
  n_samples <- length(sample_ids)

  p_by_mod <- sapply(mods_required, function(m) nrow(X_use[[m]]))
  p_total <- sum(p_by_mod)

  X_sxf <- list(
    RNAseq      = to_sample_feature(X_use[["RNAseq"]]),
    CNV         = to_sample_feature(X_use[["CNV"]]),
    Methylation = to_sample_feature(X_use[["Methylation"]]),
    miRNA       = to_sample_feature(X_use[["miRNA"]]),
    SNPs        = to_sample_feature(X_use[["SNPs"]], make_binary = TRUE)
  )

  klic_out <- NULL
  gc()
  pr <- peakRAM::peakRAM({ klic_out <<- run_klic_fixed(X_sxf) })
  elapsed_s <- as.numeric(pr$Elapsed_Time[1])
  peak_mib  <- as.numeric(pr$Peak_RAM_Used_MiB[1])

  res <- klic_out$res
  WKM <- klic_out$WKM

  clusters_df <- data.frame(
    Sample.ID = gsub("\\.", "-", sample_ids),
    Cluster_pred = as.integer(res$clustering),
    stringsAsFactors = FALSE
  )

  ari_gt <- ari_to_ground_truth(ground_truth_labels, clusters_df)

  avg_width <- NA_real_
  if (length(unique(clusters_df$Cluster_pred)) > 1L) {
    sil <- cluster::silhouette(as.integer(clusters_df$Cluster_pred), as.dist(1 - WKM))
    avg_width <- summary(sil)$avg.width
  }

  clust_file <- file.path(out_root, sprintf("%s_clusters_%dpct.tsv.gz", algorithm, pct))
  data.table::fwrite(as.data.table(clusters_df), clust_file, sep = "\t", quote = FALSE, compress = "gzip")

  perf_rows[[i]] <- data.table(
    Algorithm = algorithm,
    Feature_Centile = centile,
    Feature_Percent = pct,
    n_samples = as.integer(n_samples),
    p_total = as.integer(p_total),
    p_RNAseq = as.integer(p_by_mod[["RNAseq"]]),
    p_CNV = as.integer(p_by_mod[["CNV"]]),
    p_Methylation = as.integer(p_by_mod[["Methylation"]]),
    p_miRNA = as.integer(p_by_mod[["miRNA"]]),
    p_SNPs = as.integer(p_by_mod[["SNPs"]]),
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

  rm(input_sub, X_use, X_sxf, klic_out, res, WKM, clusters_df)
  gc()
}

perf_dt <- data.table::rbindlist(perf_rows, use.names = TRUE, fill = TRUE)
perf_file <- file.path(out_root, sprintf("%s_feature_perturbations_performance.tsv", algorithm))
data.table::fwrite(perf_dt, perf_file, sep = "\t", quote = FALSE, na = "NA")

writeLines(capture.output(sessionInfo()),
           file.path(out_root, sprintf("%s_session_info.txt", algorithm)))

message("Wrote: ", perf_file)

