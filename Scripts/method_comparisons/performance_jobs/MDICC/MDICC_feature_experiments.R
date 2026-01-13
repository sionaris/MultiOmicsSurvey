#!/usr/bin/env Rscript

# ==========================
# MDICC feature perturbations
# ==========================

packages <- c(
  "data.table","dplyr","openxlsx","peakRAM",
  "cluster","mclust","Rfast","SNFtool",
  "reticulate","Rcpp","parallel","Matrix","methods"
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
  library(reticulate)
  library(Matrix)
  library(methods)
})

wd <- Sys.getenv("SLURM_SUBMIT_DIR", getwd())
setwd(wd)

RNGversion("4.2.2")
set.seed(123)

algorithm <- "MDICC"
centiles <- c(0.10, 0.20, 0.50, 0.75, 0.90)

# Fixed hyperparameters
aff_matrix_neighbors <- 18L
cc_fixed <- 6L
k2_fixed <- 43L
k3_fixed <- 2L

input_path <- "Resources/mm_input.rds"
feat_rank_path <- "Resources/Feature_perturbations/feature_rankings.tsv.gz"
gt_path <- file.path("Resources/algorithms", paste0(algorithm, "_ground_truth_labels.xlsx"))

mdicc_dir <- file.path(getwd(), "Resources", "MDICC")
netfusion_r <- file.path(mdicc_dir, "NetworkFusion.R")
projsplx_c <- file.path(mdicc_dir, "projsplx_R.c")
projsplx_so <- file.path(mdicc_dir, "projsplx_R.so")

py_local_aff <- file.path(mdicc_dir, "LocalAffinityMatrix.py")
py_score     <- file.path(mdicc_dir, "score.py")
py_label     <- file.path(mdicc_dir, "label.py")

out_root <- Sys.getenv("BENCH_OUT_DIR", "")
if (!nzchar(out_root)) {
  job_id <- Sys.getenv("SLURM_JOB_ID", "")
  out_root <- file.path("Results","Performance","Feature_perturbations", algorithm,
                        if (nzchar(job_id)) paste0("job_", job_id) else "local")
}
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

# ---- helpers ----
ari_to_ground_truth <- function(gt_df, pred_df) {
  merged <- merge(gt_df, pred_df, by = "Sample.ID")
  if (nrow(merged) == 0) return(NA_real_)
  mclust::adjustedRandIndex(merged[["Cluster_gt"]], merged[["Cluster_pred"]])
}

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
      if (length(keep) == 0L && nrow(out[[mod]]) > 0L) keep <- rownames(out[[mod]])[1]
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

# MDICC expects samples x features
to_sample_feature <- function(mat, make_binary = FALSE) {
  x <- t(as.matrix(mat))
  if (make_binary) {
    x <- apply(x, 2, as.numeric)
    x <- as.matrix(x)
    x[x != 0] <- 1
  }
  x
}

get_mdicc_field <- function(obj, name) {
  if (methods::is(obj, "S4")) {
    if (name %in% methods::slotNames(obj)) return(methods::slot(obj, name))
    stop("MDICC returned S4 but has no slot '", name, "'. Slots: ",
         paste(methods::slotNames(obj), collapse = ", "))
  }
  if (is.list(obj) && !is.null(obj[[name]])) return(obj[[name]])
  out <- try(obj[[name]], silent = TRUE)
  if (!inherits(out, "try-error") && !is.null(out)) return(out)
  stop("Cannot extract field '", name, "' from class: ", paste(class(obj), collapse = ", "))
}

same_conda_prefix <- function(a, b) {
  a <- normalizePath(a, winslash = "/", mustWork = FALSE)
  b <- normalizePath(b, winslash = "/", mustWork = FALSE)
  # compare up to .../envs/<envname>
  chop <- function(x) {
    parts <- strsplit(x, "/", fixed = TRUE)[[1]]
    if (!"envs" %in% parts) return(dirname(dirname(x)))
    i <- match("envs", parts)
    paste(parts[seq_len(min(length(parts), i + 1L))], collapse = "/")
  }
  identical(chop(a), chop(b))
}

ensure_python <- function() {
  expected <- Sys.getenv("RETICULATE_PYTHON", Sys.getenv("MDICC_PYTHON",""))
  if (!nzchar(expected)) expected <- Sys.which("python")
  if (!nzchar(expected) || !file.exists(expected)) stop("Python not found. Set RETICULATE_PYTHON in the .sh.")
  expected <- normalizePath(expected, winslash = "/", mustWork = TRUE)

  # Only request python if it isn't initialised yet
  if (!reticulate::py_available(initialize = FALSE)) {
    reticulate::use_python(expected, required = TRUE)
  }

  cfg <- reticulate::py_config()
  message("---- reticulate diagnostics ----")
  message("reticulate::py_config() python: ", cfg$python)
  message("reticulate::py_config() version: ", cfg$version)
  message("Expected python: ", expected)
  message("--------------------------------")

  if (!same_conda_prefix(cfg$python, expected)) {
    warning("reticulate initialised a different python than expected.\n",
            "cfg$python=", cfg$python, "\nexpected=", expected,
            "\nFix by exporting RETICULATE_PYTHON in the job .sh before running R.")
  }

  needed <- c("numpy","pandas","scipy","sklearn")
  missing <- needed[!vapply(needed, reticulate::py_module_available, logical(1))]
  if (length(missing) > 0) {
    stop("Missing python modules under reticulate's python: ",
         paste(missing, collapse = ", "),
         "\n(Your conda env may have them, but reticulate is not using that python.)")
  }
  invisible(TRUE)
}

ensure_projsplx <- function() {
  if (file.exists(projsplx_so)) return(invisible(TRUE))
  if (!file.exists(projsplx_c)) stop("Missing C source: ", projsplx_c)

  owd <- getwd()
  setwd(mdicc_dir)
  on.exit(setwd(owd), add = TRUE)

  cmd <- sprintf('R CMD SHLIB "%s" -o "%s"', basename(projsplx_c), basename(projsplx_so))
  message("Compiling projsplx_R.so via: ", cmd)
  rc <- system(cmd)
  if (rc != 0 || !file.exists(projsplx_so)) {
    stop("Failed to compile projsplx_R.so. If you see Windows/COFF relocation errors, replace projsplx_R.c with the Linux source.")
  }
  invisible(TRUE)
}

fallback_testaff <- function() {
  testaff <<- function(D, K) {
    D <- as.matrix(D)
    diag(D) <- 0
    SNFtool::affinityMatrix(D, K = as.integer(K), alpha = 0.5)
  }
}

fallback_MDICClabel <- function() {
  MDICClabel <<- function(S, k) {
    S <- as.matrix(S)
    S[is.na(S)] <- 0
    diag(S) <- 0
    k <- as.integer(k)

    d <- rowSums(S)
    d[d <= 0] <- .Machine$double.eps
    DinvSqrt <- diag(1 / sqrt(d))
    L <- diag(nrow(S)) - DinvSqrt %*% S %*% DinvSqrt

    ev <- eigen(L, symmetric = TRUE)
    U <- ev$vectors[, seq_len(k), drop = FALSE]
    U <- U / sqrt(rowSums(U^2) + 1e-12)

    km <- stats::kmeans(U, centers = k, nstart = 20, iter.max = 200)
    as.integer(km$cluster - 1L)  # 0-based to match your +1 later
  }
}

load_mdicc_stack <- function() {
  if (!file.exists(netfusion_r)) stop("Missing: ", netfusion_r)
  for (ff in c(py_local_aff, py_score, py_label)) {
    if (!file.exists(ff)) stop("Missing: ", ff)
  }

  ensure_python()
  ensure_projsplx()
  dyn.load(projsplx_so)

  source(netfusion_r, local = FALSE)

  # These may or may not export what we need; keep for parity with original stack
  reticulate::source_python(py_local_aff, envir = globalenv())
  reticulate::source_python(py_score,     envir = globalenv())
  reticulate::source_python(py_label,     envir = globalenv())

  if (!exists("MDICC", mode = "function", inherits = TRUE)) {
    stop("After sourcing NetworkFusion.R, function MDICC() is missing. Your Resources/MDICC/NetworkFusion.R is not the expected MDICC implementation.")
  }

  if (!exists("testaff", mode = "function", inherits = TRUE)) {
    message("testaff() missing; defining SNFtool-based fallback affinity builder.")
    fallback_testaff()
  }

  if (!exists("MDICClabel", mode = "function", inherits = TRUE)) {
    message("MDICClabel() missing; defining spectral-clustering fallback labeler.")
    fallback_MDICClabel()
  }

  invisible(TRUE)
}

# ---- load MDICC code once ----
load_mdicc_stack()

# ---- inputs ----
input_full <- readRDS(input_path)

ground_truth_labels <- openxlsx::read.xlsx(gt_path)
if (!all(c("Sample.ID","Cluster") %in% colnames(ground_truth_labels))) {
  stop("Ground truth must have Sample.ID, Cluster")
}
ground_truth_labels$Sample.ID <- gsub("\\.", "-", ground_truth_labels$Sample.ID)
colnames(ground_truth_labels)[colnames(ground_truth_labels) == "Cluster"] <- "Cluster_gt"

feat_rank <- data.table::fread(feat_rank_path)

mods_required <- c("RNAseq","CNV","Methylation","miRNA","SNPs")
if (!all(mods_required %in% names(input_full))) {
  stop("mm_input.rds missing required modalities: ", paste(setdiff(mods_required, names(input_full)), collapse=", "))
}

perf_rows <- vector("list", length(centiles))

run_mdicc_fixed <- function(X_sf) {
  continuous <- c("RNAseq","CNV","Methylation","miRNA")
  categorical <- c("SNPs")

  input_dists <- lapply(continuous, function(m) {
    x <- as.matrix(X_sf[[m]])
    SNFtool::dist2(x, x)
  })
  names(input_dists) <- continuous
  input_dists[["SNPs"]] <- as.matrix(dist(as.matrix(X_sf[["SNPs"]]), method = "binary"))

  aff_input <- list()
  for (m in c(continuous, categorical)) {
    aff_input[[m]] <- testaff(as.matrix(input_dists[[m]]), aff_matrix_neighbors)
  }

  res <- suppressMessages(MDICC(aff_input, c = cc_fixed, k = k2_fixed))

  S <- as.matrix(get_mdicc_field(res, "S"))
  label_vec <- MDICClabel(S, k3_fixed)

  list(S = S, label = as.integer(label_vec), kernel_weights = get_mdicc_field(res, "kernel_weights"))
}

for (i in seq_along(centiles)) {
  centile <- centiles[i]
  pct <- as.integer(round(centile * 100))
  message(sprintf("[%s] %s feature subset %d%%", Sys.time(), algorithm, pct))

  top_feats <- get_top_features(feat_rank[modality %in% mods_required], centile)
  input_sub <- subset_input_features(input_full, top_feats)

  X_use <- input_sub[mods_required]
  X_use <- Filter(function(m) !is.null(m) && nrow(m) > 0L && ncol(m) > 0L, X_use)
  if (length(X_use) < 5L) stop("MDICC requires 5 modalities; got: ", paste(names(X_use), collapse=", "))

  X_use <- align_by_common_samples(X_use, ref_mod = "RNAseq")
  sample_ids <- colnames(X_use[["RNAseq"]])
  n_samples <- length(sample_ids)

  p_by_mod <- sapply(mods_required, function(m) if (!is.null(input_sub[[m]])) nrow(input_sub[[m]]) else NA_integer_)
  p_total <- sum(p_by_mod, na.rm = TRUE)

  X_sf <- list(
    RNAseq      = to_sample_feature(X_use[["RNAseq"]]),
    CNV         = to_sample_feature(X_use[["CNV"]]),
    Methylation = to_sample_feature(X_use[["Methylation"]]),
    miRNA       = to_sample_feature(X_use[["miRNA"]]),
    SNPs        = to_sample_feature(X_use[["SNPs"]], make_binary = TRUE)
  )

  mdicc_out <- NULL
  gc()
  pr <- peakRAM::peakRAM({
    mdicc_out <<- run_mdicc_fixed(X_sf)
  })
  elapsed_s <- as.numeric(pr$Elapsed_Time[1])
  peak_mib <- as.numeric(pr$Peak_RAM_Used_MiB[1])

  S <- mdicc_out$S
  labels <- mdicc_out$label

  clusters_df <- data.frame(
    Sample.ID = gsub("\\.", "-", sample_ids),
    Cluster_pred = as.integer(labels) + 1L,
    stringsAsFactors = FALSE
  )

  ari_gt <- ari_to_ground_truth(ground_truth_labels, clusters_df)

  avg_width <- NA_real_
  if (length(unique(clusters_df$Cluster_pred)) > 1L) {
    sil <- cluster::silhouette(as.integer(clusters_df$Cluster_pred), as.dist(1 - S))
    avg_width <- summary(sil)$avg.width
  }

  clust_file <- file.path(out_root, sprintf("%s_clusters_%dpct.tsv.gz", algorithm, pct))
  data.table::fwrite(as.data.table(clusters_df), clust_file, sep="\t", quote=FALSE, compress="gzip")

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
    MDICC_aff_neighbors = as.integer(aff_matrix_neighbors),
    MDICC_c = as.integer(cc_fixed),
    MDICC_k2 = as.integer(k2_fixed),
    n_clusters = as.integer(k3_fixed),
    ARI_to_Ground_Truth = ari_gt,
    Average_Silhouette_Width = avg_width,
    MDICC_Elapsed_Seconds = elapsed_s,
    MDICC_PeakRAM_MiB = peak_mib
  )

  rm(input_sub, X_use, X_sf, mdicc_out, S, labels, clusters_df)
  gc()
}

perf_dt <- data.table::rbindlist(perf_rows, use.names=TRUE, fill=TRUE)
perf_file <- file.path(out_root, sprintf("%s_feature_perturbations_performance.tsv", algorithm))
data.table::fwrite(perf_dt, perf_file, sep="\t", quote=FALSE, na="NA")

writeLines(capture.output(sessionInfo()),
           file.path(out_root, sprintf("%s_session_info.txt", algorithm)))

message("Wrote: ", perf_file)

