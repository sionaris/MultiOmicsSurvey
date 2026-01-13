#!/usr/bin/env Rscript

# ==========================================================
# iClusterBayes benchmarking
# Enforce user-writable R library (R_LIBS_USER) for installs
# and serialise any installs to avoid array-task clashes.
# ==========================================================

mode <- Sys.getenv("ICB_MODE", "local")  # local | worker | merge

user_lib <- Sys.getenv("R_LIBS_USER", "")
if (!nzchar(user_lib)) {
  user_lib <- file.path(getwd(), "R", "library")
  Sys.setenv(R_LIBS_USER = user_lib)
}
dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(user_lib, .libPaths()))

options(repos = c(CRAN = "https://cran.r-project.org"))
options(Ncpus = 1L)

bench_meta <- Sys.getenv("BENCH_META_DIR", "")
if (nzchar(bench_meta)) dir.create(bench_meta, recursive = TRUE, showWarnings = FALSE)
err_file <- if (nzchar(bench_meta)) file.path(bench_meta, "error.txt") else ""

options(error = function() {
  msg <- c("ERROR:", conditionMessage(geterrmessage()), "\nTRACEBACK:")
  tb <- capture.output(traceback(2))
  if (nzchar(err_file)) writeLines(c(msg, tb), err_file)
  quit(save = "no", status = 1)
})

with_install_lock <- function(expr) {
  lockdir <- file.path(user_lib, ".pkg_install_lockdir")
  repeat {
    if (dir.create(lockdir, showWarnings = FALSE)) break
    Sys.sleep(2)
  }
  on.exit(unlink(lockdir, recursive = TRUE, force = TRUE), add = TRUE)
  force(expr)
}

ensure_cran_pkg <- function(pkg) {
  if (requireNamespace(pkg, quietly = TRUE)) return(invisible(TRUE))
  with_install_lock({
    install.packages(pkg, lib = user_lib, repos = "https://cran.r-project.org")
  })
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Failed to install CRAN package: ", pkg)
  invisible(TRUE)
}

ensure_bioc_pkg <- function(pkg) {
  ensure_cran_pkg("BiocManager")
  if (requireNamespace(pkg, quietly = TRUE)) return(invisible(TRUE))
  with_install_lock({
    BiocManager::install(pkg, lib = user_lib, ask = FALSE, update = FALSE)
  })
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Failed to install Bioconductor package: ", pkg)
  invisible(TRUE)
}

# ---- packages ----
ensure_bioc_pkg("iClusterPlus")
packages <- c("dplyr", "data.table", "openxlsx", "Rfast", "cluster", "mclust", "peakRAM", "R.utils")
invisible(lapply(packages, ensure_cran_pkg))

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(openxlsx)
  library(Rfast)
  library(cluster)
  library(mclust)
  library(peakRAM)
  library(iClusterPlus)
})

RNGversion("4.2.2")
set.seed(123)

algorithm <- "iClusterBayes"
centiles <- c(0.10, 0.20, 0.50, 0.75, 0.90)

# ---- fixed optimal settings you provided ----
optk <- 5L
K_values <- 4L  # because number of clusters = K + 1
sdev <- 0.015
beta_var_scale <- 0.5

# ---- fixed hyperparameters ----
thin <- 3
pp_cutoff <- 0.5
n_burnin <- 1200
n_draw <- 1800
prior_gamma <- c(0.5, 0.5, 0.5, 0.5, 0.5)

cpus <- as.integer(Sys.getenv("ICB_CPUS", "1"))
if (is.na(cpus) || cpus < 1L) cpus <- 1L

data_types <- c("binomial", "gaussian", "gaussian", "gaussian", "gaussian")
mods_required <- c("SNPs", "RNAseq", "CNV", "miRNA", "Methylation")

input_path <- "Resources/mm_input.rds"

job_out <- Sys.getenv("BENCH_OUT_DIR", "")
if (nzchar(job_out)) {
  out_root <- job_out
} else {
  job_id <- Sys.getenv("SLURM_JOB_ID", "")
  out_root <- file.path("Results", "Performance", "Feature_perturbations", algorithm,
                        if (nzchar(job_id)) paste0("job_", job_id) else "local")
}
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

# ---- merge mode: aggregate task outputs into RUN_DIR/out ----
if (identical(mode, "merge")) {
  run_dir <- Sys.getenv("BENCH_RUN_DIR", "")
  if (!nzchar(run_dir)) stop("BENCH_RUN_DIR not set in merge mode")
  task_dirs <- list.dirs(run_dir, full.names = TRUE, recursive = FALSE)
  task_dirs <- task_dirs[grepl("/task_[0-9]+$", task_dirs)]
  if (length(task_dirs) == 0) stop("No task_* directories found under: ", run_dir)

  perf_files <- file.path(task_dirs, "out", sprintf("%s_feature_perturbations_performance.tsv", algorithm))
  perf_files <- perf_files[file.exists(perf_files)]

  if (length(perf_files) == 0) stop("No per-task performance files found to merge.")

  perf_list <- lapply(perf_files, function(f) fread(f))
  perf_dt <- rbindlist(perf_list, use.names = TRUE, fill = TRUE)
  if ("Feature_Percent" %in% colnames(perf_dt)) setorder(perf_dt, Feature_Percent)

  perf_file <- file.path(out_root, sprintf("%s_feature_perturbations_performance.tsv", algorithm))
  fwrite(perf_dt, perf_file, sep = "\t", quote = FALSE, na = "NA")

  # copy cluster files up (best-effort)
  cl_files <- unlist(lapply(task_dirs, function(td) {
    list.files(file.path(td, "out"), pattern = sprintf("^%s_clusters_.*pct\\.tsv\\.gz$", algorithm),
               full.names = TRUE)
  }))
  if (length(cl_files) > 0) {
    file.copy(cl_files, out_root, overwrite = TRUE)
  }

  writeLines(capture.output(sessionInfo()),
             file.path(out_root, sprintf("%s_session_info.txt", algorithm)))

  message("Wrote: ", perf_file)
  quit(save = "no", status = 0)
}

# ---- load inputs (worker/local) ----
input_full <- readRDS(input_path)

ground_truth_labels <- openxlsx::read.xlsx(
  file.path("Resources/algorithms", paste0(algorithm, "_ground_truth_labels.xlsx"))
)
if (!all(c("Sample.ID", "Cluster") %in% colnames(ground_truth_labels))) {
  stop("Ground truth file must contain columns: Sample.ID, Cluster")
}
ground_truth_labels$Sample.ID <- gsub("\\.", "-", ground_truth_labels$Sample.ID)
colnames(ground_truth_labels)[colnames(ground_truth_labels) == "Cluster"] <- "Cluster_gt"

feat_rank <- fread("Resources/Feature_perturbations/feature_rankings.tsv.gz")
if (!all(c("modality", "rank", "feature_id") %in% colnames(feat_rank))) {
  stop("feature_rankings.tsv.gz must contain columns: modality, rank, feature_id")
}

if (!all(mods_required %in% names(input_full))) {
  stop("mm_input.rds is missing required modalities for iClusterBayes: ",
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
      if (length(keep) == 0L && nrow(out[[mod]]) > 0L) keep <- rownames(out[[mod]])[1]
      out[[mod]] <- out[[mod]][keep, , drop = FALSE]
      if (nrow(out[[mod]]) == 0L) out[[mod]] <- NULL
    }
  }
  out
}

align_modalities <- function(X_list, ref_mod = NULL) {
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
  adjustedRandIndex(merged[["Cluster_gt"]], merged[["Cluster_pred"]])
}

prep_dt <- function(mat, type) {
  x <- t(as.matrix(mat))  # samples x features
  if (type == "binomial") {
    x <- apply(x, 2, as.numeric)
    x <- as.matrix(x)
    x[is.na(x)] <- 0
    x[x != 0] <- 1
    keep <- apply(x, 2, function(col) length(unique(col)) <= 2)
    x <- x[, keep, drop = FALSE]
  } else {
    x <- apply(x, 2, as.numeric)
    x <- as.matrix(x)
  }
  x
}

run_icb_fixed <- function(X_list) {
  dt1 <- prep_dt(X_list[["SNPs"]], "binomial")
  dt2 <- prep_dt(X_list[["RNAseq"]], "gaussian")
  dt3 <- prep_dt(X_list[["CNV"]], "gaussian")
  dt4 <- prep_dt(X_list[["miRNA"]], "gaussian")
  dt5 <- prep_dt(X_list[["Methylation"]], "gaussian")

  tune_results <- iClusterPlus::tune.iClusterBayes(
    dt1 = dt1, dt2 = dt2, dt3 = dt3, dt4 = dt4, dt5 = dt5,
    type = data_types,
    K = K_values,
    sdev = sdev,
    beta.var = beta_var_scale,
    n.burnin = n_burnin,
    n.draw = n_draw,
    thin = thin,
    pp.cutoff = pp_cutoff,
    prior.gamma = prior_gamma,
    cpus = cpus
  )

  fits <- tune_results$fit
  if (is.null(fits) || length(fits) < 1L) stop("tune.iClusterBayes returned no fits")
  optimal.fit <- fits[[1]]
  clusters <- optimal.fit$clusters
  if (is.null(clusters) || length(clusters) == 0L) stop("No clusters returned by iClusterBayes fit")

  list(clusters = clusters, tune_results = tune_results)
}

# ---- select which centiles to run ----
centiles_to_run <- centiles
if (identical(mode, "worker")) {
  idx <- as.integer(Sys.getenv("ICB_CENTILE_INDEX", "NA"))
  if (is.na(idx) || idx < 1L || idx > length(centiles)) {
    stop("Invalid ICB_CENTILE_INDEX: ", Sys.getenv("ICB_CENTILE_INDEX", ""))
  }
  centiles_to_run <- centiles[idx]
}

perf_rows <- vector("list", length(centiles_to_run))

for (ii in seq_along(centiles_to_run)) {
  centile <- centiles_to_run[ii]
  pct <- as.integer(round(centile * 100))
  message(sprintf("[%s] %s feature subset %d%%", Sys.time(), algorithm, pct))

  top_feats <- get_top_features(feat_rank[modality %in% mods_required], centile)
  input_sub <- subset_input_features(input_full, top_feats)

  X_use <- input_sub[mods_required]
  X_use <- Filter(function(m) !is.null(m) && nrow(m) > 0L && ncol(m) > 0L, X_use)
  if (length(X_use) < length(mods_required)) {
    stop("Some required modalities became NULL after feature subsetting at ", pct, "%")
  }

  X_use <- align_modalities(X_use, ref_mod = "RNAseq")
  sample_ids <- colnames(X_use[["RNAseq"]])
  n_samples <- length(sample_ids)

  p_by_mod <- sapply(mods_required, function(m) nrow(X_use[[m]]))
  p_total <- sum(p_by_mod, na.rm = TRUE)

  icb_out <- NULL
  gc()
  pr <- peakRAM::peakRAM({
    icb_out <<- run_icb_fixed(X_use)
  })
  icb_elapsed <- as.numeric(pr$Elapsed_Time[1])
  icb_peak_mib <- as.numeric(pr$Peak_RAM_Used_MiB[1])

  clusters_df <- data.frame(
    Sample.ID = gsub("\\.", "-", sample_ids),
    Cluster_pred = as.integer(icb_out$clusters),
    stringsAsFactors = FALSE
  )

  ari_gt <- ari_to_ground_truth(ground_truth_labels, clusters_df)
  avg_width <- NA_real_

  clust_file <- file.path(out_root, sprintf("%s_clusters_%dpct.tsv.gz", algorithm, pct))
  fwrite(as.data.table(clusters_df), clust_file, sep = "\t", quote = FALSE, compress = "gzip")

  perf_rows[[ii]] <- data.table(
    Algorithm = algorithm,
    Feature_Centile = centile,
    Feature_Percent = pct,
    n_samples = as.integer(n_samples),
    p_total = as.integer(p_total),
    p_SNPs = as.integer(p_by_mod[["SNPs"]]),
    p_RNAseq = as.integer(p_by_mod[["RNAseq"]]),
    p_CNV = as.integer(p_by_mod[["CNV"]]),
    p_miRNA = as.integer(p_by_mod[["miRNA"]]),
    p_Methylation = as.integer(p_by_mod[["Methylation"]]),
    optk = as.integer(optk),
    K_values = as.integer(K_values),
    sdev = as.numeric(sdev),
    beta_var_scale = as.numeric(beta_var_scale),
    thin = as.integer(thin),
    pp_cutoff = as.numeric(pp_cutoff),
    n_burnin = as.integer(n_burnin),
    n_draw = as.integer(n_draw),
    prior_gamma = paste(prior_gamma, collapse = ","),
    cpus = as.integer(cpus),
    ARI_to_Ground_Truth = ari_gt,
    Average_Silhouette_Width = avg_width,
    iClusterBayes_Elapsed_Seconds = icb_elapsed,
    iClusterBayes_PeakRAM_MiB = icb_peak_mib
  )

  rm(input_sub, X_use, icb_out, clusters_df)
  gc()
}

perf_dt <- rbindlist(perf_rows, use.names = TRUE, fill = TRUE)
perf_file <- file.path(out_root, sprintf("%s_feature_perturbations_performance.tsv", algorithm))
fwrite(perf_dt, perf_file, sep = "\t", quote = FALSE, na = "NA")

writeLines(capture.output(sessionInfo()),
           file.path(out_root, sprintf("%s_session_info.txt", algorithm)))

message("Wrote: ", perf_file)

