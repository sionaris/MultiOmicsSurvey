#!/usr/bin/env Rscript

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

# Robust error handler (avoid conditionMessage() on a character)
options(error = function() {
  msg <- c("ERROR:", geterrmessage(), "", "TRACEBACK:")
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

optk <- 5L
K_values <- 4L
sdev <- 0.015
beta_var_scale <- 0.5

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
subsets_path <- "Resources/Sample_perturbations/sample_subsets.tsv.gz"

job_out <- Sys.getenv("BENCH_OUT_DIR", "")
if (nzchar(job_out)) {
  out_root <- job_out
} else {
  job_id <- Sys.getenv("SLURM_JOB_ID", "")
  out_root <- file.path("Results", "Performance", "Sample_perturbations", algorithm,
                        if (nzchar(job_id)) paste0("job_", job_id) else "local")
}
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

# ---- merge mode ----
if (identical(mode, "merge")) {
  run_dir <- Sys.getenv("BENCH_RUN_DIR", "")
  if (!nzchar(run_dir)) stop("BENCH_RUN_DIR not set in merge mode")
  task_dirs <- list.dirs(run_dir, full.names = TRUE, recursive = FALSE)
  task_dirs <- task_dirs[grepl("/task_[0-9]+$", task_dirs)]
  if (length(task_dirs) == 0) stop("No task_* directories found under: ", run_dir)

  perf_files <- file.path(task_dirs, "out", sprintf("%s_sample_perturbations_performance.tsv", algorithm))
  perf_files <- perf_files[file.exists(perf_files)]

  # Skip empty/invalid perf files (e.g., created after a crash)
  perf_files <- perf_files[file.info(perf_files)$size > 0]
  if (length(perf_files) == 0) stop("No non-empty per-task performance files found to merge.")

  perf_list <- lapply(perf_files, function(f) {
    dt <- tryCatch(data.table::fread(f), error = function(e) NULL)
    if (is.null(dt) || ncol(dt) == 0) return(NULL)
    dt
  })
  perf_list <- Filter(Negate(is.null), perf_list)
  if (length(perf_list) == 0) stop("All per-task performance files were empty/invalid; nothing to merge.")

  perf_dt <- data.table::rbindlist(perf_list, use.names = TRUE, fill = TRUE)
  if (all(c("Sample_Fraction", "Replicate") %in% colnames(perf_dt))) {
    data.table::setorder(perf_dt, Sample_Fraction, Replicate)
  }

  perf_file <- file.path(out_root, sprintf("%s_sample_perturbations_performance.tsv", algorithm))
  data.table::fwrite(perf_dt, perf_file, sep = "\t", quote = FALSE, na = "NA")

  cl_files <- unlist(lapply(task_dirs, function(td) {
    list.files(file.path(td, "out"), pattern = sprintf("^%s_clusters_samples_.*pct_rep.*\\.tsv\\.gz$", algorithm),
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

# ---- worker/local inputs ----
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

if (!all(mods_required %in% names(input_full))) {
  stop("mm_input.rds is missing required modalities for iClusterBayes: ",
       paste(setdiff(mods_required, names(input_full)), collapse = ", "))
}

ref_ids <- colnames(input_full[["RNAseq"]])
sample_subsets <- sample_subsets[sample_id %in% ref_ids]
if (nrow(sample_subsets) == 0) stop("No sample IDs in sample_subsets match the input colnames.")

pairs <- unique(sample_subsets[, .(fraction, replicate)])
setorder(pairs, fraction, replicate)

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
  mclust::adjustedRandIndex(merged[["Cluster_gt"]], merged[["Cluster_pred"]])
}

prep_dt <- function(mat, type) {
  x <- t(as.matrix(mat))  # samples x features
  x <- apply(x, 2, as.numeric)
  x <- as.matrix(x)

  if (type == "binomial") {
    x[is.na(x)] <- 0
    # Binarise: presence/absence (keeps values in {0,1})
    x[x != 0] <- 1

    # IMPORTANT: iClusterBayes requires binomial columns to have exactly 2 categories.
    # Drop any SNP columns that are constant after binarisation (only 0s or only 1s).
    keep <- apply(x, 2, function(col) length(unique(col)) == 2L)
    x <- x[, keep, drop = FALSE]
  } else {
    # gaussian
    # (leave NAs as-is; iClusterPlus handles internally)
  }

  x
}

run_icb_fixed <- function(X_list) {
  dt1 <- prep_dt(X_list[["SNPs"]], "binomial")
  if (ncol(dt1) == 0L) stop("SNPs binomial matrix has 0 valid columns after filtering to exactly-2-category features.")

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

p_by_mod_full <- sapply(mods_required, function(m) nrow(input_full[[m]]))
p_total_full <- sum(p_by_mod_full, na.rm = TRUE)

# ---- choose which pairs to run ----
pairs_idx <- seq_len(nrow(pairs))
if (identical(mode, "worker")) {
  idx <- as.integer(Sys.getenv("ICB_PAIR_INDEX", "NA"))
  if (is.na(idx) || idx < 1L || idx > nrow(pairs)) {
    stop("Invalid ICB_PAIR_INDEX: ", Sys.getenv("ICB_PAIR_INDEX", ""))
  }
  pairs_idx <- idx
}

perf_rows <- vector("list", length(pairs_idx))

for (jj in seq_along(pairs_idx)) {
  i <- pairs_idx[jj]
  frac <- as.integer(pairs$fraction[i])
  repi <- as.integer(pairs$replicate[i])

  message(sprintf("[%s] %s sample subset %d%% rep %d", Sys.time(), algorithm, frac, repi))

  keep_ids <- unique(sample_subsets[fraction == frac & replicate == repi, sample_id])
  input_sub <- subset_input_cols(input_full, keep_ids)

  X_use <- input_sub[mods_required]
  X_use <- Filter(function(m) !is.null(m) && nrow(m) > 0L && ncol(m) > 0L, X_use)

  # If any required modality is missing, record a skip row (do not error)
  if (length(X_use) < length(mods_required)) {
    perf_rows[[jj]] <- data.table(
      Algorithm = algorithm,
      Sample_Fraction = frac,
      Replicate = repi,
      n_samples = as.integer(0L),
      p_total = as.integer(p_total_full),
      optk = as.integer(optk),
      ARI_to_Ground_Truth = NA_real_,
      iClusterBayes_Elapsed_Seconds = NA_real_,
      iClusterBayes_PeakRAM_MiB = NA_real_,
      iClusterBayes_Status = "SKIP_MISSING_MODALITY",
      iClusterBayes_Error = NA_character_
    )
    rm(input_sub); gc()
    next
  }

  # Align + run, but never crash the whole worker: record failure as a perf row
  run_res <- tryCatch({
    X_use_aligned <- align_modalities(X_use, ref_mod = "RNAseq")
    sample_ids <- colnames(X_use_aligned[["RNAseq"]])
    n_samples <- length(sample_ids)

    icb_out <- NULL
    gc()
    pr <- peakRAM::peakRAM({
      icb_out <<- run_icb_fixed(X_use_aligned)
    })

    elapsed_s <- as.numeric(pr$Elapsed_Time[1])
    peak_mib <- as.numeric(pr$Peak_RAM_Used_MiB[1])

    clusters_df <- data.frame(
      Sample.ID = gsub("\\.", "-", sample_ids),
      Cluster_pred = as.integer(icb_out$clusters),
      stringsAsFactors = FALSE
    )
    ari_gt <- ari_to_ground_truth(ground_truth_labels, clusters_df)

    clust_file <- file.path(out_root, sprintf("%s_clusters_samples_%dpct_rep%02d.tsv.gz", algorithm, frac, repi))
    data.table::fwrite(as.data.table(clusters_df), clust_file, sep = "\t", quote = FALSE, compress = "gzip")

    list(
      ok = TRUE,
      n_samples = n_samples,
      ari = ari_gt,
      elapsed = elapsed_s,
      peak = peak_mib,
      err = NA_character_
    )
  }, error = function(e) {
    list(ok = FALSE, n_samples = NA_integer_, ari = NA_real_, elapsed = NA_real_, peak = NA_real_,
         err = conditionMessage(e))
  })

  if (isTRUE(run_res$ok)) {
    perf_rows[[jj]] <- data.table(
      Algorithm = algorithm,
      Sample_Fraction = frac,
      Replicate = repi,
      n_samples = as.integer(run_res$n_samples),
      p_total = as.integer(p_total_full),
      p_SNPs = as.integer(p_by_mod_full[["SNPs"]]),
      p_RNAseq = as.integer(p_by_mod_full[["RNAseq"]]),
      p_CNV = as.integer(p_by_mod_full[["CNV"]]),
      p_miRNA = as.integer(p_by_mod_full[["miRNA"]]),
      p_Methylation = as.integer(p_by_mod_full[["Methylation"]]),
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
      ARI_to_Ground_Truth = as.numeric(run_res$ari),
      iClusterBayes_Elapsed_Seconds = as.numeric(run_res$elapsed),
      iClusterBayes_PeakRAM_MiB = as.numeric(run_res$peak),
      iClusterBayes_Status = "OK",
      iClusterBayes_Error = NA_character_
    )
  } else {
    # Common cause: SNP columns constant after binarisation -> no valid 2-category binomial features
    perf_rows[[jj]] <- data.table(
      Algorithm = algorithm,
      Sample_Fraction = frac,
      Replicate = repi,
      n_samples = as.integer(NA_integer_),
      p_total = as.integer(p_total_full),
      p_SNPs = as.integer(p_by_mod_full[["SNPs"]]),
      p_RNAseq = as.integer(p_by_mod_full[["RNAseq"]]),
      p_CNV = as.integer(p_by_mod_full[["CNV"]]),
      p_miRNA = as.integer(p_by_mod_full[["miRNA"]]),
      p_Methylation = as.integer(p_by_mod_full[["Methylation"]]),
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
      ARI_to_Ground_Truth = NA_real_,
      iClusterBayes_Elapsed_Seconds = NA_real_,
      iClusterBayes_PeakRAM_MiB = NA_real_,
      iClusterBayes_Status = "FAIL",
      iClusterBayes_Error = as.character(run_res$err)
    )
  }

  rm(input_sub, X_use)
  gc()
}

perf_dt <- data.table::rbindlist(perf_rows, use.names = TRUE, fill = TRUE)
perf_file <- file.path(out_root, sprintf("%s_sample_perturbations_performance.tsv", algorithm))
data.table::fwrite(perf_dt, perf_file, sep = "\t", quote = FALSE, na = "NA")

writeLines(capture.output(sessionInfo()),
           file.path(out_root, sprintf("%s_session_info.txt", algorithm)))

message("Wrote: ", perf_file)

