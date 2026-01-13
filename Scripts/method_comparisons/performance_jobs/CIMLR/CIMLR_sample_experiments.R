#!/usr/bin/env Rscript
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools", repos = "https://cran.r-project.org")
}

if (!requireNamespace("CIMLR", quietly = TRUE)) {
  devtools::install_github("danro9685/CIMLR", ref = "R")
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cran.r-project.org")
}

if (!requireNamespace("M3C", quietly = TRUE)) {
  BiocManager::install("M3C", ask = FALSE, update = FALSE)
}

packages <- c("dplyr", "data.table", "openxlsx", "Rfast", "cluster", "mclust",
              "peakRAM", "Matrix", "parallel")
invisible(lapply(packages, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg, repos = "https://cran.r-project.org")
}))


# --- Fairness: hard-cap implicit threading (BLAS/OpenMP)
Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1"
)

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(parallel)
  library(openxlsx)
  library(CIMLR)
  library(M3C)
  library(Rfast)
  library(cluster)
  library(mclust)
  library(peakRAM)
  library(Matrix)
})

RNGversion("4.2.2")
set.seed(123)

algorithm <- "CIMLR"

optk <- 2L
k_neighbors <- 15L

# --- IMPORTANT: prevent CIMLR_mod internal parallel blow-up (OOM) ---
# CIMLR_mod uses detectCores() * cores.ratio internally. Using SLURM_CPUS_PER_TASK
# here is NOT safe because detectCores() can be much larger than your allocation.
host_cores <- as.integer(parallel::detectCores(logical = FALSE))
if (is.na(host_cores) || host_cores < 1L) host_cores <- 1L
cimlr_cores_ratio <- 1 / host_cores  # ≈1 core inside CIMLR_mod() regardless of Slurm allocation
options(mc.cores = 1L)
# ---------------------------------------------------------------

binary_distance <- "binary"
nonbinary_distance <- "sqeuclidean"

# M3C params (not benchmarked)
m3c_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
m3c_iters <- as.integer(Sys.getenv("M3C_ITERS", "25"))
m3c_repsref <- as.integer(Sys.getenv("M3C_REPSREF", "100"))
m3c_repsreal <- as.integer(Sys.getenv("M3C_REPSREAL", "100"))
m3c_maxK <- 2L
m3c_objective <- "entropy"

input_path <- "Resources/mm_input.rds"
subsets_path <- "Resources/Sample_perturbations/sample_subsets.tsv.gz"
scheme_path <- "Resources/scheme.rds"
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

mode <- Sys.getenv("CIMLR_MODE", "worker")  # "worker" or "merge"

# ---------------- helpers ----------------
if (!file.exists(custom_fn_path)) stop("Missing: ", custom_fn_path)
source(custom_fn_path) # must define CIMLR_mod()

ari_to_ground_truth <- function(gt_df, pred_df) {
  merged <- merge(gt_df, pred_df, by = "Sample.ID")
  if (nrow(merged) == 0) return(NA_real_)
  adjustedRandIndex(merged[["Cluster_gt"]], merged[["Cluster_pred"]])
}

extract_m3c_scores <- function(m3c_obj, k = 2L) {
  sc <- m3c_obj$scores
  if (is.null(sc)) return(list(p = NA_real_, entropy = NA_real_, rcsi = NA_real_))
  sc <- as.data.frame(sc)
  row <- if ("K" %in% colnames(sc)) sc[sc$K == k, , drop = FALSE] else sc[k - 1L, , drop = FALSE]
  if (nrow(row) == 0) row <- sc[1, , drop = FALSE]

  p_col <- if ("NORM_P" %in% colnames(row)) "NORM_P" else if ("MONTECARLO_P" %in% colnames(row)) "MONTECARLO_P" else if ("BETA_P" %in% colnames(row)) "BETA_P" else NA_character_
  e_col <- if ("Entropy" %in% colnames(row)) "Entropy" else {
    cand <- grep("entropy", colnames(row), ignore.case = TRUE, value = TRUE)
    if (length(cand)) cand[1] else NA_character_
  }
  r_col <- if ("RCSI" %in% colnames(row)) "RCSI" else {
    cand <- grep("^RCSI", colnames(row), value = TRUE)
    if (length(cand)) cand[1] else NA_character_
  }

  list(
    p = if (!is.na(p_col)) as.numeric(row[[p_col]][1]) else NA_real_,
    entropy = if (!is.na(e_col)) as.numeric(row[[e_col]][1]) else NA_real_,
    rcsi = if (!is.na(r_col)) as.numeric(row[[r_col]][1]) else NA_real_
  )
}

extract_m3c_clusters <- function(m3c_obj, k = 2L) {
  if (!is.null(m3c_obj$realdataresults) && length(m3c_obj$realdataresults) >= k) {
    ann <- m3c_obj$realdataresults[[k]]$ordered_annotation
    if (!is.null(ann) && "consensuscluster" %in% colnames(ann)) {
      ids <- if ("ID" %in% colnames(ann)) ann$ID else rownames(ann)
      return(data.frame(Sample.ID = ids, Cluster_pred = as.integer(ann$consensuscluster)))
    }
  }
  NULL
}

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

# ---------------- merge mode ----------------
if (identical(tolower(mode), "merge")) {
  run_dir <- Sys.getenv("BENCH_RUN_DIR", dirname(out_root))
  row_files <- list.files(run_dir,
                          pattern = "^perf_row_samples_[0-9]+pct_rep[0-9]+\\.tsv$",
                          recursive = TRUE, full.names = TRUE)
  if (length(row_files) == 0) stop("No perf_row_samples_*.tsv files found under: ", run_dir)

  dt_list <- lapply(row_files, function(f) data.table::fread(f))
  perf_dt <- data.table::rbindlist(dt_list, use.names = TRUE, fill = TRUE)
  if (all(c("Sample_Fraction","Replicate") %in% names(perf_dt))) {
    data.table::setorder(perf_dt, Sample_Fraction, Replicate)
  }

  perf_file <- file.path(out_root, sprintf("%s_sample_perturbations_performance.tsv", algorithm))
  data.table::fwrite(perf_dt, perf_file, sep = "\t", quote = FALSE, na = "NA")

  writeLines(capture.output(sessionInfo()),
             file.path(out_root, sprintf("%s_session_info.txt", algorithm)))

  message("Merged and wrote: ", perf_file)
  quit(save = "no", status = 0)
}

# ---------------- worker mode ----------------
task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID", "")
if (!nzchar(task_id)) stop("Worker mode requires SLURM_ARRAY_TASK_ID (use a job array)")

# load data inside worker
input_full <- readRDS(input_path)

ground_truth_labels <- openxlsx::read.xlsx(
  file.path("Resources/algorithms", paste0(algorithm, "_ground_truth_labels.xlsx"))
)
if (!all(c("Sample.ID", "Cluster") %in% colnames(ground_truth_labels))) {
  stop("Ground truth file must contain columns: Sample.ID, Cluster")
}
ground_truth_labels$Sample.ID <- gsub("\\.", "-", ground_truth_labels$Sample.ID)
colnames(ground_truth_labels)[colnames(ground_truth_labels) == "Cluster"] <- "Cluster_gt"

scheme <- readRDS(scheme_path)
if (is.null(scheme$annCol)) stop("scheme.rds must contain $annCol")
m3c_des <- scheme$annCol
if (is.null(rownames(m3c_des))) stop("scheme$annCol must have rownames = sample IDs")

er_col <- NULL
for (nm in c("ER status", "ER_status", "ER_status ", "ER Status")) {
  if (nm %in% colnames(m3c_des)) { er_col <- nm; break }
}
m3c_des$class <- if (!is.null(er_col)) m3c_des[[er_col]] else NA
m3c_des$ID <- rownames(m3c_des)

sample_subsets <- data.table::fread(subsets_path)
need_cols <- c("fraction", "replicate", "sample_id")
if (!all(need_cols %in% colnames(sample_subsets))) {
  stop("sample_subsets.tsv.gz must contain columns: fraction, replicate, sample_id")
}

non_null_mods <- names(input_full)[!vapply(input_full, is.null, logical(1))]
ref_mod_full <- if ("RNAseq" %in% non_null_mods) "RNAseq" else non_null_mods[1]
valid_ids <- colnames(input_full[[ref_mod_full]])

sample_subsets <- sample_subsets[sample_id %in% valid_ids]
if (nrow(sample_subsets) == 0) stop("No sample IDs in sample_subsets match the input colnames.")

pairs <- unique(sample_subsets[, .(fraction, replicate)])
setorder(pairs, fraction, replicate)

idx <- as.integer(task_id)
if (is.na(idx) || idx < 1L || idx > nrow(pairs)) {
  stop("Invalid SLURM_ARRAY_TASK_ID=", task_id, " for n_pairs=", nrow(pairs))
}

frac <- as.integer(pairs$fraction[idx])
repi <- as.integer(pairs$replicate[idx])

message(sprintf("[%s] %s sample subset %d%% rep %d (task %d/%d)",
                Sys.time(), algorithm, frac, repi, idx, nrow(pairs)))

keep_ids <- unique(sample_subsets[fraction == frac & replicate == repi, sample_id])
input_sub <- subset_input_cols(input_full, keep_ids)

X_use <- input_sub[non_null_mods]
X_use <- Filter(function(m) !is.null(m) && nrow(m) > 0 && ncol(m) > 0, X_use)
if (length(X_use) == 0) stop("All modalities became empty after column subsetting for ", frac, "% rep ", repi)

ref_mod <- if ("RNAseq" %in% names(X_use)) "RNAseq" else names(X_use)[1]
sample_ids <- colnames(X_use[[ref_mod]])
n_samples <- length(sample_ids)

p_by_mod_full <- sapply(non_null_mods, function(m) if (!is.null(input_full[[m]])) nrow(input_full[[m]]) else NA_integer_)
p_total_full <- sum(p_by_mod_full, na.rm = TRUE)

binary_flags <- ifelse(names(X_use) == "SNPs", "Yes", "No")

# --- measure CIMLR_mod only ---
gc()
pr_cimlr <- peakRAM::peakRAM({
  sim <- CIMLR_mod(
    X = X_use,
    c = optk,
    k = k_neighbors,
    cores.ratio = cimlr_cores_ratio,
    binary_flags = binary_flags,
    binary_distance = binary_distance,
    nonbinary_distance = nonbinary_distance
  )
})
cimlr_elapsed <- as.numeric(pr_cimlr$Elapsed_Time[1])
cimlr_peak_mib <- as.numeric(pr_cimlr$Peak_RAM_Used_MiB[1])

Fmat <- sim$F
sim$S <- NULL; sim$LF <- NULL; sim$ydata <- NULL
rm(sim); gc()

if (is.null(Fmat)) stop("CIMLR_mod returned NULL $F")
Fmat <- as.matrix(Fmat)

if (nrow(Fmat) == n_samples && ncol(Fmat) == optk) {
  rownames(Fmat) <- sample_ids
} else if (ncol(Fmat) == n_samples && nrow(Fmat) == optk) {
  Fmat <- t(Fmat)
  rownames(Fmat) <- sample_ids
} else {
  if (nrow(Fmat) == n_samples) rownames(Fmat) <- sample_ids
}

m3c_input <- as.data.frame(t(Fmat))
rownames(m3c_input) <- paste0("t_SNE", seq_len(nrow(m3c_input)))
colnames(m3c_input) <- sample_ids

des_sub <- m3c_des[m3c_des$ID %in% sample_ids, , drop = FALSE]
if (nrow(des_sub) == 0) {
  des_sub <- data.frame(ID = sample_ids, class = NA)
  rownames(des_sub) <- des_sub$ID
}

gc()
m3c_res <- M3C::M3C(
  mydata = m3c_input,
  des = des_sub,
  cores = m3c_cores,
  iters = m3c_iters,
  repsref = m3c_repsref,
  repsreal = m3c_repsreal,
  seed = 123,
  clusteralg = "km",
  maxK = m3c_maxK,
  objective = m3c_objective,
  removeplots = TRUE,
  silent = TRUE
)

m3c_sc <- extract_m3c_scores(m3c_res, k = optk)
clusters_df <- extract_m3c_clusters(m3c_res, k = optk)
rm(m3c_res); gc()

if (is.null(clusters_df) || nrow(clusters_df) == 0) {
  km <- kmeans(Fmat, centers = optk, nstart = 200)
  clusters_df <- data.frame(Sample.ID = rownames(Fmat), Cluster_pred = as.integer(km$cluster))
}

clusters_df$Sample.ID <- gsub("\\.", "-", clusters_df$Sample.ID)
clusters_df <- clusters_df[!is.na(clusters_df$Sample.ID) & nzchar(clusters_df$Sample.ID), , drop = FALSE]

ari_gt <- ari_to_ground_truth(ground_truth_labels, clusters_df)

cl_vec <- clusters_df$Cluster_pred
names(cl_vec) <- clusters_df$Sample.ID
F_for_sil <- Fmat
rownames(F_for_sil) <- gsub("\\.", "-", rownames(F_for_sil))
common <- intersect(rownames(F_for_sil), names(cl_vec))

avg_width <- NA_real_
if (length(common) >= 5L) {
  Demb <- Rfast::Dist(F_for_sil[common, , drop = FALSE], method = "euclidean")
  sil <- silhouette(as.integer(cl_vec[common]), dist(Demb))
  avg_width <- summary(sil)$avg.width
}

clust_file <- file.path(out_root, sprintf("%s_clusters_samples_%dpct_rep%02d.tsv.gz", algorithm, frac, repi))
data.table::fwrite(data.table::as.data.table(clusters_df), clust_file, sep = "\t", quote = FALSE, compress = "gzip")

get_int1 <- function(x) if (is.null(x) || length(x) == 0L || is.na(x)) NA_integer_ else as.integer(x[1])

perf_row <- data.table(
  Algorithm = algorithm,
  Sample_Fraction = frac,
  Replicate = repi,
  n_samples = n_samples,
  p_total = p_total_full,
  p_RNAseq = get_int1(p_by_mod_full[["RNAseq"]]),
  p_CNV = get_int1(p_by_mod_full[["CNV"]]),
  p_Methylation = get_int1(p_by_mod_full[["Methylation"]]),
  p_miRNA = get_int1(p_by_mod_full[["miRNA"]]),
  p_SNPs = get_int1(p_by_mod_full[["SNPs"]]),
  n_clusters = optk,
  CIMLR_K = k_neighbors,
  CIMLR_cores_ratio = cimlr_cores_ratio,
  M3C_NORM_P = m3c_sc$p,
  M3C_Entropy = m3c_sc$entropy,
  M3C_RCSI = m3c_sc$rcsi,
  ARI_to_Ground_Truth = ari_gt,
  Average_Silhouette_Width = avg_width,
  CIMLR_Elapsed_Seconds = cimlr_elapsed,
  CIMLR_PeakRAM_MiB = cimlr_peak_mib
)

row_file <- file.path(out_root, sprintf("perf_row_samples_%dpct_rep%02d.tsv", frac, repi))
data.table::fwrite(perf_row, row_file, sep = "\t", quote = FALSE, na = "NA")

message("Wrote partial row: ", row_file)
