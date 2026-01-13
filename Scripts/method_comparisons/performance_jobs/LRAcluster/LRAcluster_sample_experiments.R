#!/usr/bin/env Rscript

# ---- package install checks ----
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cran.r-project.org")
}
if (!requireNamespace("M3C", quietly = TRUE)) {
  BiocManager::install("M3C", ask = FALSE, update = FALSE)
}

packages <- c("dplyr","data.table","openxlsx","Rfast","cluster","mclust","peakRAM","Matrix","parallel")
invisible(lapply(packages, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg, repos = "https://cran.r-project.org")
}))

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(openxlsx)
  library(Rfast)
  library(cluster)
  library(mclust)
  library(peakRAM)
  library(M3C)
})

wd <- Sys.getenv("SLURM_SUBMIT_DIR", getwd())
setwd(wd)

RNGversion("4.2.2")
set.seed(123)

algorithm <- "LRAcluster"

optr_dim <- 7L
optk <- 4L

# M3C not benchmarked
m3c_cores <- as.integer(Sys.getenv("M3C_CORES", "1"))
if (is.na(m3c_cores) || m3c_cores < 1L) m3c_cores <- 1L
m3c_iters    <- as.integer(Sys.getenv("M3C_ITERS", "100"))
m3c_repsref  <- as.integer(Sys.getenv("M3C_REPSREF", "250"))
m3c_repsreal <- as.integer(Sys.getenv("M3C_REPSREAL", "250"))
m3c_maxK     <- optk

input_path <- "Resources/mm_input.rds"
subsets_path <- "Resources/Sample_perturbations/sample_subsets.tsv.gz"
scheme_path <- "Resources/scheme.rds"
gt_path <- file.path("Resources/algorithms", paste0(algorithm, "_ground_truth_labels.xlsx"))
lra_source_path <- "Resources/full_LRAcluster_source.R"

mode <- Sys.getenv("LRA_MODE", "worker")  # worker|merge

job_out <- Sys.getenv("BENCH_OUT_DIR", "")
if (nzchar(job_out)) {
  out_root <- job_out
} else {
  job_id <- Sys.getenv("SLURM_JOB_ID", "")
  out_root <- file.path("Results", "Performance", "Sample_perturbations", algorithm,
                        if (nzchar(job_id)) paste0("job_", job_id) else "local")
}
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

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

align_modalities <- function(X_list, ref_mod = NULL) {
  ids_list <- lapply(X_list, colnames)
  common <- Reduce(intersect, ids_list)
  if (length(common) == 0L) stop("No common sample IDs across modalities.")
  if (!is.null(ref_mod) && ref_mod %in% names(X_list)) {
    common <- common[common %in% colnames(X_list[[ref_mod]])]
  }
  X_list <- lapply(X_list, function(m) m[, common, drop = FALSE])
  X_list
}

ari_to_ground_truth <- function(gt_df, pred_df) {
  merged <- merge(gt_df, pred_df, by = "Sample.ID")
  if (nrow(merged) == 0) return(NA_real_)
  mclust::adjustedRandIndex(merged[["Cluster_gt"]], merged[["Cluster_pred"]])
}

extract_m3c_clusters <- function(m3c_obj, k = 4L) {
  if (!is.null(m3c_obj$realdataresults) && length(m3c_obj$realdataresults) >= k) {
    ann <- m3c_obj$realdataresults[[k]]$ordered_annotation
    if (!is.null(ann) && "consensuscluster" %in% colnames(ann)) {
      ids <- if ("ID" %in% colnames(ann)) ann$ID else rownames(ann)
      return(data.frame(Sample.ID = ids, Cluster_pred = as.integer(ann$consensuscluster)))
    }
  }
  NULL
}

types_by_mod <- function(mod) if (identical(mod, "SNPs")) "binary" else "gaussian"
mod_order <- c("SNPs", "RNAseq", "CNV", "miRNA", "Methylation")

# ---- merge mode ----
if (identical(tolower(mode), "merge")) {
  run_dir <- Sys.getenv("BENCH_RUN_DIR", dirname(out_root))
  row_files <- list.files(run_dir, pattern = "^perf_row_samples_[0-9]+pct_rep[0-9]+\\.tsv$", recursive = TRUE, full.names = TRUE)
  if (length(row_files) == 0) stop("No perf_row_samples_*.tsv files found under: ", run_dir)

  dt_list <- lapply(row_files, data.table::fread)
  perf_dt <- data.table::rbindlist(dt_list, use.names = TRUE, fill = TRUE)
  if (all(c("Sample_Fraction","Replicate") %in% names(perf_dt))) data.table::setorder(perf_dt, Sample_Fraction, Replicate)

  perf_file <- file.path(out_root, sprintf("%s_sample_perturbations_performance.tsv", algorithm))
  data.table::fwrite(perf_dt, perf_file, sep = "\t", quote = FALSE, na = "NA")
  writeLines(capture.output(sessionInfo()),
             file.path(out_root, sprintf("%s_session_info.txt", algorithm)))
  message("Merged and wrote: ", perf_file)
  quit(save = "no", status = 0)
}

# ---- worker mode ----
task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID", "")
if (!nzchar(task_id)) stop("Worker mode requires SLURM_ARRAY_TASK_ID (use a job array)")

if (!file.exists(lra_source_path)) stop("Missing: ", lra_source_path)
source(lra_source_path)

input_full <- readRDS(input_path)

ground_truth_labels <- openxlsx::read.xlsx(gt_path)
if (!all(c("Sample.ID", "Cluster") %in% colnames(ground_truth_labels))) stop("Ground truth must have Sample.ID, Cluster")
ground_truth_labels$Sample.ID <- gsub("\\.", "-", ground_truth_labels$Sample.ID)
colnames(ground_truth_labels)[colnames(ground_truth_labels) == "Cluster"] <- "Cluster_gt"

scheme <- readRDS(scheme_path)
m3c_des <- scheme$annCol
er_col <- NULL
for (nm in c("ER status", "ER_status", "ER_status ", "ER Status")) if (nm %in% colnames(m3c_des)) { er_col <- nm; break }
m3c_des$class <- if (!is.null(er_col)) m3c_des[[er_col]] else NA
m3c_des$ID <- rownames(m3c_des)

sample_subsets <- data.table::fread(subsets_path)
need_cols <- c("fraction","replicate","sample_id")
if (!all(need_cols %in% colnames(sample_subsets))) stop("sample_subsets.tsv.gz must contain fraction, replicate, sample_id")

non_null_mods <- names(input_full)[!vapply(input_full, is.null, logical(1))]
mods_present <- intersect(mod_order, non_null_mods)
ref_mod_full <- if ("RNAseq" %in% mods_present) "RNAseq" else mods_present[1]
valid_ids <- colnames(input_full[[ref_mod_full]])

sample_subsets <- sample_subsets[sample_id %in% valid_ids]
pairs <- unique(sample_subsets[, .(fraction, replicate)])
data.table::setorder(pairs, fraction, replicate)

idx <- as.integer(task_id)
if (is.na(idx) || idx < 1L || idx > nrow(pairs)) stop("Invalid SLURM_ARRAY_TASK_ID=", task_id)

frac <- as.integer(pairs$fraction[idx])
repi <- as.integer(pairs$replicate[idx])

message(sprintf("[%s] %s sample subset %d%% rep %d (task %d/%d)", Sys.time(), algorithm, frac, repi, idx, nrow(pairs)))

keep_ids <- unique(sample_subsets[fraction == frac & replicate == repi, sample_id])
input_sub <- subset_input_cols(input_full, keep_ids)

X_use <- input_sub[mods_present]
X_use <- Filter(function(m) !is.null(m) && nrow(m) > 0L && ncol(m) > 0L, X_use)
if (length(X_use) < 2L) stop("Too few modalities after column subsetting at ", frac, "% rep ", repi)

X_use <- align_modalities(X_use, ref_mod = if ("RNAseq" %in% names(X_use)) "RNAseq" else names(X_use)[1])
ref_mod <- if ("RNAseq" %in% names(X_use)) "RNAseq" else names(X_use)[1]
sample_ids <- colnames(X_use[[ref_mod]])
n_samples <- length(sample_ids)

X_use <- X_use[intersect(mod_order, names(X_use))]
types <- vapply(names(X_use), types_by_mod, character(1))

p_by_mod_full <- sapply(c("RNAseq","CNV","Methylation","miRNA","SNPs"), function(m) if (!is.null(input_full[[m]])) nrow(input_full[[m]]) else NA_integer_)
p_total_full <- sum(p_by_mod_full, na.rm = TRUE)

# ---- measure ONLY LRAcluster() ----
res <- NULL
gc()
pr <- peakRAM::peakRAM({
  res <<- LRAcluster(data = X_use, types = types, dimension = optr_dim, names = names(X_use))
})
elapsed_s <- as.numeric(pr$Elapsed_Time[1])
peak_mib <- as.numeric(pr$Peak_RAM_Used_MiB[1])

coord <- as.matrix(res[["coordinate"]])
rm(res); gc()

m3c_input <- as.data.frame(coord)
rownames(m3c_input) <- paste0("LRA_", rownames(coord))
colnames(m3c_input) <- sample_ids

des_sub <- m3c_des[m3c_des$ID %in% sample_ids, , drop = FALSE]
if (nrow(des_sub) == 0) {
  des_sub <- data.frame(ID = sample_ids, class = NA)
  rownames(des_sub) <- des_sub$ID
}

m3c_res <- M3C::M3C(
  mydata = m3c_input, des = des_sub,
  cores = m3c_cores, iters = m3c_iters, repsref = m3c_repsref, repsreal = m3c_repsreal,
  seed = 123, clusteralg = "km", maxK = m3c_maxK,
  objective = "entropy", removeplots = TRUE, silent = TRUE
)

clusters_df <- extract_m3c_clusters(m3c_res, k = optk)
rm(m3c_res); gc()

if (is.null(clusters_df) || nrow(clusters_df) == 0) {
  km <- kmeans(t(coord), centers = optk, nstart = 100)
  clusters_df <- data.frame(Sample.ID = colnames(coord), Cluster_pred = as.integer(km$cluster))
}

clusters_df$Sample.ID <- gsub("\\.", "-", clusters_df$Sample.ID)

ari_gt <- ari_to_ground_truth(ground_truth_labels, clusters_df)

avg_width <- NA_real_
coord_t <- t(coord)
rownames(coord_t) <- gsub("\\.", "-", rownames(coord_t))
coord_t <- coord_t[clusters_df$Sample.ID, , drop = FALSE]
if (nrow(coord_t) >= optk + 1L && length(unique(clusters_df$Cluster_pred)) > 1L) {
  sil <- cluster::silhouette(as.integer(clusters_df$Cluster_pred), Rfast::Dist(coord_t, method = "euclidean"))
  avg_width <- summary(sil)$avg.width
}

clust_file <- file.path(out_root, sprintf("%s_clusters_samples_%dpct_rep%02d.tsv.gz", algorithm, frac, repi))
data.table::fwrite(data.table::as.data.table(clusters_df), clust_file, sep = "\t", quote = FALSE, compress = "gzip")

perf_row <- data.table::data.table(
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
  n_clusters = as.integer(optk),
  LRAcluster_dim = as.integer(optr_dim),
  ARI_to_Ground_Truth = ari_gt,
  Average_Silhouette_Width = avg_width,
  LRAcluster_Elapsed_Seconds = elapsed_s,
  LRAcluster_PeakRAM_MiB = peak_mib
)

row_file <- file.path(out_root, sprintf("perf_row_samples_%dpct_rep%02d.tsv", frac, repi))
data.table::fwrite(perf_row, row_file, sep = "\t", quote = FALSE, na = "NA")
message("Wrote partial row: ", row_file)

