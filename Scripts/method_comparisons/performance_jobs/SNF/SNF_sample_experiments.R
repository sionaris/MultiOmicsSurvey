#!/usr/bin/env Rscript

if (!requireNamespace("heatmap.plus", quietly = TRUE)) {
  install.packages("https://cran.r-project.org/src/contrib/Archive/heatmap.plus/heatmap.plus_1.3.tar.gz",
                   repos = NULL, type = "source")
}

if (!requireNamespace("SNFtool", quietly = TRUE)) {
  if (!requireNamespace("devtools", quietly = TRUE) || packageVersion("devtools") < "1.6") {
    install.packages("devtools")
  }
  devtools::install_github("maxconway/SNFtool")
}

packages = c("dplyr", "data.table", "openxlsx", "Rfast", "cluster", "mclust", "peakRAM")
invisible(lapply(packages, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cran.r-project.org")
  }
}))

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(openxlsx)
  library(SNFtool)
  library(Rfast)
  library(cluster)
  library(mclust)
  library(peakRAM)
})

RNGversion("4.2.2")
set.seed(123)

algorithm <- "SNF"

K_snf <- 25L
sigma_snf <- 0.5
t_snf <- 50L
n_clusters <- 2L

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

input_full <- readRDS(input_path)

ground_truth_labels <- openxlsx::read.xlsx(
  file.path("Resources/algorithms", paste0(algorithm, "_ground_truth_labels.xlsx"))
)
if (!all(c("Sample.ID", "Cluster") %in% colnames(ground_truth_labels))) {
  stop("Ground truth file must contain columns: Sample.ID, Cluster")
}
ground_truth_labels$Sample.ID <- gsub("\\.", "-", ground_truth_labels$Sample.ID)
colnames(ground_truth_labels)[colnames(ground_truth_labels) == "Cluster"] <- "Cluster_gt"

sample_subsets <- fread(subsets_path)
need_cols <- c("fraction", "replicate", "sample_id")
if (!all(need_cols %in% colnames(sample_subsets))) {
  stop("sample_subsets.tsv.gz must contain columns: fraction, replicate, sample_id")
}

non_null_mods <- names(input_full)[!vapply(input_full, is.null, logical(1))]
ref_mod <- if ("RNAseq" %in% non_null_mods) "RNAseq" else non_null_mods[1]
valid_ids <- colnames(input_full[[ref_mod]])

sample_subsets <- sample_subsets[sample_id %in% valid_ids]
if (nrow(sample_subsets) == 0) stop("No sample IDs in sample_subsets match the input colnames.")

continuous <- intersect(c("RNAseq", "CNV", "Methylation", "miRNA"), non_null_mods)

spectralClustering_eig <- function(affinity, K, type = 3) {
  d <- rowSums(affinity)
  d[d == 0] <- .Machine$double.eps
  D <- diag(d)
  L <- D - affinity
  if (type == 1) {
    NL <- L
  } else if (type == 2) {
    NL <- diag(1 / d) %*% L
  } else {
    Di <- diag(1 / sqrt(d))
    NL <- Di %*% L %*% Di
  }
  eig <- eigen(NL)
  res <- sort(abs(eig$values), index.return = TRUE)
  U <- eig$vectors[, res$ix[1:K], drop = FALSE]
  if (type == 3) {
    normalize <- function(x) x / sqrt(sum(x^2))
    U <- t(apply(U, 1, normalize))
  }
  eigDiscrete <- SNFtool:::.discretisation(U)$discrete
  labels <- apply(eigDiscrete, 1, which.max)
  U <- as.data.frame(cbind(U, labels))
  U$Sample.ID <- colnames(affinity)
  colnames(U)[(ncol(U) - 1):ncol(U)] <- c("Cluster", "Sample.ID")
  U
}

ari_to_ground_truth <- function(gt_df, pred_df) {
  merged <- merge(gt_df, pred_df, by = "Sample.ID")
  if (nrow(merged) == 0) return(NA_real_)
  adjustedRandIndex(merged[["Cluster_gt"]], merged[["Cluster_pred"]])
}

subset_input_cols <- function(input_obj, keep_ids) {
  out <- input_obj
  for (nm in names(out)) {
    if (!is.null(out[[nm]])) {
      kk <- keep_ids[keep_ids %in% colnames(out[[nm]])]
      out[[nm]] <- out[[nm]][, kk, drop = FALSE]
    }
  }
  out
}

pairs <- unique(sample_subsets[, .(fraction, replicate)])
setorder(pairs, fraction, replicate)

p_by_mod_full <- sapply(non_null_mods, function(m) nrow(input_full[[m]]))
p_total_full <- sum(p_by_mod_full, na.rm = TRUE)

perf_rows <- vector("list", nrow(pairs))

for (i in seq_len(nrow(pairs))) {
  frac <- pairs$fraction[i]
  repi <- pairs$replicate[i]

  keep_ids <- sample_subsets[fraction == frac & replicate == repi, sample_id]
  keep_ids <- unique(keep_ids)

  input_sub <- subset_input_cols(input_full, keep_ids)

  m0 <- ref_mod
  n_samples <- ncol(input_sub[[m0]])

  if (n_samples < 3L) {
    perf_rows[[i]] <- data.table(
      Algorithm = algorithm,
      Sample_Fraction = frac,
      Replicate = repi,
      n_samples = n_samples,
      p_total = p_total_full,
      SNF_K = K_snf,
      SNF_sigma = sigma_snf,
      SNF_t = t_snf,
      n_clusters = n_clusters,
      ARI_to_Ground_Truth = NA_real_,
      Average_Silhouette_Width = NA_real_,
      SNF_Elapsed_Seconds = NA_real_,
      SNF_PeakRAM_MiB = NA_real_
    )
    next
  }

  K_use <- min(K_snf, n_samples - 1L)

  message(sprintf("[%s] %s sample subset %d%% rep %d (n=%d)", Sys.time(), algorithm, frac, repi, n_samples))

  input_t <- lapply(input_sub[non_null_mods], t)

  input_dists <- list()
  for (m in intersect(continuous, non_null_mods)) {
    x <- as.matrix(input_t[[m]])
    input_dists[[m]] <- SNFtool::dist2(x, x)
  }
  if ("SNPs" %in% non_null_mods) {
    x <- as.matrix(input_t[["SNPs"]])
    input_dists[["SNPs"]] <- as.matrix(dist(x, method = "binary"))
  }

  aff <- lapply(input_dists, function(d) SNFtool::affinityMatrix(as.matrix(d), K = K_use, sigma = sigma_snf))

  gc()
  pr <- peakRAM::peakRAM({
    fusion <- SNFtool::SNF(aff, K = K_use, t = t_snf)
  })
  snf_elapsed <- as.numeric(pr$Elapsed_Time[1])
  snf_peak_mib <- as.numeric(pr$Peak_RAM_Used_MiB[1])

  res <- spectralClustering_eig(fusion, min(n_clusters, n_samples))
  col_index <- ncol(res) - 2
  sil <- silhouette(as.integer(res$Cluster),
                    dist = Rfast::Dist(res[, 1:col_index], method = "euclidean"))
  avg_width <- summary(sil)$avg.width

  clusters <- res[, c("Sample.ID", "Cluster")]
  clusters$Sample.ID <- gsub("\\.", "-", clusters$Sample.ID)
  colnames(clusters)[2] <- "Cluster_pred"

  ari_gt <- ari_to_ground_truth(ground_truth_labels, clusters)

  clust_file <- file.path(out_root, sprintf("%s_clusters_samples_%dpct_rep%02d.tsv.gz", algorithm, frac, repi))
  fwrite(as.data.table(clusters), clust_file, sep = "\t", quote = FALSE, compress = "gzip")

  perf_rows[[i]] <- data.table(
    Algorithm = algorithm,
    Sample_Fraction = frac,
    Replicate = repi,
    n_samples = n_samples,
    p_total = p_total_full,
    p_RNAseq = as.integer(p_by_mod_full[["RNAseq"]]),
    p_CNV = as.integer(p_by_mod_full[["CNV"]]),
    p_Methylation = as.integer(p_by_mod_full[["Methylation"]]),
    p_miRNA = as.integer(p_by_mod_full[["miRNA"]]),
    p_SNPs = as.integer(p_by_mod_full[["SNPs"]]),
    SNF_K = K_use,
    SNF_sigma = sigma_snf,
    SNF_t = t_snf,
    n_clusters = n_clusters,
    ARI_to_Ground_Truth = ari_gt,
    Average_Silhouette_Width = avg_width,
    SNF_Elapsed_Seconds = snf_elapsed,
    SNF_PeakRAM_MiB = snf_peak_mib
  )

  rm(input_sub, input_t, input_dists, aff, fusion, res, sil, clusters)
  gc()
}

perf_dt <- rbindlist(perf_rows, use.names = TRUE, fill = TRUE)

perf_file <- file.path(out_root, sprintf("%s_sample_perturbations_performance.tsv", algorithm))
fwrite(perf_dt, perf_file, sep = "\t", quote = FALSE, na = "NA")

writeLines(capture.output(sessionInfo()),
           file.path(out_root, sprintf("%s_session_info.txt", algorithm)))

message("Wrote: ", perf_file)

