#!/usr/bin/env Rscript

# Install required packages if not available
if (!requireNamespace("heatmap.plus", quietly = TRUE)) {
  install.packages("https://cran.r-project.org/src/contrib/Archive/heatmap.plus/heatmap.plus_1.3.tar.gz",
                   repos = NULL,
                   type = "source")
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
    install.packages(pkg,repos='https://cran.r-project.org' )
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

# Fixed centiles (fractions, not percents)
centiles <- c(0.10, 0.20, 0.50, 0.75, 0.90)

# Fixed SNF parameters
K_snf <- 25L
sigma_snf <- 0.5
t_snf <- 50L
n_clusters <- 2L

input_path <- "Resources/mm_input.rds"

# Make outputs job-safe by default (can be overridden)
job_out <- Sys.getenv("BENCH_OUT_DIR", "")
if (nzchar(job_out)) {
  out_root <- job_out
} else {
  job_id <- Sys.getenv("SLURM_JOB_ID", "")
  out_root <- file.path("Results", "Performance", "Feature_perturbations", algorithm,
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

feat_rank <- fread("Resources/Feature_perturbations/feature_rankings.tsv.gz")
if (!all(c("modality", "rank", "feature_id") %in% colnames(feat_rank))) {
  stop("feature_rankings.tsv.gz must contain columns: modality, rank, feature_id")
}

modalities_ranked <- intersect(unique(feat_rank$modality), names(input_full))
if (length(modalities_ranked) == 0) stop("No overlap between ranked modalities and input modalities.")

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
      out[[mod]] <- out[[mod]][keep, , drop = FALSE]
    }
  }
  out
}

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
  adjustedRandIndex(merged[["Cluster_gt"]], merged[["Cluster_pred"]])
}

continuous <- intersect(c("RNAseq", "CNV", "Methylation", "miRNA"), names(input_full))

perf_rows <- vector("list", length(centiles))

for (i in seq_along(centiles)) {
  centile <- centiles[i]
  pct <- as.integer(round(centile * 100))
  
  message(sprintf("[%s] %s feature subset %d%%", Sys.time(), algorithm, pct))
  
  top_feats <- get_top_features(feat_rank[modality %in% modalities_ranked], centile)
  input_sub <- subset_input_features(input_full, top_feats)
  
  p_by_mod <- sapply(modalities_ranked, function(m) if (!is.null(input_sub[[m]])) nrow(input_sub[[m]]) else NA_integer_)
  m0 <- modalities_ranked[which(!vapply(input_sub[modalities_ranked], is.null, logical(1)))[1]]
  n_samples <- ncol(input_sub[[m0]])
  p_total <- sum(p_by_mod, na.rm = TRUE)
  
  input_t <- lapply(input_sub[modalities_ranked], t)
  
  input_dists <- list()
  for (m in intersect(continuous, modalities_ranked)) {
    x <- as.matrix(input_t[[m]])
    input_dists[[m]] <- SNFtool::dist2(x, x)
  }
  if ("SNPs" %in% modalities_ranked) {
    x <- as.matrix(input_t[["SNPs"]])
    input_dists[["SNPs"]] <- as.matrix(dist(x, method = "binary"))
  }
  
  aff <- lapply(input_dists, function(d) SNFtool::affinityMatrix(as.matrix(d), K = K_snf, sigma = sigma_snf))
  
  gc()
  pr <- peakRAM::peakRAM({
    fusion <- SNFtool::SNF(aff, K = K_snf, t = t_snf)
  })
  snf_elapsed <- as.numeric(pr$Elapsed_Time[1])
  snf_peak_mb <- as.numeric(pr$Peak_RAM_Used_MiB[1])
  
  res <- spectralClustering_eig(fusion, n_clusters)
  col_index <- ncol(res) - 2
  sil <- silhouette(as.integer(res$Cluster),
                    dist = Rfast::Dist(res[, 1:col_index], method = "euclidean"))
  avg_width <- summary(sil)$avg.width
  
  clusters <- res[, c("Sample.ID", "Cluster")]
  clusters$Sample.ID <- gsub("\\.", "-", clusters$Sample.ID)
  colnames(clusters)[2] <- "Cluster_pred"
  
  ari_gt <- ari_to_ground_truth(ground_truth_labels, clusters)
  
  clust_file <- file.path(out_root, sprintf("%s_clusters_%dpct.tsv.gz", algorithm, pct))
  fwrite(as.data.table(clusters), clust_file, sep = "\t", quote = FALSE, compress = "gzip")
  
  perf_rows[[i]] <- data.table(
    Algorithm = algorithm,
    Feature_Centile = centile,
    Feature_Percent = pct,
    n_samples = n_samples,
    p_total = p_total,
    p_RNAseq = as.integer(p_by_mod[["RNAseq"]]),
    p_CNV = as.integer(p_by_mod[["CNV"]]),
    p_Methylation = as.integer(p_by_mod[["Methylation"]]),
    p_miRNA = as.integer(p_by_mod[["miRNA"]]),
    p_SNPs = as.integer(p_by_mod[["SNPs"]]),
    SNF_K = K_snf,
    SNF_sigma = sigma_snf,
    SNF_t = t_snf,
    n_clusters = n_clusters,
    ARI_to_Ground_Truth = ari_gt,
    Average_Silhouette_Width = avg_width,
    SNF_Elapsed_Seconds = snf_elapsed,
    SNF_PeakRAM_MB = snf_peak_mb
  )
  
  rm(input_sub, input_t, input_dists, aff, fusion, res, sil, clusters)
  gc()
}

perf_dt <- rbindlist(perf_rows, use.names = TRUE, fill = TRUE)

perf_file <- file.path(out_root, sprintf("%s_feature_perturbations_performance.tsv", algorithm))
fwrite(perf_dt, perf_file, sep = "\t", quote = FALSE, na = "NA")

writeLines(capture.output(sessionInfo()),
           file.path(out_root, sprintf("%s_session_info.txt", algorithm)))

message("Wrote: ", perf_file)
