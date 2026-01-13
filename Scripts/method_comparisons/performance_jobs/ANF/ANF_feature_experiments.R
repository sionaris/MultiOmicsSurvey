#!/usr/bin/env Rscript

if (!requireNamespace("heatmap.plus", quietly = TRUE)) {
  install.packages("https://cran.r-project.org/src/contrib/Archive/heatmap.plus/heatmap.plus_1.3.tar.gz",
                   repos = NULL, type = "source")
}

if (!requireNamespace("SNFtool", quietly = TRUE)) {
  if (!requireNamespace("devtools", quietly = TRUE) || packageVersion("devtools") < "1.6") {
    install.packages("devtools", repos = "https://cran.r-project.org")
  }
  devtools::install_github("maxconway/SNFtool")
}

if (!requireNamespace("ANF", quietly = TRUE)) {
  if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools", repos = "https://cran.r-project.org")
  }
  devtools::install_github("BeautyOfWeb/ANF")
}

packages <- c("dplyr", "data.table", "openxlsx", "Rfast", "cluster", "mclust", "peakRAM")
invisible(lapply(packages, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg, repos = "https://cran.r-project.org")
}))

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(openxlsx)
  library(SNFtool)
  library(ANF)
  library(Rfast)
  library(cluster)
  library(mclust)
  library(peakRAM)
})

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b

RNGversion("4.2.2")
set.seed(123)

algorithm <- "ANF"

centiles <- c(0.10, 0.20, 0.50, 0.75, 0.90)

optN <- 30L
optk <- 2L

alpha_ab <- 1/6
beta_ab <- 1/6
anf_type <- "two-step"
anf_weight <- c(1, 1, 0, 0, 0, 0, 0, 0)

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
      if (nrow(out[[mod]]) == 0L) out[[mod]] <- NULL
    }
  }
  out
}

spectralClustering_eig <- function(affinity, K, type = 3) {
  affinity <- as.matrix(affinity)
  affinity <- (affinity + t(affinity)) / 2
  diag(affinity) <- 0

  ids <- colnames(affinity) %||% rownames(affinity) %||% paste0("S", seq_len(nrow(affinity)))

  d <- rowSums(affinity)
  d[d == 0] <- .Machine$double.eps
  L <- diag(d) - affinity

  if (type == 1) {
    NL <- L
    NL <- (NL + t(NL)) / 2
    eig <- eigen(NL, symmetric = TRUE)
  } else if (type == 2) {
    NL <- diag(1 / d) %*% L
    eig <- eigen(NL)
  } else {
    Di <- diag(1 / sqrt(d))
    NL <- Di %*% L %*% Di
    NL <- (NL + t(NL)) / 2
    eig <- eigen(NL, symmetric = TRUE)
  }

  vals <- Re(eig$values)
  vecs <- Re(eig$vectors)

  ix <- order(abs(vals))[seq_len(K)]
  U <- vecs[, ix, drop = FALSE]

  if (type == 3) {
    normalize <- function(x) x / sqrt(sum(x^2))
    U <- t(apply(U, 1, normalize))
  }

  eigDiscrete <- SNFtool:::.discretisation(U)$discrete
  labels <- apply(eigDiscrete, 1, which.max)

  out <- as.data.frame(U)
  out$Cluster <- labels
  out$Sample.ID <- ids
  out
}

ari_to_ground_truth <- function(gt_df, pred_df) {
  merged <- merge(gt_df, pred_df, by = "Sample.ID")
  adjustedRandIndex(merged[["Cluster_gt"]], merged[["Cluster_pred"]])
}

continuous <- intersect(c("RNAseq", "CNV", "Methylation", "miRNA"), names(input_full))
non_null_mods_full <- names(input_full)[!vapply(input_full, is.null, logical(1))]

perf_rows <- vector("list", length(centiles))

for (i in seq_along(centiles)) {
  centile <- centiles[i]
  pct <- as.integer(round(centile * 100))
  message(sprintf("[%s] %s feature subset %d%%", Sys.time(), algorithm, pct))
  
  top_feats <- get_top_features(feat_rank[modality %in% modalities_ranked], centile)
  input_sub <- subset_input_features(input_full, top_feats)
  
  mods_use <- names(input_sub)[!vapply(input_sub, is.null, logical(1))]
  if (length(mods_use) == 0) stop("All modalities became NULL after feature subsetting at ", pct, "%")
  
  ref_mod <- if ("RNAseq" %in% mods_use) "RNAseq" else mods_use[1]
  sample_ids <- colnames(input_sub[[ref_mod]])
  if (is.null(sample_ids) || length(sample_ids) == 0) stop("Missing sample IDs (colnames) in reference modality.")
  n_samples <- length(sample_ids)
  
  K_use <- min(optN, n_samples - 1L)
  
  p_by_mod <- sapply(mods_use, function(m) nrow(input_sub[[m]]))
  p_total <- sum(p_by_mod, na.rm = TRUE)
  
  input_t <- lapply(input_sub[mods_use], function(mat) {
    tm <- t(mat)
    rownames(tm) <- sample_ids
    tm
  })
  
  input_dists <- list()
  for (m in intersect(continuous, mods_use)) {
    x <- as.matrix(input_t[[m]])
    d <- SNFtool::dist2(x, x)
    dimnames(d) <- list(sample_ids, sample_ids)
    input_dists[[m]] <- d
  }
  if ("SNPs" %in% mods_use) {
    x <- as.matrix(input_t[["SNPs"]])
    d <- as.matrix(dist(x, method = "binary"))
    dimnames(d) <- list(sample_ids, sample_ids)
    input_dists[["SNPs"]] <- d
  }
  
  aff_list <- lapply(input_dists, function(d) {
    a <- ANF::affinity_matrix(as.matrix(d), k = K_use, alpha = alpha_ab, beta = beta_ab)
    dimnames(a) <- list(sample_ids, sample_ids)
    a
  })
  
  fusion <- NULL
  gc()
  pr <- peakRAM::peakRAM({
    fusion <<- ANF::ANF(aff_list, K = K_use, weight = NULL, alpha = anf_weight,
                        type = anf_type, verbose = FALSE)
  })
  anf_elapsed <- as.numeric(pr$Elapsed_Time[1])
  anf_peak_mib <- as.numeric(pr$Peak_RAM_Used_MiB[1])
  dimnames(fusion) <- list(sample_ids, sample_ids)
  
  res <- spectralClustering_eig(fusion, optk)
  col_index <- ncol(res) - 2
  sil <- silhouette(as.integer(res$Cluster),
                    dist = Rfast::Dist(res[, 1:col_index, drop = FALSE], method = "euclidean"))
  avg_width <- summary(sil)$avg.width
  
  clusters <- res[, c("Sample.ID", "Cluster")]
  clusters$Sample.ID <- gsub("\\.", "-", clusters$Sample.ID)
  colnames(clusters)[2] <- "Cluster_pred"
  
  ari_gt <- ari_to_ground_truth(ground_truth_labels, clusters)
  
  clust_file <- file.path(out_root, sprintf("%s_clusters_%dpct.tsv.gz", algorithm, pct))
  fwrite(as.data.table(clusters), clust_file, sep = "\t", quote = FALSE, compress = "gzip")
  
  get_int1 <- function(x) if (is.null(x) || length(x) == 0L || is.na(x)) NA_integer_ else as.integer(x[1])
  
  perf_rows[[i]] <- data.table(
    Algorithm = algorithm,
    Feature_Centile = centile,
    Feature_Percent = pct,
    n_samples = as.integer(n_samples),
    p_total = as.integer(p_total),
    p_RNAseq = get_int1(p_by_mod[["RNAseq"]]),
    p_CNV = get_int1(p_by_mod[["CNV"]]),
    p_Methylation = get_int1(p_by_mod[["Methylation"]]),
    p_miRNA = get_int1(p_by_mod[["miRNA"]]),
    p_SNPs = get_int1(p_by_mod[["SNPs"]]),
    optN = as.integer(K_use),
    optk = as.integer(optk),
    alpha = as.numeric(alpha_ab),
    beta = as.numeric(beta_ab),
    weight = paste(anf_weight, collapse = ","),
    type = anf_type,
    ARI_to_Ground_Truth = ari_gt,
    Average_Silhouette_Width = avg_width,
    ANF_Elapsed_Seconds = anf_elapsed,
    ANF_PeakRAM_MiB = anf_peak_mib
  )
  
  rm(input_sub, input_t, input_dists, aff_list, fusion, res, sil, clusters)
  gc()
}

perf_dt <- rbindlist(perf_rows, use.names = TRUE, fill = TRUE)
perf_file <- file.path(out_root, sprintf("%s_feature_perturbations_performance.tsv", algorithm))
fwrite(perf_dt, perf_file, sep = "\t", quote = FALSE, na = "NA")

writeLines(capture.output(sessionInfo()),
           file.path(out_root, sprintf("%s_session_info.txt", algorithm)))

message("Wrote: ", perf_file)
