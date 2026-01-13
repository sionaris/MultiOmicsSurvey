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

if (!requireNamespace("abSNF", quietly = TRUE)) {
  if (!requireNamespace("devtools", quietly = TRUE) || packageVersion("devtools") < "1.6") {
    install.packages("devtools")
  }
  devtools::install_github("pfruan/abSNF")
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
  library(Rfast)
  library(cluster)
  library(mclust)
  library(peakRAM)
  library(abSNF)
})

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b

RNGversion("4.2.2")
set.seed(123)

algorithm <- "ab-SNF"

K_nn <- 10L
sigma_snf <- 0.5
t_snf <- 50L
optk <- 2L

input_path <- "Resources/mm_input.rds"
subsets_path <- "Resources/Sample_perturbations/sample_subsets.tsv.gz"
cosmic_path <- "Resources/COSMIC_CGC_Breast_somatic.csv"

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

COSMIC_BC_drivers <- tryCatch({
  dd <- read.csv(cosmic_path, stringsAsFactors = FALSE)
  if (!"Gene.Symbol" %in% colnames(dd)) stop("COSMIC file missing 'Gene.Symbol' column.")
  unique(as.character(dd$Gene.Symbol))
}, error = function(e) {
  stop("Failed reading COSMIC drivers: ", conditionMessage(e))
})

non_null_mods <- names(input_full)[!vapply(input_full, is.null, logical(1))]
ref_mod_full <- if ("RNAseq" %in% non_null_mods) "RNAseq" else non_null_mods[1]
valid_ids <- colnames(input_full[[ref_mod_full]])

sample_subsets <- sample_subsets[sample_id %in% valid_ids]
if (nrow(sample_subsets) == 0) stop("No sample IDs in sample_subsets match the input colnames.")

continuous <- intersect(c("RNAseq", "CNV", "Methylation", "miRNA"), non_null_mods)

mad_based_weights <- function(X) {
  v <- apply(X, 2, mad, na.rm = TRUE)
  v[!is.finite(v)] <- 0
  s <- sum(v)
  if (s <= 0) rep(0, length(v)) else v / s
}

weighted_euclidean_dist <- function(X, w, ids) {
  Xw <- sweep(X, 2, sqrt(w), `*`)
  d <- abSNF::dist2(as.matrix(Xw), as.matrix(Xw))
  dimnames(d) <- list(ids, ids)
  d
}

weighted_hamming_fraction_fast <- function(Xbin, features, drivers, driver_w = 0.8, nondriver_w = 0.2, ids) {
  X <- (as.matrix(Xbin) != 0) * 1
  w <- ifelse(features %in% drivers, driver_w, nondriver_w)
  w[!is.finite(w)] <- 0
  total_w <- sum(w)
  if (total_w <= 0) {
    d <- matrix(0, nrow = nrow(X), ncol = nrow(X))
    dimnames(d) <- list(ids, ids)
    return(d)
  }
  Xw <- sweep(X, 2, sqrt(w), `*`)
  C <- Xw %*% t(Xw)
  s <- as.numeric(X %*% w)
  d <- outer(s, s, "+") - 2 * C
  d <- pmax(d, 0)
  d <- d / total_w
  dimnames(d) <- list(ids, ids)
  d
}

spectralClustering_eig <- function(affinity, K, type = 3) {
  ids <- colnames(affinity) %||% rownames(affinity) %||% paste0("S", seq_len(nrow(affinity)))
  
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
  
  out <- as.data.frame(U)
  out$Cluster <- labels
  out$Sample.ID <- ids
  out
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
      if (ncol(out[[nm]]) == 0L) out[[nm]] <- NULL
    }
  }
  out
}

pairs <- unique(sample_subsets[, .(fraction, replicate)])
setorder(pairs, fraction, replicate)

p_by_mod_full <- sapply(non_null_mods, function(m) nrow(input_full[[m]]))
p_total_full <- sum(p_by_mod_full, na.rm = TRUE)
get_int1 <- function(x) if (is.null(x) || length(x) == 0L || is.na(x)) NA_integer_ else as.integer(x[1])

perf_rows <- vector("list", nrow(pairs))

for (i in seq_len(nrow(pairs))) {
  frac <- pairs$fraction[i]
  repi <- pairs$replicate[i]
  
  keep_ids <- unique(sample_subsets[fraction == frac & replicate == repi, sample_id])
  input_sub <- subset_input_cols(input_full, keep_ids)
  
  mods_use <- names(input_sub)[!vapply(input_sub, is.null, logical(1))]
  if (length(mods_use) == 0L) {
    perf_rows[[i]] <- data.table(
      Algorithm = algorithm, Sample_Fraction = frac, Replicate = repi,
      n_samples = 0L, p_total = as.integer(p_total_full),
      K_neighbors = NA_integer_, sigma = sigma_snf, t = t_snf, optk = optk,
      ARI_to_Ground_Truth = NA_real_, Average_Silhouette_Width = NA_real_,
      SNF_Elapsed_Seconds = NA_real_, SNF_PeakRAM_MiB = NA_real_
    )
    rm(input_sub); gc()
    next
  }
  
  ref_mod <- if ("RNAseq" %in% mods_use) "RNAseq" else mods_use[1]
  sample_ids <- colnames(input_sub[[ref_mod]])
  n_samples <- length(sample_ids)
  
  if (n_samples < 3L) {
    perf_rows[[i]] <- data.table(
      Algorithm = algorithm, Sample_Fraction = frac, Replicate = repi,
      n_samples = as.integer(n_samples), p_total = as.integer(p_total_full),
      K_neighbors = as.integer(max(1L, n_samples - 1L)), sigma = sigma_snf, t = t_snf, optk = optk,
      ARI_to_Ground_Truth = NA_real_, Average_Silhouette_Width = NA_real_,
      SNF_Elapsed_Seconds = NA_real_, SNF_PeakRAM_MiB = NA_real_
    )
    rm(input_sub); gc()
    next
  }
  
  K_use <- min(K_nn, n_samples - 1L)
  
  message(sprintf("[%s] %s sample subset %d%% rep %d (n=%d)", Sys.time(), algorithm, frac, repi, n_samples))
  
  input_t <- lapply(input_sub[mods_use], function(mat) {
    tm <- t(mat)
    rownames(tm) <- sample_ids
    tm
  })
  
  input_dists <- list()
  for (m in intersect(continuous, mods_use)) {
    X <- as.matrix(input_t[[m]])
    w <- mad_based_weights(X)
    input_dists[[m]] <- weighted_euclidean_dist(X, w, sample_ids)
  }
  if ("SNPs" %in% mods_use) {
    X <- as.matrix(input_t[["SNPs"]])
    input_dists[["SNPs"]] <- weighted_hamming_fraction_fast(
      Xbin = X,
      features = colnames(X),
      drivers = COSMIC_BC_drivers,
      driver_w = 0.8,
      nondriver_w = 0.2,
      ids = sample_ids
    )
  }
  gc()
  
  aff <- lapply(input_dists, function(d) {
    a <- abSNF::affinityMatrix(as.matrix(d), K = K_use, sigma = sigma_snf)
    dimnames(a) <- list(sample_ids, sample_ids)
    a
  })
  
  fusion <- NULL
  gc()
  pr <- peakRAM::peakRAM({
    fusion <<- abSNF::SNF(aff, K = K_use, t = t_snf)
  })
  snf_elapsed <- as.numeric(pr$Elapsed_Time[1])
  snf_peak_mib <- as.numeric(pr$Peak_RAM_Used_MiB[1])
  dimnames(fusion) <- list(sample_ids, sample_ids)
  
  res <- spectralClustering_eig(fusion, min(optk, n_samples))
  col_index <- ncol(res) - 2
  sil <- silhouette(as.integer(res$Cluster),
                    dist = Rfast::Dist(res[, 1:col_index, drop = FALSE], method = "euclidean"))
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
    n_samples = as.integer(n_samples),
    p_total = as.integer(p_total_full),
    p_RNAseq = get_int1(p_by_mod_full[["RNAseq"]]),
    p_CNV = get_int1(p_by_mod_full[["CNV"]]),
    p_Methylation = get_int1(p_by_mod_full[["Methylation"]]),
    p_miRNA = get_int1(p_by_mod_full[["miRNA"]]),
    p_SNPs = get_int1(p_by_mod_full[["SNPs"]]),
    K_neighbors = as.integer(K_use),
    sigma = as.numeric(sigma_snf),
    t = as.integer(t_snf),
    optk = as.integer(optk),
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
