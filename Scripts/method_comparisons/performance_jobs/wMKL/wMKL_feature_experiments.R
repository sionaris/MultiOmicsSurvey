#!/usr/bin/env Rscript

# ---- packages (except wMKL: must exist, otherwise stop) ----
pkgs <- c("dplyr", "data.table", "openxlsx", "Rfast", "cluster", "mclust", "peakRAM")
invisible(lapply(pkgs, function(p) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cran.r-project.org")
}))

if (!requireNamespace("wMKL", quietly = TRUE)) {
  stop('wMKL installation not found. Check installation instructions in https://github.com/sionaris/MultiOmicsSurvey/blob/main/Scripts/single_algorithm/wMKL.R (lines: 16-72)')
}

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(openxlsx)
  library(Rfast)
  library(cluster)
  library(mclust)
  library(peakRAM)
})

# ---- reproducibility ----
RNGversion("4.2.2")
set.seed(123)

algorithm <- "wMKL"

# Fixed centiles (fractions, not percents)
centiles <- c(0.10, 0.20, 0.50, 0.75, 0.90)

# Fixed wMKL run params
c_fixed <- 8L
k_fixed <- 32L

input_path <- "Resources/mm_input.rds"
cosmic_path <- "Resources/COSMIC_CGC_Breast_somatic.csv"
custom_fns_path <- "Resources/custom_performance_functions.R"

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

if (!file.exists(custom_fns_path)) stop("Missing: ", custom_fns_path)
source(custom_fns_path)

if (!exists("CIMLR.weight_mod", mode = "function")) {
  stop("CIMLR.weight_mod() not found after sourcing ", custom_fns_path)
}

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

COSMIC_BC_drivers <- tryCatch({
  dd <- read.csv(cosmic_path, stringsAsFactors = FALSE)
  if (!"Gene.Symbol" %in% colnames(dd)) stop("COSMIC file missing 'Gene.Symbol' column.")
  unique(as.character(dd$Gene.Symbol))
}, error = function(e) {
  stop("Failed reading COSMIC drivers: ", conditionMessage(e))
})

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

# ---- weights (features in rows, samples in columns) ----
mad_based_weights <- function(X) {
  # MAD per feature (row-wise), then normalise
  v <- apply(X, 1, mad, na.rm = TRUE)
  v[!is.finite(v)] <- 0
  if (sum(v) == 0) {
    w <- rep(1, length(v))
    w <- w / sum(w)
    names(w) <- rownames(X)
  } else {
    w <- v / sum(v)
    names(w) <- rownames(X)
  }
  w
}

weighted_mutations <- function(features, drivers, driver_weight = 0.8, nondriver_weight = 0.2) {
  w <- ifelse(features %in% drivers, driver_weight, nondriver_weight)
  w[!is.finite(w)] <- 0
  if (sum(w) == 0) w <- rep(1, length(w))
  w <- w / sum(w)
  names(w) <- features
  w
}

ari_to_ground_truth <- function(gt_df, pred_df) {
  merged <- merge(gt_df, pred_df, by = "Sample.ID")
  if (nrow(merged) == 0) return(NA_real_)
  adjustedRandIndex(merged[["Cluster_gt"]], merged[["Cluster_pred"]])
}

# ---- silhouette eigenspace extractor (as per your snippet) ----
extract_eigenspace <- function(affinity, K, type = 3) {
  d <- rowSums(affinity)
  d[d == 0] <- .Machine$double.eps
  D <- diag(d)
  L <- D - affinity
  if (type == 1) {
    NL <- L
  } else if (type == 2) {
    Di <- diag(1 / d)
    NL <- Di %*% L
  } else {
    Di <- diag(1 / sqrt(d))
    NL <- Di %*% L %*% Di
  }
  eig <- eigen(NL)
  res <- sort(abs(eig$values), index.return = TRUE)
  U <- eig$vectors[, res$ix[1:K], drop = FALSE]
  normalize <- function(x) x / sqrt(sum(x^2))
  if (type == 3) {
    U <- t(apply(U, 1, normalize))
  }
  eigDiscrete <- wMKL:::.discretisation(U)
  eigDiscrete <- eigDiscrete$discrete
  labels <- apply(eigDiscrete, 1, which.max)
  list(labels = labels, U = U, eigDiscrete = eigDiscrete)
}

safe_avg_sil <- function(labels, U) {
  cl <- as.integer(labels)
  if (length(unique(cl)) < 2L || nrow(U) < 3L) return(NA_real_)
  sil <- silhouette(cl, dist = Rfast::Dist(as.matrix(U), method = "euclidean"))
  summary(sil)$avg.width
}

continuous <- intersect(c("RNAseq", "CNV", "Methylation", "miRNA"), names(input_full))

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

  p_by_mod <- sapply(mods_use, function(m) nrow(input_sub[[m]]))
  p_total <- sum(p_by_mod, na.rm = TRUE)

  # weights per modality
  cont_weights <- lapply(intersect(continuous, mods_use), function(m) mad_based_weights(as.matrix(input_sub[[m]])))
  names(cont_weights) <- intersect(continuous, mods_use)

  weights <- cont_weights
  if ("SNPs" %in% mods_use) {
    bin_w <- weighted_mutations(
      features = rownames(input_sub[["SNPs"]]),
      drivers = COSMIC_BC_drivers,
      driver_weight = 0.8,
      nondriver_weight = 0.2
    )
    weights[["SNPs"]] <- bin_w
  }

  # align weights to input order used in CIMLR.weight_mod
  X_use <- input_sub[mods_use]
  weights <- weights[names(X_use)]

  # methods aligned to X_use order
  methods <- vapply(names(X_use), function(m) if (m == "SNPs") "binary" else "sqeuclidean", character(1))

  if (n_samples < c_fixed) {
    # cannot form 8 clusters meaningfully; record NA like your small-n guards elsewhere
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
      wMKL_c = as.integer(c_fixed),
      wMKL_k = as.integer(k_fixed),
      ARI_to_Ground_Truth = NA_real_,
      Average_Silhouette_Width = NA_real_,
      wMKL_Elapsed_Seconds = NA_real_,
      wMKL_PeakRAM_MiB = NA_real_
    )
    rm(input_sub, X_use, weights, methods); gc()
    next
  }

  cluster_results <- NULL
  gc()
  pr <- peakRAM::peakRAM({
    cluster_results <<- CIMLR.weight_mod(
      X = X_use,
      c = c_fixed,
      cores.ratio = 0,
      k = k_fixed,
      weight = weights,
      methods = methods
    )
  })
  elapsed <- as.numeric(pr$Elapsed_Time[1])
  peak_mib <- as.numeric(pr$Peak_RAM_Used_MiB[1])

  # assignments
  clusters <- data.frame(
    Sample.ID = gsub("\\.", "-", sample_ids),
    Cluster_pred = as.integer(cluster_results$y_spectral),
    stringsAsFactors = FALSE
  )

  ari_gt <- ari_to_ground_truth(ground_truth_labels, clusters)

  # eigenspace for silhouette
  U <- extract_eigenspace(affinity = cluster_results$S, K = c_fixed)$U
  rownames(U) <- sample_ids
  colnames(U) <- paste0("eig", seq_len(c_fixed))
  avg_width <- safe_avg_sil(clusters$Cluster_pred, U)

  clust_file <- file.path(out_root, sprintf("%s_clusters_%dpct.tsv.gz", algorithm, pct))
  fwrite(as.data.table(clusters), clust_file, sep = "\t", quote = FALSE, compress = "gzip")

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
    wMKL_c = as.integer(c_fixed),
    wMKL_k = as.integer(k_fixed),
    ARI_to_Ground_Truth = ari_gt,
    Average_Silhouette_Width = avg_width,
    wMKL_Elapsed_Seconds = elapsed,
    wMKL_PeakRAM_MiB = peak_mib
  )

  rm(input_sub, X_use, weights, methods, cluster_results, clusters, U)
  gc()
}

perf_dt <- rbindlist(perf_rows, use.names = TRUE, fill = TRUE)
perf_file <- file.path(out_root, sprintf("%s_feature_perturbations_performance.tsv", algorithm))
fwrite(perf_dt, perf_file, sep = "\t", quote = FALSE, na = "NA")

writeLines(capture.output(sessionInfo()),
           file.path(out_root, sprintf("%s_session_info.txt", algorithm)))

message("Wrote: ", perf_file)

