#!/usr/bin/env Rscript

# ---- packages ----
packages <- c(
  "dplyr", "data.table", "openxlsx", "peakRAM",
  "cluster", "mclust", "fastcluster", "vegan", "coca"
)

invisible(lapply(packages, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cran.r-project.org")
  }
}))

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(openxlsx)
  library(peakRAM)
  library(cluster)
  library(mclust)
  library(fastcluster)
  library(vegan)
  library(coca)
})

RNGversion("4.2.2")
set.seed(123)

algorithm <- "COCA"
centiles <- c(0.10, 0.20, 0.50, 0.75, 0.90)

# COCA params
optk <- 2L
maxK <- 10L
methods <- "hclust"

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

# --- helpers ---
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

ari_to_ground_truth <- function(gt_df, pred_df) {
  merged <- merge(gt_df, pred_df, by = "Sample.ID")
  if (nrow(merged) == 0) return(NA_real_)
  adjustedRandIndex(merged[["Cluster_gt"]], merged[["Cluster_pred"]])
}

# SNPs -> binary, others -> euclidean
coca_distances_for <- function(mod_names) {
  vapply(mod_names, function(m) if (identical(m, "SNPs")) "binary" else "euclidean", character(1))
}

# buildMOC expects N x P (samples x features), so align by rownames (samples)
align_by_samples <- function(X_list) {
  ids_list <- lapply(X_list, rownames)
  common <- Reduce(intersect, ids_list)
  if (length(common) == 0) stop("No common sample IDs across modalities after transpose.")
  common <- common[common %in% ids_list[[1]]]
  lapply(X_list, function(m) m[common, , drop = FALSE])
}

# Transpose mm_input matrices: (features x samples) -> (samples x features)
to_samples_x_features <- function(mat) {
  stopifnot(!is.null(rownames(mat)), !is.null(colnames(mat)))
  tm <- t(as.matrix(mat))
  # tm rownames become original colnames (sample IDs)
  # tm colnames become original rownames (feature IDs)
  tm
}

preferred_order <- c("SNPs", "RNAseq", "CNV", "Methylation", "miRNA")
non_null_mods_full <- names(input_full)[!vapply(input_full, is.null, logical(1))]
modalities_ranked <- intersect(unique(feat_rank$modality), non_null_mods_full)
if (length(modalities_ranked) == 0) stop("No overlap between ranked modalities and input modalities.")
modalities_ranked <- unique(c(intersect(preferred_order, modalities_ranked),
                              setdiff(modalities_ranked, preferred_order)))

perf_rows <- vector("list", length(centiles))

for (i in seq_along(centiles)) {
  centile <- centiles[i]
  pct <- as.integer(round(centile * 100))
  message(sprintf("[%s] %s feature subset %d%%", Sys.time(), algorithm, pct))

  top_feats <- get_top_features(feat_rank[modality %in% modalities_ranked], centile)
  input_sub <- subset_input_features(input_full, top_feats)

  # Drop empty modalities
  X_use_raw <- input_sub[modalities_ranked]
  X_use_raw <- Filter(function(m) !is.null(m) && nrow(m) > 0L && ncol(m) > 0L, X_use_raw)
  if (length(X_use_raw) == 0) stop("All modalities empty after feature subsetting at ", pct, "%")

  # counts in original orientation (features are rows)
  p_by_mod <- sapply(modalities_ranked, function(m) if (!is.null(input_sub[[m]])) nrow(input_sub[[m]]) else NA_integer_)
  p_total <- sum(p_by_mod, na.rm = TRUE)

  # transpose for COCA
  X_use <- lapply(X_use_raw, to_samples_x_features)
  X_use <- align_by_samples(X_use)

  mod_names <- names(X_use)
  distances <- coca_distances_for(mod_names)
  M_use <- length(X_use)

  sample_ids <- rownames(X_use[[1]])
  n_samples <- length(sample_ids)

  if (n_samples < 3L) {
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
      n_clusters = as.integer(optk),
      COCA_M = as.integer(M_use),
      COCA_maxK = as.integer(maxK),
      COCA_method = methods,
      COCA_distances = paste(distances, collapse = ","),
      ARI_to_Ground_Truth = NA_real_,
      Average_Silhouette_Width = NA_real_,
      COCA_Elapsed_Seconds = NA_real_,
      COCA_PeakRAM_MiB = NA_real_
    )
    rm(input_sub, X_use_raw, X_use); gc()
    next
  }

  # ---- measure ONLY buildMOC() ----
  COCA_moc <- NULL
  gc()
  pr <- peakRAM::peakRAM({
    COCA_moc <<- coca::buildMOC(
      data = X_use,
      M = M_use,
      maxK = maxK,
      methods = methods,
      distances = distances
    )
  })
  coca_elapsed <- as.numeric(pr$Elapsed_Time[1])
  coca_peak_mib <- as.numeric(pr$Peak_RAM_Used_MiB[1])

  moc <- as.matrix(COCA_moc$moc)
  rm(COCA_moc); gc()

  if (is.null(rownames(moc)) && nrow(moc) == n_samples) rownames(moc) <- sample_ids

  # ---- downstream clustering (NOT benchmarked) ----
  vgd <- vegan::vegdist(moc, method = "jaccard")
  hcs <- fastcluster::hclust(vgd, method = "ward.D")
  labels <- cutree(hcs, k = optk)

  clusters <- data.frame(
    Sample.ID = rownames(moc),
    Cluster_pred = as.integer(labels),
    stringsAsFactors = FALSE
  )
  clusters$Sample.ID <- gsub("\\.", "-", clusters$Sample.ID)

  avg_width <- NA_real_
  if (length(unique(labels)) > 1L) {
    sil <- cluster::silhouette(labels, vgd)
    avg_width <- mean(sil[, "sil_width"])
  }

  ari_gt <- ari_to_ground_truth(ground_truth_labels, clusters)

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
    n_clusters = as.integer(optk),
    COCA_M = as.integer(M_use),
    COCA_maxK = as.integer(maxK),
    COCA_method = methods,
    COCA_distances = paste(distances, collapse = ","),
    ARI_to_Ground_Truth = ari_gt,
    Average_Silhouette_Width = avg_width,
    COCA_Elapsed_Seconds = coca_elapsed,
    COCA_PeakRAM_MiB = coca_peak_mib
  )

  rm(input_sub, X_use_raw, X_use, moc, vgd, hcs, labels, clusters)
  gc()
}

perf_dt <- rbindlist(perf_rows, use.names = TRUE, fill = TRUE)
perf_file <- file.path(out_root, sprintf("%s_feature_perturbations_performance.tsv", algorithm))
fwrite(perf_dt, perf_file, sep = "\t", quote = FALSE, na = "NA")

writeLines(capture.output(sessionInfo()),
           file.path(out_root, sprintf("%s_session_info.txt", algorithm)))

message("Wrote: ", perf_file)

