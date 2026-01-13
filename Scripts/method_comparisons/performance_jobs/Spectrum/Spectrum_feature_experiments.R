#!/usr/bin/env Rscript

# ---- packages ----
packages <- c("dplyr", "data.table", "openxlsx", "Rfast", "cluster", "mclust", "peakRAM", "Spectrum")
invisible(lapply(packages, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cran.r-project.org")
  }
}))

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(openxlsx)
  library(Rfast)
  library(cluster)
  library(mclust)
  library(peakRAM)
  library(Spectrum)
})

RNGversion("4.2.2")
set.seed(123)

algorithm <- "Spectrum"

# Fixed centiles (fractions, not percents)
centiles <- c(0.10, 0.20, 0.50, 0.75, 0.90)

# Spectrum parameters (fixed)
spec_method <- 2L
spec_diffusion <- TRUE
spec_kerneltype <- "density"
spec_maxk <- 10L
spec_NN <- 3L
spec_NN2 <- 7L
spec_frac <- 2
spec_thresh <- 7L
spec_tunekernel <- TRUE
spec_clusteralg <- "km"
spec_diffusion_iters <- 5L
spec_KNNs_p <- 20L
spec_fontsize <- 18
spec_dotsize <- 1.25

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

ari_to_ground_truth <- function(gt_df, pred_df) {
  merged <- merge(gt_df, pred_df, by = "Sample.ID")
  if (nrow(merged) == 0) return(NA_real_)
  adjustedRandIndex(merged[["Cluster_gt"]], merged[["Cluster_pred"]])
}

ref_modality <- function(input_list, prefer = "RNAseq") {
  non_null <- names(input_list)[!vapply(input_list, is.null, logical(1))]
  if (length(non_null) == 0) stop("All modalities are NULL.")
  if (prefer %in% non_null) return(prefer)
  non_null[1]
}

safe_silhouette_avg <- function(assignments, embeddings) {
  cl <- as.integer(assignments)
  if (length(unique(cl)) < 2L || nrow(embeddings) < 3L) return(NA_real_)
  sil <- silhouette(cl, dist = Rfast::Dist(as.matrix(embeddings), method = "euclidean"))
  summary(sil)$avg.width
}

# ---- main loop ----
perf_rows <- vector("list", length(centiles))

for (i in seq_along(centiles)) {
  centile <- centiles[i]
  pct <- as.integer(round(centile * 100))

  message(sprintf("[%s] %s feature subset %d%%", Sys.time(), algorithm, pct))

  top_feats <- get_top_features(feat_rank[modality %in% modalities_ranked], centile)
  input_sub <- subset_input_features(input_full, top_feats)

  # counts
  p_by_mod <- sapply(modalities_ranked, function(m) if (!is.null(input_sub[[m]])) nrow(input_sub[[m]]) else NA_integer_)
  m0 <- ref_modality(input_sub)
  n_samples <- ncol(input_sub[[m0]])
  p_total <- sum(p_by_mod, na.rm = TRUE)

  # prepare for Spectrum: list of data.frames, samples as columns, features as rows
  data_for_spectrum <- lapply(input_sub[modalities_ranked], function(x) {
    if (is.null(x)) return(NULL)
    as.data.frame(as.matrix(x))
  })
  data_for_spectrum <- data_for_spectrum[!vapply(data_for_spectrum, is.null, logical(1))]

  gc()
  pr <- peakRAM::peakRAM({
    spec_res <- Spectrum::Spectrum(
      data_for_spectrum,
      method = spec_method,
      silent = TRUE,
      showres = FALSE,
      diffusion = spec_diffusion,
      kerneltype = spec_kerneltype,
      maxk = spec_maxk,
      NN = spec_NN,
      NN2 = spec_NN2,
      frac = spec_frac,
      thresh = spec_thresh,
      tunekernel = spec_tunekernel,
      clusteralg = spec_clusteralg,
      diffusion_iters = spec_diffusion_iters,
      KNNs_p = spec_KNNs_p,
      fontsize = spec_fontsize,
      dotsize = spec_dotsize,
      missing = FALSE
    )
  })
  spec_elapsed <- as.numeric(pr$Elapsed_Time[1])
  spec_peak_mib <- as.numeric(pr$Peak_RAM_Used_MiB[1])

  K_found <- as.integer(spec_res$K)
  embeddings <- as.data.frame(spec_res$eigensystem$vectors)[, seq_len(K_found), drop = FALSE]
  avg_width <- safe_silhouette_avg(spec_res$assignments, embeddings)

  # clusters (use reference modality colnames as sample IDs)
  sample_ids <- colnames(input_sub[[m0]])
  clusters <- data.frame(
    Sample.ID = gsub("\\.", "-", sample_ids),
    Cluster_pred = as.integer(spec_res$assignments),
    stringsAsFactors = FALSE
  )

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
    Spectrum_method = spec_method,
    Spectrum_kerneltype = spec_kerneltype,
    Spectrum_diffusion = spec_diffusion,
    Spectrum_maxk = spec_maxk,
    Spectrum_NN = spec_NN,
    Spectrum_NN2 = spec_NN2,
    Spectrum_frac = spec_frac,
    Spectrum_thresh = spec_thresh,
    Spectrum_tunekernel = spec_tunekernel,
    Spectrum_clusteralg = spec_clusteralg,
    Spectrum_diffusion_iters = spec_diffusion_iters,
    Spectrum_KNNs_p = spec_KNNs_p,
    Spectrum_K_found = K_found,
    ARI_to_Ground_Truth = ari_gt,
    Average_Silhouette_Width = avg_width,
    Spectrum_Elapsed_Seconds = spec_elapsed,
    Spectrum_PeakRAM_MiB = spec_peak_mib
  )

  rm(input_sub, data_for_spectrum, spec_res, embeddings, clusters)
  gc()
}

perf_dt <- rbindlist(perf_rows, use.names = TRUE, fill = TRUE)

perf_file <- file.path(out_root, sprintf("%s_feature_perturbations_performance.tsv", algorithm))
fwrite(perf_dt, perf_file, sep = "\t", quote = FALSE, na = "NA")

writeLines(capture.output(sessionInfo()),
           file.path(out_root, sprintf("%s_session_info.txt", algorithm)))

message("Wrote: ", perf_file)

