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
if (length(non_null_mods) == 0) stop("All modalities in mm_input.rds are NULL.")

ref_mod <- if ("RNAseq" %in% non_null_mods) "RNAseq" else non_null_mods[1]
valid_ids <- colnames(input_full[[ref_mod]])

sample_subsets <- sample_subsets[sample_id %in% valid_ids]
if (nrow(sample_subsets) == 0) stop("No sample IDs in sample_subsets match the input colnames.")

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

safe_silhouette_avg <- function(assignments, embeddings) {
  cl <- as.integer(assignments)
  if (length(unique(cl)) < 2L || nrow(embeddings) < 3L) return(NA_real_)
  sil <- silhouette(cl, dist = Rfast::Dist(as.matrix(embeddings), method = "euclidean"))
  summary(sil)$avg.width
}

pairs <- unique(sample_subsets[, .(fraction, replicate)])
setorder(pairs, fraction, replicate)

p_by_mod_full <- sapply(non_null_mods, function(m) nrow(input_full[[m]]))
p_total_full <- sum(p_by_mod_full, na.rm = TRUE)

perf_rows <- vector("list", nrow(pairs))

for (i in seq_len(nrow(pairs))) {
  frac <- pairs$fraction[i]
  repi <- pairs$replicate[i]

  keep_ids <- unique(sample_subsets[fraction == frac & replicate == repi, sample_id])

  input_sub <- subset_input_cols(input_full, keep_ids)
  n_samples <- ncol(input_sub[[ref_mod]])

  if (n_samples < 3L) {
    perf_rows[[i]] <- data.table(
      Algorithm = algorithm,
      Sample_Fraction = frac,
      Replicate = repi,
      n_samples = n_samples,
      p_total = p_total_full,
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
      Spectrum_K_found = NA_integer_,
      ARI_to_Ground_Truth = NA_real_,
      Average_Silhouette_Width = NA_real_,
      Spectrum_Elapsed_Seconds = NA_real_,
      Spectrum_PeakRAM_MiB = NA_real_
    )
    next
  }

  message(sprintf("[%s] %s sample subset %d%% rep %d (n=%d)", Sys.time(), algorithm, frac, repi, n_samples))

  data_for_spectrum <- lapply(input_sub[non_null_mods], function(x) as.data.frame(as.matrix(x)))

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

  sample_ids <- colnames(input_sub[[ref_mod]])
  clusters <- data.frame(
    Sample.ID = gsub("\\.", "-", sample_ids),
    Cluster_pred = as.integer(spec_res$assignments),
    stringsAsFactors = FALSE
  )

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

perf_file <- file.path(out_root, sprintf("%s_sample_perturbations_performance.tsv", algorithm))
fwrite(perf_dt, perf_file, sep = "\t", quote = FALSE, na = "NA")

writeLines(capture.output(sessionInfo()),
           file.path(out_root, sprintf("%s_session_info.txt", algorithm)))

message("Wrote: ", perf_file)

