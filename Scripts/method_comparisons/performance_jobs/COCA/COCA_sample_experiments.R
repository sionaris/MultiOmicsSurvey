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
optk <- 2L
maxK <- 10L
methods <- "hclust"

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
preferred_order <- c("SNPs", "RNAseq", "CNV", "Methylation", "miRNA")
mods_use_full <- unique(c(intersect(preferred_order, non_null_mods),
                          setdiff(non_null_mods, preferred_order)))

ref_mod_full <- if ("RNAseq" %in% mods_use_full) "RNAseq" else mods_use_full[1]
valid_ids <- colnames(input_full[[ref_mod_full]])
sample_subsets <- sample_subsets[sample_id %in% valid_ids]
if (nrow(sample_subsets) == 0) stop("No sample IDs in sample_subsets match the input colnames.")

# --- helpers ---
ari_to_ground_truth <- function(gt_df, pred_df) {
  merged <- merge(gt_df, pred_df, by = "Sample.ID")
  if (nrow(merged) == 0) return(NA_real_)
  adjustedRandIndex(merged[["Cluster_gt"]], merged[["Cluster_pred"]])
}

coca_distances_for <- function(mod_names) {
  vapply(mod_names, function(m) if (identical(m, "SNPs")) "binary" else "euclidean", character(1))
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

to_samples_x_features <- function(mat) {
  stopifnot(!is.null(rownames(mat)), !is.null(colnames(mat)))
  t(as.matrix(mat))
}

align_by_samples <- function(X_list) {
  ids_list <- lapply(X_list, rownames)
  common <- Reduce(intersect, ids_list)
  if (length(common) == 0) stop("No common sample IDs across modalities after transpose.")
  common <- common[common %in% ids_list[[1]]]
  lapply(X_list, function(m) m[common, , drop = FALSE])
}

pairs <- unique(sample_subsets[, .(fraction, replicate)])
setorder(pairs, fraction, replicate)

p_by_mod_full <- sapply(mods_use_full, function(m) if (!is.null(input_full[[m]])) nrow(input_full[[m]]) else NA_integer_)
p_total_full <- sum(p_by_mod_full, na.rm = TRUE)
get_int1 <- function(x) if (is.null(x) || length(x) == 0L || is.na(x)) NA_integer_ else as.integer(x[1])

perf_rows <- vector("list", nrow(pairs))

for (i in seq_len(nrow(pairs))) {
  frac <- as.integer(pairs$fraction[i])
  repi <- as.integer(pairs$replicate[i])

  message(sprintf("[%s] %s sample subset %d%% rep %d", Sys.time(), algorithm, frac, repi))

  keep_ids <- unique(sample_subsets[fraction == frac & replicate == repi, sample_id])
  input_sub <- subset_input_cols(input_full, keep_ids)

  X_use_raw <- input_sub[mods_use_full]
  X_use_raw <- Filter(function(m) !is.null(m) && nrow(m) > 0L && ncol(m) > 0L, X_use_raw)
  if (length(X_use_raw) == 0) {
    perf_rows[[i]] <- data.table(
      Algorithm = algorithm,
      Sample_Fraction = frac,
      Replicate = repi,
      n_samples = 0L,
      p_total = as.integer(p_total_full),
      p_RNAseq = get_int1(p_by_mod_full[["RNAseq"]]),
      p_CNV = get_int1(p_by_mod_full[["CNV"]]),
      p_Methylation = get_int1(p_by_mod_full[["Methylation"]]),
      p_miRNA = get_int1(p_by_mod_full[["miRNA"]]),
      p_SNPs = get_int1(p_by_mod_full[["SNPs"]]),
      n_clusters = as.integer(optk),
      COCA_M = NA_integer_,
      COCA_maxK = as.integer(maxK),
      COCA_method = methods,
      COCA_distances = NA_character_,
      ARI_to_Ground_Truth = NA_real_,
      Average_Silhouette_Width = NA_real_,
      COCA_Elapsed_Seconds = NA_real_,
      COCA_PeakRAM_MiB = NA_real_
    )
    rm(input_sub); gc()
    next
  }

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
      Sample_Fraction = frac,
      Replicate = repi,
      n_samples = as.integer(n_samples),
      p_total = as.integer(p_total_full),
      p_RNAseq = get_int1(p_by_mod_full[["RNAseq"]]),
      p_CNV = get_int1(p_by_mod_full[["CNV"]]),
      p_Methylation = get_int1(p_by_mod_full[["Methylation"]]),
      p_miRNA = get_int1(p_by_mod_full[["miRNA"]]),
      p_SNPs = get_int1(p_by_mod_full[["SNPs"]]),
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
perf_file <- file.path(out_root, sprintf("%s_sample_perturbations_performance.tsv", algorithm))
fwrite(perf_dt, perf_file, sep = "\t", quote = FALSE, na = "NA")

writeLines(capture.output(sessionInfo()),
           file.path(out_root, sprintf("%s_session_info.txt", algorithm)))

message("Wrote: ", perf_file)

