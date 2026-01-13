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

# Fixed wMKL run params
c_fixed <- 8L
k_fixed <- 32L

input_path <- "Resources/mm_input.rds"
subsets_path <- "Resources/Sample_perturbations/sample_subsets.tsv.gz"
cosmic_path <- "Resources/COSMIC_CGC_Breast_somatic.csv"
custom_fns_path <- "Resources/custom_performance_functions.R"

# Make outputs job-safe by default (can be overridden)
job_out <- Sys.getenv("BENCH_OUT_DIR", "")
if (nzchar(job_out)) {
  out_root <- job_out
} else {
  job_id <- Sys.getenv("SLURM_JOB_ID", "")
  out_root <- file.path("Results", "Performance", "Sample_perturbations", algorithm,
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
if (length(non_null_mods) == 0) stop("All modalities in mm_input.rds are NULL.")

ref_mod_full <- if ("RNAseq" %in% non_null_mods) "RNAseq" else non_null_mods[1]
valid_ids <- colnames(input_full[[ref_mod_full]])

sample_subsets <- sample_subsets[sample_id %in% valid_ids]
if (nrow(sample_subsets) == 0) stop("No sample IDs in sample_subsets match the input colnames.")

# ---- weights (features in rows, samples in columns) ----
mad_based_weights <- function(X) {
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

pairs <- unique(sample_subsets[, .(fraction, replicate)])
setorder(pairs, fraction, replicate)

continuous <- intersect(c("RNAseq", "CNV", "Methylation", "miRNA"), non_null_mods)

p_by_mod_full <- sapply(non_null_mods, function(m) nrow(input_full[[m]]))
p_total_full <- sum(p_by_mod_full, na.rm = TRUE)

perf_rows <- vector("list", nrow(pairs))

for (i in seq_len(nrow(pairs))) {
  frac <- pairs$fraction[i]
  repi <- pairs$replicate[i]

  keep_ids <- unique(sample_subsets[fraction == frac & replicate == repi, sample_id])
  input_sub <- subset_input_cols(input_full, keep_ids)

  mods_use <- names(input_sub)[!vapply(input_sub, is.null, logical(1))]
  if (length(mods_use) == 0L) {
    perf_rows[[i]] <- data.table(
      Algorithm = algorithm,
      Sample_Fraction = frac,
      Replicate = repi,
      n_samples = 0L,
      p_total = as.integer(p_total_full),
      wMKL_c = as.integer(c_fixed),
      wMKL_k = as.integer(k_fixed),
      ARI_to_Ground_Truth = NA_real_,
      Average_Silhouette_Width = NA_real_,
      wMKL_Elapsed_Seconds = NA_real_,
      wMKL_PeakRAM_MiB = NA_real_
    )
    rm(input_sub); gc()
    next
  }

  ref_mod <- if ("RNAseq" %in% mods_use) "RNAseq" else mods_use[1]
  sample_ids <- colnames(input_sub[[ref_mod]])
  n_samples <- length(sample_ids)

  message(sprintf("[%s] %s sample subset %d%% rep %d (n=%d)", Sys.time(), algorithm, frac, repi, n_samples))

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

  X_use <- input_sub[mods_use]
  weights <- weights[names(X_use)]
  methods <- vapply(names(X_use), function(m) if (m == "SNPs") "binary" else "sqeuclidean", character(1))

  if (n_samples < c_fixed) {
    perf_rows[[i]] <- data.table(
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

  clusters <- data.frame(
    Sample.ID = gsub("\\.", "-", sample_ids),
    Cluster_pred = as.integer(cluster_results$y_spectral),
    stringsAsFactors = FALSE
  )

  ari_gt <- ari_to_ground_truth(ground_truth_labels, clusters)

  U <- extract_eigenspace(affinity = cluster_results$S, K = c_fixed)$U
  rownames(U) <- sample_ids
  colnames(U) <- paste0("eig", seq_len(c_fixed))
  avg_width <- safe_avg_sil(clusters$Cluster_pred, U)

  clust_file <- file.path(out_root, sprintf("%s_clusters_samples_%dpct_rep%02d.tsv.gz", algorithm, frac, repi))
  fwrite(as.data.table(clusters), clust_file, sep = "\t", quote = FALSE, compress = "gzip")

  perf_rows[[i]] <- data.table(
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
perf_file <- file.path(out_root, sprintf("%s_sample_perturbations_performance.tsv", algorithm))
fwrite(perf_dt, perf_file, sep = "\t", quote = FALSE, na = "NA")

writeLines(capture.output(sessionInfo()),
           file.path(out_root, sprintf("%s_session_info.txt", algorithm)))

message("Wrote: ", perf_file)

