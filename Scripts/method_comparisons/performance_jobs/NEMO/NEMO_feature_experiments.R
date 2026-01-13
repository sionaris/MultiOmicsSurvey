#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

req <- c("data.table", "openxlsx", "Rfast", "cluster", "mclust", "peakRAM", "SNFtool", "NEMO")
missing <- req[!vapply(req, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop(
    "Missing R packages: ", paste(missing, collapse = ", "), "\n",
    "Install once in your R library (do not install inside jobs).\n",
    "For NEMO: devtools::install_github('Shamir-Lab/NEMO/NEMO')"
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(openxlsx)
  library(Rfast)
  library(cluster)
  library(mclust)
  library(peakRAM)
  library(SNFtool)
  library(NEMO)
})

RNGversion("4.2.2")
set.seed(123)

algorithm <- "NEMO"

# Fixed centiles (fractions, not percents)
centiles <- c(0.10, 0.20, 0.50, 0.75, 0.90)

# Fixed NEMO params (optimal run)
num_neighbors <- 10L
sigma_nemo <- 0.5
n_clusters <- 2L

input_path <- "Resources/mm_input.rds"
rank_path  <- "Resources/Feature_perturbations/feature_rankings.tsv.gz"

custom_fun_path <- "Resources/custom_performance_functions.R"
if (!file.exists(custom_fun_path)) stop("Missing custom performance functions: ", custom_fun_path)
source(custom_fun_path)
if (!exists("nemo.affinity.graph_mod")) stop("nemo.affinity.graph_mod not found after sourcing: ", custom_fun_path)

# Output root (job-safe)
job_out <- Sys.getenv("BENCH_OUT_DIR", "")
if (nzchar(job_out)) {
  out_root <- job_out
} else {
  job_id <- Sys.getenv("SLURM_JOB_ID", "")
  out_root <- file.path("Results", "Performance", "Feature_perturbations", algorithm,
                        if (nzchar(job_id)) paste0("job_", job_id) else "local")
}
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_path)) stop("Missing input RDS: ", input_path)
if (!file.exists(rank_path)) stop("Missing feature rankings: ", rank_path)

input_full <- readRDS(input_path)

gt_path <- file.path("Resources/algorithms", paste0(algorithm, "_ground_truth_labels.xlsx"))
if (!file.exists(gt_path)) stop("Missing ground truth file: ", gt_path)

ground_truth_labels <- openxlsx::read.xlsx(gt_path)
if (!all(c("Sample.ID", "Cluster") %in% colnames(ground_truth_labels))) {
  stop("Ground truth file must contain columns: Sample.ID, Cluster")
}
ground_truth_labels$Sample.ID <- gsub("\\.", "-", ground_truth_labels$Sample.ID)
colnames(ground_truth_labels)[colnames(ground_truth_labels) == "Cluster"] <- "Cluster_gt"

feat_rank <- data.table::fread(rank_path)
if (!all(c("modality", "rank", "feature_id") %in% colnames(feat_rank))) {
  stop("feature_rankings.tsv.gz must contain columns: modality, rank, feature_id")
}

modalities_ranked <- intersect(unique(feat_rank$modality), names(input_full))
if (length(modalities_ranked) == 0) stop("No overlap between ranked modalities and input modalities.")

# Keep a stable modality order (helps binary_flags alignment)
preferred_order <- c("RNAseq", "CNV", "Methylation", "miRNA", "SNPs")
modalities <- preferred_order[preferred_order %in% modalities_ranked]
if (length(modalities) == 0) modalities <- modalities_ranked

get_top_features <- function(feat_rank_dt, centile) {
  mods <- unique(feat_rank_dt$modality)
  out <- vector("list", length(mods)); names(out) <- mods
  for (mod in mods) {
    mod_dt <- feat_rank_dt[modality == mod][order(rank)]
    if (nrow(mod_dt) == 0) stop("No ranked features for modality: ", mod)
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
      if (length(keep) == 0) stop("After feature subsetting, 0 features kept for modality: ", mod)
      out[[mod]] <- out[[mod]][keep, , drop = FALSE]
    }
  }
  out
}

align_common_samples <- function(input_obj, mods) {
  ids_list <- lapply(mods, function(m) colnames(input_obj[[m]]))
  common <- Reduce(intersect, ids_list)
  if (length(common) < 3L) stop("Too few common samples across modalities: ", length(common))
  ref <- if ("RNAseq" %in% mods) "RNAseq" else mods[1]
  common <- common[common %in% colnames(input_obj[[ref]])]
  for (m in mods) input_obj[[m]] <- input_obj[[m]][, common, drop = FALSE]
  list(obj = input_obj, sample_ids = common)
}

# Required spectral clustering (your adapted version)
spectralClustering_eig <- function (affinity, K, type = 3) {
  d = rowSums(affinity)
  d[d == 0] = .Machine$double.eps
  D = diag(d)
  L = D - affinity
  if (type == 1) {
    NL = L
  } else if (type == 2) {
    Di = diag(1/d)
    NL = Di %*% L
  } else if (type == 3) {
    Di = diag(1/sqrt(d))
    NL = Di %*% L %*% Di
  }
  eig = eigen(NL)
  res = sort(abs(eig$values), index.return = TRUE)
  U = eig$vectors[, res$ix[1:K]]
  normalize <- function(x) x/sqrt(sum(x^2))
  if (type == 3) {
    U = t(apply(U, 1, normalize))
  }
  eigDiscrete = SNFtool:::.discretisation(U)
  eigDiscrete = eigDiscrete$discrete
  labels = apply(eigDiscrete, 1, which.max)
  U = as.data.frame(cbind(U, labels))
  U$Sample.ID = colnames(affinity)
  colnames(U)[(ncol(U)-1):ncol(U)] = c("Cluster", "Sample.ID")
  return(U)
}

ari_to_ground_truth <- function(gt_df, pred_df) {
  merged <- merge(gt_df, pred_df, by = "Sample.ID")
  if (nrow(merged) == 0) return(NA_real_)
  mclust::adjustedRandIndex(merged[["Cluster_gt"]], merged[["Cluster_pred"]])
}

perf_rows <- vector("list", length(centiles))

for (i in seq_along(centiles)) {
  centile <- centiles[i]
  pct <- as.integer(round(centile * 100))

  message(sprintf("[%s] %s feature subset %d%%", Sys.time(), algorithm, pct))

  top_feats <- get_top_features(feat_rank[modality %in% modalities], centile)
  input_sub <- subset_input_features(input_full, top_feats)

  # Align samples across modalities + keep deterministic ordering
  aligned <- align_common_samples(input_sub, modalities)
  input_sub <- aligned$obj
  sample_ids <- aligned$sample_ids
  n_samples <- length(sample_ids)

  p_by_mod <- sapply(modalities, function(m) nrow(input_sub[[m]]))
  p_total <- sum(p_by_mod, na.rm = TRUE)

  # Build binary_flags aligned to modalities
  binary_flags <- ifelse(modalities == "SNPs", "Yes", "No")

  # PeakRAM around affinity computation (dominant cost)
  aff_mat <- NULL
  pr <- peakRAM::peakRAM({
    aff_mat <<- nemo.affinity.graph_mod(
      input_sub,
      k = num_neighbors,
      sigma = sigma_nemo,
      binary_flags = binary_flags,
      binary_distance = "binary"
    )
  })
  nemo_elapsed <- as.numeric(pr$Elapsed_Time[1])
  nemo_peak_mib <- as.numeric(pr$Peak_RAM_Used_MiB[1])

  if (is.null(aff_mat)) stop("nemo.affinity.graph_mod returned NULL.")
  if (is.null(colnames(aff_mat))) colnames(aff_mat) <- sample_ids
  if (is.null(rownames(aff_mat))) rownames(aff_mat) <- sample_ids

  # Cluster
  res <- spectralClustering_eig(aff_mat, K = n_clusters)

  # Silhouette (on eigenvectors)
  col_index <- ncol(res) - 2
  sil <- cluster::silhouette(
    as.integer(res$Cluster),
    dist = Rfast::Dist(res[, 1:col_index], method = "euclidean")
  )
  avg_width <- summary(sil)$avg.width

  clusters <- res[, c("Sample.ID", "Cluster")]
  clusters$Sample.ID <- gsub("\\.", "-", clusters$Sample.ID)
  colnames(clusters)[2] <- "Cluster_pred"

  ari_gt <- ari_to_ground_truth(ground_truth_labels, clusters)

  clust_file <- file.path(out_root, sprintf("%s_clusters_%dpct.tsv.gz", algorithm, pct))
  data.table::fwrite(as.data.table(clusters), clust_file, sep = "\t", quote = FALSE, compress = "gzip")

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
    NEMO_num_neighbors = as.integer(num_neighbors),
    NEMO_sigma = as.numeric(sigma_nemo),
    n_clusters = as.integer(n_clusters),
    ARI_to_Ground_Truth = as.numeric(ari_gt),
    Average_Silhouette_Width = as.numeric(avg_width),
    NEMO_Elapsed_Seconds = as.numeric(nemo_elapsed),
    NEMO_PeakRAM_MiB = as.numeric(nemo_peak_mib)
  )

  rm(input_sub, aligned, aff_mat, res, sil, clusters)
  gc()
}

perf_dt <- data.table::rbindlist(perf_rows, use.names = TRUE, fill = TRUE)
perf_file <- file.path(out_root, sprintf("%s_feature_perturbations_performance.tsv", algorithm))
data.table::fwrite(perf_dt, perf_file, sep = "\t", quote = FALSE, na = "NA")

writeLines(capture.output(sessionInfo()),
           file.path(out_root, sprintf("%s_session_info.txt", algorithm)))

message("Wrote: ", perf_file)

