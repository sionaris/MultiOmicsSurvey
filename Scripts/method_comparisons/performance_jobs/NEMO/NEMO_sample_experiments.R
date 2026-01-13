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

# Fixed NEMO params (optimal run)
num_neighbors <- 10L
sigma_nemo <- 0.5
n_clusters <- 2L

input_path  <- "Resources/mm_input.rds"
subsets_path <- "Resources/Sample_perturbations/sample_subsets.tsv.gz"

custom_fun_path <- "Resources/custom_performance_functions.R"
if (!file.exists(custom_fun_path)) stop("Missing custom performance functions: ", custom_fun_path)
source(custom_fun_path)
if (!exists("nemo.affinity.graph_mod")) stop("nemo.affinity.graph_mod not found after sourcing: ", custom_fun_path)

job_out <- Sys.getenv("BENCH_OUT_DIR", "")
if (nzchar(job_out)) {
  out_root <- job_out
} else {
  job_id <- Sys.getenv("SLURM_JOB_ID", "")
  out_root <- file.path("Results", "Performance", "Sample_perturbations", algorithm,
                        if (nzchar(job_id)) paste0("job_", job_id) else "local")
}
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_path)) stop("Missing input RDS: ", input_path)
if (!file.exists(subsets_path)) stop("Missing sample subsets: ", subsets_path)

input_full <- readRDS(input_path)

gt_path <- file.path("Resources/algorithms", paste0(algorithm, "_ground_truth_labels.xlsx"))
if (!file.exists(gt_path)) stop("Missing ground truth file: ", gt_path)

ground_truth_labels <- openxlsx::read.xlsx(gt_path)
if (!all(c("Sample.ID", "Cluster") %in% colnames(ground_truth_labels))) {
  stop("Ground truth file must contain columns: Sample.ID, Cluster")
}
ground_truth_labels$Sample.ID <- gsub("\\.", "-", ground_truth_labels$Sample.ID)
colnames(ground_truth_labels)[colnames(ground_truth_labels) == "Cluster"] <- "Cluster_gt"

sample_subsets <- data.table::fread(subsets_path)
need_cols <- c("fraction", "replicate", "sample_id")
if (!all(need_cols %in% colnames(sample_subsets))) {
  stop("sample_subsets.tsv.gz must contain columns: fraction, replicate, sample_id")
}

non_null_mods <- names(input_full)[!vapply(input_full, is.null, logical(1))]
preferred_order <- c("RNAseq", "CNV", "Methylation", "miRNA", "SNPs")
modalities <- preferred_order[preferred_order %in% non_null_mods]
if (length(modalities) == 0) modalities <- non_null_mods

ref_mod <- if ("RNAseq" %in% modalities) "RNAseq" else modalities[1]
valid_ids <- colnames(input_full[[ref_mod]])
sample_subsets <- sample_subsets[sample_id %in% valid_ids]
if (nrow(sample_subsets) == 0) stop("No sample IDs in sample_subsets match the input colnames.")

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

align_common_samples <- function(input_obj, mods) {
  ids_list <- lapply(mods, function(m) colnames(input_obj[[m]]))
  common <- Reduce(intersect, ids_list)
  if (length(common) < 3L) return(NULL)
  ref <- if ("RNAseq" %in% mods) "RNAseq" else mods[1]
  common <- common[common %in% colnames(input_obj[[ref]])]
  for (m in mods) input_obj[[m]] <- input_obj[[m]][, common, drop = FALSE]
  list(obj = input_obj, sample_ids = common)
}

pairs <- unique(sample_subsets[, .(fraction, replicate)])
setorder(pairs, fraction, replicate)

p_by_mod_full <- sapply(modalities, function(m) nrow(input_full[[m]]))
p_total_full <- sum(p_by_mod_full, na.rm = TRUE)

perf_rows <- vector("list", nrow(pairs))

for (i in seq_len(nrow(pairs))) {
  frac <- as.integer(pairs$fraction[i])
  repi <- as.integer(pairs$replicate[i])

  keep_ids <- unique(sample_subsets[fraction == frac & replicate == repi, sample_id])

  input_sub <- subset_input_cols(input_full[modalities], keep_ids)
  aligned <- align_common_samples(input_sub, modalities)

  if (is.null(aligned)) {
    perf_rows[[i]] <- data.table(
      Algorithm = algorithm,
      Sample_Fraction = frac,
      Replicate = repi,
      n_samples = as.integer(0),
      p_total = as.integer(p_total_full),
      NEMO_num_neighbors = as.integer(num_neighbors),
      NEMO_sigma = as.numeric(sigma_nemo),
      n_clusters = as.integer(n_clusters),
      ARI_to_Ground_Truth = NA_real_,
      Average_Silhouette_Width = NA_real_,
      NEMO_Elapsed_Seconds = NA_real_,
      NEMO_PeakRAM_MiB = NA_real_
    )
    next
  }

  input_sub <- aligned$obj
  sample_ids <- aligned$sample_ids
  n_samples <- length(sample_ids)

  message(sprintf("[%s] %s sample subset %d%% rep %d (n=%d)", Sys.time(), algorithm, frac, repi, n_samples))

  binary_flags <- ifelse(modalities == "SNPs", "Yes", "No")

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

  if (is.null(aff_mat)) stop("nemo.affinity.graph_mod returned NULL for frac=", frac, " rep=", repi)
  if (is.null(colnames(aff_mat))) colnames(aff_mat) <- sample_ids
  if (is.null(rownames(aff_mat))) rownames(aff_mat) <- sample_ids

  res <- spectralClustering_eig(aff_mat, K = min(n_clusters, n_samples))

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

  clust_file <- file.path(out_root, sprintf("%s_clusters_samples_%dpct_rep%02d.tsv.gz", algorithm, frac, repi))
  data.table::fwrite(as.data.table(clusters), clust_file, sep = "\t", quote = FALSE, compress = "gzip")

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
    NEMO_num_neighbors = as.integer(num_neighbors),
    NEMO_sigma = as.numeric(sigma_nemo),
    n_clusters = as.integer(min(n_clusters, n_samples)),
    ARI_to_Ground_Truth = as.numeric(ari_gt),
    Average_Silhouette_Width = as.numeric(avg_width),
    NEMO_Elapsed_Seconds = as.numeric(nemo_elapsed),
    NEMO_PeakRAM_MiB = as.numeric(nemo_peak_mib)
  )

  rm(input_sub, aligned, aff_mat, res, sil, clusters)
  gc()
}

perf_dt <- data.table::rbindlist(perf_rows, use.names = TRUE, fill = TRUE)
perf_file <- file.path(out_root, sprintf("%s_sample_perturbations_performance.tsv", algorithm))
data.table::fwrite(perf_dt, perf_file, sep = "\t", quote = FALSE, na = "NA")

writeLines(capture.output(sessionInfo()),
           file.path(out_root, sprintf("%s_session_info.txt", algorithm)))

message("Wrote: ", perf_file)

