#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)

req <- c("data.table","openxlsx","SNFtool","peakRAM","cluster","mclust","Rfast","kernlab",
         "foreach", "doParallel", "iterators")
missing <- req[!vapply(req, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Missing R packages: ", paste(missing, collapse = ", "),
       "\nInstall once in your R library; do not install inside jobs.")
}

suppressPackageStartupMessages({
  library(data.table)
  library(openxlsx)
  library(SNFtool)
  library(peakRAM)
  library(cluster)
  library(mclust)
  library(Rfast)
  library(kernlab)
  library(foreach)
  library(doParallel)
})

foreach::registerDoSEQ()
options(mc.cores = 1)

RNGversion("4.2.2")
set.seed(123)

`%||%` <- function(x,y) if (!is.null(x) && length(x) && !is.na(x)) x else y

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list()
  i <- 1
  while (i <= length(args)) {
    if (startsWith(args[i], "--")) {
      key <- sub("^--", "", args[i])
      val <- if (i + 1 <= length(args)) args[i + 1] else ""
      out[[key]] <- val
      i <- i + 2
    } else i <- i + 1
  }
  out
}
a <- parse_args()

task_index <- suppressWarnings(as.integer(a[["task_index"]]))
if (is.na(task_index) || task_index < 1L) stop("Invalid --task_index")

out_dir <- Sys.getenv("BENCH_OUT_DIR", "")
meta_dir <- Sys.getenv("BENCH_META_DIR", "")
run_dir <- Sys.getenv("BENCH_RUN_DIR", "")
if (!nzchar(out_dir) || !nzchar(meta_dir) || !nzchar(run_dir)) stop("Missing BENCH_* env vars.")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(meta_dir, recursive = TRUE, showWarnings = FALSE)

algorithm <- "RWR-F"

# Fixed optimal params (RWR-F)
aff_sigma <- 0.5
aff_num_neighbors <- 50L
n_clusters <- 2L
iteration_max <- 1000L
gamma_fixed <- 0.7

# Source RWR-F code
rwr_src <- Sys.getenv("RWR_SOURCE", "Resources/RWR-F_source.R")
if (!file.exists(rwr_src)) stop("Missing RWR source file: ", rwr_src)
source(rwr_src)

input_path <- Sys.getenv("MM_INPUT_RDS", "Resources/mm_input.rds")
subsets_path <- Sys.getenv("SAMPLE_SUBSETS", "Resources/Sample_perturbations/sample_subsets.tsv.gz")
if (!file.exists(input_path)) stop("Missing input RDS: ", input_path)
if (!file.exists(subsets_path)) stop("Missing sample_subsets: ", subsets_path)

# Ground truth (override if needed)
gt_override <- Sys.getenv("GROUND_TRUTH_XLSX", "")
gt_candidates <- c(
  gt_override,
  file.path("Resources/algorithms", paste0(algorithm, "_ground_truth_labels.xlsx")),
  file.path("Resources/algorithms", "RWRF_ground_truth_labels.xlsx"),
  file.path("Resources/algorithms", "RWR-F_ground_truth_labels.xlsx")
)
gt_candidates <- gt_candidates[nzchar(gt_candidates)]
gt_path <- gt_candidates[file.exists(gt_candidates)][1]
if (is.na(gt_path) || !nzchar(gt_path)) {
  stop("Missing ground truth xlsx. Tried:\n- ", paste(gt_candidates, collapse = "\n- "))
}

input_full <- readRDS(input_path)
required_mods <- intersect(c("RNAseq","CNV","Methylation","miRNA","SNPs"), names(input_full))
if (length(required_mods) < 2) stop("Too few modalities in mm_input.rds for RWR-F.")

# Determine valid ids from reference modality
ref_mod <- if ("RNAseq" %in% required_mods) "RNAseq" else required_mods[1]
valid_ids <- colnames(input_full[[ref_mod]])

sample_subsets <- data.table::fread(subsets_path)
need_cols <- c("fraction","replicate","sample_id")
if (!all(need_cols %in% names(sample_subsets))) stop("sample_subsets must have fraction, replicate, sample_id")

sample_subsets <- sample_subsets[sample_id %in% valid_ids]
pairs <- unique(sample_subsets[, .(fraction = as.integer(fraction), replicate = as.integer(replicate))])
setorder(pairs, fraction, replicate)

if (task_index > nrow(pairs)) stop("task_index out of range: ", task_index, " > ", nrow(pairs))

frac <- pairs$fraction[task_index]
repi <- pairs$replicate[task_index]

keep_ids <- unique(sample_subsets[fraction == frac & replicate == repi, sample_id])
if (length(keep_ids) < 3L) stop("Too few samples in subset for task_index=", task_index)

# Subset columns across modalities and align intersection
input_sub <- input_full
for (m in required_mods) {
  mat <- input_sub[[m]]
  kk <- intersect(colnames(mat), keep_ids)
  input_sub[[m]] <- mat[, kk, drop = FALSE]
}
ids_list <- lapply(required_mods, function(m) colnames(input_sub[[m]]))
common <- Reduce(intersect, ids_list)
if (length(common) < 3L) stop("Too few common samples across modalities after subsetting.")
if ("RNAseq" %in% required_mods) common <- common[common %in% colnames(input_sub[["RNAseq"]])]
for (m in required_mods) input_sub[[m]] <- input_sub[[m]][, common, drop = FALSE]

sample_ids <- gsub("\\.", "-", common)
n_samples <- length(sample_ids)
p_by_mod <- sapply(required_mods, function(m) nrow(input_sub[[m]]))
p_total <- as.integer(sum(p_by_mod))

# Distance list
continuous <- intersect(c("RNAseq","CNV","Methylation","miRNA"), required_mods)

distL <- list()
for (m in continuous) {
  x <- t(as.matrix(input_sub[[m]]))
  d <- SNFtool::dist2(x, x)
  dimnames(d) <- list(sample_ids, sample_ids)
  distL[[m]] <- d
}
if ("SNPs" %in% required_mods) {
  x <- t(as.matrix(input_sub[["SNPs"]]))
  d <- as.matrix(dist(x, method = "binary"))
  dimnames(d) <- list(sample_ids, sample_ids)
  distL[["SNPs"]] <- d
}

# Affinity matrices (NOT measured by peakRAM)
affinityL <- lapply(distL, function(d) {
  a <- SNFtool::affinityMatrix(as.matrix(d), K = aff_num_neighbors, sigma = aff_sigma)
  dimnames(a) <- dimnames(d)
  a
})

gc()

# --- peakRAM ONLY on fusion ---
fused_rwrf <- NULL
pr <- peakRAM::peakRAM({
  fused_rwrf <<- RWR_fusion(sim_list = affinityL, iteration_max = iteration_max, gama = gamma_fixed)
})

peak_mib <- as.numeric(pr$Peak_RAM_Used_MiB[1])
elapsed_s <- as.numeric(pr$Elapsed_Time[1])

# Clustering
specc_with_embedding <- function(Kmat, centers, iterations = 200) {
  X <- as.matrix(Kmat)
  d <- 1 / sqrt(rowSums(X))
  L <- d * X %*% diag(d)
  xi <- eigen(L)$vectors[, 1:centers, drop = FALSE]
  yi <- xi / sqrt(rowSums(xi^2))
  km <- kmeans(yi, centers = centers, iter.max = iterations)
  list(cluster = km$cluster, embedding = yi)
}

km <- fused_rwrf
dimnames(km) <- list(sample_ids, sample_ids)
res <- specc_with_embedding(km, centers = n_clusters)

sil <- silhouette(as.integer(res$cluster), dist = Rfast::Dist(res$embedding, method = "euclidean"))
avg_width <- summary(sil)$avg.width

clusters <- data.table(
  Sample.ID = sample_ids,
  Cluster_pred = as.integer(res$cluster)
)

# ARI
gt <- openxlsx::read.xlsx(gt_path)
if (!all(c("Sample.ID","Cluster") %in% colnames(gt))) stop("Ground truth must have Sample.ID, Cluster")
gt$Sample.ID <- gsub("\\.", "-", gt$Sample.ID)
setDT(gt); setnames(gt, "Cluster", "Cluster_gt")

merged <- merge(gt, clusters, by = "Sample.ID")
ari_gt <- if (nrow(merged) > 0) mclust::adjustedRandIndex(merged$Cluster_gt, merged$Cluster_pred) else NA_real_

# Outputs
data.table::fwrite(clusters, file.path(out_dir, "clusters.tsv.gz"), sep = "\t", quote = FALSE, compress = "gzip")

prep_meta <- data.table(
  Algorithm = algorithm,
  Mode = "sample",
  Sample_Percent = as.integer(frac),
  Replicate = as.integer(repi),
  n_samples = as.integer(n_samples),
  p_total = as.integer(p_total)
)
for (m in c("RNAseq","CNV","Methylation","miRNA","SNPs")) prep_meta[[paste0("p_", m)]] <- as.integer(p_by_mod[m] %||% NA_integer_)
data.table::fwrite(prep_meta, file.path(meta_dir, "prep_meta.tsv"), sep = "\t", quote = FALSE)

perf_row <- data.table(
  Algorithm = algorithm,
  Sample_Percent = as.integer(frac),
  Replicate = as.integer(repi),
  n_samples = as.integer(n_samples),
  p_total = as.integer(p_total),
  Affinity_K = as.integer(aff_num_neighbors),
  Affinity_sigma = as.numeric(aff_sigma),
  RWR_iteration_max = as.integer(iteration_max),
  RWR_gamma = as.numeric(gamma_fixed),
  n_clusters = as.integer(n_clusters),
  ARI_to_Ground_Truth = as.numeric(ari_gt),
  Average_Silhouette_Width = as.numeric(avg_width),
  PeakRAM_MiB_Fusion = as.numeric(peak_mib),
  PeakRAM_Elapsed_Seconds_Fusion = as.numeric(elapsed_s)
)
data.table::fwrite(perf_row, file.path(out_dir, "perf_row.tsv"), sep = "\t", quote = FALSE, na = "NA")

writeLines(capture.output(sessionInfo()), file.path(meta_dir, "sessionInfo.txt"))

message(sprintf("[%s] %s sample %d%% rep %d done | peakMiB=%.1f elapsed=%.2f",
                Sys.time(), algorithm, frac, repi, peak_mib, elapsed_s))

