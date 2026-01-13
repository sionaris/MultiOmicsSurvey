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
  out <- list(); i <- 1
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

centile <- suppressWarnings(as.numeric(a[["centile"]]))
if (is.na(centile) || centile <= 0 || centile > 1) stop("Invalid --centile")

out_dir <- Sys.getenv("BENCH_OUT_DIR", "")
meta_dir <- Sys.getenv("BENCH_META_DIR", "")
run_dir <- Sys.getenv("BENCH_RUN_DIR", "")
if (!nzchar(out_dir) || !nzchar(meta_dir) || !nzchar(run_dir)) stop("Missing BENCH_* env vars from wrapper.")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(meta_dir, recursive = TRUE, showWarnings = FALSE)

algorithm <- "RWR-NF"
pct <- as.integer(round(centile * 100))

# Fixed optimal params (RWR-NF)
aff_sigma <- 0.5
aff_num_neighbors <- 10L
n_clusters <- 2L
iteration_max <- 1000L

gamma_fixed <- 0.7
neighbor_num_fixed <- 10L
alpha_fixed <- 0.9
beta_fixed <- 0.9

# Source RWR code (same source file)
rwr_src <- Sys.getenv("RWR_SOURCE", "Resources/RWR-F_source.R")
if (!file.exists(rwr_src)) stop("Missing RWR source file: ", rwr_src)
source(rwr_src)

input_path <- Sys.getenv("MM_INPUT_RDS", "Resources/mm_input.rds")
rank_path  <- Sys.getenv("FEATURE_RANKINGS", "Resources/Feature_perturbations/feature_rankings.tsv.gz")
if (!file.exists(input_path)) stop("Missing input RDS: ", input_path)
if (!file.exists(rank_path)) stop("Missing feature rankings: ", rank_path)

gt_override <- Sys.getenv("GROUND_TRUTH_XLSX", "")
gt_candidates <- c(
  gt_override,
  file.path("Resources/algorithms", paste0(algorithm, "_ground_truth_labels.xlsx")),
  file.path("Resources/algorithms", "RWRNF_ground_truth_labels.xlsx"),
  file.path("Resources/algorithms", "RWR-NF_ground_truth_labels.xlsx")
)
gt_candidates <- gt_candidates[nzchar(gt_candidates)]
gt_path <- gt_candidates[file.exists(gt_candidates)][1]
if (is.na(gt_path) || !nzchar(gt_path)) {
  stop("Missing ground truth xlsx. Tried:\n- ", paste(gt_candidates, collapse = "\n- "))
}

input_full <- readRDS(input_path)

required_mods <- intersect(c("RNAseq","CNV","Methylation","miRNA","SNPs"), names(input_full))
if (length(required_mods) < 2) stop("Too few modalities in mm_input.rds for RWR-NF.")

feat_rank <- data.table::fread(rank_path)
need_cols <- c("modality","rank","feature_id")
if (!all(need_cols %in% names(feat_rank))) stop("feature_rankings must have columns: modality, rank, feature_id")
feat_rank <- feat_rank[modality %in% required_mods][order(modality, rank)]
if (nrow(feat_rank) == 0) stop("No rankings for required modalities.")

get_top_features <- function(dt, cent) {
  out <- vector("list", length(required_mods)); names(out) <- required_mods
  for (m in required_mods) {
    d <- dt[modality == m][order(rank)]
    if (nrow(d) == 0) stop("No ranked features for modality: ", m)
    top_n <- max(1L, ceiling(cent * nrow(d)))
    out[[m]] <- d$feature_id[seq_len(top_n)]
  }
  out
}

subset_input_features <- function(input_obj, top_features) {
  out <- input_obj
  for (m in names(top_features)) {
    mat <- out[[m]]
    keep <- intersect(top_features[[m]], rownames(mat))
    if (length(keep) == 0L) stop("0 features kept for modality: ", m)
    out[[m]] <- mat[keep, , drop = FALSE]
  }
  out
}

top_feats <- get_top_features(feat_rank, centile)
input_sub <- subset_input_features(input_full, top_feats)

ids_list <- lapply(required_mods, function(m) colnames(input_sub[[m]]))
common <- Reduce(intersect, ids_list)
if (length(common) < 3L) stop("Too few common samples after feature subsetting.")
if ("RNAseq" %in% required_mods) common <- common[common %in% colnames(input_sub[["RNAseq"]])]
for (m in required_mods) input_sub[[m]] <- input_sub[[m]][, common, drop = FALSE]

sample_ids <- gsub("\\.", "-", common)
n_samples <- length(sample_ids)
p_by_mod <- sapply(required_mods, function(m) nrow(input_sub[[m]]))
p_total <- as.integer(sum(p_by_mod))

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

affinityL <- lapply(distL, function(d) {
  a <- SNFtool::affinityMatrix(as.matrix(d), K = aff_num_neighbors, sigma = aff_sigma)
  dimnames(a) <- dimnames(d)
  a
})

gc()

fused_rwrnf <- NULL
pr <- peakRAM::peakRAM({
  fused_rwrnf <<- RWR_fusion_neighbor(sim_list = affinityL, iteration_max = iteration_max, gama = gamma_fixed,
                            neighbor_num  = neighbor_num_fixed, alpha = alpha_fixed, beta = beta_fixed)
})

peak_mib <- as.numeric(pr$Peak_RAM_Used_MiB[1])
elapsed_s <- as.numeric(pr$Elapsed_Time[1])

specc_with_embedding <- function(Kmat, centers, iterations = 200) {
  X <- as.matrix(Kmat)
  d <- 1 / sqrt(rowSums(X))
  L <- d * X %*% diag(d)
  xi <- eigen(L)$vectors[, 1:centers, drop = FALSE]
  yi <- xi / sqrt(rowSums(xi^2))
  km <- kmeans(yi, centers = centers, iter.max = iterations)
  list(cluster = km$cluster, embedding = yi)
}

km <- fused_rwrnf
dimnames(km) <- list(sample_ids, sample_ids)
res <- specc_with_embedding(km, centers = n_clusters)

sil <- silhouette(as.integer(res$cluster), dist = Rfast::Dist(res$embedding, method = "euclidean"))
avg_width <- summary(sil)$avg.width

clusters <- data.table(Sample.ID = sample_ids, Cluster_pred = as.integer(res$cluster))

gt <- openxlsx::read.xlsx(gt_path)
if (!all(c("Sample.ID","Cluster") %in% colnames(gt))) stop("Ground truth must have Sample.ID, Cluster")
gt$Sample.ID <- gsub("\\.", "-", gt$Sample.ID)
setDT(gt); setnames(gt, "Cluster", "Cluster_gt")

merged <- merge(gt, clusters, by = "Sample.ID")
ari_gt <- if (nrow(merged) > 0) mclust::adjustedRandIndex(merged$Cluster_gt, merged$Cluster_pred) else NA_real_

data.table::fwrite(clusters, file.path(out_dir, "clusters.tsv.gz"),
                   sep = "\t", quote = FALSE, compress = "gzip")

prep_meta <- data.table(
  Algorithm = algorithm,
  Mode = "feature",
  Feature_Centile = centile,
  Feature_Percent = pct,
  n_samples = as.integer(n_samples),
  p_total = as.integer(p_total)
)
for (m in c("RNAseq","CNV","Methylation","miRNA","SNPs"))
  prep_meta[[paste0("p_", m)]] <- as.integer(p_by_mod[m] %||% NA_integer_)
data.table::fwrite(prep_meta, file.path(meta_dir, "prep_meta.tsv"), sep = "\t", quote = FALSE)

perf_row <- data.table(
  Algorithm = algorithm,
  Feature_Centile = centile,
  Feature_Percent = pct,
  n_samples = as.integer(n_samples),
  p_total = as.integer(p_total),
  Affinity_K = as.integer(aff_num_neighbors),
  Affinity_sigma = as.numeric(aff_sigma),
  RWR_iteration_max = as.integer(iteration_max),
  RWR_gamma = as.numeric(gamma_fixed),
  RWRNF_neighbor_num = as.integer(neighbor_num_fixed),
  RWRNF_alpha = as.numeric(alpha_fixed),
  RWRNF_beta = as.numeric(beta_fixed),
  n_clusters = as.integer(n_clusters),
  ARI_to_Ground_Truth = as.numeric(ari_gt),
  Average_Silhouette_Width = as.numeric(avg_width),
  PeakRAM_MiB_FusionNeighbor = as.numeric(peak_mib),
  PeakRAM_Elapsed_Seconds_FusionNeighbor = as.numeric(elapsed_s)
)
data.table::fwrite(perf_row, file.path(out_dir, "perf_row.tsv"), sep = "\t", quote = FALSE, na = "NA")

writeLines(capture.output(sessionInfo()), file.path(meta_dir, "sessionInfo.txt"))

message(sprintf("[%s] %s feature %d%% done | peakMiB=%.1f elapsed=%.2f",
                Sys.time(), algorithm, pct, peak_mib, elapsed_s))

