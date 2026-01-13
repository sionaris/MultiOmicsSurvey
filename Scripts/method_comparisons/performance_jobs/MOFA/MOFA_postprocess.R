#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)

`%||%` <- function(x, y) if (!is.null(x) && nzchar(x)) x else y

req <- c("data.table", "jsonlite", "openxlsx", "mclust", "cluster", "Rfast")
missing <- req[!vapply(req, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Missing R packages: ", paste(missing, collapse = ", "),
       "\nInstall once in your R library; do not install inside jobs.")
}
if (!requireNamespace("MOFA2", quietly = TRUE)) {
  stop("Missing R package: MOFA2\nInstall with BiocManager::install('MOFA2')")
}
if (!requireNamespace("M3C", quietly = TRUE)) {
  stop("Missing R package: M3C\nInstall from CRAN (or your standard source).")
}

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
  library(openxlsx)
  library(mclust)
  library(cluster)
  library(Rfast)
  library(MOFA2)
  library(M3C)
})

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
mode <- a[["mode"]] %||% ""
prep_meta <- a[["prep_meta"]] %||% ""
sample_ids_file <- a[["sample_ids"]] %||% ""
mofa_hdf5 <- a[["mofa_hdf5"]] %||% ""
metrics_json <- a[["metrics_json"]] %||% ""
timev_file <- a[["timev_file"]] %||% ""
out_clusters <- a[["out_clusters"]] %||% ""
out_perf_row <- a[["out_perf_row"]] %||% ""

if (!nzchar(mode) || !(mode %in% c("feature", "sample"))) stop("Invalid --mode")
for (p in c(prep_meta, sample_ids_file, mofa_hdf5, metrics_json, out_clusters, out_perf_row)) {
  if (!nzchar(p)) stop("Missing required path argument.")
}
if (!file.exists(prep_meta)) stop("Missing prep_meta: ", prep_meta)
if (!file.exists(sample_ids_file)) stop("Missing sample_ids: ", sample_ids_file)
if (!file.exists(mofa_hdf5)) stop("Missing mofa_hdf5: ", mofa_hdf5)
if (!file.exists(metrics_json)) stop("Missing metrics_json: ", metrics_json)
if (nzchar(timev_file) && !file.exists(timev_file)) stop("Missing timev_file: ", timev_file)

algorithm <- "MOFA"

# Parse /usr/bin/time -v max RSS (kbytes)
parse_max_rss_kb <- function(path) {
  if (!nzchar(path) || !file.exists(path)) return(NA_real_)
  ln <- readLines(path, warn = FALSE)
  hit <- ln[grepl("^\\s*Maximum resident set size \\(kbytes\\):", ln)]
  if (!length(hit)) return(NA_real_)
  suppressWarnings(as.numeric(trimws(sub(".*:\\s*", "", hit[1]))))
}
max_rss_kb <- parse_max_rss_kb(timev_file)
max_rss_mib <- if (is.finite(max_rss_kb)) max_rss_kb / 1024 else NA_real_

prep <- data.table::fread(prep_meta)
met <- jsonlite::fromJSON(metrics_json)

build_s <- met$build_seconds %||% NA_real_
run_s <- met$run_seconds %||% NA_real_
buildrun_s <- met$build_plus_run_seconds %||% NA_real_
peak_gpu_mib <- met$peak_gpu_mem_mib %||% NA_real_

# Ground truth (expected to exist like your other methods)
gt_path <- Sys.getenv("GROUND_TRUTH_XLSX",
                      file.path("Resources/algorithms", paste0(algorithm, "_ground_truth_labels.xlsx")))
if (!file.exists(gt_path)) stop("Missing ground truth xlsx: ", gt_path)

gt <- openxlsx::read.xlsx(gt_path)
if (!all(c("Sample.ID", "Cluster") %in% colnames(gt))) stop("Ground truth must have Sample.ID, Cluster")
gt$Sample.ID <- gsub("\\.", "-", gt$Sample.ID)
colnames(gt)[colnames(gt) == "Cluster"] <- "Cluster_gt"

# Load model + embeddings
final_MOFA <- MOFA2::load_model(mofa_hdf5)
emb <- final_MOFA@expectations$Z$group0
emb <- as.matrix(emb)
colnames(emb) <- paste0("MOFA", seq_len(ncol(emb)))

# Ensure rownames are real Sample.IDs: fall back to sample_ids file order
sample_ids <- data.table::fread(sample_ids_file)[["Sample.ID"]]
sample_ids <- gsub("\\.", "-", sample_ids)

if (is.null(rownames(emb)) || anyNA(rownames(emb)) || !all(rownames(emb) %in% sample_ids)) {
  if (nrow(emb) != length(sample_ids)) stop("Embedding rows do not match sample_ids length.")
  rownames(emb) <- sample_ids
}

# M3C clustering (fixed optk=3, as requested)
optk <- 3L
m3c_des <- data.frame(ID = rownames(emb), class = rep("all", nrow(emb)), stringsAsFactors = FALSE)
m3c_input <- as.data.frame(t(emb))  # features (MOFA factors) in rows

RNGversion("4.2.2"); set.seed(123)
m3c_res <- M3C::M3C(
  m3c_input, des = m3c_des,
  iters = 100, repsref = 250, repsreal = 250,
  seed = 123, fsize = 18, lthick = 2, dotsize = 1.25,
  clusteralg = "km", maxK = 10
)

clusters_df <- data.frame(
  Sample.ID = rownames(m3c_res[["realdataresults"]][[optk]][["ordered_annotation"]]),
  Cluster_pred = as.integer(m3c_res[["realdataresults"]][[optk]][["ordered_annotation"]][["consensuscluster"]]),
  stringsAsFactors = FALSE
)
clusters_df$Sample.ID <- gsub("\\.", "-", clusters_df$Sample.ID)

# ARI
merged <- merge(gt, clusters_df, by = "Sample.ID")
ari_gt <- if (nrow(merged) > 0) mclust::adjustedRandIndex(merged$Cluster_gt, merged$Cluster_pred) else NA_real_

# Silhouette on embeddings
sil <- NA_real_
if (length(unique(clusters_df$Cluster_pred)) > 1L) {
  emb_sub <- emb[clusters_df$Sample.ID, , drop = FALSE]
  sil_obj <- cluster::silhouette(as.integer(clusters_df$Cluster_pred), Rfast::Dist(emb_sub))
  sil <- summary(sil_obj)$avg.width
}

# Write clusters
dir.create(dirname(out_clusters), recursive = TRUE, showWarnings = FALSE)
data.table::fwrite(as.data.table(clusters_df), out_clusters, sep = "\t", quote = FALSE, compress = "gzip")

# Perf row
dir.create(dirname(out_perf_row), recursive = TRUE, showWarnings = FALSE)

n_samples <- nrow(emb)
n_clusters <- length(unique(clusters_df$Cluster_pred))

row <- NULL
if (mode == "feature") {
  row <- data.table(
    Algorithm = algorithm,
    Feature_Centile = as.numeric(prep$Feature_Centile[1]),
    Feature_Percent = as.integer(prep$Feature_Percent[1]),
    n_samples = as.integer(n_samples),
    p_total = as.integer(prep$p_total[1]),
    p_RNAseq = as.integer(prep$p_RNAseq[1]),
    p_CNV = as.integer(prep$p_CNV[1]),
    p_SNPs = as.integer(prep$p_SNPs[1]),
    p_miRNA = as.integer(prep$p_miRNA[1]),
    p_Methylation = as.integer(prep$p_Methylation[1]),
    n_clusters = as.integer(n_clusters),
    ARI_to_Ground_Truth = ari_gt,
    Average_Silhouette_Width = sil,
    MOFA_Build_Seconds = as.numeric(build_s),
    MOFA_Run_Seconds = as.numeric(run_s),
    MOFA_BuildPlusRun_Seconds = as.numeric(buildrun_s),
    MOFA_PeakRSS_MiB = as.numeric(max_rss_mib),
    MOFA_PeakGPU_MiB = as.numeric(peak_gpu_mib)
  )
} else {
  row <- data.table(
    Algorithm = algorithm,
    Sample_Percent = as.integer(prep$Sample_Percent[1]),
    Replicate = as.integer(prep$Replicate[1]),
    n_samples = as.integer(n_samples),
    p_total = as.integer(prep$p_total[1]),
    p_RNAseq = as.integer(prep$p_RNAseq[1]),
    p_CNV = as.integer(prep$p_CNV[1]),
    p_SNPs = as.integer(prep$p_SNPs[1]),
    p_miRNA = as.integer(prep$p_miRNA[1]),
    p_Methylation = as.integer(prep$p_Methylation[1]),
    n_clusters = as.integer(n_clusters),
    ARI_to_Ground_Truth = ari_gt,
    Average_Silhouette_Width = sil,
    MOFA_Build_Seconds = as.numeric(build_s),
    MOFA_Run_Seconds = as.numeric(run_s),
    MOFA_BuildPlusRun_Seconds = as.numeric(buildrun_s),
    MOFA_PeakRSS_MiB = as.numeric(max_rss_mib),
    MOFA_PeakGPU_MiB = as.numeric(peak_gpu_mib)
  )
}

data.table::fwrite(row, out_perf_row, sep = "\t", quote = FALSE, na = "NA")
message("Wrote: ", out_clusters)
message("Wrote: ", out_perf_row)

