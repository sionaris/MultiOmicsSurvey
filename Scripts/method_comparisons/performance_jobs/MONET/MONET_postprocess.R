#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)

`%||%` <- function(x, y) if (!is.null(x) && nzchar(x)) x else y

req <- c("data.table", "openxlsx", "cluster", "mclust", "jsonlite")
missing <- req[!vapply(req, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Missing R packages: ", paste(missing, collapse = ", "),
       "\nInstall once in your R library; do not install inside jobs.")
}

suppressPackageStartupMessages({
  library(data.table)
})

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
    } else {
      i <- i + 1
    }
  }
  out
}

a <- parse_args()

mode <- a[["mode"]] %||% ""
prep_meta <- a[["prep_meta"]] %||% ""
metrics_json <- a[["metrics_json"]] %||% ""
timev_file <- a[["timev_file"]] %||% ""
clustering_tsv <- a[["clustering_tsv"]] %||% ""
avg_adj_tsv_gz <- a[["avg_adj_tsv_gz"]] %||% ""
out_clusters <- a[["out_clusters"]] %||% ""
out_perf_row <- a[["out_perf_row"]] %||% ""

if (!nzchar(mode) || !(mode %in% c("feature", "sample"))) stop("Invalid --mode (feature|sample)")
for (p in c(prep_meta, metrics_json, clustering_tsv, avg_adj_tsv_gz, out_clusters, out_perf_row)) {
  if (!nzchar(p)) stop("Missing required path argument.")
}
if (!file.exists(prep_meta)) stop("Missing prep_meta: ", prep_meta)
if (!file.exists(metrics_json)) stop("Missing metrics_json: ", metrics_json)
if (!file.exists(clustering_tsv)) stop("Missing clustering_tsv: ", clustering_tsv)
if (!file.exists(avg_adj_tsv_gz)) stop("Missing avg_adj_tsv_gz: ", avg_adj_tsv_gz)
if (nzchar(timev_file) && !file.exists(timev_file)) stop("Missing timev_file: ", timev_file)

algorithm <- "MONET"

gt_path <- Sys.getenv("GROUND_TRUTH_XLSX", file.path("Resources/algorithms", paste0(algorithm, "_ground_truth_labels.xlsx")))
if (!file.exists(gt_path)) stop("Missing ground truth xlsx: ", gt_path)

prep <- data.table::fread(prep_meta)
met <- jsonlite::fromJSON(metrics_json)

# Parse time -v max RSS
parse_max_rss_kb <- function(path) {
  if (!nzchar(path) || !file.exists(path)) return(NA_real_)
  ln <- readLines(path, warn = FALSE)
  # Avoid \s (not a valid escape in base R strings); use POSIX [[:space:]]
  hit <- ln[grepl("^[[:space:]]*Maximum resident set size \\(kbytes\\):", ln)]
  if (!length(hit)) return(NA_real_)
  val <- suppressWarnings(as.numeric(trimws(sub(".*:[[:space:]]*", "", hit[1]))))
  val
}
max_rss_kb <- parse_max_rss_kb(timev_file)
max_rss_mib <- if (is.finite(max_rss_kb)) max_rss_kb / 1024 else NA_real_

# Load clustering
clust <- data.table::fread(clustering_tsv)
if (!all(c("Sample.ID", "Cluster") %in% names(clust))) {
  stop("clustering_tsv must contain columns: Sample.ID, Cluster")
}
clust[, Sample.ID := gsub("\\.", "-", Sample.ID)]
setnames(clust, "Cluster", "Cluster_pred")

# Ground truth
gt <- openxlsx::read.xlsx(gt_path)
if (!all(c("Sample.ID", "Cluster") %in% colnames(gt))) {
  stop("Ground truth xlsx must contain columns: Sample.ID, Cluster")
}
gt$Sample.ID <- gsub("\\.", "-", gt$Sample.ID)
colnames(gt)[colnames(gt) == "Cluster"] <- "Cluster_gt"

# ARI
merged <- merge(gt, clust, by = "Sample.ID")
ari_gt <- if (nrow(merged) > 0) mclust::adjustedRandIndex(merged$Cluster_gt, merged$Cluster_pred) else NA_real_

# Load avg adjacency (gzipped wide TSV)
adj_dt <- data.table::fread(avg_adj_tsv_gz)
id_col <- names(adj_dt)[1]
ids_raw <- adj_dt[[id_col]]
ids <- gsub("\\.", "-", ids_raw)

adj_mat <- as.matrix(adj_dt[, -1, with = FALSE])
rownames(adj_mat) <- ids
colnames(adj_mat) <- gsub("\\.", "-", colnames(adj_mat))

# Align ordering to clustering
common <- intersect(rownames(adj_mat), clust$Sample.ID)
if (length(common) < 3L) stop("Too few common samples between adjacency and clustering.")
adj_mat <- adj_mat[common, common, drop = FALSE]
clust_sub <- clust[Sample.ID %in% common]
setkey(clust_sub, Sample.ID)
clust_sub <- clust_sub[common]

labs <- as.integer(clust_sub$Cluster_pred)

# Silhouette on (1 - avg_adjacency)/2
dmatrix <- (1 - adj_mat) / 2
dmatrix[dmatrix < 0] <- 0
diag(dmatrix) <- 0

avg_sil <- NA_real_
if (length(unique(labs)) > 1L) {
  sil <- cluster::silhouette(labs, dmatrix = dmatrix)
  avg_sil <- summary(sil)$avg.width
}

# Output clusters (TSV.GZ)
dir.create(dirname(out_clusters), recursive = TRUE, showWarnings = FALSE)
data.table::fwrite(
  clust[, .(Sample.ID, Cluster_pred)],
  out_clusters, sep = "\t", quote = FALSE, compress = "gzip"
)

# Performance row
dir.create(dirname(out_perf_row), recursive = TRUE, showWarnings = FALSE)

main_loop_elapsed <- met$main_loop_elapsed_seconds %||% NA_real_
n_modules <- met$n_modules %||% NA_integer_
n_samples <- met$n_samples %||% nrow(prep)

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
    p_Methylation = as.integer(prep$p_Methylation[1]),
    p_miRNA = as.integer(prep$p_miRNA[1]),
    p_SNPs = as.integer(prep$p_SNPs[1]),
    n_clusters = as.integer(length(unique(labs))),
    n_modules = as.integer(n_modules),
    ARI_to_Ground_Truth = ari_gt,
    Average_Silhouette_Width = avg_sil,
    MONET_MainLoop_Elapsed_Seconds = as.numeric(main_loop_elapsed),
    MONET_PeakRSS_MiB = as.numeric(max_rss_mib)
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
    p_Methylation = as.integer(prep$p_Methylation[1]),
    p_miRNA = as.integer(prep$p_miRNA[1]),
    p_SNPs = as.integer(prep$p_SNPs[1]),
    n_clusters = as.integer(length(unique(labs))),
    n_modules = as.integer(n_modules),
    ARI_to_Ground_Truth = ari_gt,
    Average_Silhouette_Width = avg_sil,
    MONET_MainLoop_Elapsed_Seconds = as.numeric(main_loop_elapsed),
    MONET_PeakRSS_MiB = as.numeric(max_rss_mib)
  )
}

data.table::fwrite(row, out_perf_row, sep = "\t", quote = FALSE, na = "NA")

