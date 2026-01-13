#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)

req <- c("data.table")
missing <- req[!vapply(req, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "))

suppressPackageStartupMessages(library(data.table))

`%||%` <- function(x, y) if (!is.null(x) && nzchar(x)) x else y

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
out_dir <- a[["out_dir"]] %||% ""
prep_meta <- a[["prep_meta"]] %||% ""
sample_ids_file <- a[["sample_ids"]] %||% ""

if (is.na(centile) || !is.finite(centile) || centile <= 0 || centile > 1) stop("Invalid --centile")
if (!nzchar(out_dir)) stop("Missing --out_dir")
if (!nzchar(prep_meta)) stop("Missing --prep_meta")
if (!nzchar(sample_ids_file)) stop("Missing --sample_ids")

RNGversion("4.2.2"); set.seed(123)

input_path <- Sys.getenv("MM_INPUT_RDS", "Resources/mm_input.rds")
rank_path  <- Sys.getenv("FEATURE_RANKINGS", "Resources/Feature_perturbations/feature_rankings.tsv.gz")

if (!file.exists(input_path)) stop("Missing input RDS: ", input_path)
if (!file.exists(rank_path)) stop("Missing feature rankings: ", rank_path)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(prep_meta), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(sample_ids_file), recursive = TRUE, showWarnings = FALSE)

mods <- c("RNAseq", "CNV", "SNPs", "miRNA", "Methylation")
input_full <- readRDS(input_path)
if (!all(mods %in% names(input_full))) {
  stop("mm_input.rds missing required modalities: ",
       paste(setdiff(mods, names(input_full)), collapse = ", "))
}

feat_rank <- data.table::fread(rank_path)
need_cols <- c("modality", "rank", "feature_id")
if (!all(need_cols %in% names(feat_rank))) stop("feature_rankings missing required columns")
feat_rank <- feat_rank[modality %in% mods]
if (nrow(feat_rank) == 0) stop("No ranked features for required modalities")

get_top_features <- function(dt, cent) {
  out <- vector("list", length(mods)); names(out) <- mods
  for (m in mods) {
    d <- dt[modality == m][order(rank)]
    if (nrow(d) == 0) stop("No ranked features for modality: ", m)
    top_n <- max(1L, ceiling(cent * nrow(d)))
    out[[m]] <- d$feature_id[seq_len(top_n)]
  }
  out
}

top_feats <- get_top_features(feat_rank, centile)

# subset rows by top features + align common samples
input_sub <- input_full
for (m in mods) {
  mat <- input_sub[[m]]
  keep <- intersect(top_feats[[m]], rownames(mat))
  if (!length(keep)) stop("0 kept features for modality: ", m)
  input_sub[[m]] <- mat[keep, , drop = FALSE]
}

ids_list <- lapply(mods, function(m) colnames(input_sub[[m]]))
common <- Reduce(intersect, ids_list)
if (!length(common)) stop("No common samples across modalities after feature subsetting.")
common <- common[common %in% colnames(input_sub[["RNAseq"]])]

for (m in mods) input_sub[[m]] <- input_sub[[m]][, common, drop = FALSE]
sample_ids <- common

# clean SNPs to {0,1}
if ("SNPs" %in% mods) {
  snp <- as.matrix(input_sub[["SNPs"]])
  snp <- suppressWarnings(apply(snp, 2, as.numeric))
  snp[is.na(snp)] <- 0
  snp[snp != 0] <- 1
  input_sub[["SNPs"]] <- snp
  rownames(input_sub[["SNPs"]]) <- rownames(snp)
  colnames(input_sub[["SNPs"]]) <- sample_ids
}

# write CSVs (features in rows, samples in columns)
for (m in mods) {
  mat <- as.matrix(input_sub[[m]])
  mat <- suppressWarnings(apply(mat, 2, as.numeric))
  mat[is.na(mat)] <- 0
  mat <- as.matrix(mat)
  rownames(mat) <- rownames(input_sub[[m]])
  colnames(mat) <- sample_ids
  dt <- data.table::as.data.table(mat)
  dt <- data.table::data.table(Feature = rownames(mat), dt)
  data.table::fwrite(dt, file.path(out_dir, paste0(m, ".csv")), sep = ",", quote = FALSE)
}

pct <- as.integer(round(centile * 100))
p_by_mod <- sapply(mods, function(m) nrow(input_sub[[m]]))

meta_dt <- data.table(
  Algorithm = "MOFA",
  Mode = "feature",
  Feature_Centile = centile,
  Feature_Percent = pct,
  n_samples = as.integer(length(sample_ids)),
  p_RNAseq = as.integer(p_by_mod[["RNAseq"]]),
  p_CNV = as.integer(p_by_mod[["CNV"]]),
  p_SNPs = as.integer(p_by_mod[["SNPs"]]),
  p_miRNA = as.integer(p_by_mod[["miRNA"]]),
  p_Methylation = as.integer(p_by_mod[["Methylation"]]),
  p_total = as.integer(sum(p_by_mod))
)

data.table::fwrite(meta_dt, prep_meta, sep = "\t", quote = FALSE)
data.table::fwrite(data.table(Sample.ID = sample_ids), sample_ids_file, sep = "\t", quote = FALSE)

message(sprintf("[%s] MOFA feature subset %d%%: wrote inputs for n=%d samples",
                Sys.time(), pct, length(sample_ids)))

