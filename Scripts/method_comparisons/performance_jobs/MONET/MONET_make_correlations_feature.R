#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)

# MONET correlation-prep (feature perturbations)
# - mm_input.rds has features in rows, samples in columns.
# - HiClimR::fastCor correlates columns, so fastCor(features x samples) -> sample x sample correlation.

req <- c("data.table", "HiClimR")
missing <- req[!vapply(req, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Missing R packages: ", paste(missing, collapse = ", "),
       "\nInstall once in your R library; do not install inside jobs.")
}

suppressPackageStartupMessages({
  library(data.table)
})

`%||%` <- function(x, y) if (!is.null(x) && nzchar(x)) x else y

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

centile <- suppressWarnings(as.numeric(a[["centile"]]))
corr_dir <- a[["corr_dir"]] %||% ""
prep_meta <- a[["prep_meta"]] %||% ""

if (is.na(centile) || !is.finite(centile) || centile <= 0 || centile > 1) {
  stop("Invalid or missing --centile (expected (0,1]). Got: ", a[["centile"]])
}
if (!nzchar(corr_dir)) stop("Missing --corr_dir")
if (!nzchar(prep_meta)) stop("Missing --prep_meta")

RNGversion("4.2.2")
set.seed(123)

algorithm <- "MONET"

input_path <- Sys.getenv("MM_INPUT_RDS", "Resources/mm_input.rds")
rank_path  <- Sys.getenv("FEATURE_RANKINGS", "Resources/Feature_perturbations/feature_rankings.tsv.gz")

if (!file.exists(input_path)) stop("Missing input RDS: ", input_path)
if (!file.exists(rank_path)) stop("Missing feature rankings: ", rank_path)

dir.create(corr_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(prep_meta), recursive = TRUE, showWarnings = FALSE)

required_mods <- c("RNAseq", "CNV", "Methylation", "miRNA", "SNPs")

input_full <- readRDS(input_path)
if (!all(required_mods %in% names(input_full))) {
  stop("mm_input.rds missing required modalities: ",
       paste(setdiff(required_mods, names(input_full)), collapse = ", "))
}

feat_rank <- data.table::fread(rank_path)
need_cols <- c("modality", "rank", "feature_id")
if (!all(need_cols %in% names(feat_rank))) {
  stop("feature_rankings.tsv.gz must contain columns: modality, rank, feature_id")
}

feat_rank <- feat_rank[modality %in% required_mods]
if (nrow(feat_rank) == 0) stop("feature_rankings has no rows for required modalities.")

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
    if (length(keep) == 0L) stop("After feature subsetting, modality has 0 kept rows: ", m)
    out[[m]] <- mat[keep, , drop = FALSE]
  }
  out
}

align_common_samples <- function(input_obj, mods) {
  ids_list <- lapply(mods, function(m) colnames(input_obj[[m]]))
  common <- Reduce(intersect, ids_list)
  if (length(common) == 0L) stop("No common sample IDs across modalities.")
  # Keep consistent order: RNAseq first if possible
  common <- common[common %in% colnames(input_obj[["RNAseq"]])]
  for (m in mods) input_obj[[m]] <- input_obj[[m]][, common, drop = FALSE]
  list(obj = input_obj, sample_ids = common)
}

top_feats <- get_top_features(feat_rank, centile)
input_sub <- subset_input_features(input_full, top_feats)

aligned <- align_common_samples(input_sub, required_mods)
input_sub <- aligned$obj
sample_ids <- aligned$sample_ids
n_samples <- length(sample_ids)

p_by_mod <- sapply(required_mods, function(m) nrow(input_sub[[m]]))
pct <- as.integer(round(centile * 100))

message(sprintf("[%s] %s feature subset %d%%: computing correlations (n=%d)",
                Sys.time(), algorithm, pct, n_samples))

for (m in required_mods) {
  mat <- as.matrix(input_sub[[m]])
  cm <- HiClimR::fastCor(mat)
  dimnames(cm) <- list(sample_ids, sample_ids)
  diag(cm) <- 0
  write.csv(cm, file.path(corr_dir, paste0(m, ".csv")), quote = FALSE)
}

meta_dt <- data.table(
  Algorithm = algorithm,
  Mode = "feature",
  Feature_Centile = centile,
  Feature_Percent = pct,
  n_samples = as.integer(n_samples),
  p_RNAseq = as.integer(p_by_mod[["RNAseq"]]),
  p_CNV = as.integer(p_by_mod[["CNV"]]),
  p_Methylation = as.integer(p_by_mod[["Methylation"]]),
  p_miRNA = as.integer(p_by_mod[["miRNA"]]),
  p_SNPs = as.integer(p_by_mod[["SNPs"]]),
  p_total = as.integer(sum(p_by_mod))
)

data.table::fwrite(meta_dt, prep_meta, sep = "\t", quote = FALSE)
