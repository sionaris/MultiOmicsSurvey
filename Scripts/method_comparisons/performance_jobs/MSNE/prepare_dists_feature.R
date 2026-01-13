#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)

req <- c("data.table", "SNFtool")
missing <- req[!vapply(req, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing R packages: ", paste(missing, collapse=", "),
                          "\nInstall once in your R library; do not install inside jobs.")

suppressPackageStartupMessages({
  library(data.table)
  library(SNFtool)
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
centile <- as.numeric(a[["centile"]])
out_dir <- a[["out_dir"]]

if (!is.finite(centile) || centile <= 0 || centile > 1) stop("Invalid --centile (expected (0,1]).")
if (!nzchar(out_dir)) stop("Missing --out_dir")

RNGversion("4.2.2"); set.seed(123)

input_path <- Sys.getenv("MM_INPUT_RDS", "Resources/TCGA/mm_input.rds")
rank_path  <- Sys.getenv("FEATURE_RANKINGS", "Resources/Performance/Feature_perturbations/feature_rankings.tsv.gz")

if (!file.exists(input_path)) stop("Missing input RDS: ", input_path)
if (!file.exists(rank_path))  stop("Missing feature rankings: ", rank_path)

required_mods <- c("RNAseq", "CNV", "Methylation", "miRNA", "SNPs")
input_full <- readRDS(input_path)
if (!all(required_mods %in% names(input_full))) {
  stop("mm_input.rds missing modalities: ", paste(setdiff(required_mods, names(input_full)), collapse=", "))
}

feat_rank <- data.table::fread(rank_path)
need_cols <- c("modality","rank","feature_id")
if (!all(need_cols %in% names(feat_rank))) stop("feature_rankings must have: modality, rank, feature_id")

feat_rank <- feat_rank[modality %in% required_mods]

get_top <- function(dt, cent) {
  out <- vector("list", length(required_mods)); names(out) <- required_mods
  for (m in required_mods) {
    d <- dt[modality == m][order(rank)]
    top_n <- max(1L, ceiling(cent * nrow(d)))
    out[[m]] <- d$feature_id[seq_len(top_n)]
  }
  out
}

top_feats <- get_top(feat_rank, centile)

# subset features (rows)
input_sub <- input_full
for (m in required_mods) {
  mat <- input_sub[[m]]
  keep <- intersect(top_feats[[m]], rownames(mat))
  input_sub[[m]] <- mat[keep, , drop=FALSE]
  if (nrow(input_sub[[m]]) == 0L) stop("0 features after subsetting for modality: ", m)
}

# align common samples (columns)
ids_list <- lapply(required_mods, function(m) colnames(input_sub[[m]]))
common <- Reduce(intersect, ids_list)
if (length(common) == 0L) stop("No common sample IDs across modalities.")
for (m in required_mods) input_sub[[m]] <- input_sub[[m]][, common, drop=FALSE]

# compute distances (samples x samples)
dir.create(file.path(out_dir, "MO_Dists"), recursive=TRUE, showWarnings=FALSE)

cont <- c("RNAseq","CNV","Methylation","miRNA")
for (m in cont) {
  X <- t(as.matrix(input_sub[[m]]))               # samples x features
  D <- SNFtool::dist2(X, X)                       # squared Euclidean
  dimnames(D) <- list(common, common)
  write.csv(D, file.path(out_dir, "MO_Dists", paste0(m, ".csv")), quote=FALSE)
}

# SNPs binary distance
X <- t(as.matrix(input_sub[["SNPs"]]))
D <- as.matrix(dist(X, method="binary"))
dimnames(D) <- list(common, common)
write.csv(D, file.path(out_dir, "MO_Dists", "SNPs.csv"), quote=FALSE)

meta <- data.table(
  Algorithm="MSNE",
  Mode="feature",
  Feature_Centile=centile,
  Feature_Percent=as.integer(round(centile*100)),
  n_samples=length(common)
)
data.table::fwrite(meta, file.path(out_dir, "prep_meta.tsv"), sep="\t", quote=FALSE)

