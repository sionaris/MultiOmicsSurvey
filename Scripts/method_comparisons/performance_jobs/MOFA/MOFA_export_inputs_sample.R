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
task_index <- suppressWarnings(as.integer(a[["task_index"]]))
out_dir <- a[["out_dir"]] %||% ""
prep_meta <- a[["prep_meta"]] %||% ""
sample_ids_file <- a[["sample_ids"]] %||% ""

if (is.na(task_index) || task_index < 1L) stop("Invalid --task_index")
if (!nzchar(out_dir)) stop("Missing --out_dir")
if (!nzchar(prep_meta)) stop("Missing --prep_meta")
if (!nzchar(sample_ids_file)) stop("Missing --sample_ids")

RNGversion("4.2.2"); set.seed(123)

input_path <- Sys.getenv("MM_INPUT_RDS", "Resources/mm_input.rds")
subsets_path <- Sys.getenv("SAMPLE_SUBSETS", "Resources/Sample_perturbations/sample_subsets.tsv.gz")

if (!file.exists(input_path)) stop("Missing input RDS: ", input_path)
if (!file.exists(subsets_path)) stop("Missing sample subsets: ", subsets_path)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(prep_meta), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(sample_ids_file), recursive = TRUE, showWarnings = FALSE)

mods <- c("RNAseq", "CNV", "SNPs", "miRNA", "Methylation")
input_full <- readRDS(input_path)
if (!all(mods %in% names(input_full))) {
  stop("mm_input.rds missing required modalities: ",
       paste(setdiff(mods, names(input_full)), collapse = ", "))
}

ref_ids <- colnames(input_full[["RNAseq"]])

sample_subsets <- data.table::fread(subsets_path)
need_cols <- c("fraction", "replicate", "sample_id")
if (!all(need_cols %in% colnames(sample_subsets))) stop("sample_subsets missing required columns")

sample_subsets <- sample_subsets[sample_id %in% ref_ids]
pairs <- unique(sample_subsets[, .(fraction = as.integer(fraction), replicate = as.integer(replicate))])
setorder(pairs, fraction, replicate)

if (task_index > nrow(pairs)) stop("task_index out of range")

frac <- pairs$fraction[task_index]
repi <- pairs$replicate[task_index]
keep_ids <- unique(sample_subsets[fraction == frac & replicate == repi, sample_id])
if (length(keep_ids) < 3L) stop("Too few samples in subset")

# subset columns across modalities
input_sub <- input_full
for (m in mods) {
  mat <- input_sub[[m]]
  common <- intersect(colnames(mat), keep_ids)
  input_sub[[m]] <- mat[, common, drop = FALSE]
}

# align common samples across all modalities
ids_list <- lapply(mods, function(m) colnames(input_sub[[m]]))
common <- Reduce(intersect, ids_list)
if (length(common) < 3L) stop("Too few common samples across modalities after subsetting.")
common <- common[common %in% colnames(input_sub[["RNAseq"]])]
for (m in mods) input_sub[[m]] <- input_sub[[m]][, common, drop = FALSE]
sample_ids <- common

# clean SNPs to {0,1}
snp <- as.matrix(input_sub[["SNPs"]])
snp <- suppressWarnings(apply(snp, 2, as.numeric))
snp[is.na(snp)] <- 0
snp[snp != 0] <- 1
input_sub[["SNPs"]] <- snp
rownames(input_sub[["SNPs"]]) <- rownames(snp)
colnames(input_sub[["SNPs"]]) <- sample_ids

# write CSVs
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

p_by_mod <- sapply(mods, function(m) nrow(input_sub[[m]]))

meta_dt <- data.table(
  Algorithm = "MOFA",
  Mode = "sample",
  Sample_Percent = as.integer(frac),
  Replicate = as.integer(repi),
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

message(sprintf("[%s] MOFA sample subset %d%% rep %d: wrote inputs for n=%d samples [task %d/%d]",
                Sys.time(), frac, repi, length(sample_ids), task_index, nrow(pairs)))

