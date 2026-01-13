#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)

# MONET correlation-prep (sample perturbations)
# Uses Resources/Sample_perturbations/sample_subsets.tsv.gz for *shared* sample subsets across methods.

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

task_index <- suppressWarnings(as.integer(a[["task_index"]]))
if (is.na(task_index)) task_index <- suppressWarnings(as.integer(Sys.getenv("MONET_TASK_INDEX", "")))
if (is.na(task_index)) task_index <- suppressWarnings(as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "")))

corr_dir <- a[["corr_dir"]] %||% ""
prep_meta <- a[["prep_meta"]] %||% ""

if (is.na(task_index) || task_index < 1L) stop("Missing/invalid --task_index (or MONET_TASK_INDEX/SLURM_ARRAY_TASK_ID).")
if (!nzchar(corr_dir)) stop("Missing --corr_dir")
if (!nzchar(prep_meta)) stop("Missing --prep_meta")

RNGversion("4.2.2")
set.seed(123)

algorithm <- "MONET"

input_path <- Sys.getenv("MM_INPUT_RDS", "Resources/mm_input.rds")
subsets_path <- Sys.getenv("SAMPLE_SUBSETS", "Resources/Sample_perturbations/sample_subsets.tsv.gz")

if (!file.exists(input_path)) stop("Missing input RDS: ", input_path)
if (!file.exists(subsets_path)) stop("Missing sample subsets: ", subsets_path)

dir.create(corr_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(prep_meta), recursive = TRUE, showWarnings = FALSE)

required_mods <- c("RNAseq", "CNV", "Methylation", "miRNA", "SNPs")

input_full <- readRDS(input_path)
if (!all(required_mods %in% names(input_full))) {
  stop("mm_input.rds missing required modalities: ",
       paste(setdiff(required_mods, names(input_full)), collapse = ", "))
}

# Determine valid ids from a reference modality
ref_mod_full <- if ("RNAseq" %in% required_mods) "RNAseq" else required_mods[1]
valid_ids <- colnames(input_full[[ref_mod_full]])

sample_subsets <- data.table::fread(subsets_path)
need_cols <- c("fraction","replicate","sample_id")
if (!all(need_cols %in% colnames(sample_subsets))) stop("sample_subsets.tsv.gz must contain fraction, replicate, sample_id")

sample_subsets <- sample_subsets[sample_id %in% valid_ids]

pairs <- unique(sample_subsets[, .(fraction = as.integer(fraction), replicate = as.integer(replicate))])
setorder(pairs, fraction, replicate)

if (task_index > nrow(pairs)) stop("task_index out of range: ", task_index, " > ", nrow(pairs))

frac <- pairs$fraction[task_index]
repi <- pairs$replicate[task_index]

keep_ids <- unique(sample_subsets[fraction == frac & replicate == repi, sample_id])
if (length(keep_ids) < 3L) stop("Too few samples in subset: fraction=", frac, " replicate=", repi)

# Subset columns across modalities
input_sub <- input_full
for (m in required_mods) {
  mat <- input_sub[[m]]
  common <- intersect(colnames(mat), keep_ids)
  input_sub[[m]] <- mat[, common, drop = FALSE]
}

# Align common samples across modalities and keep a consistent order
ids_list <- lapply(required_mods, function(m) colnames(input_sub[[m]]))
common <- Reduce(intersect, ids_list)
if (length(common) < 3L) stop("No (or too few) common sample IDs across modalities after subsetting.")
common <- common[common %in% colnames(input_sub[["RNAseq"]])]

for (m in required_mods) input_sub[[m]] <- input_sub[[m]][, common, drop = FALSE]

sample_ids <- common
n_samples <- length(sample_ids)

p_by_mod <- sapply(required_mods, function(m) nrow(input_sub[[m]]))

message(sprintf("[%s] %s sample subset %d%% rep %d: computing correlations (n=%d) [task %d/%d]",
                Sys.time(), algorithm, frac, repi, n_samples, task_index, nrow(pairs)))

for (m in required_mods) {
  mat <- as.matrix(input_sub[[m]])
  cm <- HiClimR::fastCor(mat)
  dimnames(cm) <- list(sample_ids, sample_ids)
  diag(cm) <- 0
  write.csv(cm, file.path(corr_dir, paste0(m, ".csv")), quote = FALSE)
}

meta_dt <- data.table(
  Algorithm = algorithm,
  Mode = "sample",
  Sample_Percent = as.integer(frac),
  Replicate = as.integer(repi),
  n_samples = as.integer(n_samples),
  p_RNAseq = as.integer(p_by_mod[["RNAseq"]]),
  p_CNV = as.integer(p_by_mod[["CNV"]]),
  p_Methylation = as.integer(p_by_mod[["Methylation"]]),
  p_miRNA = as.integer(p_by_mod[["miRNA"]]),
  p_SNPs = as.integer(p_by_mod[["SNPs"]]),
  p_total = as.integer(sum(p_by_mod))
)

data.table::fwrite(meta_dt, prep_meta, sep = "\t", quote = FALSE)
