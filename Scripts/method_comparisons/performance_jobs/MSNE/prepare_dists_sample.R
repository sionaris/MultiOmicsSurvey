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
frac <- as.numeric(a[["frac"]])
rep  <- as.integer(a[["rep"]])
out_dir <- a[["out_dir"]]

if (!is.finite(frac) || frac <= 0 || frac > 1) stop("Invalid --frac (expected (0,1]).")
if (!is.finite(rep) || rep < 1) stop("Invalid --rep (expected >=1).")
if (!nzchar(out_dir)) stop("Missing --out_dir")

RNGversion("4.2.2"); set.seed(123 + rep)

input_path <- Sys.getenv("MM_INPUT_RDS", "Resources/TCGA/mm_input.rds")
if (!file.exists(input_path)) stop("Missing input RDS: ", input_path)

required_mods <- c("RNAseq","CNV","Methylation","miRNA","SNPs")
input_full <- readRDS(input_path)

# align common samples first
ids_list <- lapply(required_mods, function(m) colnames(input_full[[m]]))
common <- Reduce(intersect, ids_list)
if (length(common) == 0L) stop("No common sample IDs across modalities.")

n_keep <- max(2L, as.integer(ceiling(frac * length(common))))
keep_ids <- sort(sample(common, n_keep, replace=FALSE))

input_sub <- input_full
for (m in required_mods) input_sub[[m]] <- input_sub[[m]][, keep_ids, drop=FALSE]

dir.create(file.path(out_dir, "MO_Dists"), recursive=TRUE, showWarnings=FALSE)

cont <- c("RNAseq","CNV","Methylation","miRNA")
for (m in cont) {
  X <- t(as.matrix(input_sub[[m]]))
  D <- SNFtool::dist2(X, X)
  dimnames(D) <- list(keep_ids, keep_ids)
  write.csv(D, file.path(out_dir, "MO_Dists", paste0(m, ".csv")), quote=FALSE)
}

X <- t(as.matrix(input_sub[["SNPs"]]))
D <- as.matrix(dist(X, method="binary"))
dimnames(D) <- list(keep_ids, keep_ids)
write.csv(D, file.path(out_dir, "MO_Dists", "SNPs.csv"), quote=FALSE)

meta <- data.table(
  Algorithm="MSNE",
  Mode="sample",
  Sample_Frac=frac,
  Sample_Percent=as.integer(round(frac*100)),
  Replicate=rep,
  n_samples=length(keep_ids)
)
data.table::fwrite(meta, file.path(out_dir, "prep_meta.tsv"), sep="\t", quote=FALSE)

