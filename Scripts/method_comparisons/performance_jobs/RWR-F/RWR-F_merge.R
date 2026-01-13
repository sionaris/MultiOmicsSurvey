#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)

req <- c("data.table")
missing <- req[!vapply(req, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing R packages: ", paste(missing, collapse=", "))

suppressPackageStartupMessages(library(data.table))

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
mode <- a[["mode"]]
if (!nzchar(mode) || !(mode %in% c("feature","sample"))) stop("Use --mode feature|sample")

run_dir <- Sys.getenv("BENCH_RUN_DIR", "")
if (!nzchar(run_dir)) stop("Missing BENCH_RUN_DIR in env")

out_dir <- file.path(run_dir, "out")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

rows <- list.files(run_dir, pattern = "perf_row\\.tsv$", recursive = TRUE, full.names = TRUE)
if (!length(rows)) stop("No perf_row.tsv found under: ", run_dir)

dt <- rbindlist(lapply(rows, fread), use.names = TRUE, fill = TRUE)

# Sort
if (mode == "feature" && all(c("Feature_Percent") %in% names(dt))) {
  setorder(dt, Algorithm, Feature_Percent)
} else if (mode == "sample" && all(c("Sample_Percent","Replicate") %in% names(dt))) {
  setorder(dt, Algorithm, Sample_Percent, Replicate)
} else {
  setorder(dt, Algorithm)
}

outfile <- file.path(out_dir, sprintf("RWR-F_%s_perturbations_performance.tsv", mode))
fwrite(dt, outfile, sep = "\t", quote = FALSE, na = "NA")

message("Wrote: ", outfile)

