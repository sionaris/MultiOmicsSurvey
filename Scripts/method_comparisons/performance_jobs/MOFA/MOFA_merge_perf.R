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
run_dir <- a[["run_dir"]] %||% ""
mode <- a[["mode"]] %||% ""
out_file <- a[["out_file"]] %||% ""

if (!nzchar(run_dir) || !dir.exists(run_dir)) stop("Invalid --run_dir")
if (!nzchar(out_file)) stop("Missing --out_file")
if (!nzchar(mode) || !(mode %in% c("feature", "sample"))) stop("Invalid --mode")

task_dirs <- list.dirs(run_dir, full.names = TRUE, recursive = FALSE)
task_dirs <- task_dirs[grepl("/task_[0-9]+$", task_dirs)]
if (!length(task_dirs)) stop("No task_* dirs under: ", run_dir)

perf_files <- file.path(task_dirs, "out", "mofa_perf_row.tsv")
perf_files <- perf_files[file.exists(perf_files)]
if (!length(perf_files)) stop("No mofa_perf_row.tsv found to merge under: ", run_dir)

dt_list <- lapply(perf_files, fread)
merged <- rbindlist(dt_list, use.names = TRUE, fill = TRUE)

if (mode == "feature" && all(c("Feature_Percent") %in% names(merged))) setorder(merged, Feature_Percent)
if (mode == "sample" && all(c("Sample_Percent", "Replicate") %in% names(merged))) setorder(merged, Sample_Percent, Replicate)

dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
fwrite(merged, out_file, sep = "\t", quote = FALSE, na = "NA")
message("Wrote: ", out_file)

