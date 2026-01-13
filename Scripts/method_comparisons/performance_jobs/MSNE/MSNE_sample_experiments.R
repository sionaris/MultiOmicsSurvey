#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
step <- if (length(args) >= 1) args[[1]] else ""
if (!step %in% c("prepare", "post", "merge")) {
  stop("Usage: Rscript Jobs/MSNE/MSNE_sample_experiments.R <prepare|post|merge>")
}

required_pkgs <- c("data.table", "openxlsx", "SNFtool", "mclust", "cluster", "Rfast")
missing <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  stop("Missing R packages: ", paste(missing, collapse = ", "),
       "\nInstall once in your R library; do not install inside jobs.")
}

suppressPackageStartupMessages({
  library(data.table)
  library(openxlsx)
  library(SNFtool)
  library(mclust)
  library(cluster)
  library(Rfast)
})

algorithm <- "MSNE"

bench_out <- Sys.getenv("BENCH_OUT_DIR", "")
bench_meta <- Sys.getenv("BENCH_META_DIR", "")
bench_run <- Sys.getenv("BENCH_RUN_DIR", "")
if (!nzchar(bench_out) && step != "merge") stop("BENCH_OUT_DIR not set")
if (!nzchar(bench_meta) && step != "merge") stop("BENCH_META_DIR not set")
if (!nzchar(bench_run) && step == "merge") stop("BENCH_RUN_DIR not set")

input_path <- Sys.getenv("MSNE_INPUT_RDS", "Resources/mm_input.rds")
subsets_path <- Sys.getenv("MSNE_SAMPLE_SUBSETS", "Resources/Sample_perturbations/sample_subsets.tsv.gz")
gt_path <- Sys.getenv("MSNE_GROUND_TRUTH", file.path("Resources/algorithms", "MSNE_ground_truth_labels.xlsx"))

parse_timev <- function(path) {
  if (!file.exists(path)) return(list(elapsed_s = NA_real_, maxrss_kb = NA_real_))
  x <- readLines(path, warn = FALSE)

  get_val <- function(prefix) {
    ln <- x[grepl(paste0("^", prefix), x)]
    if (length(ln) == 0) return(NA_character_)
    sub(paste0("^", prefix, "\\s*:\\s*"), "", ln[[1]])
  }

  elapsed_txt <- get_val("Elapsed \\(wall clock\\) time \\(h:mm:ss or m:ss\\)")
  maxrss_txt <- get_val("Maximum resident set size \\(kbytes\\)")

  to_seconds <- function(s) {
    if (is.na(s)) return(NA_real_)
    s <- trimws(s)
    parts <- strsplit(s, ":", fixed = TRUE)[[1]]
    parts <- as.numeric(parts)
    if (any(is.na(parts))) return(NA_real_)
    if (length(parts) == 2) return(parts[1] * 60 + parts[2])
    if (length(parts) == 3) return(parts[1] * 3600 + parts[2] * 60 + parts[3])
    NA_real_
  }

  list(
    elapsed_s = to_seconds(elapsed_txt),
    maxrss_kb = suppressWarnings(as.numeric(maxrss_txt))
  )
}

subset_input_cols <- function(input_obj, keep_ids) {
  out <- input_obj
  for (nm in names(out)) {
    if (!is.null(out[[nm]])) {
      kk <- keep_ids[keep_ids %in% colnames(out[[nm]])]
      out[[nm]] <- out[[nm]][, kk, drop = FALSE]
      if (ncol(out[[nm]]) == 0L || nrow(out[[nm]]) == 0L) out[[nm]] <- NULL
    }
  }
  out
}

if (step == "prepare") {
  RNGversion("4.2.2"); set.seed(123)

  if (!file.exists(input_path)) stop("Missing input: ", input_path)
  if (!file.exists(subsets_path)) stop("Missing sample subsets: ", subsets_path)

  input_full <- readRDS(input_path)
  sample_subsets <- data.table::fread(subsets_path)
  need_cols <- c("fraction", "replicate", "sample_id")
  if (!all(need_cols %in% colnames(sample_subsets))) {
    stop("sample_subsets.tsv.gz must contain columns: fraction, replicate, sample_id")
  }

  non_null_mods <- names(input_full)[!vapply(input_full, is.null, logical(1))]
  ref_mod <- if ("RNAseq" %in% non_null_mods) "RNAseq" else non_null_mods[1]
  valid_ids <- colnames(input_full[[ref_mod]])
  sample_subsets <- sample_subsets[sample_id %in% valid_ids]
  if (nrow(sample_subsets) == 0) stop("No sample IDs in sample_subsets match input colnames")

  pairs <- unique(sample_subsets[, .(fraction, replicate)])
  setorder(pairs, fraction, replicate)

  idx <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
  if (idx < 1 || idx > nrow(pairs)) stop("Array index out of range")

  frac <- as.integer(pairs$fraction[idx])
  repi <- as.integer(pairs$replicate[idx])

  keep_ids <- unique(sample_subsets[fraction == frac & replicate == repi, sample_id])
  input_sub <- subset_input_cols(input_full, keep_ids)

  mods_req <- c("RNAseq", "CNV", "Methylation", "miRNA", "SNPs")
  mods_use <- intersect(mods_req, names(input_sub)[!vapply(input_sub, is.null, logical(1))])
  if (length(mods_use) < 5) stop("Missing required modalities after sample subsetting")

  # align to common samples across modalities
  ids_common <- Reduce(intersect, lapply(input_sub[mods_use], colnames))
  if (length(ids_common) < 3) stop("Too few common samples after alignment")
  for (m in mods_use) input_sub[[m]] <- input_sub[[m]][, ids_common, drop = FALSE]

  out_dists <- file.path(bench_out, "MO_Dists")
  dir.create(out_dists, recursive = TRUE, showWarnings = FALSE)

  continuous <- c("RNAseq", "CNV", "Methylation", "miRNA")
  for (m in continuous) {
    x <- t(as.matrix(input_sub[[m]]))
    d <- SNFtool::dist2(x, x)
    dimnames(d) <- list(ids_common, ids_common)
    data.table::fwrite(data.table::as.data.table(d, keep.rownames = "Sample.ID"),
                       file.path(out_dists, paste0(m, ".csv")), sep = ",", quote = FALSE)
  }

  x <- t(as.matrix(input_sub[["SNPs"]]))
  x <- apply(x, 2, as.numeric)
  d <- as.matrix(dist(x, method = "binary"))
  dimnames(d) <- list(ids_common, ids_common)
  data.table::fwrite(data.table::as.data.table(d, keep.rownames = "Sample.ID"),
                     file.path(out_dists, "SNPs.csv"), sep = ",", quote = FALSE)

  # record task meta
  data.table::fwrite(
    data.table(Algorithm=algorithm, Sample_Fraction=frac, Replicate=repi, n_samples=length(ids_common)),
    file.path(bench_meta, "subset_meta.tsv"),
    sep="\t", quote=FALSE
  )
  writeLines(ids_common, file.path(bench_meta, "sample_ids.txt"))
  writeLines(capture.output(sessionInfo()), file.path(bench_meta, sprintf("%s_session_info_prepare.txt", algorithm)))

  message("Prepared distance matrices for ", algorithm, " ", frac, "% rep ", repi)
}

if (step == "post") {
  RNGversion("4.2.2"); set.seed(123)

  if (!file.exists(gt_path)) stop("Missing ground-truth labels: ", gt_path)
  gt <- openxlsx::read.xlsx(gt_path)
  if (!all(c("Sample.ID", "Cluster") %in% colnames(gt))) {
    stop("Ground truth file must contain columns: Sample.ID, Cluster")
  }
  gt$Sample.ID <- gsub("\\.", "-", gt$Sample.ID)
  colnames(gt)[colnames(gt) == "Cluster"] <- "Cluster_gt"

  # recover frac/rep from subset_meta
  sm <- data.table::fread(file.path(bench_meta, "subset_meta.tsv"))
  frac <- as.integer(sm$Sample_Fraction[1])
  repi <- as.integer(sm$Replicate[1])

  emb_path <- file.path(bench_out, "msne_embeddings.csv")
  clu_path <- file.path(bench_out, "msne_clusters.csv")
  if (!file.exists(emb_path) || !file.exists(clu_path)) {
    stop("Missing MSNE outputs (expected msne_embeddings.csv and msne_clusters.csv in BENCH_OUT_DIR)")
  }

  emb <- data.table::fread(emb_path)
  emb_ids <- emb[[1]]
  emb_mat <- as.matrix(emb[, -1, with = FALSE])
  rownames(emb_mat) <- emb_ids

  clu <- data.table::fread(clu_path)
  clu_ids <- clu[[1]]
  labs <- as.integer(clu[[2]])
  names(labs) <- clu_ids

  common <- intersect(rownames(emb_mat), names(labs))
  emb_mat <- emb_mat[common, , drop = FALSE]
  labs <- labs[common]

  clusters_df <- data.frame(
    Sample.ID = gsub("\\.", "-", common),
    Cluster_pred = as.integer(labs),
    stringsAsFactors = FALSE
  )

  merged <- merge(gt, clusters_df, by = "Sample.ID")
  ari_gt <- if (nrow(merged) == 0) NA_real_ else mclust::adjustedRandIndex(merged$Cluster_gt, merged$Cluster_pred)

  avg_width <- NA_real_
  if (length(unique(labs)) > 1L && nrow(emb_mat) > 2L) {
    sil <- cluster::silhouette(as.integer(labs),
                              dist = Rfast::Dist(emb_mat, method = "euclidean"))
    avg_width <- summary(sil)$avg.width
  }

  tv <- parse_timev(file.path(bench_meta, "timev_msne.txt"))
  elapsed_s <- tv$elapsed_s
  peak_mib <- if (is.na(tv$maxrss_kb)) NA_real_ else tv$maxrss_kb / 1024

  clust_file <- file.path(bench_out, sprintf("%s_clusters_samples_%dpct_rep%02d.tsv.gz", algorithm, frac, repi))
  data.table::fwrite(data.table::as.data.table(clusters_df),
                     clust_file, sep = "\t", quote = FALSE, compress = "gzip")

  perf <- data.table(
    Algorithm = algorithm,
    Sample_Fraction = frac,
    Replicate = repi,
    n_samples = as.integer(length(common)),
    n_clusters = 5L,
    MSNE_k = 50L,
    MSNE_num_walks = 200L,
    MSNE_embed_size = 50L,
    MSNE_window_size = 5L,
    MSNE_walk_length = 40L,
    MSNE_workers = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "10")),
    ARI_to_Ground_Truth = ari_gt,
    Average_Silhouette_Width = avg_width,
    MSNE_Elapsed_Seconds = elapsed_s,
    MSNE_PeakRAM_MiB = peak_mib
  )

  row_path <- file.path(bench_out, sprintf("%s_perf_row_samples_%dpct_rep%02d.tsv", algorithm, frac, repi))
  data.table::fwrite(perf, row_path, sep = "\t", quote = FALSE, na = "NA")

  writeLines(capture.output(sessionInfo()), file.path(bench_meta, sprintf("%s_session_info_post.txt", algorithm)))
  message("Postprocessed ", algorithm, " ", frac, "% rep ", repi)
}

if (step == "merge") {
  out_root <- file.path(bench_run, "out")
  meta_root <- file.path(bench_run, "meta")
  dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
  dir.create(meta_root, recursive = TRUE, showWarnings = FALSE)

  rows <- list.files(bench_run, pattern = "MSNE_perf_row_samples_.*\\.tsv$", recursive = TRUE, full.names = TRUE)
  if (length(rows) == 0) stop("No perf rows found under: ", bench_run)

  dt <- data.table::rbindlist(lapply(rows, data.table::fread), use.names = TRUE, fill = TRUE)
  if (all(c("Sample_Fraction", "Replicate") %in% colnames(dt))) data.table::setorder(dt, Sample_Fraction, Replicate)

  perf_file <- file.path(out_root, sprintf("%s_sample_perturbations_performance.tsv", algorithm))
  data.table::fwrite(dt, perf_file, sep = "\t", quote = FALSE, na = "NA")

  writeLines(capture.output(sessionInfo()), file.path(meta_root, sprintf("%s_session_info_merge.txt", algorithm)))
  message("Wrote: ", perf_file)
}

