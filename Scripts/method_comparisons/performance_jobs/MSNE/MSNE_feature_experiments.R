#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
step <- if (length(args) >= 1) args[[1]] else ""
if (!step %in% c("prepare", "post", "merge")) {
  stop("Usage: Rscript Jobs/MSNE/MSNE_feature_experiments.R <prepare|post|merge>")
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

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && nzchar(a)) a else b

algorithm <- "MSNE"

bench_out <- Sys.getenv("BENCH_OUT_DIR", "")
bench_meta <- Sys.getenv("BENCH_META_DIR", "")
bench_run <- Sys.getenv("BENCH_RUN_DIR", "")
if (!nzchar(bench_out) && step != "merge") stop("BENCH_OUT_DIR not set")
if (!nzchar(bench_meta) && step != "merge") stop("BENCH_META_DIR not set")
if (!nzchar(bench_run) && step == "merge") stop("BENCH_RUN_DIR not set")

input_path <- Sys.getenv("MSNE_INPUT_RDS", "Resources/mm_input.rds")
feat_rank_path <- Sys.getenv("MSNE_FEAT_RANK", "Resources/Feature_perturbations/feature_rankings.tsv.gz")
gt_path <- Sys.getenv("MSNE_GROUND_TRUTH", file.path("Resources/algorithms", "MSNE_ground_truth_labels.xlsx"))

centile <- as.numeric(Sys.getenv("MSNE_CENTILE", "NA"))
pct <- as.integer(round(centile * 100))

parse_timev <- function(path) {
  # GNU time -v -o file => lines include:
  # "Elapsed (wall clock) time (h:mm:ss or m:ss): 0:26.53"
  # "Maximum resident set size (kbytes): 10836"
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
    # formats: "m:ss.xx" or "h:mm:ss"
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

if (step == "prepare") {
  RNGversion("4.2.2"); set.seed(123)

  if (!file.exists(input_path)) stop("Missing input: ", input_path)
  if (!file.exists(feat_rank_path)) stop("Missing feature rankings: ", feat_rank_path)

  input_full <- readRDS(input_path)

  feat_rank <- data.table::fread(feat_rank_path)
  need_cols <- c("modality", "rank", "feature_id")
  if (!all(need_cols %in% colnames(feat_rank))) {
    stop("feature_rankings.tsv.gz must contain columns: modality, rank, feature_id")
  }

  modalities <- intersect(unique(feat_rank$modality), names(input_full))
  if (length(modalities) == 0) stop("No overlap between ranked modalities and mm_input.rds")

  get_top_features <- function(dt, cval) {
    out <- vector("list", length(modalities)); names(out) <- modalities
    for (m in modalities) {
      md <- dt[modality == m][order(rank)]
      top_n <- max(1L, ceiling(cval * nrow(md)))
      out[[m]] <- md$feature_id[seq_len(top_n)]
    }
    out
  }

  subset_input <- function(inp, top_feats) {
    out <- inp
    for (m in names(top_feats)) {
      if (!is.null(out[[m]])) {
        keep <- intersect(top_feats[[m]], rownames(out[[m]]))
        out[[m]] <- out[[m]][keep, , drop = FALSE]
        if (nrow(out[[m]]) == 0L) out[[m]] <- NULL
      }
    }
    out
  }

  top_feats <- get_top_features(feat_rank, centile)
  input_sub <- subset_input(input_full, top_feats)

  mods_use <- names(input_sub)[!vapply(input_sub, is.null, logical(1))]
  if (length(mods_use) == 0) stop("All modalities became NULL after subsetting at ", pct, "%")

  ref_mod <- if ("RNAseq" %in% mods_use) "RNAseq" else mods_use[1]
  sample_ids <- colnames(input_sub[[ref_mod]])
  if (length(sample_ids) == 0) stop("Missing sample IDs in reference modality")

  out_dists <- file.path(bench_out, "MO_Dists")
  dir.create(out_dists, recursive = TRUE, showWarnings = FALSE)

  continuous <- intersect(c("RNAseq", "CNV", "Methylation", "miRNA"), mods_use)

  for (m in continuous) {
    x <- t(as.matrix(input_sub[[m]]))  # samples x features
    d <- SNFtool::dist2(x, x)
    dimnames(d) <- list(sample_ids, sample_ids)
    fp <- file.path(out_dists, paste0(m, ".csv"))
    data.table::fwrite(data.table::as.data.table(d, keep.rownames = "Sample.ID"),
                       fp, sep = ",", quote = FALSE)
  }

  if ("SNPs" %in% mods_use) {
    x <- t(as.matrix(input_sub[["SNPs"]]))
    x <- apply(x, 2, as.numeric)
    d <- as.matrix(dist(x, method = "binary"))
    dimnames(d) <- list(sample_ids, sample_ids)
    fp <- file.path(out_dists, "SNPs.csv")
    data.table::fwrite(data.table::as.data.table(d, keep.rownames = "Sample.ID"),
                       fp, sep = ",", quote = FALSE)
  }

  p_by_mod <- sapply(mods_use, function(m) nrow(input_sub[[m]]))
  p_total <- sum(p_by_mod, na.rm = TRUE)

  meta <- data.table(
    Algorithm = algorithm,
    Feature_Centile = centile,
    Feature_Percent = pct,
    n_samples = as.integer(length(sample_ids)),
    p_total = as.integer(p_total),
    p_RNAseq = as.integer(p_by_mod[["RNAseq"]] %||% NA_integer_),
    p_CNV = as.integer(p_by_mod[["CNV"]] %||% NA_integer_),
    p_Methylation = as.integer(p_by_mod[["Methylation"]] %||% NA_integer_),
    p_miRNA = as.integer(p_by_mod[["miRNA"]] %||% NA_integer_),
    p_SNPs = as.integer(p_by_mod[["SNPs"]] %||% NA_integer_)
  )
  data.table::fwrite(meta, file.path(bench_meta, sprintf("%s_prepare_meta_%dpct.tsv", algorithm, pct)),
                     sep = "\t", quote = FALSE)

  writeLines(sample_ids, file.path(bench_meta, "sample_ids.txt"))
  writeLines(capture.output(sessionInfo()), file.path(bench_meta, sprintf("%s_session_info_prepare.txt", algorithm)))

  message("Prepared distance matrices for ", algorithm, " ", pct, "%")
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

  # clusters output
  clusters_df <- data.frame(
    Sample.ID = gsub("\\.", "-", common),
    Cluster_pred = as.integer(labs),
    stringsAsFactors = FALSE
  )

  # ARI
  merged <- merge(gt, clusters_df, by = "Sample.ID")
  ari_gt <- if (nrow(merged) == 0) NA_real_ else mclust::adjustedRandIndex(merged$Cluster_gt, merged$Cluster_pred)

  # silhouette on embeddings
  avg_width <- NA_real_
  if (length(unique(labs)) > 1L && nrow(emb_mat) > 2L) {
    sil <- cluster::silhouette(as.integer(labs),
                              dist = Rfast::Dist(emb_mat, method = "euclidean"))
    avg_width <- summary(sil)$avg.width
  }

  # time/mem from GNU time
  timev_path <- file.path(bench_meta, "timev_msne.txt")
  tv <- parse_timev(timev_path)
  elapsed_s <- tv$elapsed_s
  peak_mib <- if (is.na(tv$maxrss_kb)) NA_real_ else tv$maxrss_kb / 1024

  # write clusters
  clust_file <- file.path(bench_out, sprintf("%s_clusters_%dpct.tsv.gz", algorithm, pct))
  data.table::fwrite(data.table::as.data.table(clusters_df),
                     clust_file, sep = "\t", quote = FALSE, compress = "gzip")

  perf <- data.table(
    Algorithm = algorithm,
    Feature_Centile = centile,
    Feature_Percent = pct,
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

  row_path <- file.path(bench_out, sprintf("%s_perf_row_%dpct.tsv", algorithm, pct))
  data.table::fwrite(perf, row_path, sep = "\t", quote = FALSE, na = "NA")

  writeLines(capture.output(sessionInfo()), file.path(bench_meta, sprintf("%s_session_info_post.txt", algorithm)))

  message("Postprocessed ", algorithm, " ", pct, "%")
}

if (step == "merge") {
  out_root <- file.path(bench_run, "out")
  meta_root <- file.path(bench_run, "meta")
  dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
  dir.create(meta_root, recursive = TRUE, showWarnings = FALSE)

  rows <- list.files(bench_run, pattern = "MSNE_perf_row_.*pct\\.tsv$", recursive = TRUE, full.names = TRUE)
  if (length(rows) == 0) stop("No perf rows found under: ", bench_run)

  dt <- data.table::rbindlist(lapply(rows, data.table::fread), use.names = TRUE, fill = TRUE)
  if ("Feature_Percent" %in% colnames(dt)) data.table::setorder(dt, Feature_Percent)

  perf_file <- file.path(out_root, sprintf("%s_feature_perturbations_performance.tsv", algorithm))
  data.table::fwrite(dt, perf_file, sep = "\t", quote = FALSE, na = "NA")

  writeLines(capture.output(sessionInfo()), file.path(meta_root, sprintf("%s_session_info_merge.txt", algorithm)))
  message("Wrote: ", perf_file)
}

