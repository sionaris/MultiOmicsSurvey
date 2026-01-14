# Scripts/automated_scripts/benchmark_plotting_functions.R

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(cowplot)
  library(ggrepel)
})

if (!requireNamespace("rcartocolor", quietly = TRUE)) {
  stop("Package 'rcartocolor' is required (for carto_pal). Install with install.packages('rcartocolor').")
}

`%||%` <- function(x, y) if (!is.null(x) && length(x) && nzchar(as.character(x[1]))) x else y
dir_exists <- function(p) is.character(p) && length(p) == 1L && dir.exists(p)

theme_benchmark <- function(base_size = 10, legend = FALSE) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line.x = element_line(color = "black"),
      axis.line.y = element_line(color = "black"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.title = element_text(face = "bold"),
      legend.position = if (isTRUE(legend)) "right" else "none",
      plot.margin = margin(5.5, 70, 5.5, 5.5, "pt")
    )
}

sanitise_title_for_filename <- function(title) {
  x <- tolower(title)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

guess_benchmark_dirs <- function(root = "Results/Performance_benchmarks") {
  feature_dir <- file.path(root, "Feature_perturbations")
  
  sample_dir1 <- file.path(root, "Sample_perturbations")
  sample_dir2 <- file.path(root, "Smaple_perutrbations")
  sample_dir <- if (dir_exists(sample_dir1)) sample_dir1 else sample_dir2
  
  if (!dir_exists(feature_dir)) stop("Feature directory not found: ", feature_dir)
  if (!dir_exists(sample_dir)) {
    stop("Sample directory not found: ", sample_dir,
         " (checked both Sample_perturbations and Smaple_perutrbations)")
  }
  
  list(feature_dir = feature_dir, sample_dir = sample_dir)
}

file_nonempty <- function(paths) {
  if (length(paths) == 0) return(logical(0))
  paths <- as.character(paths)
  ok <- file.exists(paths)
  out <- rep(FALSE, length(paths))
  if (!any(ok)) return(out)
  fi <- file.info(paths[ok])
  sz <- fi$size
  out[ok] <- is.finite(sz) & !is.na(sz) & (sz > 0)
  out
}

# -----------------------------
# Categories (your mapping)
# -----------------------------
similarity_network_methods <- c("ab-SNF", "ANF", "MDICC", "MSNE", "NEMO", "RWR-F", "RWR-NF", "SNF", "Spectrum")
multiple_kernel_learning   <- c("CIMLR", "KLIC", "wMKL") # "mixKernel"
matrix_factorization       <- c("MFA", "MOFA", "LRAcluster")
graph_methods              <- c("MONET") #, "PAMOGK")
bayesian_methods           <- c("iClusterBayes")
cc_ensemble                <- c("COCA")

category_colors <- c(
  "Similarity Network" = rcartocolor::carto_pal("Bold", n = 12)[1],
  "Multiple Kernel Learning" = rcartocolor::carto_pal("Bold", n = 12)[2],
  "Matrix Factorization" = rcartocolor::carto_pal("Antique", n = 12)[5],
  "Graph-based Methods" = rcartocolor::carto_pal("Bold", n = 12)[4],
  "Bayesian" = rcartocolor::carto_pal("Bold", n = 12)[11],
  "Consensus/Ensemble Clustering" = rcartocolor::carto_pal("Bold", n = 12)[9]
)

method_categories <- c(
  setNames(rep("Similarity Network", length(similarity_network_methods)), similarity_network_methods),
  setNames(rep("Multiple Kernel Learning", length(multiple_kernel_learning)), multiple_kernel_learning),
  setNames(rep("Matrix Factorization", length(matrix_factorization)), matrix_factorization),
  setNames(rep("Graph-based Methods", length(graph_methods)), graph_methods),
  setNames(rep("Bayesian", length(bayesian_methods)), bayesian_methods),
  setNames(rep("Consensus/Ensemble Clustering", length(cc_ensemble)), cc_ensemble)
)

add_categories <- function(dt, method_categories_map = method_categories) {
  dt[, Category := unname(method_categories_map[as.character(Algorithm)])]
  dt[is.na(Category) | Category == "", Category := "Matrix Factorization"]
  dt
}

# -----------------------------
# Robust normalisation helpers (ROW-WISE)
# -----------------------------
normalise_percent_vec <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  out <- x
  
  # fraction -> percent
  idx <- !is.na(out) & out <= 1.01
  out[idx] <- out[idx] * 100
  
  # basis points -> percent
  idx <- !is.na(out) & out > 110
  out[idx] <- out[idx] / 100
  
  out
}

coalesce_numeric_cols <- function(dt, cols) {
  cols <- intersect(cols, names(dt))
  if (length(cols) == 0) return(rep(NA_real_, nrow(dt)))
  vecs <- lapply(cols, function(cn) suppressWarnings(as.numeric(dt[[cn]])))
  do.call(data.table::fcoalesce, vecs)
}

# -----------------------------
# Subset size standardisation (ROW-WISE)
# -----------------------------
standardise_subset_cols <- function(dt, mode = c("feature", "sample")) {
  mode <- match.arg(mode)
  
  if (mode == "feature") {
    fp1 <- if ("Feature_Percent" %in% names(dt)) normalise_percent_vec(dt[["Feature_Percent"]]) else rep(NA_real_, nrow(dt))
    fp2 <- if ("Feature_Centile" %in% names(dt)) normalise_percent_vec(dt[["Feature_Centile"]]) else rep(NA_real_, nrow(dt))
    fp3 <- if ("Feature_Percentile" %in% names(dt)) normalise_percent_vec(100 * suppressWarnings(as.numeric(dt[["Feature_Percentile"]]))) else rep(NA_real_, nrow(dt))
    fp <- data.table::fcoalesce(fp1, fp2, fp3)
    
    dt[, Subset_Percent := as.integer(round(fp))]
    dt[, Subset_Absolute := if ("p_total" %in% names(dt)) as.integer(p_total) else NA_integer_]
  } else {
    sp1 <- if ("Sample_Percent" %in% names(dt)) normalise_percent_vec(dt[["Sample_Percent"]]) else rep(NA_real_, nrow(dt))
    sp2 <- if ("Sample_Fraction" %in% names(dt)) normalise_percent_vec(dt[["Sample_Fraction"]]) else rep(NA_real_, nrow(dt))
    sp <- data.table::fcoalesce(sp1, sp2)
    
    dt[, Subset_Percent := as.integer(round(sp))]
    dt[, Subset_Absolute := if ("n_samples" %in% names(dt)) as.integer(n_samples) else NA_integer_]
  }
  
  dt
}

# -----------------------------
# Runtime + memory standardisation (covers your variants)
# -----------------------------
standardise_time_memory <- function(dt) {
  nm <- names(dt)
  
  # TIME candidates
  time_cols <- nm[
    grepl("Elapsed_Seconds", nm, ignore.case = TRUE) |
      grepl("MainLoop_Elapsed_Seconds", nm, ignore.case = TRUE) |
      grepl("Elapsed_Time", nm, ignore.case = TRUE)
  ]
  dt[, Time_s := coalesce_numeric_cols(dt, time_cols)]
  
  # MEMORY candidates
  mem_cols <- nm[
    grepl("PeakRAM", nm, ignore.case = TRUE) |
      grepl("PeakRSS", nm, ignore.case = TRUE) |
      grepl("MaxRSS", nm, ignore.case = TRUE) |
      grepl("Peak.*MiB", nm, ignore.case = TRUE) |
      grepl("Peak.*MB", nm, ignore.case = TRUE)
  ]
  mem_cols <- mem_cols[!grepl("Elapsed", mem_cols, ignore.case = TRUE)]
  
  mb_cols <- mem_cols[grepl("(_|\\b)MB(\\b|_)", mem_cols) & !grepl("MiB", mem_cols, ignore.case = TRUE)]
  mib_cols <- setdiff(mem_cols, mb_cols)
  
  mem_mib <- if (length(mib_cols) > 0) coalesce_numeric_cols(dt, mib_cols) else rep(NA_real_, nrow(dt))
  mem_mb  <- if (length(mb_cols)  > 0) coalesce_numeric_cols(dt, mb_cols)  else rep(NA_real_, nrow(dt))
  
  dt[, Peak_MiB := data.table::fcoalesce(mem_mib, mem_mb * (1000 / 1024))]
  
  dt
}

safe_median <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (all(is.na(x))) return(NA_real_)
  stats::median(x, na.rm = TRUE)
}

aggregate_by_subset <- function(dt) {
  dt[, .(
    ARI_to_Ground_Truth = safe_median(ARI_to_Ground_Truth),
    Average_Silhouette_Width = safe_median(Average_Silhouette_Width),
    Time_s = safe_median(Time_s),
    Peak_MiB = safe_median(Peak_MiB),
    n_rows = .N
  ), by = .(Algorithm, Mode, Category, Subset_Percent, Subset_Absolute)]
}

# -----------------------------
# File discovery: out/ dirs and a perturbations TSV
# -----------------------------
find_out_dirs <- function(method_dir) {
  if (!dir.exists(method_dir)) return(character(0))
  dd <- list.dirs(method_dir, recursive = TRUE, full.names = TRUE)
  unique(dd[basename(dd) == "out"])
}

pick_perf_file <- function(out_dirs, method_name, mode = c("feature", "sample")) {
  mode <- match.arg(mode)
  if (length(out_dirs) == 0) return(character(0))
  
  candidates <- character(0)
  
  for (od in out_dirs) {
    if (!dir.exists(od)) next
    
    # MONET special-case
    if (toupper(method_name) == "MONET") {
      mf <- file.path(od, "monet_perf_rows.tsv")
      if (file.exists(mf) && file_nonempty(mf)) return(mf)
    }
    
    files <- list.files(od, full.names = TRUE)
    files <- files[grepl("\\.tsv$", files, ignore.case = TRUE)]
    files <- files[file_nonempty(files)]
    if (length(files) == 0) next
    
    pert <- files[grepl("perturbations", basename(files), ignore.case = TRUE)]
    if (length(pert) > 0) {
      if (mode == "feature") {
        ms <- pert[grepl("feature", basename(pert), ignore.case = TRUE)]
      } else {
        ms <- pert[grepl("sample", basename(pert), ignore.case = TRUE)]
      }
      if (length(ms) > 0) pert <- ms
      sz <- file.info(pert)$size
      candidates <- c(candidates, pert[order(-sz)][1])
      next
    }
    
    pr <- files[grepl("perf_row", basename(files), ignore.case = TRUE)]
    if (length(pr) > 0) candidates <- c(candidates, pr)
  }
  
  candidates <- unique(candidates)
  if (length(candidates) == 0) return(character(0))
  
  sz <- file.info(candidates)$size
  candidates[order(-sz, nchar(candidates))][1]
}

# -----------------------------
# Data ingestion
# -----------------------------
read_mode_perf <- function(mode = c("feature", "sample"),
                           root = "Results/Performance_benchmarks",
                           method_categories_map = method_categories) {
  mode <- match.arg(mode)
  dirs <- guess_benchmark_dirs(root)
  base_dir <- if (mode == "feature") dirs$feature_dir else dirs$sample_dir
  
  method_dirs <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)
  method_dirs <- method_dirs[normalizePath(method_dirs, winslash = "/") != normalizePath(base_dir, winslash = "/")]
  if (length(method_dirs) == 0) stop("No method subdirectories found under: ", base_dir)
  
  all_rows <- list()
  
  for (mdir in method_dirs) {
    method_name <- basename(mdir)
    out_dirs <- find_out_dirs(mdir)
    pf <- pick_perf_file(out_dirs, method_name, mode = mode)
    if (length(pf) == 0 || !file.exists(pf) || !file_nonempty(pf)) next
    
    dt <- tryCatch(data.table::fread(pf), error = function(e) NULL)
    if (is.null(dt) || nrow(dt) == 0 || ncol(dt) == 0) next
    
    if (!("Algorithm" %in% names(dt))) dt[, Algorithm := method_name]
    dt[, Mode := mode]
    dt[, Source_File := pf]
    
    all_rows[[length(all_rows) + 1L]] <- dt
  }
  
  if (length(all_rows) == 0) stop("No readable performance TSVs found for mode=", mode)
  
  out <- data.table::rbindlist(all_rows, use.names = TRUE, fill = TRUE)
  
  out <- add_categories(out, method_categories_map = method_categories_map)
  out <- standardise_subset_cols(out, mode = mode)
  out <- standardise_time_memory(out)
  
  out <- aggregate_by_subset(out)
  out
}

# -----------------------------
# Plotting
# -----------------------------
add_end_labels <- function(p, dt, x_col, y_col, label_col = "Algorithm") {
  dt2 <- dt %>%
    as.data.frame() %>%
    dplyr::filter(!is.na(.data[[x_col]]), !is.na(.data[[y_col]])) %>%
    dplyr::group_by(.data[[label_col]]) %>%
    dplyr::slice_max(order_by = .data[[x_col]], n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()
  
  if (nrow(dt2) == 0) return(p)
  
  p + ggrepel::geom_text_repel(
    data = dt2,
    aes(x = .data[[x_col]], y = .data[[y_col]], label = .data[[label_col]]),
    direction = "y",
    hjust = 0,
    nudge_x = 0.8,
    segment.size = 0.2,
    min.segment.length = 0,
    size = 3,
    max.overlaps = Inf
  ) +
    coord_cartesian(clip = "off")
}

plot_metric_lines <- function(dt,
                              mode = c("feature", "sample"),
                              metric = c("ARI_to_Ground_Truth", "Peak_MiB", "Time_s"),
                              x_axis = c("percent", "absolute"),
                              title = NULL,
                              ylab = NULL,
                              category_colors_map = category_colors,
                              legend = FALSE) {
  mode <- match.arg(mode)
  metric <- match.arg(metric)
  x_axis <- match.arg(x_axis)
  
  x_col <- if (x_axis == "percent") "Subset_Percent" else "Subset_Absolute"
  x_lab <- if (mode == "feature") {
    if (x_axis == "percent") "Feature subset size (%)" else "Number of features (p)"
  } else {
    if (x_axis == "percent") "Sample subset size (%)" else "Number of samples (n)"
  }
  
  plot_dt <- dt %>%
    as.data.frame() %>%
    mutate(
      Algorithm = as.factor(Algorithm),
      Category = factor(Category, levels = names(category_colors_map))
    ) %>%
    filter(!is.na(.data[[x_col]]), !is.na(.data[[metric]]))
  
  p <- ggplot(plot_dt, aes(x = .data[[x_col]], y = .data[[metric]], group = Algorithm, color = Category))
  
  if (x_axis == "percent") {
    p <- p +
      geom_vline(xintercept = seq(0, 100, 10), linetype = "dotted", linewidth = 0.2, alpha = 0.7) +
      scale_x_continuous(limits = c(-5, 110), breaks = seq(0, 100, 10))
  }
  
  p <- p +
    geom_line(linewidth = 0.55, alpha = 0.9) +
    geom_point(size = 1.1) +
    scale_color_manual(values = category_colors_map, drop = FALSE) +
    labs(
      title = title %||% "",
      x = x_lab,
      y = ylab %||% metric
    ) +
    theme_benchmark(base_size = 10, legend = legend)
  
  add_end_labels(p, plot_dt, x_col = x_col, y_col = metric, label_col = "Algorithm")
}

save_plot_png_pdf <- function(p, title, out_dir,
                              dpi = 700,
                              width_px = 3840, height_px = 2160) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  fn <- paste0("benchmarks_", sanitise_title_for_filename(title))
  
  ggsave(
    filename = file.path(out_dir, paste0(fn, ".png")),
    plot = p,
    dpi = dpi,
    width = width_px, height = height_px, units = "px"
  )
  ggsave(
    filename = file.path(out_dir, paste0(fn, ".pdf")),
    plot = p,
    dpi = dpi,
    width = width_px, height = height_px, units = "px"
  )
}

# -----------------------------
# Empirical scaling models (TIME)
# -----------------------------

o_notation <- function(var, exponent, tol = 0.15) {
  if (!is.finite(exponent)) return(NA_character_)
  k <- round(exponent)
  if (abs(exponent - k) < tol && k >= 0 && k <= 8) {
    if (k == 0) return("O(1)")
    if (k == 1) return(sprintf("O(%s)", var))
    return(sprintf("O(%s^%d)", var, k))
  }
  sprintf("O(%s^%.2f)", var, exponent)
}

fit_powerlaw <- function(x, y) {
  x <- suppressWarnings(as.numeric(x))
  y <- suppressWarnings(as.numeric(y))
  keep <- is.finite(x) & is.finite(y) & x > 0 & y > 0
  x <- x[keep]; y <- y[keep]
  if (length(x) < 4 || length(unique(x)) < 4) return(NULL)
  
  fit <- stats::lm(log(y) ~ log(x))
  sm <- summary(fit)
  ci <- tryCatch(confint(fit), error = function(e) NULL)
  
  slope <- unname(coef(fit)[2])
  intercept <- unname(coef(fit)[1])
  
  list(
    intercept = intercept,
    slope = slope,
    r2 = unname(sm$r.squared),
    n = length(x),
    slope_lo = if (!is.null(ci)) unname(ci["log(x)", 1]) else NA_real_,
    slope_hi = if (!is.null(ci)) unname(ci["log(x)", 2]) else NA_real_
  )
}

# Decide whether Subset_Absolute is already the subset size, or whether we should derive it
# from percent via (max base) * percent/100. Do this PER ALGORITHM.
choose_effective_size <- function(base, pct) {
  base <- suppressWarnings(as.numeric(base))
  pct  <- suppressWarnings(as.numeric(pct))
  
  cand1 <- base
  full  <- suppressWarnings(max(base, na.rm = TRUE))
  cand2 <- if (is.finite(full)) full * pct / 100 else rep(NA_real_, length(base))
  
  # If cand1 barely varies, it's probably "full size" repeated -> use cand2
  u1 <- unique(cand1[is.finite(cand1)])
  if (length(u1) < 4) return(cand2)
  
  # Otherwise pick whichever tracks percent better (Spearman)
  c1 <- suppressWarnings(cor(cand1, pct, method = "spearman", use = "complete.obs"))
  c2 <- suppressWarnings(cor(cand2, pct, method = "spearman", use = "complete.obs"))
  
  if (is.finite(c2) && (!is.finite(c1) || abs(c2) > abs(c1))) cand2 else cand1
}

fit_time_scaling_models <- function(root = "Results/Performance_benchmarks",
                                    method_categories_map = method_categories,
                                    min_points = 4) {
  # Use the same processed/aggregated tables your plotting uses
  feat <- read_mode_perf("feature", root = root, method_categories_map = method_categories_map)
  samp <- read_mode_perf("sample",  root = root, method_categories_map = method_categories_map)
  
  # Effective sizes per algorithm
  feat[, p_eff := choose_effective_size(Subset_Absolute, Subset_Percent), by = Algorithm]
  samp[, n_eff := choose_effective_size(Subset_Absolute, Subset_Percent), by = Algorithm]
  
  res_feat <- feat[, {
    fit <- fit_powerlaw(p_eff, Time_s)
    if (is.null(fit) || fit$n < min_points) {
      NULL
    } else {
      data.table(
        Mode = "feature",
        Algorithm = .BY$Algorithm,
        Category = .BY$Category,
        Predictor = "p",
        Intercept = fit$intercept,
        Exponent = fit$slope,
        Exponent_Lo = fit$slope_lo,
        Exponent_Hi = fit$slope_hi,
        R2 = fit$r2,
        N_points = fit$n,
        O_notation = o_notation("p", fit$slope)
      )
    }
  }, by = .(Algorithm, Category)]
  
  res_samp <- samp[, {
    fit <- fit_powerlaw(n_eff, Time_s)
    if (is.null(fit) || fit$n < min_points) {
      NULL
    } else {
      data.table(
        Mode = "sample",
        Algorithm = .BY$Algorithm,
        Category = .BY$Category,
        Predictor = "n",
        Intercept = fit$intercept,
        Exponent = fit$slope,
        Exponent_Lo = fit$slope_lo,
        Exponent_Hi = fit$slope_hi,
        R2 = fit$r2,
        N_points = fit$n,
        O_notation = o_notation("n", fit$slope)
      )
    }
  }, by = .(Algorithm, Category)]
  
  data.table::rbindlist(list(res_feat, res_samp), use.names = TRUE, fill = TRUE)
}

write_time_scaling_table <- function(root = "Results/Performance_benchmarks",
                                     out_dir = "Results/Comparisons",
                                     filename = "scaling_time_models.csv",
                                     min_points = 4,
                                     method_categories_map = method_categories) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  tab <- fit_time_scaling_models(
    root = root,
    method_categories_map = method_categories_map,
    min_points = min_points
  )
  out_path <- file.path(out_dir, filename)
  data.table::fwrite(tab, out_path)
  invisible(tab)
}

# -----------------------------
# Public API
# -----------------------------
make_benchmark_plots <- function(
    root = "Results/Performance_benchmarks",
    out_dir = "Results/Comparisons",
    x_axis = c("percent", "absolute"),
    save = TRUE,
    dpi = 700,
    width_px = 3840,
    height_px = 2160,
    method_categories_map = method_categories,
    category_colors_map = category_colors,
    legend = FALSE
) {
  x_axis <- match.arg(x_axis)
  
  feat <- read_mode_perf("feature", root = root, method_categories_map = method_categories_map)
  samp <- read_mode_perf("sample", root = root, method_categories_map = method_categories_map)
  
  t_feat_ari  <- "Agreement with ground-truth clustering under feature perturbations"
  t_feat_mem  <- "Peak memory under feature perturbations"
  t_feat_time <- "Runtime under feature perturbations"
  
  t_samp_ari  <- "Agreement with ground-truth clustering under sample perturbations"
  t_samp_mem  <- "Peak memory under sample perturbations"
  t_samp_time <- "Runtime under sample perturbations"
  
  p_feat_ari <- plot_metric_lines(feat, "feature", "ARI_to_Ground_Truth", x_axis, t_feat_ari,
                                  "ARI (subset vs ground truth)", category_colors_map, legend)
  p_feat_mem <- plot_metric_lines(feat, "feature", "Peak_MiB", x_axis, t_feat_mem,
                                  "Peak memory (MiB)", category_colors_map, legend)
  p_feat_time <- plot_metric_lines(feat, "feature", "Time_s", x_axis, t_feat_time,
                                   "Runtime (seconds)", category_colors_map, legend)
  
  p_samp_ari <- plot_metric_lines(samp, "sample", "ARI_to_Ground_Truth", x_axis, t_samp_ari,
                                  "ARI (subset vs ground truth)", category_colors_map, legend)
  p_samp_mem <- plot_metric_lines(samp, "sample", "Peak_MiB", x_axis, t_samp_mem,
                                  "Peak memory (MiB)", category_colors_map, legend)
  p_samp_time <- plot_metric_lines(samp, "sample", "Time_s", x_axis, t_samp_time,
                                   "Runtime (seconds)", category_colors_map, legend)
  
  plots <- list(
    feature_ARI = list(plot = p_feat_ari, title = t_feat_ari),
    feature_memory = list(plot = p_feat_mem, title = t_feat_mem),
    feature_time = list(plot = p_feat_time, title = t_feat_time),
    sample_ARI = list(plot = p_samp_ari, title = t_samp_ari),
    sample_memory = list(plot = p_samp_mem, title = t_samp_mem),
    sample_time = list(plot = p_samp_time, title = t_samp_time)
  )
  
  if (isTRUE(save)) {
    for (nm in names(plots)) {
      save_plot_png_pdf(
        plots[[nm]]$plot,
        plots[[nm]]$title,
        out_dir = out_dir,
        dpi = dpi,
        width_px = width_px,
        height_px = height_px
      )
    }
  }
  
  invisible(plots)
}
