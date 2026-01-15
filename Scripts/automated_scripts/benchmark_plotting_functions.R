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
# Empirical scaling models (TIME) -- robust to missing absolute sizes (e.g., MSNE feature mode)
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

estimate_full_size_global <- function(dt, abs_col = "Subset_Absolute", pct_col = "Subset_Percent") {
  abs <- suppressWarnings(as.numeric(dt[[abs_col]]))
  pct <- suppressWarnings(as.numeric(dt[[pct_col]]))
  keep <- is.finite(abs) & is.finite(pct) & pct > 0 & pct <= 100 & abs > 0
  if (!any(keep)) return(NA_real_)
  # full ~= abs / (pct/100)
  full_est <- abs[keep] / (pct[keep] / 100)
  stats::median(full_est[is.finite(full_est)], na.rm = TRUE)
}

build_effective_size <- function(dt_group, full_global, abs_col = "Subset_Absolute", pct_col = "Subset_Percent") {
  abs <- suppressWarnings(as.numeric(dt_group[[abs_col]]))
  pct <- suppressWarnings(as.numeric(dt_group[[pct_col]]))
  
  abs_finite <- abs[is.finite(abs)]
  # If abs varies meaningfully (>=4 distinct), treat it as effective size already
  if (length(unique(abs_finite)) >= 4) {
    return(list(size = abs, source = "recorded_absolute"))
  }
  
  # Otherwise reconstruct from percent using global full size
  if (is.finite(full_global)) {
    size <- full_global * (pct / 100)
    return(list(size = size, source = "reconstructed_from_percent"))
  }
  
  # Last resort: cannot build effective size
  list(size = rep(NA_real_, nrow(dt_group)), source = "unavailable")
}

fit_time_scaling_models <- function(root = "Results/Performance_benchmarks",
                                    method_categories_map = method_categories,
                                    min_points = 4) {
  feat <- read_mode_perf("feature", root = root, method_categories_map = method_categories_map)
  samp <- read_mode_perf("sample",  root = root, method_categories_map = method_categories_map)
  
  # Global "full" sizes inferred from methods that have Subset_Absolute + Subset_Percent
  p_full_global <- estimate_full_size_global(feat, abs_col = "Subset_Absolute", pct_col = "Subset_Percent")
  n_full_global <- estimate_full_size_global(samp, abs_col = "Subset_Absolute", pct_col = "Subset_Percent")
  
  # Compute effective sizes per algorithm (plus record how it was derived)
  feat[, c("p_eff", "p_size_source") := {
    tmp <- build_effective_size(.SD, full_global = p_full_global,
                                abs_col = "Subset_Absolute", pct_col = "Subset_Percent")
    list(tmp$size, tmp$source)
  }, by = Algorithm]
  
  samp[, c("n_eff", "n_size_source") := {
    tmp <- build_effective_size(.SD, full_global = n_full_global,
                                abs_col = "Subset_Absolute", pct_col = "Subset_Percent")
    list(tmp$size, tmp$source)
  }, by = Algorithm]
  
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
        Size_Source = p_size_source[1],
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
        Size_Source = n_size_source[1],
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
# Stability at 90% sample perturbations
# -----------------------------

# ---- small utilities ----
if (!exists("%||%")) `%||%` <- function(a, b) if (!is.null(a)) a else b

if (!exists("file_nonempty")) {
  file_nonempty <- function(x) {
    x <- as.character(x)
    ok <- file.exists(x)
    ok[ok] <- file.info(x[ok])$size > 0
    ok
  }
}

if (!exists("guess_benchmark_dirs")) {
  guess_benchmark_dirs <- function(root) {
    if (!dir.exists(root)) stop("Root directory not found: ", root)
    top <- list.dirs(root, recursive = FALSE, full.names = TRUE)
    looks_like_method_root <- any(vapply(top, function(p) {
      kids <- list.dirs(p, recursive = FALSE, full.names = FALSE)
      any(grepl("^(job|array)_", kids, ignore.case = TRUE))
    }, logical(1)))
    if (looks_like_method_root) return(list(sample_dir = root))
    cand <- top[grepl("sample", basename(top), ignore.case = TRUE)]
    if (length(cand) >= 1) return(list(sample_dir = cand[1]))
    list(sample_dir = root)
  }
}

# ---- categories (uses your vectors if present; else uses existing method_categories if already defined) ----
build_method_categories <- function() {
  groups <- list(
    "Similarity Network"            = "similarity_network_methods",
    "Multiple Kernel Learning"      = "multiple_kernel_learning",
    "Matrix Factorization"          = "matrix_factorization",
    "Graph-based Methods"           = "graph_methods",
    "Bayesian"                      = "bayesian",
    "Consensus/Ensemble Clustering" = "cc_ensemble"
  )
  
  out <- character(0)
  for (pretty in names(groups)) {
    var <- groups[[pretty]]
    if (!exists(var, inherits = TRUE)) next
    v <- get(var, inherits = TRUE)
    v <- as.character(v)
    v <- v[nzchar(v)]
    if (!length(v)) next
    out <- c(out, stats::setNames(rep(pretty, length(v)), v))
  }
  out
}

get_method_category <- function(alg) {
  mc <- NULL
  if (exists("method_categories", inherits = TRUE)) mc <- get("method_categories", inherits = TRUE)
  if (is.null(mc) || length(mc) == 0) mc <- build_method_categories()
  
  if (is.environment(mc)) return(mc[[alg]] %||% NA_character_)
  if (is.list(mc)) return(mc[[alg]] %||% NA_character_)
  
  if (is.atomic(mc) && !is.null(names(mc))) {
    v <- unname(mc[alg])
    if (length(v) == 1 && !is.na(v)) return(as.character(v))
    return(NA_character_)
  }
  
  if (is.data.frame(mc) && all(c("Algorithm", "Category") %in% names(mc))) {
    hit <- mc$Category[match(alg, mc$Algorithm)]
    return(if (is.na(hit)) NA_character_ else as.character(hit))
  }
  
  NA_character_
}

# ---- percent inference ----
normalize_percent_value <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (!is.finite(x)) return(NA_real_)
  if (x <= 1.2) return(100 * x)         # 0.9 -> 90
  if (x > 100) return(x / 100)          # 9000 -> 90
  x
}

task_id_to_pct <- function(task_id) {
  task_id <- suppressWarnings(as.integer(task_id))
  if (!is.finite(task_id)) return(NA_real_)
  if (task_id >= 1L  && task_id <= 10L) return(10)
  if (task_id >= 11L && task_id <= 20L) return(20)
  if (task_id >= 21L && task_id <= 30L) return(50)
  if (task_id >= 31L && task_id <= 40L) return(70)
  if (task_id >= 41L && task_id <= 50L) return(90)
  NA_real_
}

infer_task_id_from_path <- function(path) {
  path <- as.character(path)[1]
  m <- regmatches(path, regexpr("(?i)task_([0-9]+)", path, perl = TRUE))
  if (!length(m) || !nzchar(m)) return(NA_integer_)
  suppressWarnings(as.integer(sub("(?i)task_", "", m, perl = TRUE)))
}

infer_pct_from_task <- function(path) task_id_to_pct(infer_task_id_from_path(path))

infer_pct_from_filename <- function(path) {
  path <- as.character(path)[1]
  bn <- basename(path)
  
  m1 <- regmatches(bn, regexpr("[0-9]+(?=pct)", bn, perl = TRUE, ignore.case = TRUE))
  if (length(m1) && nzchar(m1)) return(normalize_percent_value(as.numeric(m1)))
  
  m2 <- regmatches(path, regexpr("[0-9]+(?=pct)", path, perl = TRUE, ignore.case = TRUE))
  if (length(m2) && nzchar(m2)) return(normalize_percent_value(as.numeric(m2)))
  
  pct_task <- infer_pct_from_task(path)
  if (is.finite(pct_task)) return(pct_task)
  
  NA_real_
}

# ---- locating cluster files across layouts ----
find_cluster_files_in_out <- function(out_dir) {
  if (!dir.exists(out_dir)) return(character(0))
  list.files(
    out_dir,
    pattern = "clusters.*\\.tsv(\\.gz)?$",
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = TRUE
  )
}

find_out_dirs <- function(method_dir) {
  if (!dir.exists(method_dir)) return(character(0))
  method_dir <- normalizePath(method_dir, winslash = "/", mustWork = FALSE)
  
  children <- list.dirs(method_dir, recursive = FALSE, full.names = TRUE)
  run_roots <- unique(c(method_dir, children))
  run_roots <- run_roots[
    vapply(run_roots, function(p) {
      p == method_dir || grepl("^(job|array|run)_", basename(p), ignore.case = TRUE)
    }, logical(1))
  ]
  
  out_dirs <- character(0)
  
  for (rr in run_roots) {
    main_out <- file.path(rr, "out")
    if (dir.exists(main_out) && length(find_cluster_files_in_out(main_out)) > 0) {
      out_dirs <- c(out_dirs, main_out)
      next
    }
    
    task_dirs <- list.dirs(rr, recursive = FALSE, full.names = TRUE)
    task_dirs <- task_dirs[grepl("^task_[0-9]+$", basename(task_dirs), ignore.case = TRUE)]
    task_outs <- file.path(task_dirs, "out")
    task_outs <- task_outs[dir.exists(task_outs)]
    out_dirs <- c(out_dirs, task_outs)
  }
  
  out_dirs <- unique(normalizePath(out_dirs, winslash = "/", mustWork = FALSE))
  out_dirs[vapply(out_dirs, function(od) length(find_cluster_files_in_out(od)) > 0, logical(1))]
}

infer_outdir_percents_from_perf <- function(out_dir) {
  if (!dir.exists(out_dir)) return(numeric(0))
  
  cand <- character(0)
  if (file.exists(file.path(out_dir, "monet_perf_rows.tsv"))) {
    cand <- file.path(out_dir, "monet_perf_rows.tsv")
  } else {
    f <- list.files(out_dir, full.names = TRUE, recursive = FALSE)
    f <- f[grepl("\\.tsv$", f, ignore.case = TRUE)]
    f <- f[file_nonempty(f)]
    f <- f[grepl("perturbations|perf_row|perf_rows|performance", basename(f), ignore.case = TRUE)]
    if (length(f) > 0) cand <- f[1]
  }
  
  if (length(cand) == 0 || !file.exists(cand) || !file_nonempty(cand)) {
    pct_task <- infer_pct_from_task(out_dir)
    if (is.finite(pct_task)) return(unique(pct_task))
    return(numeric(0))
  }
  
  dt <- tryCatch(data.table::fread(cand), error = function(e) NULL)
  if (is.null(dt) || nrow(dt) == 0) {
    pct_task <- infer_pct_from_task(out_dir)
    if (is.finite(pct_task)) return(unique(pct_task))
    return(numeric(0))
  }
  
  nm <- names(dt)
  pct <- numeric(0)
  
  if ("Sample_Percent" %in% nm) {
    pct <- dt[["Sample_Percent"]]
  } else if ("Sample_Fraction" %in% nm) {
    pct <- dt[["Sample_Fraction"]]
  } else if ("Sample_Centile" %in% nm) {
    pct <- 100 * dt[["Sample_Centile"]]
  } else {
    pct_cols <- nm[grepl("percent|pct|fraction|centile", nm, ignore.case = TRUE)]
    if (length(pct_cols) >= 1) pct <- dt[[pct_cols[1]]]
  }
  
  pct <- vapply(pct, normalize_percent_value, numeric(1))
  pct <- pct[is.finite(pct)]
  if (length(pct)) return(unique(pct))
  
  pct_task <- infer_pct_from_task(out_dir)
  if (is.finite(pct_task)) return(unique(pct_task))
  
  numeric(0)
}

collect_cluster_files_for_percent <- function(method_dir, target_pct = 90, tol = 0.5) {
  out_dirs <- find_out_dirs(method_dir)
  if (length(out_dirs) == 0) return(character(0))
  
  hits <- character(0)
  
  for (od in out_dirs) {
    cf <- find_cluster_files_in_out(od)
    if (length(cf) == 0) next
    
    pct_from_name <- vapply(cf, infer_pct_from_filename, numeric(1))
    keep <- is.finite(pct_from_name) & abs(pct_from_name - target_pct) <= tol
    
    if (any(keep)) {
      hits <- c(hits, cf[keep])
      next
    }
    
    od_pcts <- infer_outdir_percents_from_perf(od)
    if (any(abs(od_pcts - target_pct) <= tol)) hits <- c(hits, cf)
  }
  
  unique(hits)
}

# ---- reading clusters + ARI ----
read_clusters_any <- function(path) {
  if (!is.character(path) || length(path) != 1L || !file.exists(path)) return(NULL)
  dt <- tryCatch(data.table::fread(path), error = function(e) NULL)
  if (is.null(dt) || nrow(dt) == 0) return(NULL)
  
  nm <- names(dt)
  
  if (!("Sample.ID" %in% nm)) {
    sid_cand <- intersect(nm, c("SampleID", "sample_id", "sample", "ID", "Id", "id"))
    if (length(sid_cand) >= 1) data.table::setnames(dt, sid_cand[1], "Sample.ID")
    else data.table::setnames(dt, nm[1], "Sample.ID")
  }
  
  if (!("Cluster" %in% names(dt))) {
    clcand <- intersect(names(dt), c(
      "Cluster_pred", "cluster", "ClusterPred", "consensuscluster", "ConsensusCluster",
      "ClusterLabel", "label", "Label", "pred", "Pred"
    ))
    if (length(clcand) >= 1) data.table::setnames(dt, clcand[1], "Cluster")
  }
  
  if (!all(c("Sample.ID", "Cluster") %in% names(dt))) return(NULL)
  
  dt[, .(
    Sample.ID = gsub("\\.", "-", as.character(Sample.ID)),
    Cluster   = suppressWarnings(as.integer(Cluster))
  )]
}

ari_between_clusters <- function(dt_a, dt_b, min_common = 3L) {
  if (is.null(dt_a) || is.null(dt_b)) return(list(ari = NA_real_, overlap = 0L))
  common <- intersect(dt_a$Sample.ID, dt_b$Sample.ID)
  ov <- length(common)
  if (ov < min_common) return(list(ari = NA_real_, overlap = ov))
  
  a <- dt_a[match(common, dt_a$Sample.ID), Cluster]
  b <- dt_b[match(common, dt_b$Sample.ID), Cluster]
  
  if (!requireNamespace("mclust", quietly = TRUE)) {
    stop("Package 'mclust' required. Install with install.packages('mclust').")
  }
  
  list(ari = mclust::adjustedRandIndex(a, b), overlap = ov)
}

pairwise_ari_table <- function(cluster_list, min_common = 3L) {
  k <- length(cluster_list)
  if (k < 2) return(data.table::data.table())
  
  res <- vector("list", k * (k - 1L) / 2L)
  idx <- 0L
  for (i in seq_len(k - 1L)) {
    for (j in (i + 1L):k) {
      tmp <- ari_between_clusters(cluster_list[[i]], cluster_list[[j]], min_common = min_common)
      idx <- idx + 1L
      res[[idx]] <- data.table::data.table(rep_i = i, rep_j = j, ARI = tmp$ari, Overlap = tmp$overlap)
    }
  }
  data.table::rbindlist(res, use.names = TRUE, fill = TRUE)
}

# ---- main ----
compute_sample_stability_90pct <- function(root = "Results/Performance_benchmarks",
                                           out_dir = "Results/Comparisons",
                                           pct = 90,
                                           min_common_sample = 20L,
                                           save = TRUE,
                                           dpi = 700,
                                           width_px = 3840,
                                           height_px = 2160) {
  if (!requireNamespace("data.table", quietly = TRUE)) stop("Install data.table")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Install ggplot2")
  
  dirs <- guess_benchmark_dirs(root)
  base_dir <- dirs$sample_dir
  if (!dir.exists(base_dir)) stop("Not found: ", base_dir)
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  method_dirs <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)
  method_dirs <- method_dirs[normalizePath(method_dirs, winslash = "/") != normalizePath(base_dir, winslash = "/")]
  
  pair_rows <- list()
  sum_rows  <- list()
  skip_rows <- list()
  idx <- 0L
  
  for (mdir in method_dirs) {
    alg  <- basename(mdir)
    catg <- get_method_category(alg)
    
    cfiles <- collect_cluster_files_for_percent(mdir, target_pct = pct, tol = 0.5)
    cfiles <- cfiles[file.exists(cfiles)]
    
    if (length(cfiles) < 2) {
      skip_rows[[length(skip_rows) + 1L]] <- data.table::data.table(
        Algorithm = alg, Category = catg,
        Reason = if (length(cfiles) == 0) "No cluster files matched target percent"
        else "Only one cluster file matched target percent"
      )
      next
    }
    
    cl <- lapply(cfiles, read_clusters_any)
    cl <- cl[!vapply(cl, is.null, logical(1))]
    if (length(cl) < 2) {
      skip_rows[[length(skip_rows) + 1L]] <- data.table::data.table(
        Algorithm = alg, Category = catg, Reason = "Cluster files unreadable or missing required columns"
      )
      next
    }
    
    pw <- pairwise_ari_table(cl, min_common = min_common_sample)
    pw <- pw[is.finite(ARI)]
    if (nrow(pw) == 0) {
      skip_rows[[length(skip_rows) + 1L]] <- data.table::data.table(
        Algorithm = alg, Category = catg, Reason = "No valid ARIs (likely overlap < min_common_sample)"
      )
      next
    }
    
    idx <- idx + 1L
    pair_rows[[idx]] <- cbind(data.table::data.table(Algorithm = alg, Category = catg, Subset_Percent = pct), pw)
    
    sum_rows[[idx]] <- data.table::data.table(
      Algorithm        = alg,
      Category         = catg,
      Subset_Percent   = pct,
      Stability_Median = stats::median(pw$ARI, na.rm = TRUE),
      Stability_Q25    = stats::quantile(pw$ARI, 0.25, na.rm = TRUE, names = FALSE),
      Stability_Q75    = stats::quantile(pw$ARI, 0.75, na.rm = TRUE, names = FALSE),
      N_clusterings    = length(cl),
      N_pairs          = nrow(pw),
      Overlap_Median   = stats::median(pw$Overlap, na.rm = TRUE),
      Overlap_Q25      = stats::quantile(pw$Overlap, 0.25, na.rm = TRUE, names = FALSE),
      Overlap_Q75      = stats::quantile(pw$Overlap, 0.75, na.rm = TRUE, names = FALSE)
    )
  }
  
  pairs <- data.table::rbindlist(pair_rows, use.names = TRUE, fill = TRUE)
  summ  <- data.table::rbindlist(sum_rows,  use.names = TRUE, fill = TRUE)
  skips <- data.table::rbindlist(skip_rows, use.names = TRUE, fill = TRUE)
  
  if (is.null(pairs) || ncol(pairs) == 0) {
    pairs <- data.table::data.table(
      Algorithm = character(), Category = character(), Subset_Percent = numeric(),
      rep_i = integer(), rep_j = integer(), ARI = numeric(), Overlap = integer()
    )
  }
  if (is.null(summ) || ncol(summ) == 0) {
    summ <- data.table::data.table(
      Algorithm = character(), Category = character(), Subset_Percent = numeric(),
      Stability_Median = numeric(), Stability_Q25 = numeric(), Stability_Q75 = numeric(),
      N_clusterings = integer(), N_pairs = integer(),
      Overlap_Median = numeric(), Overlap_Q25 = numeric(), Overlap_Q75 = numeric()
    )
  }
  if (is.null(skips) || ncol(skips) == 0) {
    skips <- data.table::data.table(
      Algorithm = character(), Category = character(), Reason = character()
    )
  }
  
  data.table::fwrite(
    if (nrow(pairs)) pairs[order(Category, Algorithm, rep_i, rep_j)] else pairs,
    file.path(out_dir, "stability_pairs_90pct.csv")
  )
  data.table::fwrite(
    if (nrow(summ)) summ[order(Category, Algorithm)] else summ,
    file.path(out_dir, "stability_summary_90pct.csv")
  )
  data.table::fwrite(
    if (nrow(skips)) skips[order(Category, Algorithm)] else skips,
    file.path(out_dir, "stability_skipped_90pct.csv")
  )
  
  if (nrow(pairs) == 0) {
    warning("No pairwise ARIs computed at ", pct, "%.")
    return(invisible(list(pairs = pairs, summary = summ, skipped = skips, plot = NULL)))
  }
  
  ord <- summ[order(Category, -Stability_Median)]$Algorithm
  pairs[, Algorithm := factor(Algorithm, levels = ord)]
  
  p <- ggplot2::ggplot(pairs, ggplot2::aes(x = Algorithm, y = ARI, fill = Category)) +
    ggplot2::geom_violin(trim = TRUE, alpha = 0.9, linewidth = 0.2) +
    ggplot2::geom_boxplot(width = 0.18, outlier.size = 0.4, linewidth = 0.25, alpha = 0.9) +
    ggplot2::scale_y_continuous(limits = c(-0.1, 1.1)) +
    ggplot2::labs(
      title    = sprintf("Replicate stability at %d%% sample perturbation", pct),
      # subtitle = sprintf("Pairwise ARI across replicate clusterings; min overlap = %d", min_common_sample),
      x = "Method", y = "Pairwise ARI"
    )
  
  if (exists("category_colors", inherits = TRUE)) {
    p <- p + ggplot2::scale_fill_manual(values = get("category_colors", inherits = TRUE))
  }
  
  if (exists("theme_benchmark", inherits = TRUE)) {
    p <- p + theme_benchmark(base_size = 10, legend = TRUE) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position = "right",
        legend.title = element_text(face = "bold")
      )
  } else {
    p <- p + ggplot2::theme_bw(base_size = 10) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position = "right",
        legend.title = element_text(face = "bold")
      )
  }
  
  if (isTRUE(save)) {
    ggplot2::ggsave(file.path(out_dir, "stability_90pct_violin.png"),
                    plot = p, dpi = dpi, width = width_px, height = height_px, units = "px")
    ggplot2::ggsave(file.path(out_dir, "stability_90pct_violin.pdf"),
                    plot = p, dpi = dpi, width = width_px, height = height_px, units = "px")
  }
  
  invisible(list(pairs = pairs, summary = summ, skipped = skips, plot = p))
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
