################################################################################
# iClusterBayes Parallel Processing with Additional Robust Metrics
#
# This script:
#   1) Iterates over combinations of hyperparameters (sdev, beta.var.scale).
#   2) For each combination, it reads the corresponding .rds results file,
#      generates plots, and produces two PDF reports:
#         - summary_report.pdf   (original summary information)
#         - robust_report.pdf    (new robust metrics and final score)
#   3) Combines all summary_report.pdf into a single PDF and similarly for all
#      robust_report.pdf into another single PDF.
#   4) Stores robust metrics in a list object for programmatic inspection.
################################################################################

library(ggplot2)
library(gridExtra)
library(doParallel)
library(foreach)
library(rmarkdown)
library(pdftools)
library(jsonlite)

################################################################################
# Helper Functions
################################################################################

# Prepares a global data frame for Z.ar across all K values, for plotting
prepare_global_data <- function(data_list, K_values, var_name) {
  df <- do.call(rbind, lapply(seq_along(data_list), function(i) {
    data.frame(value = data_list[[i]], K = K_values[i] + 1)
  }))
  df$K <- as.factor(df$K)
  df$variable <- var_name
  return(df)
}

# Prepares a modality-specific data frame for Beta/Gamma acceptance across K
prepare_modality_data <- function(data_list, K_values, var_name) {
  do.call(rbind, lapply(seq_along(data_list), function(i) {
    data.frame(
      value    = unlist(data_list[[i]]),
      modality = rep(modalities, sapply(data_list[[i]], length)),
      K        = K_values[i] + 1,
      variable = var_name
    )
  }))
}

# Captures the standard summary statistics as text output
get_summary_statistics_text <- function(fit_list, K_values) {
  summary_text <- capture.output({
    for (K in K_values) {
      cat("\n--- Summary for K =", K, "---\n")
      
      if (!is.null(fit_list[[K]]$Z.ar)) {
        cat("Z.ar:\n")
        print(summary(fit_list[[K]]$Z.ar))
      }
      
      if (!is.null(fit_list[[K]]$beta.ar)) {
        cat("\nbeta.ar:\n")
        for (m in seq_along(fit_list[[K]]$beta.ar)) {
          cat("  Modality", m, ":\n")
          print(summary(fit_list[[K]]$beta.ar[[m]]))
        }
      }
      
      if (!is.null(fit_list[[K]]$gamma.ar)) {
        cat("\ngamma.ar:\n")
        for (m in seq_along(fit_list[[K]]$gamma.ar)) {
          cat("  Modality", m, ":\n")
          print(summary(fit_list[[K]]$gamma.ar[[m]]))
        }
      }
    }
  })
  return(summary_text)
}

# Creates a global density plot for Z.ar
create_global_plot <- function(data, title) {
  ggplot(data, aes(x = value, color = K, fill = K)) +
    geom_density(alpha = 0.3, linewidth = 0.2) +
    labs(title = title, x = "Value", y = "Density") +
    theme_minimal(base_size = 7) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "bottom"
    )
}

# Creates modality-specific density plots for Beta/Gamma acceptance
create_modality_plot <- function(data, K, title) {
  ggplot(data[data$K == K+1, ], aes(x = value, color = modality, fill = modality)) +
    geom_density(alpha = 0.3, linewidth = 0.2) +
    labs(title = paste(title, "- K =", K+1), x = "Value", y = "Density") +
    scale_color_manual(values = modality_colours) +
    scale_fill_manual(values = modality_colours) +
    theme_minimal(base_size = 7) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "bottom"
    )
}

################################################################################
# New: Compute Additional Robust Metrics and Score
################################################################################

compute_robust_metrics_and_score <- function(tune_results, K_values) {
  metrics_list <- list()
  final_score  <- 0
  
  for (K in K_values) {
    z_ar     <- tune_results$fit[[K]]$Z.ar
    beta_ar  <- tune_results$fit[[K]]$beta.ar
    gamma_ar <- tune_results$fit[[K]]$gamma.ar
    
    if (is.null(z_ar)) next
    
    median_z   <- median(z_ar, na.rm = TRUE)
    all_beta   <- unlist(beta_ar)
    median_beta <- median(all_beta, na.rm = TRUE)
    all_gamma  <- unlist(gamma_ar)
    median_gamma <- median(all_gamma, na.rm = TRUE)
    
    # Example scoring: sum of distances from 0.45 for each median
    # Lower is "better" in this example
    score_k <- (abs(median_z - 0.45) +
                  abs(median_beta - 0.45) +
                  abs(median_gamma - 0.45))
    
    metrics_list[[K]] <- list(
      K               = K,
      median_z_ar     = median_z,
      median_beta_ar  = median_beta,
      median_gamma_ar = median_gamma,
      score_k         = score_k
    )
    final_score <- final_score + score_k
  }
  
  return(list(
    per_K_metrics = metrics_list,
    final_score   = final_score
  ))
}

################################################################################
# Define the Parameter Grid (Real Set)
################################################################################

sdev_values <- c(0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.05)
beta_values <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.8, 1.0)

dirs <- c()
for (sd in sdev_values) {
  for (bv in beta_values) {
    dir_str <- paste0("sdev_", sd, "_beta_", bv)
    dirs <- c(dirs, dir_str)
  }
}

################################################################################
# Parallel Execution Setup
################################################################################

num_cores <- 4  # Adjust as appropriate
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Vectors to store paths for summary and robust PDFs
pdf_paths_summary <- vector("character", length(dirs))
pdf_paths_robust  <- vector("character", length(dirs))

# Will store robust scores for each combo
robust_scores_list <- list()

################################################################################
# Main Loop
################################################################################

results_output <- foreach(i = seq_along(dirs), .combine = "c") %dopar% {
  library(ggplot2)
  library(gridExtra)
  library(rmarkdown)
  library(jsonlite)
  
  tokens <- strsplit(dirs[i], "_")[[1]]
  sdev_val <- as.numeric(tokens[2])
  beta_val <- as.numeric(tokens[4])
  
  results_dir <- file.path("Results/single_algorithm/iClusterBayes", dirs[i])
  if (!dir.exists(results_dir)) {
    dir.create(results_dir, recursive = TRUE)
  }
  
  # Attempt to read the .rds file
  rds_file <- file.path("Resources/HPC output/iCB_HPC", paste0("tune_results_", dirs[i], ".rds"))
  if (!file.exists(rds_file)) {
    return(paste0(dirs[i], "|NA|NA|", toJSON(list())))
  }
  
  tune_results <- readRDS(rds_file)
  
  K_values <- 1:9
  modalities <- c("SNPs", "RNAseq", "CNV", "miRNA", "Methylation")
  modality_colours <- c("red", "blue", "green", "purple", "orange")
  
  # Extract acceptance arrays
  Z_ar_data    <- lapply(K_values, function(K) tune_results$fit[[K]]$Z.ar)
  beta_ar_data <- lapply(K_values, function(K) tune_results$fit[[K]]$beta.ar)
  gamma_ar_data<- lapply(K_values, function(K) tune_results$fit[[K]]$gamma.ar)
  
  # Prepare data frames
  Z_ar_df     <- prepare_global_data(Z_ar_data, K_values, "Z.ar")
  beta_ar_df  <- prepare_modality_data(beta_ar_data, K_values, "beta.ar")
  gamma_ar_df <- prepare_modality_data(gamma_ar_data, K_values, "gamma.ar")
  
  # Generate plots
  plot_Z_ar <- create_global_plot(Z_ar_df, "Overlaid Distributions of Z.ar")
  beta_ar_plots <- lapply(K_values, function(K) create_modality_plot(beta_ar_df, K, "beta.ar Distributions"))
  gamma_ar_plots<- lapply(K_values, function(K) create_modality_plot(gamma_ar_df, K, "gamma.ar Distributions"))
  
  # Save plots (individual PDF directory)
  ggsave(
    filename = file.path(getwd(), results_dir, "Z_ar_plot.png"),
    plot = plot_Z_ar,
    dpi = 400, width = 6, height = 4
  )
  for (k in K_values) {
    ggsave(
      filename = file.path(getwd(), results_dir, paste0("beta_ar_K_", k+1, ".png")),
      plot = beta_ar_plots[[k]],
      dpi = 400, width = 6, height = 4
    )
    ggsave(
      filename = file.path(getwd(), results_dir, paste0("gamma_ar_K_", k+1, ".png")),
      plot = gamma_ar_plots[[k]],
      dpi = 400, width = 6, height = 4
    )
  }
  
  # Capture standard summary stats
  summary_stats_text <- get_summary_statistics_text(tune_results$fit, K_values)
  summary_stats_str  <- paste(summary_stats_text, collapse = "\n")
  
  # Produce summary_report
  rmd_file_summary <- file.path(getwd(), results_dir, "summary_report.Rmd")
  rmd_content_summary <- c(
    "---",
    "title: \"iClusterBayes Summary Report\"",
    "author: \"Automated Report\"",
    "output: pdf_document",
    "---",
    "",
    "# Hyperparameters",
    "```{r}",
    paste0("cat('sdev = ", sdev_val, ", beta.var.scale = ", beta_val, "\\n')"),
    "```",
    "",
    "# Key Summary Statistics",
    "```{r}",
    "cat(paste0('```\\n', '", summary_stats_str, "', '\\n```'))",
    "```",
    "",
    "# Conclusion",
    "This brief report includes acceptance rate summaries, posterior probabilities, ",
    "and basic statistics for each K (1–9). The `.png` files show Z.ar, beta.ar, ",
    "and gamma.ar distributions."
  )
  
  writeLines(rmd_content_summary, con = rmd_file_summary)
  pdf_output_path_summary <- file.path(getwd(), results_dir, "summary_report.pdf")
  
  rmarkdown::render(
    input       = rmd_file_summary,
    output_file = pdf_output_path_summary,
    quiet       = TRUE
  )
  
  # Compute robust metrics
  robust_metrics <- compute_robust_metrics_and_score(tune_results, K_values)
  
  # Produce robust_report
  robust_text <- capture.output({
    cat("## Per-K Metrics\n")
    for (k in K_values) {
      if (!is.null(robust_metrics$per_K_metrics[[k]])) {
        mk <- robust_metrics$per_K_metrics[[k]]
        cat("\n--- K =", k, "---\n")
        cat("Median Z.ar:       ", mk$median_z_ar,    "\n")
        cat("Median Beta.ar:    ", mk$median_beta_ar, "\n")
        cat("Median Gamma.ar:   ", mk$median_gamma_ar,"\n")
        cat("Score segment:     ", mk$score_k,        "\n")
      }
    }
    cat("\n## Final Score\n")
    cat("Final Score:", robust_metrics$final_score, "\n")
  })
  robust_text_str <- paste(robust_text, collapse = "\n")
  
  rmd_file_robust <- file.path(getwd(), results_dir, "robust_report.Rmd")
  rmd_content_robust <- c(
    "---",
    "title: \"Robust Metrics Report\"",
    "author: \"Automated Report\"",
    "output: pdf_document",
    "---",
    "",
    "# Hyperparameters",
    "```{r}",
    paste0("cat('sdev = ", sdev_val, ", beta.var.scale = ", beta_val, "\\n')"),
    "```",
    "",
    "# Robust Metrics & Score",
    "```{r}",
    "cat(paste0('```\\n', '", robust_text_str, "', '\\n```'))",
    "```"
  )
  writeLines(rmd_content_robust, con = rmd_file_robust)
  pdf_output_path_robust <- file.path(getwd(), results_dir, "robust_report.pdf")
  
  rmarkdown::render(
    input       = rmd_file_robust,
    output_file = pdf_output_path_robust,
    quiet       = TRUE
  )
  
  # Convert robust metrics to JSON
  robust_obj <- list(
    combo_name = dirs[i],
    sdev_val   = sdev_val,
    beta_val   = beta_val,
    metrics    = robust_metrics
  )
  robust_json_str <- toJSON(robust_obj, auto_unbox = TRUE)
  
  # Return a combined string so we can parse later
  paste0(dirs[i], "|",
         pdf_output_path_summary, "|",
         pdf_output_path_robust, "|",
         robust_json_str)
}

stopCluster(cl)

################################################################################
# Post-processing: parse the returned data for each combo
################################################################################

valid_lines <- grep("\\|", results_output, value = TRUE)
pdf_summary <- character(length(valid_lines))
pdf_robust  <- character(length(valid_lines))

index <- 1
for (line in valid_lines) {
  parts <- strsplit(line, "\\|")[[1]]
  if (length(parts) < 4) next
  
  combo_name    <- parts[1]
  summary_pdf   <- parts[2]
  robust_pdf    <- parts[3]
  robust_json   <- parts[4]
  
  pdf_summary[index] <- summary_pdf
  pdf_robust[index]  <- robust_pdf
  
  robust_data <- fromJSON(robust_json)
  robust_scores_list[[combo_name]] <- robust_data
  
  index <- index + 1
}

pdf_summary <- pdf_summary[file.exists(pdf_summary)]
pdf_robust  <- pdf_robust[file.exists(pdf_robust)]

# Create final combined PDF for summary
if (length(pdf_summary) > 1) {
  combined_summary_path <- file.path(
    getwd(), "Results", "single_algorithm", "iClusterBayes",
    "Global_Combined_Report.pdf"
  )
  pdf_combine(pdf_summary, output = combined_summary_path)
  message("Global combined summary PDF created at: ", combined_summary_path)
} else {
  message("No or only one summary PDF found. Skipping global summary concatenation.")
}

# Create final combined PDF for robust
if (length(pdf_robust) > 1) {
  combined_robust_path <- file.path(
    getwd(), "Results", "single_algorithm", "iClusterBayes",
    "Global_Robust_Report.pdf"
  )
  pdf_combine(pdf_robust, output = combined_robust_path)
  message("Global combined robust PDF created at: ", combined_robust_path)
} else {
  message("No or only one robust PDF found. Skipping global robust concatenation.")
}

################################################################################
# Example of printing or using the final robust_scores_list
################################################################################

cat("\n--- Final Scores for Each Combo ---\n")
for (nm in names(robust_scores_list)) {
  fs <- robust_scores_list[[nm]]$metrics$final_score
  cat(nm, " => final_score =", fs, "\n")
}

