# Load necessary libraries
library(ggplot2)
library(gridExtra)
library(doParallel)
library(foreach)
library(rmarkdown)

# Define useful functions

# Helper function for Z.ar data preparation
prepare_global_data <- function(data_list, K_values, var_name) {
  data <- do.call(rbind, lapply(seq_along(data_list), function(i) {
    data.frame(value = data_list[[i]], K = K_values[i] + 1)
  }))
  data$K <- as.factor(data$K)
  data$variable <- var_name
  return(data)
}

# Helper function for modality-specific data preparation
prepare_modality_data <- function(data_list, K_values, var_name) {
  do.call(rbind, lapply(seq_along(data_list), function(i) {
    df <- data.frame(
      value = unlist(data_list[[i]]),
      modality = rep(modalities, sapply(data_list[[i]], length)),
      K = K_values[i] + 1,
      variable = var_name
    )
    return(df)
  }))
}

# Modified function to capture summary statistics as text instead of printing
get_summary_statistics_text <- function(fit_list, K_values) {
  # Use capture.output to store the print output as a character vector
  summary_text <- capture.output({
    for (K in K_values) {
      cat("\n--- Summary for K =", K, "---\n")
      
      # Extract Z.ar
      if (!is.null(fit_list[[K]]$Z.ar)) {
        Z_ar_stats <- summary(fit_list[[K]]$Z.ar)
        cat("Z.ar (Location: fit[[", K, "]]$Z.ar):\n", sep = "")
        print(Z_ar_stats)
      }
      
      # Extract beta.ar (for all modalities)
      if (!is.null(fit_list[[K]]$beta.ar)) {
        cat("\nbeta.ar (Location: fit[[", K, "]]$beta.ar):\n", sep = "")
        lapply(seq_along(fit_list[[K]]$beta.ar), function(modality) {
          stats <- summary(fit_list[[K]]$beta.ar[[modality]])
          cat("  Modality", modality, ":\n")
          print(stats)
        })
      }
      
      # Extract gamma.ar (for all modalities)
      if (!is.null(fit_list[[K]]$gamma.ar)) {
        cat("\ngamma.ar (Location: fit[[", K, "]]$gamma.ar):\n", sep = "")
        lapply(seq_along(fit_list[[K]]$gamma.ar), function(modality) {
          stats <- summary(fit_list[[K]]$gamma.ar[[modality]])
          cat("  Modality", modality, ":\n")
          print(stats)
        })
      }
    }
  })
  
  return(summary_text)
}

# Function to create a global density plot (Z.ar)
create_global_plot <- function(data, title) {
  ggplot(data, aes(x = value, color = K, fill = K)) +
    geom_density(alpha = 0.3, linewidth = 0.2) +
    labs(title = title, x = "Value", y = "Density", color = "K", fill = "K") +
    theme(
      plot.title = element_text(hjust = 0.5, size = 5.5, face = "bold"),
      axis.title = element_text(face = "bold", size = 4.5),
      axis.text = element_text(size = 3.5),
      axis.line = element_line(linewidth = 0.2),
      axis.ticks = element_line(linewidth = 0.1),
      legend.title = element_text(size = 4, face = "bold"),
      legend.text = element_text(size = 3.5),
      legend.key.size = unit(0.3, "lines"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_rect(fill = "white")
    )
}

# Function to create modality-specific density plots for beta.ar and gamma.ar
create_modality_plot <- function(data, K, title) {
  ggplot(data[data$K == K+1, ], aes(x = value, color = modality, fill = modality)) +
    geom_density(alpha = 0.3, linewidth = 0.2) +
    scale_color_manual(values = modality_colours) +
    scale_fill_manual(values = modality_colours) +
    labs(title = paste(title, "- K =", K+1), x = "Value", y = "Density", color = "Modality", fill = "Modality") +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, size = 5.5, face = "bold"),
      axis.title = element_text(face = "bold", size = 4.5),
      axis.text = element_text(size = 3.5),
      axis.line = element_line(linewidth = 0.2),
      axis.ticks = element_line(linewidth = 0.1),
      legend.title = element_text(size = 4, face = "bold"),
      legend.text = element_text(size = 3.5),
      legend.key.size = unit(0.3, "lines"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_rect(fill = "white")
    )
}


# Create list of output directories
dirs = c("sdev_0.01_beta_0.3",
         "sdev_0.01_beta_0.5",
         "sdev_0.01_beta_0.8",
         "sdev_0.05_beta_0.3",
         "sdev_0.05_beta_0.5",
         "sdev_0.05_beta_0.8",
         "sdev_0.025_beta_0.3",
         "sdev_0.025_beta_0.5",
         "sdev_0.025_beta_0.8")

# Set up parallel backend
num_cores <- 4 # adjust this based on the number of cores available
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Parallelize the loop using foreach
foreach(i = seq_along(dirs)) %dopar% {
  library(ggplot2)
  library(gridExtra)
  library(rmarkdown) # load inside foreach
  
  sdev = as.numeric(strsplit(dirs[i], split = "_")[[1]][2])
  beta.var.scale = as.numeric(strsplit(dirs[i], split = "_")[[1]][4])
  
  if (!dir.exists(paste0("Results/single_algorithm/iClusterBayes/", dirs[i]))) {
    dir.create(paste0("Results/single_algorithm/iClusterBayes/", dirs[i]), recursive = TRUE)
  }
  
  # Read in the tune results
  tune_results = readRDS(paste0("Resources/HPC output/tune_results_", dirs[i], ".rds"))
  
  # Directory to save plots and Rmd
  prelim_output_dir <- paste0("Results/single_algorithm/iClusterBayes/", dirs[i])
  dir.create(prelim_output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Simulated K values and modalities
  K_values <- 1:9
  modalities <- c("SNPs", "RNAseq", "CNV", "miRNA", "Methylation")
  modality_colours <- c("red", "blue", "green", "purple", "orange")
  
  hyperparameter_values <- list(
    n.burnin = 1000,
    n.draw = 1200,
    thin = 1,
    pp_cutoff = 0.5,
    prior_gamma = rep(0.1, length(modalities)),
    K_values = paste(K_values, "/ Examined number of clusters:", K_values + 1),
    sdev = sdev,
    beta.var.scale = beta.var.scale
  )
  
  # Extract Z.ar, beta.ar, and gamma.ar data
  Z_ar_data <- lapply(K_values, function(K) tune_results$fit[[K]]$Z.ar)
  beta_ar_data <- lapply(K_values, function(K) tune_results$fit[[K]]$beta.ar)
  gamma_ar_data <- lapply(K_values, function(K) tune_results$fit[[K]]$gamma.ar)
  
  # Prepare Z.ar data
  Z_ar_df <- prepare_global_data(Z_ar_data, K_values, "Z.ar")
  
  # Prepare beta.ar and gamma.ar data for modality-specific plots
  beta_ar_df <- prepare_modality_data(beta_ar_data, K_values, "beta.ar")
  gamma_ar_df <- prepare_modality_data(gamma_ar_data, K_values, "gamma.ar")
  
  # Generate the Z.ar plot
  plot_Z_ar <- create_global_plot(Z_ar_df, "Overlaid Distributions of Z.ar")
  
  # Generate the beta.ar and gamma.ar plots for each K
  beta_ar_plots <- lapply(K_values, function(K) {
    create_modality_plot(beta_ar_df, K, "beta.ar Distributions")
  })
  gamma_ar_plots <- lapply(K_values, function(K) {
    create_modality_plot(gamma_ar_df, K, "gamma.ar Distributions")
  })
  
  # Save Z.ar plot
  ggsave(
    filename = file.path(prelim_output_dir, "Z_ar_plot.png"),
    plot = plot_Z_ar,
    dpi = 700,
    width = 1920,
    height = 1080,
    units = "px"
  )
  
  # Save beta.ar and gamma.ar plots for each K
  for (K in K_values) {
    beta_plot_file <- file.path(prelim_output_dir, paste0("beta_ar_K_", K+1, ".png"))
    gamma_plot_file <- file.path(prelim_output_dir, paste0("gamma_ar_K_", K+1, ".png"))
    
    # Save beta.ar plot
    ggsave(
      filename = beta_plot_file,
      plot = beta_ar_plots[[K]],
      dpi = 700,
      width = 1920,
      height = 1080,
      units = "px"
    )
    
    # Save gamma.ar plot
    ggsave(
      filename = gamma_plot_file,
      plot = gamma_ar_plots[[K]],
      dpi = 700,
      width = 1920,
      height = 1080,
      units = "px"
    )
  }
  
  # Capture summary statistics
  summary_stats_text <- get_summary_statistics_text(tune_results$fit, K_values)
  summary_stats_str <- paste(summary_stats_text, collapse = "\n")
  
  # Create an R Markdown file for the PDF report
  rmd_file <- file.path(prelim_output_dir, "summary_report.Rmd")
  
  # Create Rmd content
  rmd_content <- c(
    "---",
    "title: \"iClusterBayes Summary Report\"",
    "author: \"Automated Report\"",
    "output: pdf_document",
    "---",
    "",
    "# Hyperparameters",
    "```{r, echo=FALSE}",
    "hyperparameter_values",
    "```",
    "",
    "# Summary Statistics",
    "```{r, echo=FALSE, results='asis'}",
    "cat(paste0(\"```\n\", summary_stats_str, \"\n```\"))",
    "```"
  )
  
  # Write the Rmd file
  writeLines(rmd_content, con = rmd_file)
  
  # Render the PDF report
  rmarkdown::render(
    input = rmd_file,
    output_file = file.path("summary_report.pdf"),
    quiet = TRUE
  )
}

# Stop the cluster
stopCluster(cl)
