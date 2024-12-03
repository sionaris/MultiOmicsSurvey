# Load necessary libraries
library(ggplot2)
library(gridExtra)

# Read in the tune results
tune_results = readRDS("Resources/HPC output/tune_results_first_run.rds")

# Directory to save plots
prelim_output_dir <- "Results/single_algorithm/iClusterBayes/preliminary_first_run"
dir.create(prelim_output_dir, recursive = TRUE, showWarnings = FALSE)

# Simulated K values and modalities
K_values <- 1:9
modalities <- c("SNPs", "RNAseq", "CNV", "miRNA", "Methylation")
modality_colours <- c("red", "blue", "green", "purple", "orange")

# Extract Z.ar, beta.ar, and gamma.ar data
Z_ar_data <- lapply(K_values, function(K) tune_results$fit[[K]]$Z.ar)
beta_ar_data <- lapply(K_values, function(K) tune_results$fit[[K]]$beta.ar)
gamma_ar_data <- lapply(K_values, function(K) tune_results$fit[[K]]$gamma.ar)

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
    data.frame(
      value = unlist(data_list[[i]]),  # Flatten the list of modality-specific values
      modality = rep(modalities, sapply(data_list[[i]], length)),  # Assign modalities dynamically
      K = K_values[i] + 1,
      variable = var_name
    )
  }))
}

# Prepare Z.ar data
Z_ar_df <- prepare_global_data(Z_ar_data, K_values, "Z.ar")

# Prepare beta.ar and gamma.ar data for modality-specific plots
beta_ar_df <- prepare_modality_data(beta_ar_data, K_values, "beta.ar")
gamma_ar_df <- prepare_modality_data(gamma_ar_data, K_values, "gamma.ar")

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

# Function to compute and print summary statistics
print_summary_statistics <- function(fit_list, K_values) {
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
}

# Call the function to print summary statistics
print_summary_statistics(tune_results$fit, K_values)

