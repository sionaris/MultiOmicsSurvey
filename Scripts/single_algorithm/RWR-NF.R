# Import data from gitignored "Resources/BRCA complete/" folder #####

# This RDS object was produced using the Scripts/MOVICS/MOVICS_baseline.R script
input = readRDS("Resources/TCGA/mm_input.rds")

# Setup environment variables for markdown #####

# Ensure reproducibility
RNGversion("4.2.2")
set.seed(123)

# Load custom helper functions
source("Scripts/automated_scripts/custom_functions.R")
source("Scripts/automated_scripts/modified_MOVICS_functions.R")

# Preamble
home = getwd()
algorithm = "RWR-NF"
alg_feature_pref = "cols" # Where does the algorithm expect the features to be
citation = fetch_citation(algorithm = algorithm)
data_source = "TCGA" # e.g. TCGA, TCGA-transNEO, transNEO-PARTNER
data_types = "RNAseq-CNV-Methylation-miRNA-SNPs" # e.g. RNAseq, RNAseq-CNV-miRNA
evaluation_source = "transNEO" # e.g. PARTNER, transNEO-PARTNER 
title = paste0("Results from ", algorithm)
subtitle = paste0("<b>Train</b>: ", data_source, " ", data_types, 
                  " | <b>Evaluation</b>: ", evaluation_source)
in_a_nutshell = fetch_in_a_nutshell(algorithm = algorithm)
ground_truth_labels = openxlsx::read.xlsx("Results/MOVICS_baseline/MOVICS_TCGA_RNAseq-CNV-Methylation-miRNA-SNPs_eval_on_transNEO_clusterings.xlsx")
ground_truth_k = 2 # optk from MOVICS
optk_boolean = "TRUE" # either TRUE or FALSE. Answers whether the algorithm suggests an optimal k
optk_text = ifelse(optk_boolean == TRUE,
                   "<u>suggests</u> an estimate of the optimal number of multi-omic clusters $k$",
                   "<u>does not suggest</u> an optimal number of multi-omic clusters $k$")

# Detailed description of the algorithm
description = paste(readLines(paste0("Resources/algorithm_descriptions/", algorithm,
                                     "_description.Rmd")),
                    collapse = "\n") # File path to .Rmd file within Resources/algorithm_descriptions

# Create algorithm directory if it doesn't exist
if (!dir.exists(paste0(home, "/Results/single_algorithm/", 
                       algorithm))) {
  dir.create(paste0(home, "/Results/single_algorithm/", 
                    algorithm))
}

# Preprocessing flags and code #####
library(stringr)
library(dplyr)

# All preprocessing for this input has already been performed using the 
# Scripts/MOVICS/MOVICS_baseline.R script

# Extract the names of the modalities that will be used
modalities = unlist(strsplit(data_types, "-"))

# Replace with TRUE wherever features are in rows
features_in_rows = rep(TRUE, length(modalities))

# Run algorithm #####

# Check what kind of arrangement the algorithm requires as input 
# (i.e. features in rows or columns)?
if (alg_feature_pref == "rows") {
  rogue_indices = which(features_in_rows == FALSE)
  if (length(rogue_indices >= 1)) {
    for (index in rogue_indices) {
      cols = colnames(input[[index]])
      rows = rownames(input[[index]])
      input[[index]] = t(input[[index]])
      rownames(input[[index]]) = cols # transpose names
      colnames(input[[index]]) = rows # transpose names
      rm(rows, cols)
    }
  }
} else if (alg_feature_pref == "cols") {
  rogue_indices = which(features_in_rows == TRUE)
  if (length(rogue_indices >= 1)) {
    for (index in rogue_indices) {
      cols = colnames(input[[index]])
      rows = rownames(input[[index]])
      input[[index]] = t(input[[index]])
      rownames(input[[index]]) = cols # transpose names
      colnames(input[[index]]) = rows # transpose names
      rm(rows, cols)
    }
  }
}
rm(rogue_indices, index); gc()

# Import clinical data for the TCGA samples of interest
clinical_data = openxlsx::read.xlsx("Resources/TCGA/clinical_data.xlsx")

# Source RWR-F code
source("Scripts/single_algorithm/RWR-F_HPC/RWR-F_source.R")
library(SNFtool)
library(clValid)       # Provides the 'dunn' index
library(kernlab)       # For spectral clustering (specc)
library(foreach)
library(doParallel)

# Hyperparameters and other parameters ###
k_range = 2:10
neighbor_step = 5
num_neighbors_range = seq(10, 50, neighbor_step)
sigma_step = 0.1
sigma_range = seq(0.3, 0.8, sigma_step)

# RWR-F hyperparameters
RWRF_gamma_fixed <- 0.7
RWRNF_gamma_fixed <- 0.7
RWRNF_num_neighbors_fixed <- 10
RWRNF_alpha_fixed <- 0.9
RWRNF_beta_fixed <- 0.9
RWR_iteration_max <- 1000

# Create distances list
# Modality types
continuous = c("RNAseq", "CNV", "Methylation", "miRNA")
categorical = c("SNPs")

# Calculate the pair-wise distance (Euclidean for continuous modalities)
distL = lapply(input[continuous], function(x) {
  x = as.matrix(x)
  x = dist2(x, x)
})

# Binary for SNPs (see ?dist for details)
distL[["SNPs"]] = as.matrix(dist(as.matrix(input$SNPs),
                                 as.matrix(input$SNPs),
                                 method = "binary"))
gc()

# Run the Scripts/single_algorithm/RWR-F_HPC/*.R and *.sh scripts at an HPC
saveRDS(distL, "Resources/RWR-F_input_dists.rds")

# Job scripts were generated with the Scripts/automated_scripts/RWR-F_scripts_generation.R
# Fusions results from the HPC are imported here
Fusions = list()
results_indices = grep(".rds", list.files("Resources/HPC output/RWR-F_HPC/"))

for (filename in list.files("Resources/HPC output/RWR-F_HPC/")[results_indices]) {
  NN_val = strsplit(filename, "_")[[1]][2]
  sigma_val = strsplit(filename, "_")[[1]][4]
  Fusions[[paste0("NN = ", NN_val, ", sigma = ", sigma_val)]] = readRDS(paste0("Resources/HPC output/RWR-F_HPC/",
                                                                               filename))
  Fusions[[paste0("NN = ", NN_val, ", sigma = ", sigma_val)]] = Fusions[[paste0("NN = ", NN_val, ", sigma = ", sigma_val)]]$fused_rwrnf
}
rm(NN_val, sigma_val); gc()

# Function to compute both Frobenius norm and Pearson correlation between matrices
compute_matrix_similarity <- function(matrices) {
  num_matrices <- length(matrices)
  similarity_frobenius <- matrix(0, nrow = num_matrices, ncol = num_matrices)
  similarity_pearson <- matrix(0, nrow = num_matrices, ncol = num_matrices)
  
  for (i in 1:num_matrices) {
    for (j in 1:num_matrices) {
      if (i != j) {
        similarity_frobenius[i, j] <- frobenius_norm(matrices[[i]], matrices[[j]])
        similarity_pearson[i, j] <- pearson_correlation(matrices[[i]], matrices[[j]])
      }
    }
  }
  
  # Set row names and column names
  
  rownames(similarity_frobenius) <- colnames(similarity_frobenius) <- 
    rownames(similarity_pearson) <- colnames(similarity_pearson) <- names(matrices)
  
  return(list(Frobenius = similarity_frobenius, Pearson = similarity_pearson))
}

# Check similarities for a given nn
nn_similarities = list()
for (nn in num_neighbors_range) {
  indices = grepl(paste0("NN = ", nn), names(Fusions))
  matrices = Fusions[indices]
  similarity_results <- compute_matrix_similarity(matrices)
  
  # Modify row names and column names for the similarity matrices
  matrix_names <- substr(names(matrices), 9, 20)
  rownames(similarity_results$Frobenius) <- matrix_names
  colnames(similarity_results$Frobenius) <- matrix_names
  rownames(similarity_results$Pearson) <- matrix_names
  colnames(similarity_results$Pearson) <- matrix_names
  
  nn_similarities[[paste0("NN = ", nn)]] <- similarity_results
}

print(nn_similarities)

# Check similarities for a given sigma
sigma_similarities = list()
for (sigma in sigma_range) {
  indices = grepl(paste0("sigma = ", sigma), names(Fusions))
  matrices = Fusions[indices]
  similarity_results <- compute_matrix_similarity(matrices)
  
  # Modify row names and column names for the similarity matrices
  matrix_names <- substr(names(matrices), 0, 7)
  rownames(similarity_results$Frobenius) <- matrix_names
  colnames(similarity_results$Frobenius) <- matrix_names
  rownames(similarity_results$Pearson) <- matrix_names
  colnames(similarity_results$Pearson) <- matrix_names
  
  sigma_similarities[[paste0("sigma = ", sigma)]] <- similarity_results
}

print(sigma_similarities)
rm(indices, matrices); gc()

# All similarities
all_similarities = compute_matrix_similarity(Fusions)

# Overall tests ###
# Get summary statistics
# Initialize nn_summary and sigma_summary
nn_summary <- data.frame(
  nn = integer(),        
  mean = numeric(),      
  median = numeric(),    
  sd = numeric()         
)

sigma_summary <- data.frame(
  sigma = numeric(),     
  mean = numeric(),      
  median = numeric(),   
  sd = numeric()        
)

# Calculate nn_summary
for (nn in names(nn_similarities)) {
  pearson_matrix <- nn_similarities[[nn]]$Pearson
  pearson_values <- pearson_matrix[lower.tri(pearson_matrix, diag = FALSE)]
  
  # Compute mean, median, and standard deviation
  mean_value <- mean(pearson_values)
  median_value <- median(pearson_values)
  sd_value <- sd(pearson_values)
  
  # Append the results to nn_summary
  nn_summary <- rbind(nn_summary, data.frame(
    nn = as.numeric(sub("NN = ", "", nn)),
    mean = mean_value,
    median = median_value,
    sd = sd_value
  ))
}

# Calculate sigma_summary
for (sigma in names(sigma_similarities)) {
  pearson_matrix <- sigma_similarities[[sigma]]$Pearson
  pearson_values <- pearson_matrix[lower.tri(pearson_matrix, diag = FALSE)]
  
  # Compute mean, median, and standard deviation
  mean_value <- mean(pearson_values)
  median_value <- median(pearson_values)
  sd_value <- sd(pearson_values)
  
  # Append the results to sigma_summary
  sigma_summary <- rbind(sigma_summary, data.frame(
    sigma = as.numeric(sub("sigma = ", "", sigma)),
    mean = mean_value,
    median = median_value,
    sd = sd_value
  ))
}

# Display the summaries
print(nn_summary)
print(sigma_summary)

# Parametric ##
# Reshape data for ANOVA
nn_reshape <- reshape_Pearson_for_tests(nn_similarities)
sigma_reshape <- reshape_Pearson_for_tests(sigma_similarities)

# Perform ANOVA for nn
anova_nn <- aov(Value ~ Factor, data = nn_reshape)
summary(anova_nn)

# Perform ANOVA for sigma
anova_sigma <- aov(Value ~ Factor, data = sigma_reshape)
summary(anova_sigma)

# Conclusion
if (summary(anova_nn)[[1]][["Pr(>F)"]][1] < 0.05) {
  conclusion1 = "Overall, the choice of sigma significantly affects the results for a given ``nn``."
  cat(conclusion1)
} else {
  conclusion1 = "Overall, the choice of sigma does not significantly affect the results for a given ``nn``."
  cat(conclusion1)
}

if (summary(anova_sigma)[[1]][["Pr(>F)"]][1] < 0.05) {
  conclusion2 = "Overall, the choice of nn significantly affects the results for a given ``sigma``."
  cat(conclusion2)
} else {
  conclusion2 = "Overall, the choice of nn does not significantly affect the results for a given ``sigma``."
  cat(conclusion2)
}

# Determine overall effect
if (summary(anova_nn)[[1]][["Pr(>F)"]][1] < 0.05 &&
    summary(anova_sigma)[[1]][["Pr(>F)"]][1] < 0.05) {
  nn_mean_diff <- max(nn_summary$mean) - min(nn_summary$mean)
  sigma_mean_diff <- max(sigma_summary$mean) - min(sigma_summary$mean)
  
  # Means
  if (nn_mean_diff > sigma_mean_diff) {
    conclusion3 = paste0("The ``nn`` effect is stronger than the ``sigma`` effect based on mean",
                         " Pearson similarities (", nn_mean_diff, " vs. ", sigma_mean_diff, ").")
    cat(conclusion3)
  } else if (nn_mean_diff < sigma_mean_diff) {
    conclusion3 = paste0("The ``sigma`` effect is stronger than the ``nn`` effect based on mean",
                         " Pearson similarities (", sigma_mean_diff, " vs. ", nn_mean_diff, ").")
    cat(conclusion3)
  }
  
  # Standard deviations
  nn_sd_diff <- max(nn_summary$sd) - min(nn_summary$sd)
  sigma_sd_diff <- max(sigma_summary$sd) - min(sigma_summary$sd)
  
  if (nn_sd_diff > sigma_sd_diff) {
    conclusion4 = paste0("The ``nn`` effect is stronger than the sigma effect based on", 
                         " the standard deviation of Pearson similarities (",
                         nn_sd_diff, " vs. ", sigma_sd_diff, ").")
    cat(conclusion4)
  } else if (nn_sd_diff < sigma_sd_diff) {
    conclusion4 = paste0("The ``sigma`` effect is stronger than the nn effect based on",  
                         " the standard deviation of Pearson similarities (",
                         sigma_sd_diff, " vs. ", nn_sd_diff, ").")
    cat(conclusion4)
  }
} else {
  conclusion3 = "The ``nn`` effect and ``sigma`` effects are practically equal based on mean Pearson similarities."
  cat(conclusion3)
}

if (exists("conclusion4")) {
  conclusion = paste(conclusion1, conclusion2, conclusion3, conclusion4)
  rm(conclusion1, conclusion2, conclusion3, conclusion4)
  sig_status = TRUE
} else {
  conclusion = paste(conclusion1, conclusion2, conclusion3)
  rm(conclusion1, conclusion2, conclusion3)
  sig_status = FALSE
}

# Non-parametric ##
# Reshape data for Kruskal-Wallis Test
nn_reshape <- reshape_Pearson_for_tests(nn_similarities)  
sigma_reshape <- reshape_Pearson_for_tests(sigma_similarities)  

# Perform Kruskal-Wallis test for nn
kruskal_nn <- kruskal.test(Value ~ Factor, data = nn_reshape)
print(kruskal_nn)

# Perform Kruskal-Wallis test for sigma
kruskal_sigma <- kruskal.test(Value ~ Factor, data = sigma_reshape)
print(kruskal_sigma)

# Perform pairwise Wilcoxon tests if Kruskal-Wallis is significant
if (kruskal_nn$p.value < 0.05) {
  pairwise_nn <- pairwise.wilcox.test(nn_reshape$Value, nn_reshape$Factor, p.adjust.method = "bonferroni")
  print(pairwise_nn)
}

if (kruskal_sigma$p.value < 0.05) {
  pairwise_sigma <- pairwise.wilcox.test(sigma_reshape$Value, sigma_reshape$Factor, p.adjust.method = "bonferroni")
  print(pairwise_sigma)
}

# Initialize conclusion variables to avoid undefined errors
np_conclusion1 <- NULL
np_conclusion2 <- NULL
np_conclusion3 <- NULL
np_conclusion4 <- NULL

# np_conclusion
if (kruskal_nn$p.value < 0.05) {
  np_conclusion1 <- "Overall, the choice of sigma significantly affects the results for a given `nn`."
  cat(np_conclusion1, "\n")
  
  # Additional comparisons for mean or median if needed
  nn_median_diff <- max(nn_summary$median) - min(nn_summary$median)
  sigma_median_diff <- max(sigma_summary$median) - min(sigma_summary$median)
  
  if (nn_median_diff > sigma_median_diff) {
    np_conclusion3 <- paste0("The `nn` effect is stronger than the `sigma` effect based on median",
                             " Pearson similarities (", nn_median_diff, " vs. ", sigma_median_diff, ").")
    cat(np_conclusion3, "\n")
  } else if (nn_median_diff < sigma_median_diff) {
    np_conclusion3 <- paste0("The `sigma` effect is stronger than the `nn` effect based on median",
                             " Pearson similarities (", sigma_median_diff, " vs. ", nn_median_diff, ").")
    cat(np_conclusion3, "\n")
  } else {
    np_conclusion3 <- "The `nn` effect and `sigma` effects are practically equal based on median Pearson similarities."
    cat(np_conclusion3, "\n")
  }
} else {
  np_conclusion1 <- "Overall, the choice of sigma does not significantly affect the results for a given `nn`."
  cat(np_conclusion1, "\n")
}

if (kruskal_sigma$p.value < 0.05) {
  np_conclusion2 <- "Overall, the choice of nn significantly affects the results for a given `sigma`."
  cat(np_conclusion2, "\n")
  
  # Additional comparisons for standard deviations if needed
  nn_sd_diff <- max(nn_summary$sd) - min(nn_summary$sd)
  sigma_sd_diff <- max(sigma_summary$sd) - min(sigma_summary$sd)
  
  if (nn_sd_diff > sigma_sd_diff) {
    np_conclusion4 <- paste0("The `nn` effect is stronger than the sigma effect based on", 
                             " the standard deviation of Pearson similarities (",
                             nn_sd_diff, " vs. ", sigma_sd_diff, ").")
    cat(np_conclusion4, "\n")
  } else if (nn_sd_diff < sigma_sd_diff) {
    np_conclusion4 <- paste0("The `sigma` effect is stronger than the nn effect based on",  
                             " the standard deviation of Pearson similarities (",
                             sigma_sd_diff, " vs. ", nn_sd_diff, ").")
    cat(np_conclusion4, "\n")
  } else {
    np_conclusion4 <- "The `nn` effect and `sigma` effects are practically equal based on the standard deviation of Pearson similarities."
    cat(np_conclusion4, "\n")
  }
} else {
  np_conclusion2 <- "Overall, the choice of nn does not significantly affect the results for a given `sigma`."
  cat(np_conclusion2, "\n")
}

# Consolidate all np_conclusions
np_conclusion <- paste(c(np_conclusion1, np_conclusion2, np_conclusion3, np_conclusion4)[!sapply(c(np_conclusion1, np_conclusion2, np_conclusion3, np_conclusion4), is.null)], collapse = " ")
np_sig_status <- (kruskal_nn$p.value < 0.05) | (kruskal_sigma$p.value < 0.05)
cat(np_conclusion)

# Handle sig_status_final
if (np_sig_status == FALSE && sig_status == FALSE) {
  sig_status_final = FALSE
} else {
  sig_status_final = TRUE
}

# If no significant differences are shown between/across hyperparameters then pick median values
if (!sig_status_final){
  optN = median(num_neighbors_range) # 30
  optSigma = median(sigma_range) # 0.5
}

# We will therefore pick the median sigma, rounded down to 0.5, a value we actually used in affinities.
optSigma = 0.5

# The choice of number of neighbors affects the final matrix more than sigma, which
# also plays a role, according only to parametric tests

# Pearson and Frobenius histograms for different nn AND sigma = 0.5
Fusions_filt = Fusions[which(grepl(paste0("sigma = ", optSigma), names(Fusions)))]
Pearson_hist_matrix = compute_matrix_similarity(Fusions_filt)$Pearson
Pearson_values <- Pearson_hist_matrix[lower.tri(Pearson_hist_matrix, diag = FALSE)]
mean_Pearson_value <- mean(Pearson_values)
median_Pearson_value <- median(Pearson_values)
sd_Pearson_value <- sd(Pearson_values)

Frobenius_hist_matrix = compute_matrix_similarity(Fusions_filt)$Frobenius
Frobenius_values <- Frobenius_hist_matrix[lower.tri(Frobenius_hist_matrix, diag = FALSE)]
mean_Frobenius_value <- mean(Frobenius_values)
median_Frobenius_value <- median(Frobenius_values)
sd_Frobenius_value <- sd(Frobenius_values)

# Plot histogram of Pearson values
library(ggplot2)
ggplot(data = data.frame(Pearson_values), aes(x = Pearson_values)) +
  geom_histogram(breaks = seq(0, 1.07, length.out = length(Pearson_values)),
                 fill = "skyblue", color = "lightblue", size = 0.15) +
  stat_density(aes(color = "Density"), geom = "line", linewidth = 0.4) +
  geom_vline(aes(xintercept = mean_Pearson_value, color = "Mean"), linewidth = 0.2) + 
  geom_vline(aes(xintercept = median_Pearson_value, color = "Median"), linewidth = 0.2) + 
  geom_vline(aes(xintercept = mean_Pearson_value - sd_Pearson_value, color = "Mean - SD"), 
             linetype = "dashed", linewidth = 0.2) + 
  geom_vline(aes(xintercept = mean_Pearson_value + sd_Pearson_value, color = "Mean + SD"), 
             linetype = "dashed", linewidth = 0.2) +
  scale_color_manual(name = "Lines", values = c("Mean" = "red", "Median" = "orange", 
                                                "Mean - SD" = "grey25", "Mean + SD" = "grey25",
                                                "Density" = "darkblue")) +
  labs(title = expression(bold(paste("Histogram of Pearson values between affinity matrices for varying NN and ",
                                     sigma, " = 0.5"))), 
       x = "Affinity Matrix Pearson Values", y = "Frequency") +
  scale_x_continuous(name = "Affinity Matrix Pearson Values", limits = c(0, 1.07),
                     breaks = seq(0, 1.07, 0.1), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(panel.background = element_blank(),
        axis.line = element_line(linewidth = 0.25),
        plot.title = element_text(face = "bold", size = 6.3),
        axis.title = element_text(face = "bold", size = 5.8),
        axis.text = element_text(size = 5),
        axis.ticks = element_line(linewidth = 0.2),
        legend.text = element_text(size = 4.5),
        legend.title = element_text(size = 5, face = "bold"),
        legend.key.spacing.y = unit(1, "mm"),
        legend.key.size = unit(0.25, "cm"),
        legend.box.background = element_rect(color = "black"))
ggsave(filename = paste0(algorithm, "_matrix_Pearson_similarity_histogram.pdf"),
       path = paste0(home, 
                     "/Results/single_algorithm/", algorithm, "/Supplement"), 
       width = 2880, height = 1820, device = 'pdf', units = "px",
       dpi = 700)
dev.off()

# Plot histogram of Frobenius values
ggplot(data = data.frame(Frobenius_values), aes(x = Frobenius_values)) +
  geom_histogram(breaks = seq(0, 2.7, length.out = length(Frobenius_values)),
                 fill = "skyblue", color = "lightblue", size = 0.15) +
  stat_density(aes(color = "Density"), geom = "line", linewidth = 0.4) +
  geom_vline(aes(xintercept = mean_Frobenius_value, color = "Mean"), linewidth = 0.2) + 
  geom_vline(aes(xintercept = median_Frobenius_value, color = "Median"), linewidth = 0.2) + 
  geom_vline(aes(xintercept = mean_Frobenius_value - sd_Frobenius_value, color = "Mean - SD"), 
             linetype = "dashed", linewidth = 0.2) + 
  geom_vline(aes(xintercept = mean_Frobenius_value + sd_Frobenius_value, color = "Mean + SD"), 
             linetype = "dashed", linewidth = 0.2) +
  scale_color_manual(name = "Lines", values = c("Mean" = "red", "Median" = "orange", 
                                                "Mean - SD" = "grey25", "Mean + SD" = "grey25",
                                                "Density" = "darkblue")) +
  labs(title = expression(bold(paste("Histogram of Frobenius values between affinity matrices for varying NN and ",
                                     sigma, " = 0.5"))), 
       x = "Affinity Matrix Frobenius Values", y = "Frequency") +
  scale_x_continuous(name = "Affinity Matrix Frobenius Values", limits = c(0, 2.7),
                     breaks = seq(0, 2.7, 0.25), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(panel.background = element_blank(),
        axis.line = element_line(linewidth = 0.25),
        plot.title = element_text(face = "bold", size = 6.3),
        axis.title = element_text(face = "bold", size = 5.8),
        axis.text = element_text(size = 5),
        axis.ticks = element_line(linewidth = 0.2),
        legend.text = element_text(size = 4.5),
        legend.title = element_text(size = 5, face = "bold"),
        legend.key.spacing.y = unit(1, "mm"),
        legend.key.size = unit(0.25, "cm"),
        legend.box.background = element_rect(color = "black"))
ggsave(filename = paste0(algorithm, "_matrix_Frobenius_similarity_histogram.pdf"),
       path = paste0(home, 
                     "/Results/single_algorithm/", algorithm, "/Supplement"), 
       width = 2880, height = 1820, device = 'pdf', units = "px",
       dpi = 700)
dev.off()

# Spectral clustering using kernlab ###
# Slightly modify the kernlab method to return embeddings as well
setMethod("specc", signature(x = "kernelMatrix"),
          function(x, centers, nystrom.red = FALSE, iterations = 200, ...) {
            m <- nrow(x)
            if (missing(centers))
              stop("centers must be a number or a matrix")
            if (length(centers) == 1) {
              nc <- centers
              if (m < centers)
                stop("more cluster centers than data points.")
            } else {
              nc <- dim(centers)[2]
            }
            
            if (dim(x)[1] != dim(x)[2]) {
              nystrom.red <- TRUE
              if (dim(x)[1] < dim(x)[2])
                x <- t(x)
              m <- nrow(x)
              n <- ncol(x)
            }
            
            if (nystrom.red == TRUE) {
              A <- x[1:n, ]
              B <- x[-(1:n), ]
              d1 <- colSums(rbind(A, B))
              d2 <- rowSums(B) + drop(matrix(colSums(B), 1) %*% .ginv(A) %*% t(B))
              dhat <- sqrt(1/c(d1, d2))
              
              A <- A * (dhat[1:n] %*% t(dhat[1:n]))
              B <- B * (dhat[(n+1):m] %*% t(dhat[1:n]))
              
              Asi <- .sqrtm(.ginv(A))
              Q <- A + Asi %*% crossprod(B) %*% Asi
              tmpres <- svd(Q)
              U <- tmpres$u
              L <- tmpres$d
              
              V <- rbind(A, B) %*% Asi %*% U %*% .ginv(sqrt(diag(L)))
              yi <- matrix(0, m, nc)
              
              for (i in 1:nc)  # Compute the normalized embedding
                yi[, i] <- V[, i] / sqrt(sum(V[, i]^2))
              
              res <- kmeans(yi, centers, iterations)
            } else {
              d <- 1/sqrt(rowSums(x))
              l <- d * x %*% diag(d)
              xi <- eigen(l)$vectors[, 1:nc]
              yi <- xi / sqrt(rowSums(xi^2))
              res <- kmeans(yi, centers, iterations)
            }
            
            # Instead of returning a new "specc" object with just cluster assignments,
            # return a list containing both the clusters and the computed embeddings.
            return(list(cluster = res$cluster,
                        embedding = yi,
                        size = res$size,
                        centers = matrix(0),      # Placeholder for centers
                        withinss = c(0),          # Placeholder for withinss
                        info = "Kernel Matrix used as input."))
          })

# Set RNG version and seed for reproducibility
RNGversion("4.2.2")
set.seed(123)
clusterings <- vector("list", length(Fusions_filt))
names(clusterings) <- names(Fusions_filt)

# Outer loop: iterate over the elements in Fusions_filt
for (i in seq_along(Fusions_filt)) {
  fusion <- Fusions_filt[[i]]
  class(fusion) <- "kernelMatrix"
  
  for(k in k_range) {
    clusterings[[names(Fusions_filt)[i]]][[paste0("k = ", k)]] = specc(fusion, centers = k)@.Data
  }
}

# Name the results and add silhouettes
sil_ranks = list()
for (i in 1:length(clusterings)) {
  for (j in 1:length(clusterings[[i]])) {
    clusterings[[i]][[j]]$silhouette = silhouette(as.integer(clusterings[[i]][[j]]$cluster),
                                                  dist = Rfast::Dist(clusterings[[i]][[j]]$embedding, 
                                                                     method = "euclidean"))
    clusterings[[i]][[j]]$avg_width = summary(clusterings[[i]][[j]]$silhouette)$avg.width
    sil_ranks[[paste0(names(clusterings)[i], ", ", names(clusterings[[i]])[j])]] = clusterings[[i]][[j]]$avg_width
  }
}

# Best clustering
best_sil = sil_ranks[which.max(unlist(sil_ranks))]
optNN = as.numeric(substr(strsplit(names(best_sil), ", ")[[1]][1], 6, 7))
optSigma = as.numeric(substr(strsplit(names(best_sil), ", ")[[1]][2], 9, 11))
optk = as.numeric(substr(strsplit(names(best_sil), ", ")[[1]][3], 5, 
                         nchar(strsplit(names(best_sil), ", ")[[1]][1])))

best_clustering = clusterings[[paste0("NN = ", optNN, ", sigma = ", optSigma)]][[paste0("k = ", optk)]]

RWRNF_clusters = as.data.frame(list(Cluster = best_clustering$cluster, 
                                   Sample.ID = colnames(Fusions_filt[[paste0("NN = ", optNN, ", sigma = ", optSigma)]])))

RWRNF_clusters$Sample.ID = gsub("\\.", "-", RWRNF_clusters$Sample.ID)
rownames(RWRNF_clusters) = RWRNF_clusters$Sample.ID

# Feature ranking
RWRNF_feature_ranks = list()
binary_flags = c(TRUE, FALSE, FALSE, FALSE, FALSE)
for (i in 1:length(input)) {
  RWRNF_feature_ranks[[i]] = rankFeaturesByNMI_parallely(data = list(input[[i]]), 
                                                        W = final_affinity_matrix,
                                                        ncores = 8,
                                                        binary = binary_flags[i],
                                                        nn = optN,
                                                        sigma = optSigma)
  cat("Done with", names(input)[i], "\n")
}
names(RWRNF_feature_ranks) = names(input)

rank_sum_nas = sapply(RWRNF_feature_ranks, function(x) sum(is.na(x$NMI_ranks)))
names(rank_sum_nas) = names(RWRNF_feature_ranks)

feature_ranks_text1 = paste0("Feature ranks could not be calculated for some features in the ",
                             paste(names(rank_sum_nas)[which(rank_sum_nas > 0)], collapse = ", "), 
                             ifelse(length(which(rank_sum_nas > 0)) > 1, " modalities (", " modality ("),
                             paste(rank_sum_nas[which(rank_sum_nas > 0)], collapse = ", "),
                             ").")
rm(rank_sum_nas); gc()

# Combine all modalities into a single data frame
for (i in 1:length(RWRNF_feature_ranks)) {
  names(RWRNF_feature_ranks[[i]]$NMI_scores) = colnames(input[[i]])
}

# Prepare data for plotting
# Combine all modalities into a single data frame
feature_data <- do.call(rbind, lapply(names(RWRNF_feature_ranks), function(modality) {
  data.frame(
    Feature = paste(modality, seq_along(RWRNF_feature_ranks[[modality]]$NMI_scores), sep = "_"),
    Modality = modality,
    NMI_Scores = RWRNF_feature_ranks[[modality]]$NMI_scores
  )
}))

# Remove features with NA NMI scores
feature_data <- feature_data[!is.na(feature_data$NMI_Scores), ]

# Calculate global ranks based on NMI scores
feature_data$Global_Rank <- rank(-feature_data$NMI_Scores, ties.method = "first")

# Filter the top 1000 features
top_features <- feature_data[order(feature_data$Global_Rank), ][1:1000, ]

# Ensure the data is ordered by rank for proper plotting
top_features <- top_features[order(top_features$Global_Rank), ]

# Create the bar plot
ggplot(top_features, aes(x = NMI_Scores, y = Global_Rank, fill = Modality)) +
  geom_bar(
    stat = "identity",
    orientation = "y",
    width = 1,  # Bars fully adjacent with no gaps
    alpha = 0.85,  # Apply transparency
    color = NA  # Removes outlines completely
  ) +
  scale_fill_manual(
    name = "Modality",
    values = c(
      "RNAseq" = "#1B9E77",
      "miRNA" = "#7570B3",
      "Methylation" = "deeppink4",
      "CNV" = "#E7298A",
      "SNPs" = "#66A61E"
    )
  ) +
  scale_x_continuous(
    breaks = seq(0, 1, by = 0.1),  
    limits = c(-0.01, 1),  # Expand lower limit slightly
    expand = c(0, 0)  # Remove extra padding on the x-axis
  ) +
  scale_y_reverse(
    breaks = c(1, seq(100, 1000, by = 100)),  # Y-axis reversed
    limits = c(1001, 0),  # Expand upper limit slightly
    expand = c(0, 0)  # Remove extra padding on the y-axis
  ) +
  labs(
    title = "Top 1000 Features Ranked by NMI with Subtypes",
    x = "NMI Score",
    y = "Global Rank"
  ) +
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 5),  # Smaller y-axis labels
    axis.text.x = element_text(size = 5),  # Larger x-axis labels
    axis.title = element_text(face = "bold", size = 6),
    axis.line.y = element_line(linewidth = 0),
    axis.line.x = element_line(linewidth = 0.1),
    axis.ticks = element_line(linewidth = 0.05),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 7),  # Bold and centrally aligned title
    legend.title = element_text(face = "bold", size = 5),
    legend.text = element_text(size = 5),
    legend.key.size = unit(0.35, "cm"),
    panel.grid.major.y = element_blank(),  # Remove horizontal grid lines
    panel.grid.major.x = element_line(color = "gray90")  # Keep vertical grid lines
  ) +
  guides(color = "none", alpha = "none")

ggsave(filename = paste0(algorithm, "_top_", nrow(top_features), "_feature_ranks.png"),
       path = paste0(home, 
                     "/Results/single_algorithm/", algorithm, "/Supplement"), 
       width = 1920*1.2, height = 1920*1.5, device = 'png', units = "px",
       dpi = 700)
dev.off()

feature_ranks_text2 = paste0("In the top ", nrow(top_features), " features, ",
                             paste(names(table(top_features$Modality)), collapse = ", "),
                             " features are found (",
                             paste(table(top_features$Modality), collapse = ", "),
                             ", respectively).")

feature_ranks_text = paste(feature_ranks_text1, feature_ranks_text2, collapse = " ")
rm(feature_ranks_text1, feature_ranks_text2); gc()

# Main results ###
# Examine cluster similarity to MOVICS by measuring NMI and ARI indices #####
# (Jaccard may be misleading)

# Calculate ARI and NMI
library(mclust)
library(clue)

ARI_to_MOVICS = calculate_ari_index(cluster_df1 = ground_truth_labels,
                                    cluster_df2 = RWRNF_clusters,
                                    sample_col = "Sample.ID",
                                    clust_col = "Cluster",
                                    suffixes = c("_MOVICS_CS",
                                                 paste0("_", algorithm)))

NMI_to_MOVICS = calculate_nmi_index(cluster_df1 = ground_truth_labels,
                                    cluster_df2 = RWRNF_clusters,
                                    sample_col = "Sample.ID",
                                    clust_col = "Cluster",
                                    suffixes = c("_MOVICS_CS",
                                                 paste0("_", algorithm)))

# Low statistics when compared to the MOVICS. Results very different

# MOVICS-like analysis #####
library(MOVICS)
library(ComplexHeatmap)

# Import coloring scheme
scheme = readRDS("Resources/scheme.rds")
annCol = scheme$annCol
annColors = scheme$annColors
cluster_colors = scheme$clust.colors
col.list = scheme$col.list
var2comp = scheme$var2comp %>%
  dplyr::select(-`Consensus Subtype`) %>%
  mutate(Sample.ID = rownames(.)) %>%
  inner_join(RWRNF_clusters, by = "Sample.ID") %>%
  tibble::column_to_rownames(var = "Sample.ID") %>%
  mutate(RWRNF = paste0(algorithm, Cluster)) %>%
  dplyr::select(RWRNF, everything()) %>%
  dplyr::select(-Cluster)
rm(scheme); gc()

# Silhouette
getSilhouette_ggplot(sil      = best_clustering$silhouette,
                     fig.path = paste0(home, "/Results/single_algorithm/", algorithm),
                     fig.name = "Silhouette",
                     height   = 5.5,
                     width    = 5.5,
                     axis_label_size = 12,
                     axis_label_font = "bold",
                     text_size = 1.5,
                     title_size = 16,
                     algorithm = algorithm,
                     save_plot = TRUE)
dev.off()

# Heatmap prep
plotdata <- lapply(lapply(input, as.matrix), 
                   function(mat) mat[, colSums(mat != 0) > 0])
plotdata <- lapply(plotdata, t)
plotdata = getStdiz(
  data = plotdata,
  halfwidth = c(NA, 3, 3, 3, 3), # No halfwidth for SNPs
  centerFlag = c(F, F, F, F, F),
  scaleFlag = c(F, F, F, F, F)
)

plot_object = list(clust.res = RWRNF_clusters %>%
                     dplyr::rename(samID = Sample.ID, clust = Cluster))

# Export consensus clustering object
clust = as.data.frame(plot_object$clust.res)
colnames(clust) = c("Sample.ID", "Cluster")
clust$Cluster = paste0(algorithm, clust$Cluster)
openxlsx::write.xlsx(clust, paste0(home, "/Results/single_algorithm/", algorithm, "/", 
                                   algorithm, "_", data_source, "_",
                                   data_types, "_eval_on_", evaluation_source,
                                   "_clusterings.xlsx"))

# Save environment
save.image(paste0(home, "/Results/single_algorithm/", 
                  algorithm, "/", algorithm, "_", data_source, "_",
                  data_types, "_eval_on_", evaluation_source,
                  "_env.RData"))