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
optN = as.numeric(substr(strsplit(names(best_sil), ", ")[[1]][1], 6, 7))
optSigma = as.numeric(substr(strsplit(names(best_sil), ", ")[[1]][2], 9, 11))
optk = as.numeric(substr(strsplit(names(best_sil), ", ")[[1]][3], 5, 
                         nchar(strsplit(names(best_sil), ", ")[[1]][1])))

best_clustering = clusterings[[paste0("NN = ", optN, ", sigma = ", optSigma)]][[paste0("k = ", optk)]]

RWRNF_clusters = as.data.frame(list(Cluster = best_clustering$cluster, 
                                   Sample.ID = colnames(Fusions_filt[[paste0("NN = ", optNN, ", sigma = ", optSigma)]])))

RWRNF_clusters$Sample.ID = gsub("\\.", "-", RWRNF_clusters$Sample.ID)
rownames(RWRNF_clusters) = RWRNF_clusters$Sample.ID

# Feature ranking
final_affinity_matrix = Fusions_filt[[paste0("NN = ", optN, ", sigma = ", optSigma)]]
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
clust = as.data.frame(plot_object$clust.res) %>%
  dplyr::select(Sample.ID = samID, Cluster = clust)
clust$Cluster = paste0(algorithm, clust$Cluster)
openxlsx::write.xlsx(clust, paste0(home, "/Results/single_algorithm/", algorithm, "/", 
                                   algorithm, "_", data_source, "_",
                                   data_types, "_eval_on_", evaluation_source,
                                   "_clusterings.xlsx"))

# Order features
feature_orders = readRDS("Resources/TCGA/mm_feature_orders.rds")
for (i in 1:length(plotdata)) {
  plotdata[[i]] = plotdata[[i]][feature_orders[[names(plotdata)[i]]], , drop = FALSE]
}

getMoHeatmap_single_algorithm2(algorithm_name = algorithm,
                               data          = plotdata,
                               row.title     = names(plotdata),
                               is.binary     = c(T,F,F,F,F), 
                               legend.name   = c("SNPs",
                                                 "Standardized RNAseq norm. counts",
                                                 "Standardized CNV",
                                                 "Standardized miRNA norm. counts",
                                                 "Standardized Methylation M-values"
                               ),
                               cluster_rows = rep(F, length(plotdata)),
                               cluster_cols = rep(F, length(plotdata)),
                               show.col.dend = rep(F, length(plotdata)),
                               show.colnames = FALSE,
                               show.row.dend = rep(F, length(plotdata)),
                               show.rownames = rep(F, length(plotdata)),
                               clust.res     = plot_object$clust.res, # consensusMOIC-like results
                               # clust.dist.row = c("manhattan", rep("euclidean", 4)),
                               # clust.method.row = rep("ward.D", length(plotdata)),
                               annRow        = NULL, # no selected features
                               color         = col.list,
                               annCol        = annCol, # annotation for samples
                               annColors     = annColors, # annotation color
                               width         = 20, # width of each subheatmap
                               height        = 10, # height of each subheatmap
                               fig.path      = paste0(home, "/Results/single_algorithm/", algorithm),
                               fig.name      = paste0("default_", algorithm, "_Comprehensive_heatmap"))
dev.off()
gc()

# Clinical variables ###
# Remove unknown levels for statistical tests
var2comp_nonas = var2comp
for (i in 1:ncol(var2comp)) {
  nas = which(var2comp[, i] == "Unknown")
  var2comp_nonas[nas, i] = NA
  empties = which(var2comp[, i] == "")
  var2comp_nonas[empties, i] = NA
}
rm(nas, empties); gc()

# Statistical comparisons
clin_comp = compClinvar_single_algorithm(algorithm_name = algorithm,
                                         moic.res = plot_object,
                                         var2comp = var2comp_nonas,
                                         strata = "RWRNF",
                                         factorVars = c("vital_status", "race_list", "ethnicity",
                                                        "history_of_neoadjuvant_treatment",
                                                        "primary_lymph_node_presentation_assessment",
                                                        "histological_type", "menopause_status",
                                                        "breast_carcinoma_progesterone_receptor_status",
                                                        "breast_carcinoma_estrogen_receptor_status",
                                                        "lab_proc_her2_neu_immunohistochemistry_receptor_status",
                                                        "distant_metastasis_present_ind2",
                                                        "stage_event_pathologic_stage"),
                                         nonnormalVars = c("days_to_birth", "days_to_death",
                                                           "days_to_last_known_alive", 
                                                           "days_to_last_followup",
                                                           "age_at_initial_pathologic_diagnosis",
                                                           "er_level_cell_percentage_category",
                                                           "progesterone_receptor_level_cell_percent_category",
                                                           "number_of_lymphnodes_positive_by_ihc",
                                                           "number_of_lymphnodes_positive_by_he"),
                                         includeNA = FALSE,
                                         doWord = TRUE,
                                         tab.name = "Summary_of_clinical_variables",
                                         res.path = paste0(home, "/Results/single_algorithm/", algorithm, "/"),
                                         output_pdf = TRUE,
                                         pdf_level_col_width = c("7em", "10em"),
                                         pdf_count_col_width = "10em",
                                         pdf_pval_col_width = "3em",
                                         pdf_test_col_width = "8em",
                                         pdf_tab_font_size = 9)

clin_ordinal_comp = compClinvar_ordinal_single_algorithm(algorithm_name = algorithm,
                                                         moic.res = plot_object,
                                                         var2comp = var2comp_nonas %>%
                                                           dplyr::select(number_of_lymphnodes_positive_by_ihc,
                                                                         number_of_lymphnodes_positive_by_he,
                                                                         RWRNF),
                                                         strata = "RWRNF",
                                                         ordinalVars = c("number_of_lymphnodes_positive_by_ihc",
                                                                         "number_of_lymphnodes_positive_by_he"),
                                                         includeNA = FALSE,
                                                         tab.name = "Summary of ordinal clinical variables",
                                                         res.path = paste0(home, "/Results/single_algorithm/", algorithm, "/"),
                                                         output_pdf = TRUE,
                                                         pdf_template_loc = paste0(home, "/Scripts/automated_scripts/clincomp_template.Rmd"),
                                                         pdf_level_col_width = c("7em", "10em"),
                                                         pdf_count_col_width = "10em",
                                                         pdf_pval_col_width = "3em",
                                                         pdf_test_col_width = "8em",
                                                         pdf_tab_font_size = 9)

# Oncoprint ###
oncoprint <- compMut_single_algorithm(algorithm_name = algorithm,
                                      moic.res  = plot_object,
                                      mut.matrix   = plotdata$SNPs, # binary somatic mutation matrix
                                      doWord       = TRUE, # generate table in .docx format
                                      doPlot       = TRUE, # draw OncoPrint
                                      freq.cutoff  = 0.05, # keep those genes that mutated in at least 5% of samples
                                      p.adj.cutoff = 0.05, # keep those genes with adjusted p value < 0.05 to draw OncoPrint
                                      innerclust   = TRUE, # perform clustering within each subtype
                                      annCol       = annCol, # same annotation for heatmap
                                      annColors    = annColors, # same annotation color for heatmap
                                      width        = 12, 
                                      height       = 6,
                                      fig.name     = paste0(algorithm, "_", data_source, "_",
                                                            data_types, "_eval_on_", evaluation_source,
                                                            "_oncoprint"),
                                      tab.name     = "Independent test between subtype and mutation",
                                      fig.path     = paste0(home, "/Results/single_algorithm/", algorithm),
                                      res.path     = paste0(home, "/Results/single_algorithm/", algorithm))

# Drug sensitivity comparison ###
drug_sensitivity <- compDrugsen_single_algorithm(algorithm_name = algorithm,
                                                 moic.res    = plot_object,
                                                 norm.expr   = plotdata$RNAseq,
                                                 drugs       = c("Cisplatin", "Paclitaxel", "Lapatinib",
                                                                 "Doxorubicin", "5-Fluorouracil",
                                                                 "Sorafenib"), # a vector of names of drug in GDSC
                                                 tissueType  = "breast", # choose specific tissue type to construct ridge regression model
                                                 test.method = "nonparametric", # statistical testing method
                                                 prefix      = "Violin_plot_of_IC50",
                                                 seed = 123,
                                                 fig.path = paste0(home, "/Results/single_algorithm/", algorithm))

# Agreement with other subtypes ###
subtype_agreement <- compAgree_single_algorithm(algorithm_name = algorithm,
                                                moic.res  = plot_object,
                                                subt2comp = annCol[, c("ER status", "PR status",
                                                                       "HER2 status", "Metastasis", "Stage")],
                                                doPlot    = TRUE,
                                                box.width = 0.2,
                                                fig.name  = "Classification_agreement",
                                                fig.path  = paste0(home, "/Results/single_algorithm/", algorithm),
                                                width     = 12)
dev.off()

# DGEA ###
dgea = runDEA_mod(dea.method = "limma", # we use normalized data as input
                  expr = plotdata$RNAseq,
                  moic.res = plot_object,
                  prefix = "dgea_",
                  sort.p = TRUE,
                  overwt = TRUE,
                  verbose = TRUE,
                  res.path = paste0(home, "/Results/single_algorithm/", algorithm),
                  algorithm = algorithm)

# # Identify unique subtype biomarkers
# # 1. Up-regulated markers
dgea.marker.up <- runMarker_single_algorithm(algorithm_name = algorithm,
                                             moic.res = plot_object,
                                             dea.method    = "limma", # name of DEA method
                                             prefix        = "dgea_", # MUST be the same of argument in runDEA()
                                             dat.path      = paste0(home, "/Results/single_algorithm/", algorithm), # path of DEA files
                                             res.path      = paste0(home, "/Results/single_algorithm/", algorithm), # path to save marker files
                                             p.cutoff      = 0.05, # p cutoff to identify significant DEGs
                                             p.adj.cutoff  = 0.05, # padj cutoff to identify significant DEGs
                                             dirct         = "up", # direction of dysregulation in expression
                                             n.marker      = 100, # number of biomarkers for each subtype
                                             doplot        = TRUE, # generate diagonal heatmap
                                             norm.expr     = plotdata$RNAseq, # use normalized expression as heatmap input
                                             annCol        = annCol, # sample annotation in heatmap
                                             annColors     = annColors, # colors for sample annotation
                                             show_rownames = TRUE, # show no rownames (biomarker name)
                                             centerFlag = F,
                                             scaleFlag = F,
                                             halfwidth = 3,
                                             fig.name      = "upregulated_biomarkers_heatmap",
                                             fig.path = paste0(home, "/Results/single_algorithm/", algorithm),
                                             width = 14,
                                             height = 12,
                                             fontsize_row = 3,
                                             name = "normalized RNA-seq")
dev.off()

# # 2. Down-regulated markers
dgea.marker.down <- runMarker_single_algorithm(algorithm_name = algorithm,
                                               moic.res = plot_object,
                                               dea.method    = "limma", # name of DEA method
                                               prefix        = "dgea_", # MUST be the same of argument in runDEA()
                                               dat.path      = paste0(home, "/Results/single_algorithm/", algorithm), # path of DEA files
                                               res.path      = paste0(home, "/Results/single_algorithm/", algorithm), # path to save marker files
                                               p.cutoff      = 0.05, # p cutoff to identify significant DEGs
                                               p.adj.cutoff  = 0.05, # padj cutoff to identify significant DEGs
                                               dirct         = "down", # direction of dysregulation in expression
                                               n.marker      = 100, # number of biomarkers for each subtype
                                               doplot        = TRUE, # generate diagonal heatmap
                                               norm.expr     = plotdata$RNAseq, # use normalized expression as heatmap input
                                               annCol        = annCol, # sample annotation in heatmap
                                               annColors     = annColors, # colors for sample annotation
                                               show_rownames = TRUE, # show no rownames (biomarker name)
                                               centerFlag = F,
                                               scaleFlag = F,
                                               halfwidth = 3,
                                               fig.name      = "downregulated_biomarkers_heatmap",
                                               fig.path = paste0(home, "/Results/single_algorithm/", algorithm),
                                               width = 14,
                                               height = 12,
                                               fontsize_row = 3,
                                               name = "normalized RNA-seq")
dev.off()

# DMEA ###
dmea = runDEA_mod(dea.method = "limma", # we use normalized data as input
                  expr = plotdata$Methylation,
                  moic.res = plot_object,
                  prefix = "dmea_",
                  sort.p = TRUE,
                  overwt = TRUE,
                  verbose = TRUE,
                  res.path = paste0(home, "/Results/single_algorithm/", algorithm),
                  algorithm = algorithm)

# # Identify unique subtype biomarkers
# # 1. Up-regulated markers
methyl.marker.up <- runMarker_single_algorithm(algorithm_name = algorithm,
                                               moic.res = plot_object,
                                               dea.method    = "limma", # name of DEA method
                                               prefix        = "dmea_", # MUST be the same of argument in runDEA()
                                               dat.path      = paste0(home, "/Results/single_algorithm/", algorithm), # path of DEA files
                                               res.path      = paste0(home, "/Results/single_algorithm/", algorithm), # path to save marker files
                                               p.cutoff      = 0.05, # p cutoff to identify significant DEGs
                                               p.adj.cutoff  = 0.05, # padj cutoff to identify significant DEGs
                                               dirct         = "up", # direction of dysregulation in expression
                                               n.marker      = 100, # number of biomarkers for each subtype
                                               doplot        = TRUE, # generate diagonal heatmap
                                               norm.expr     = plotdata$Methylation, # use normalized expression as heatmap input
                                               annCol        = annCol, # sample annotation in heatmap
                                               annColors     = annColors, # colors for sample annotation
                                               show_rownames = TRUE, # show no rownames (biomarker name)
                                               centerFlag = F,
                                               scaleFlag = F,
                                               halfwidth = 3,
                                               fig.name      = "hypermethylated_biomarkers_heatmap",
                                               fig.path = paste0(home, "/Results/single_algorithm/", algorithm),
                                               width = 14,
                                               height = 12,
                                               fontsize_row = 0, # 3 default
                                               name = "normalized Methylation")
dev.off()

# # 2. Down-regulated markers
methyl.marker.down <- runMarker_single_algorithm(algorithm_name = algorithm,
                                                 moic.res = plot_object,
                                                 dea.method    = "limma", # name of DEA method
                                                 prefix        = "dmea_", # MUST be the same of argument in runDEA()
                                                 dat.path      = paste0(home, "/Results/single_algorithm/", algorithm), # path of DEA files
                                                 res.path      = paste0(home, "/Results/single_algorithm/", algorithm), # path to save marker files
                                                 p.cutoff      = 0.05, # p cutoff to identify significant DEGs
                                                 p.adj.cutoff  = 0.05, # padj cutoff to identify significant DEGs
                                                 dirct         = "down", # direction of dysregulation in expression
                                                 n.marker      = 100, # number of biomarkers for each subtype
                                                 doplot        = TRUE, # generate diagonal heatmap
                                                 norm.expr     = plotdata$Methylation, # use normalized expression as heatmap input
                                                 annCol        = annCol, # sample annotation in heatmap
                                                 annColors     = annColors, # colors for sample annotation
                                                 show_rownames = TRUE, # show no rownames (biomarker name)
                                                 centerFlag = F,
                                                 scaleFlag = F,
                                                 halfwidth = 3,
                                                 fig.name      = "hypomethylated_biomarkers_heatmap",
                                                 fig.path = paste0(home, "/Results/single_algorithm/", algorithm),
                                                 width = 14,
                                                 height = 12,
                                                 fontsize_row = 0, # 3 default
                                                 name = "normalized Methylation")
dev.off()

# DmiREA ###
dmiRea = runDEA_mod(dea.method = "limma", # we use normalized data as input
                    expr = plotdata$miRNA,
                    moic.res = plot_object,
                    prefix = "dmiRea_",
                    sort.p = TRUE,
                    overwt = TRUE,
                    verbose = TRUE,
                    res.path = paste0(home, "/Results/single_algorithm/", algorithm),
                    algorithm = algorithm)

# # Identify unique subtype biomarkers
# # 1. Up-regulated markers
miRNA.marker.up <- runMarker_single_algorithm(algorithm_name = algorithm,
                                              moic.res = plot_object,
                                              dea.method    = "limma", # name of DEA method
                                              prefix        = "dmiRea_", # MUST be the same of argument in runDEA()
                                              dat.path      = paste0(home, "/Results/single_algorithm/", algorithm), # path of DEA files
                                              res.path      = paste0(home, "/Results/single_algorithm/", algorithm), # path to save marker files
                                              p.cutoff      = 0.05, # p cutoff to identify significant DEGs
                                              p.adj.cutoff  = 0.05, # padj cutoff to identify significant DEGs
                                              dirct         = "up", # direction of dysregulation in expression
                                              n.marker      = 100, # number of biomarkers for each subtype
                                              doplot        = TRUE, # generate diagonal heatmap
                                              norm.expr     = plotdata$miRNA, # use normalized expression as heatmap input
                                              annCol        = annCol, # sample annotation in heatmap
                                              annColors     = annColors, # colors for sample annotation
                                              show_rownames = TRUE, # show no rownames (biomarker name)
                                              centerFlag = F,
                                              scaleFlag = F,
                                              halfwidth = 3,
                                              fig.name      = "upregulated_miRNA_biomarkers_heatmap",
                                              fig.path = paste0(home, "/Results/single_algorithm/", algorithm),
                                              width = 14,
                                              height = 12,
                                              fontsize_row = 0, # 3 default
                                              name = "normalized miRNA")
dev.off()

# # 2. Down-regulated markers
miRNA.marker.down <- runMarker_single_algorithm(algorithm_name = algorithm,
                                                moic.res = plot_object,
                                                dea.method    = "limma", # name of DEA method
                                                prefix        = "dmiRea_", # MUST be the same of argument in runDEA()
                                                dat.path      = paste0(home, "/Results/single_algorithm/", algorithm), # path of DEA files
                                                res.path      = paste0(home, "/Results/single_algorithm/", algorithm), # path to save marker files
                                                p.cutoff      = 0.05, # p cutoff to identify significant DEGs
                                                p.adj.cutoff  = 0.05, # padj cutoff to identify significant DEGs
                                                dirct         = "down", # direction of dysregulation in expression
                                                n.marker      = 100, # number of biomarkers for each subtype
                                                doplot        = TRUE, # generate diagonal heatmap
                                                norm.expr     = plotdata$miRNA, # use normalized expression as heatmap input
                                                annCol        = annCol, # sample annotation in heatmap
                                                annColors     = annColors, # colors for sample annotation
                                                show_rownames = TRUE, # show no rownames (biomarker name)
                                                centerFlag = F,
                                                scaleFlag = F,
                                                halfwidth = 3,
                                                fig.name      = "downregulated_miRNA_biomarkers_heatmap",
                                                fig.path = paste0(home, "/Results/single_algorithm/", algorithm),
                                                width = 14,
                                                height = 12,
                                                fontsize_row = 0, # 3 default
                                                name = "normalized miRNA")
dev.off()

# GSEA ###
# Load MSigDb file
MSIGDB.FILE <- paste0(home, "/Resources/Pathways/GO-BP_c5.go.bp.v2024.1.Hs.symbols.gmt")

# GSEA up-regulated
RNGversion("4.2.2")
set.seed(123)
gsea.up <- runGSEA_mod_4.4_single_algorithm(algorithm_name = algorithm,
                                            moic.res     = plot_object,
                                            dea.method   = "limma", # name of DEA method
                                            prefix       = "dgea_", # MUST be the same of argument in runDEA()
                                            dat.path      = paste0(home, "/Results/single_algorithm/", algorithm), # path of DEA files
                                            res.path      = paste0(home, "/Results/single_algorithm/", algorithm), # path to save marker files
                                            msigdb.path  = MSIGDB.FILE, # MUST be the ABSOLUTE path of msigdb file
                                            norm.expr    = plotdata$RNAseq, # use normalized expression to calculate enrichment score
                                            dirct        = "up", # direction of dysregulation in pathway
                                            n.path       = 20,
                                            p.cutoff     = 0.05, # p cutoff to identify significant pathways
                                            p.adj.cutoff = 0.05, # padj cutoff to identify significant pathways
                                            gsva.method  = "gsva", # method to calculate single sample enrichment score
                                            name         = "GSVA scores", # name for colorbar
                                            norm.method  = "mean", # normalization method to calculate subtype-specific enrichment score
                                            fig.name     = "upregulated_pathway_heatmap",
                                            nPerm = 10000,
                                            minGSSize = 5,
                                            maxGSSize = 500,
                                            fig.path = paste0(home, "/Results/single_algorithm/", algorithm),
                                            width = 14, height = 12)

# GSEA down-regulated
RNGversion("4.2.2")
set.seed(123)
gsea.down <- runGSEA_mod_4.4_single_algorithm(algorithm_name = algorithm,
                                              moic.res     = plot_object,
                                              dea.method   = "limma", # name of DEA method
                                              prefix       = "dgea_", # MUST be the same of argument in runDEA()
                                              dat.path      = paste0(home, "/Results/single_algorithm/", algorithm), # path of DEA files
                                              res.path      = paste0(home, "/Results/single_algorithm/", algorithm), # path to save marker files
                                              msigdb.path  = MSIGDB.FILE, # MUST be the ABSOLUTE path of msigdb file
                                              norm.expr    = plotdata$RNAseq, # use normalized expression to calculate enrichment score
                                              dirct        = "down", # direction of dysregulation in pathway
                                              n.path       = 20,
                                              p.cutoff     = 0.05, # p cutoff to identify significant pathways
                                              p.adj.cutoff = 0.05, # padj cutoff to identify significant pathways
                                              gsva.method  = "gsva", # method to calculate single sample enrichment score
                                              name         = "GSVA scores", # name for colorbar
                                              norm.method  = "mean", # normalization method to calculate subtype-specific enrichment score
                                              fig.name     = "downregulated_pathway_heatmap",
                                              nPerm = 10000,
                                              minGSSize = 5,
                                              maxGSSize = 500,
                                              fig.path = paste0(home, "/Results/single_algorithm/", algorithm),
                                              width = 14, height = 12)

# Gene set variation analysis #####
# locate ABSOLUTE path of gene set file
GSET.FILE <- paste0(home, "/Resources/Pathways/gene_sets_of_interest.gmt")

RNGversion("4.2.2")
set.seed(123)
gsva.res = runGSVA_mod_4.4_single_algorithm(algorithm_name = algorithm,
                                            moic.res      = plot_object,
                                            norm.expr     = plotdata$RNAseq,
                                            gset.gmt.path = GSET.FILE, # ABSOLUTE path of gene set file
                                            gsva.method   = "gsva", # method to calculate single sample enrichment score
                                            annCol        = annCol,
                                            annColors     = annColors,
                                            fig.path      = paste0(home, "/Results/single_algorithm/", algorithm),
                                            fig.name      = "gene_sets_of_interest_heatmap",
                                            centerFlag    = F,
                                            scaleFlag     = F,
                                            distance      = 'euclidean',
                                            linkage       = 'average',
                                            show_rownames = TRUE,
                                            show_colnames = FALSE,
                                            height        = 8,
                                            width         = 12,
                                            name          = "GSVA scores")
dev.off()

# Hierarchical clustering of pathways
library(pathfindR)
library(fastcluster)

# Get unique pathways for each subtype
GSEAfiles_up <- sort(dir(paste0(home, "/Results/single_algorithm/", algorithm), 
                         pattern = "unique_upexpr_pathway.txt$"))
GSEAfiles_down <- sort(dir(paste0(home, "/Results/single_algorithm/", algorithm), 
                           pattern = "unique_downexpr_pathway.txt$"))

unique_upexpr_pathways = list()
for (i in 1:length(gsea.up$gsea.list)) {
  unique_upexpr_pathways[[i]] = data.table::fread(paste0(paste0(home, "/Results/single_algorithm/", algorithm), 
                                                         "/", GSEAfiles_up[i]),
                                                  header = TRUE, sep = "\t")
}

unique_downexpr_pathways = list()
for (i in 1:length(gsea.down$gsea.list)) {
  unique_downexpr_pathways[[i]] = data.table::fread(paste0(paste0(home, "/Results/single_algorithm/", algorithm), 
                                                           "/", GSEAfiles_down[i]),
                                                    header = TRUE, sep = "\t")
}

names(unique_downexpr_pathways) = names(unique_upexpr_pathways) = names(gsea.up$gsea.list)

# Filter GSEA input
gsea.up_unique = gsea.up
for (i in 1:length(unique_upexpr_pathways)) {
  unq = unique_upexpr_pathways[[i]]$V1
  gsea.up_unique$gsea.list[[i]]@result = gsea.up_unique$gsea.list[[i]]@result[gsea.up_unique$gsea.list[[i]]@result$ID %in%
                                                                                unq, ]
}

gsea.down_unique = gsea.down
for (i in 1:length(unique_downexpr_pathways)) {
  unq = unique_downexpr_pathways[[i]]$V1
  gsea.down_unique$gsea.list[[i]]@result = gsea.down_unique$gsea.list[[i]]@result[gsea.down_unique$gsea.list[[i]]@result$ID %in%
                                                                                    unq, ]
}

rm(unq); gc()

hclust_input_up = prepare_gsea_output_for_hclust(gsea_output = gsea.up_unique, 
                                                 dgea_output_name_style = "dgea_", 
                                                 dea.method = "limma", 
                                                 mo.method = "",
                                                 dat.path = paste0(home, "/Results/single_algorithm/", algorithm), 
                                                 dgea_padj_cutoff = 0.05,
                                                 logfc_cutoff = 0,
                                                 pathway_padj_cutoff = 0.05)

hclust_input_down = prepare_gsea_output_for_hclust(gsea_output = gsea.down_unique, 
                                                   dgea_output_name_style = "dgea_", 
                                                   dea.method = "limma", 
                                                   mo.method = "",
                                                   dat.path = paste0(home, "/Results/single_algorithm/", algorithm), 
                                                   dgea_padj_cutoff = 0.05,
                                                   logfc_cutoff = 0,
                                                   pathway_padj_cutoff = 0.05)

hclust_input = c(hclust_input_up, hclust_input_down)
names(hclust_input) = c(paste0(rep("up_", length(hclust_input_up)), 
                               names(hclust_input_up)),
                        paste0(rep("down_", length(hclust_input_down)), 
                               names(hclust_input_down)))

rm(hclust_input_down, hclust_input_up); gc()

# Load doParallel if not already loaded
library(parallel)
library(foreach)
library(doParallel)

# Set up the number of cores to use: minimum of length(hclust_input) or 5
cl <- makeCluster(min(length(hclust_input), 5))
registerDoParallel(cl)

# Use foreach with parallel processing
timestamp()
hclust_output <- foreach(i = 1:length(hclust_input), .packages = c("pathfindR", "fastcluster")) %dopar% {
  RNGversion("4.2.2")
  set.seed(123)
  source("Scripts/automated_scripts/fast_pathfindR_hclust.R")
  
  # Perform clustering, handle errors
  result <- cluster_enriched_terms_fast(hclust_input[[i]], method = "hierarchical", plot_clusters_graph = FALSE,
                                        use_description = FALSE, use_active_snw_genes = FALSE)
  if (is.character(result) && result == "hclust impossible") {
    return("hclust impossible")
  } else {
    return(result)
  }
}

timestamp() # ~2.5 mins
stopCluster(cl)
gc()
names(hclust_output) <- names(hclust_input)

# Are there any null sets?
which(hclust_output == "hclust impossible")

# Export
library(openxlsx)
wb = createWorkbook()
for (j in 1:length(hclust_output)) {
  addWorksheet(wb, names(hclust_output)[j])
  writeData(wb, names(hclust_output)[j], hclust_output[[j]])
}
saveWorkbook(wb, file = paste0(home, "/Results/single_algorithm/", algorithm, "/", 
                               algorithm, "_representative_pathways.xlsx"),
             overwrite = TRUE); rm(wb)

# Plot pathway heatmaps
hclust_pathway_plots_up = plot_pathway_heatmaps(gsea.lists = hclust_output[grepl("up", names(hclust_output))], 
                                                norm.expr = plotdata$RNAseq, 
                                                present_clusters = c("RWR-NF1", "RWR-NF2"),
                                                representative = TRUE, moic.res = plot_object,
                                                subtype_prefix = algorithm, n.path = 20, msigdb.path = MSIGDB.FILE,
                                                norm.method = "mean", dirct = "up",
                                                fig.name = "upregulated_pathway_heatmap",
                                                name = "GSVA scores",
                                                fig.path = paste0(home, "/Results/single_algorithm/", algorithm), 
                                                width = 15, height = 10, gsva.method = "gsva")

hclust_pathway_plots_down = plot_pathway_heatmaps(gsea.lists = hclust_output[grepl("down", names(hclust_output))], 
                                                  norm.expr = plotdata$RNAseq, 
                                                  present_clusters = c("RWR-NF1", "RWR-NF2"),
                                                  representative = TRUE, moic.res = plot_object,
                                                  subtype_prefix = algorithm, n.path = 20, msigdb.path = MSIGDB.FILE,
                                                  norm.method = "mean", dirct = "down",
                                                  fig.name = "downregulated_pathway_heatmap",
                                                  name = "GSVA scores",
                                                  fig.path = paste0(home, "/Results/single_algorithm/", algorithm), 
                                                  width = 15, height = 10, gsva.method = "gsva")

# Fraction Genome Altered ###
fga_df = readRDS("Resources/TCGA/fga_df.rds"); gc()

fga.RWRNF <- compFGA_optimized(moic.res     = plot_object,
                               segment      = fga_df,
                               iscopynumber = TRUE, 
                               test.method  = "nonparametric", # statistical testing method (Wilcoxon with asymptotic approximation. Consider Kruskall Wallis?)
                               fig.path     = paste0(home, "/Results/single_algorithm/", algorithm),
                               fig.name     = paste0("FGA_barplot_", algorithm),
                               prefix = algorithm,
                               width = 16,
                               ga_column = "ga", # genome altered column
                               clust.col = cluster_colors,
                               title = paste0(algorithm, " FGA plot: simple criteria"))

fga.RWRNF.COSMIC <- compFGA_optimized(moic.res     = plot_object,
                                      segment      = fga_df,
                                      iscopynumber = TRUE, 
                                      test.method  = "nonparametric", # statistical testing method (Wilcoxon with asymptotic approximation. Consider Kruskall Wallis?)
                                      fig.path     = paste0(home, "/Results/single_algorithm/", algorithm),
                                      fig.name     = paste0("COSMIC_criteria_FGA_barplot_", algorithm),
                                      prefix = algorithm,
                                      width = 16,
                                      ga_column = "COSMIC_ga", # genome altered column
                                      clust.col = cluster_colors,
                                      title = paste0(algorithm, " FGA plot: COSMIC criteria"))

rm(fga_df); gc()

# Evaluation #####
# Run Nearest Template Prediction in transNEO cohort ###
# Load transNEO data
transNEO_mm_inputs = readRDS("Resources/transNEO/transNEO_multimodal_inputs.rds")
transcr = readRDS("Resources/transNEO/log2.norm.counts.plus1_transNEO.rds")

# get as many templates as possible
dgea.marker.up_1000 <- runMarker_single_algorithm_no_export(algorithm_name = algorithm,
                                                            moic.res = plot_object,
                                                            n.marker = 1000,
                                                            dea.method    = "limma", # name of DEA method
                                                            prefix        = "dgea_", # MUST be the same of argument in runDEA()
                                                            dat.path      = paste0(home, "/Results/single_algorithm/", algorithm), # path of DEA files
                                                            p.cutoff      = 0.05, # p cutoff to identify significant DEGs
                                                            p.adj.cutoff  = 0.05, # padj cutoff to identify significant DEGs
                                                            norm.expr = plotdata$RNAseq,
                                                            dirct         = "up" # direction of dysregulation in expression
)

# 2. Down-regulated markers
dgea.marker.down_1000 <- runMarker_single_algorithm_no_export(algorithm_name = algorithm,
                                                              moic.res = plot_object,
                                                              n.marker = 1000,
                                                              dea.method    = "limma", # name of DEA method
                                                              prefix        = "dgea_", # MUST be the same of argument in runDEA()
                                                              dat.path      = paste0(home, "/Results/single_algorithm/", algorithm), # path of DEA files
                                                              p.cutoff      = 0.05, # p cutoff to identify significant DEGs
                                                              p.adj.cutoff  = 0.05, # padj cutoff to identify significant DEGs
                                                              norm.expr = plotdata$RNAseq,
                                                              dirct         = "down" # direction of dysregulation in expression
)

# Up-regulated expression features
RNGversion("4.2.2")
timestamp()
transNEO_ntp_expr_up = runNTP(
  expr = transcr,
  templates = dgea.marker.up_1000$templates,
  scaleFlag = TRUE,
  centerFlag = TRUE,
  nPerm = 10000,
  seed = 123,
  distance = "cosine", # default
  doPlot = TRUE,
  height = 8,
  width = 12,
  fig.path = paste0(home, "/Results/single_algorithm/", algorithm),
  fig.name = "ntp_expr_up_heatmap_transNEO")
timestamp() # 4 min

# down-regulated
RNGversion("4.2.2")
timestamp()
transNEO_ntp_expr_down = runNTP(
  expr = transcr,
  templates = dgea.marker.down_1000$templates,
  scaleFlag = TRUE,
  centerFlag = TRUE,
  nPerm = 10000,
  seed = 123,
  distance = "cosine", # default
  doPlot = TRUE,
  height = 8,
  width = 12,
  fig.path = paste0(home, "/Results/single_algorithm/", algorithm),
  fig.name = "ntp_expr_down_heatmap_transNEO")
timestamp() # 2.5 min

# Check concordance
expr_conc = as.data.frame(transNEO_ntp_expr_down$clust.res) %>%
  dplyr::rename(clust_down = clust) %>%
  inner_join(as.data.frame(transNEO_ntp_expr_up$clust.res) %>%
               dplyr::rename(clust_up = clust),
             by = "samID")

# This is counter-intuitive but due to opposite directions of deregulation this
# is how it works (perhaps this was expected)
expr_conc$agreement = ifelse(expr_conc$clust_down!=expr_conc$clust_up, "Yes", "No")
paste("Agremeent of NTP subtypes with respect to expression data from the external cohort is: ",
      length(which(expr_conc$agreement == "Yes"))/nrow(expr_conc)*100, "% (", nrow(expr_conc),
      " samples).")

# Compare clinical variables of interest across clusters
transNEO_var2comp = transNEO_mm_inputs$`Full pheno` %>%
  dplyr::select(LN.status.at.diagnosis, ER.status, HER2.status,
                Grade.pre.NAT, pCR.RD, Age, T.stage, PAM50, iC10,
                NAT.regimen, Chemo.cycles,
                aHER2.cycles, RCB.score, STAT1.gsva,
                GGI.gsva, ESC.gsva, TMB, HRD.sum, Donor.ID) %>%
  inner_join(expr_conc %>% dplyr::select(Donor.ID = samID, RWRNF = clust_up),
             by = "Donor.ID")
rownames(transNEO_var2comp) = transNEO_var2comp$Donor.ID
transNEO_var2comp = transNEO_var2comp %>% dplyr::select(-Donor.ID)

# Convert to factors
transNEO_var2comp$LN.status.at.diagnosis = factor(transNEO_var2comp$LN.status.at.diagnosis,
                                                  levels = c("NEG", "POS"),
                                                  labels = c("Negative", "Positive"))
transNEO_var2comp$ER.status = factor(transNEO_var2comp$ER.status,
                                     levels = c("NEG", "POS"),
                                     labels = c("Negative", "Positive"))
transNEO_var2comp$HER2.status = factor(transNEO_var2comp$HER2.status,
                                       levels = c("NEG", "POS"),
                                       labels = c("Negative", "Positive"))
transNEO_var2comp$Grade.pre.NAT = factor(transNEO_var2comp$Grade.pre.NAT,
                                         levels = c(1, 2, 3, 4),
                                         labels = c("Grade 1", "Grade 2", "Grade 3", "Grade 4"))
transNEO_var2comp$pCR.RD = factor(transNEO_var2comp$pCR.RD,
                                  levels = c("pCR", "RD"),
                                  labels = c("pCR", "Residual Disease"))
transNEO_var2comp$PAM50 = factor(transNEO_var2comp$PAM50,
                                 levels = c("Basal", "Her2", "LumB", "LumA", "Normal", "Unk"),
                                 labels = c("Basal-like", "HER2+", "Luminal B", "Luminal A",
                                            "Normal-like", "Unknown"))
transNEO_var2comp$iC10 = factor(transNEO_var2comp$iC10,
                                levels = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
                                labels = paste("iC", seq(1, 10, 1), sep = ""))

# Remove unknown levels for statistical tests
transNEO_var2comp_nonas = transNEO_var2comp
for (i in 1:ncol(transNEO_var2comp)) {
  nas = which(transNEO_var2comp[, i] == "Unknown")
  transNEO_var2comp_nonas[nas, i] = NA
  empties = which(transNEO_var2comp[, i] == "")
  transNEO_var2comp_nonas[empties, i] = NA
}
rm(nas, empties); gc()

transNEO_clincomp = compClinvar_single_algorithm(algorithm_name = algorithm,
                                                 moic.res = transNEO_ntp_expr_up,
                                                 var2comp = transNEO_var2comp_nonas,
                                                 strata = "RWRNF",
                                                 factorVars = c("ER.status", "HER2.status",
                                                                "NAT.regimen", 
                                                                "pCR.RD", "LN.status.at.diagnosis"),
                                                 nonnormalVars = c("Age",
                                                                   "RCB.score", "STAT1.gsva", "GGI.gsva",
                                                                   "ESC.gsva", "TMB", "HRD.sum",
                                                                   "Grade.pre.NAT", "Chemo.cycle", "aHER2.cycles"),
                                                 includeNA = FALSE,
                                                 doWord = TRUE,
                                                 tab.name = "transNEO_Summary_of_clinical_variables",
                                                 res.path = paste0(home, "/Results/single_algorithm/", algorithm, "/"),
                                                 output_pdf = TRUE,
                                                 pdf_level_col_width = c("7em", "10em"),
                                                 pdf_count_col_width = "10em",
                                                 pdf_pval_col_width = "3em",
                                                 pdf_test_col_width = "8em",
                                                 pdf_tab_font_size = 9)

transNEO_ntp_expr_up_ord = transNEO_ntp_expr_up
transNEO_ntp_expr_up_ord$clust.res$clust = gsub(algorithm, "", transNEO_ntp_expr_up_ord$clust.res$clust)
transNEO_ordinal_clincomp = compClinvar_ordinal_single_algorithm(algorithm_name = algorithm,
                                                                 moic.res = transNEO_ntp_expr_up_ord,
                                                                 var2comp = transNEO_var2comp_nonas %>%
                                                                   dplyr::select(Grade.pre.NAT, 
                                                                                 Chemo.cycles, 
                                                                                 aHER2.cycles,
                                                                                 RWRNF),
                                                                 strata = "RWRNF",
                                                                 ordinalVars = c("Grade.pre.NAT",
                                                                                 "Chemo.cycles",
                                                                                 "aHER2.cycles"),
                                                                 includeNA = FALSE,
                                                                 tab.name = "transNEO Summary of ordinal clinical variables",
                                                                 res.path = paste0(home, "/Results/single_algorithm/", algorithm, "/"),
                                                                 output_pdf = TRUE,
                                                                 pdf_template_loc = paste0(home, "/Scripts/automated_scripts/clincomp_template.Rmd"),
                                                                 pdf_level_col_width = c("7em", "10em"),
                                                                 pdf_count_col_width = "10em",
                                                                 pdf_pval_col_width = "3em",
                                                                 pdf_test_col_width = "8em",
                                                                 pdf_tab_font_size = 9)

# Run PAM ###
RNGversion("4.2.2")
set.seed(123)
transNEO_pam = runPAM_single_algorithm(algorithm_name = algorithm,
                                       train.expr = plotdata$RNAseq,
                                       moic.res   = plot_object,
                                       test.expr  = transcr)

# Check consistency across methods

# Get predictions for TCGA (discovery cohort)
RNGversion("4.2.2")
set.seed(123)
TCGA.ntp.pred = runNTP(expr = plotdata$RNAseq[, plot_object$clust.res$samID],
                       templates = dgea.marker.up_1000$templates, distance = "cosine",
                       doPlot = F, nPerm = 10000)

TCGA.pam.pred = runPAM_single_algorithm(algorithm_name = algorithm,
                                        train.expr = plotdata$RNAseq[, plot_object$clust.res$samID],
                                        moic.res = plot_object,
                                        test.expr = plotdata$RNAseq[, plot_object$clust.res$samID])

# consensus TCGA vs NTP TCGA
runKappa_single_algorithm(algorithm_name = algorithm,
                          subt1 = plot_object$clust.res$clust,
                          subt2 = gsub(algorithm, "", TCGA.ntp.pred$clust.res$clust),
                          subt1.lab = algorithm,
                          subt2.lab = "NTP TCGA",
                          height = 8,
                          width = 8,
                          fig.path = paste0(home, "/Results/single_algorithm/", algorithm),
                          fig.name = paste0("kappa_", algorithm, "_vs_NTP_TCGA"))

# consensus TCGA vs PAM TCGA
runKappa_single_algorithm(algorithm_name = algorithm,
                          subt1 = plot_object$clust.res$clust,
                          subt2 = gsub(algorithm, "", TCGA.pam.pred$clust.res$clust),
                          subt1.lab = algorithm,
                          subt2.lab = "PAM TCGA",
                          height = 8,
                          width = 8,
                          fig.path = paste0(home, "/Results/single_algorithm/", algorithm),
                          fig.name = paste0("kappa_", algorithm, "_vs_PAM_TCGA"))

# NTP transNEO vs PAM transNEO
runKappa_single_algorithm(algorithm_name = algorithm,
                          subt1 = as.numeric(gsub(algorithm, "",
                                                  transNEO_ntp_expr_up$clust.res$clust)),
                          subt2 = as.numeric(transNEO_pam$clust.res$clust),
                          subt1.lab = "transNEO NTP",
                          subt2.lab = "transNEO PAM",
                          height = 8,
                          width = 8,
                          fig.path = paste0(home, "/Results/single_algorithm/", algorithm),
                          fig.name = "kappa_NTP_vs_PAM_transNEO")

# Supplementary results #####

# Create subdirectory for supplementary plots
if (!dir.exists(paste0(home, "/Results/single_algorithm/", algorithm, "/Supplement"))) {
  dir.create(paste0(home, "/Results/single_algorithm/", algorithm, "/Supplement"))
}

# Setup for heatmaps
colors_heatmap = rev(colorRampPalette(viridisLite::magma(10))(255))
cluster_colors_heatmap = c("#2EC4B6", "#E71D36")
clust_annot_pheno = annCol %>% mutate(Sample.ID = rownames(.)) %>%
  inner_join(clust, by = "Sample.ID") %>%
  dplyr::rename(RWRNF = Cluster, samID = "Sample.ID")
rownames(clust_annot_pheno) = clust_annot_pheno$samID
afh_colnames = colnames(annCol)

# Prepare affinity matrices
aff_CNV = normalize_affinity_matrix(
  SNFtool::affinityMatrix(
    SNFtool::dist2(input$CNV,
                   input$CNV),
    K = optN, sigma = optSigma))
colnames(aff_CNV) = rownames(aff_CNV) = rownames(input$CNV)

aff_rna = normalize_affinity_matrix(
  SNFtool::affinityMatrix(
    SNFtool::dist2(input$RNAseq,
                   input$RNAseq),
    K = optN, sigma = optSigma))
colnames(aff_rna) = rownames(aff_rna) = rownames(input$RNA)

aff_miRNA = normalize_affinity_matrix(
  SNFtool::affinityMatrix(
    SNFtool::dist2(input$miRNA,
                   input$miRNA),
    K = optN, sigma = optSigma))
colnames(aff_miRNA) = rownames(aff_miRNA) = rownames(input$miRNA)

aff_Methyl = normalize_affinity_matrix(
  SNFtool::affinityMatrix(
    SNFtool::dist2(input$Methylation,
                   input$Methylation),
    K = optN, sigma = optSigma))
colnames(aff_Methyl) = rownames(aff_Methyl) = rownames(input$Methylation)

aff_SNPs = normalize_affinity_matrix(
  SNFtool::affinityMatrix(
    as.matrix(dist(as.matrix(input$SNPs),
                   as.matrix(input$SNPs),
                   method = "binary")),
    K = optN, sigma = optSigma))
colnames(aff_SNPs) = rownames(aff_SNPs) = rownames(input$SNPs)

aff_final = final_affinity_matrix

# CNV
create_MO_heatmap(matrix = aff_CNV, algorithm = algorithm, 
                  need.diag.zero = TRUE, 
                  clust_annot_pheno = clust_annot_pheno %>% dplyr::rename(`RWR-NF` = RWRNF),
                  afh_colnames = afh_colnames, 
                  colors = colors_heatmap,
                  annColors = annColors,
                  heatmap_title = "CNV first affinity heatmap",
                  cluster_colors = cluster_colors_heatmap,
                  legend_title = "Normalized affinity",
                  cluster_cols_flag = FALSE,
                  cluster_rows_flag = FALSE,
                  splits_flag = TRUE,
                  output_file_name = paste0(home, "/Results/single_algorithm/", algorithm, 
                                            "/Supplement/aff_CNV_heatmap.png"))

# RNAseq
create_MO_heatmap(matrix = aff_rna, algorithm = algorithm, 
                  need.diag.zero = TRUE, 
                  clust_annot_pheno = clust_annot_pheno %>% dplyr::rename(`RWR-NF` = RWRNF),
                  afh_colnames = afh_colnames, 
                  colors = colors_heatmap,
                  annColors = annColors,
                  heatmap_title = "RNAseq first affinity heatmap",
                  cluster_colors = cluster_colors_heatmap,
                  legend_title = "Normalized affinity",
                  cluster_cols_flag = FALSE,
                  cluster_rows_flag = FALSE,
                  splits_flag = TRUE,
                  output_file_name = paste0(home, "/Results/single_algorithm/", algorithm, 
                                            "/Supplement/aff_RNAseq_heatmap.png"))

# miRNA
create_MO_heatmap(matrix = aff_miRNA, algorithm = algorithm, 
                  need.diag.zero = TRUE, 
                  clust_annot_pheno = clust_annot_pheno %>% dplyr::rename(`RWR-NF` = RWRNF),
                  afh_colnames = afh_colnames, 
                  colors = colors_heatmap,
                  annColors = annColors,
                  heatmap_title = "miRNA first affinity heatmap",
                  cluster_colors = cluster_colors_heatmap,
                  legend_title = "Normalized affinity",
                  cluster_cols_flag = FALSE,
                  cluster_rows_flag = FALSE,
                  splits_flag = TRUE,
                  output_file_name = paste0(home, "/Results/single_algorithm/", algorithm,
                                            "/Supplement/aff_miRNA_heatmap.png"))

# Methylation
create_MO_heatmap(matrix = aff_Methyl, algorithm = algorithm, 
                  need.diag.zero = TRUE, 
                  clust_annot_pheno = clust_annot_pheno %>% dplyr::rename(`RWR-NF` = RWRNF),
                  afh_colnames = afh_colnames, 
                  colors = colors_heatmap,
                  annColors = annColors,
                  heatmap_title = "Methylation first affinity heatmap",
                  cluster_colors = cluster_colors_heatmap,
                  legend_title = "Normalized affinity",
                  cluster_cols_flag = FALSE,
                  cluster_rows_flag = FALSE,
                  splits_flag = TRUE,
                  output_file_name = paste0(home, "/Results/single_algorithm/", algorithm, 
                                            "/Supplement/aff_Methylation_heatmap.png"))

# SNPs
create_MO_heatmap(matrix = aff_SNPs, algorithm = algorithm, 
                  need.diag.zero = TRUE, 
                  clust_annot_pheno = clust_annot_pheno %>% dplyr::rename(`RWR-NF` = RWRNF),
                  afh_colnames = afh_colnames, 
                  colors = colors_heatmap,
                  annColors = annColors,
                  heatmap_title = "SNPs first affinity heatmap",
                  cluster_colors = cluster_colors_heatmap,
                  legend_title = "Normalized affinity",
                  cluster_cols_flag = FALSE,
                  cluster_rows_flag = FALSE,
                  splits_flag = TRUE,
                  output_file_name = paste0(home, "/Results/single_algorithm/", algorithm,
                                            "/Supplement/aff_SNPs_heatmap.png"))

# Final affinity matrix
create_MO_heatmap(matrix = final_affinity_matrix, algorithm = algorithm, 
                  need.diag.zero = TRUE, 
                  clust_annot_pheno = clust_annot_pheno %>% dplyr::rename(`RWR-NF` = RWRNF),
                  afh_colnames = afh_colnames, 
                  colors = colors_heatmap,
                  annColors = annColors,
                  heatmap_title = "Final affinity heatmap",
                  cluster_colors = cluster_colors_heatmap,
                  legend_title = "Normalized affinity",
                  cluster_cols_flag = FALSE,
                  cluster_rows_flag = FALSE,
                  splits_flag = TRUE,
                  output_file_name = paste0(home, "/Results/single_algorithm/", algorithm,
                                            "/Supplement/aff_final_affinity_heatmap.png"))

# Final affinity matrix with clustered rows and columns
create_MO_heatmap(matrix = final_affinity_matrix, algorithm = algorithm, 
                  need.diag.zero = TRUE, 
                  clust_annot_pheno = clust_annot_pheno %>% dplyr::rename(`RWR-NF` = RWRNF),
                  afh_colnames = afh_colnames, 
                  colors = colors_heatmap,
                  annColors = annColors,
                  heatmap_title = "Final affinity heatmap",
                  cluster_colors = cluster_colors_heatmap,
                  legend_title = "Normalized affinity",
                  cluster_cols_flag = TRUE,
                  cluster_rows_flag = TRUE,
                  splits_flag = FALSE,
                  output_file_name = paste0(home, "/Results/single_algorithm/", algorithm,
                                            "/Supplement/hclust_aff_final_affinity_heatmap.png"))

# PCA ###
# CNV
pca_from_sim_matrix(sim_matrix = aff_CNV, algorithm = algorithm, 
                    clust_res = clust_annot_pheno %>% dplyr::select(samID, `RWR-NF` = RWRNF),
                    cluster_colors = cluster_colors_heatmap, 
                    output_path = paste0(home, "/Results/single_algorithm/", algorithm, "/Supplement"), 
                    title_add = "CNV")

# RNAseq
pca_from_sim_matrix(sim_matrix = aff_rna, algorithm = algorithm, 
                    clust_res = clust_annot_pheno %>% dplyr::select(samID, `RWR-NF` = RWRNF),
                    cluster_colors = cluster_colors_heatmap, 
                    output_path = paste0(home, "/Results/single_algorithm/", algorithm, "/Supplement"), 
                    title_add = "RNAseq")

# miRNA
pca_from_sim_matrix(sim_matrix = aff_miRNA, algorithm = algorithm, 
                    clust_res = clust_annot_pheno %>% dplyr::select(samID, `RWR-NF` = RWRNF),
                    cluster_colors = cluster_colors_heatmap, 
                    output_path = paste0(home, "/Results/single_algorithm/", algorithm, "/Supplement"), 
                    title_add = "miRNA")

# Methylation
pca_from_sim_matrix(sim_matrix = aff_Methyl, algorithm = algorithm, 
                    clust_res = clust_annot_pheno %>% dplyr::select(samID, `RWR-NF` = RWRNF),
                    cluster_colors = cluster_colors_heatmap, 
                    output_path = paste0(home, "/Results/single_algorithm/", algorithm, "/Supplement"), 
                    title_add = "Methylation")

# SNPs
pca_from_sim_matrix(sim_matrix = aff_SNPs, algorithm = algorithm, 
                    clust_res = clust_annot_pheno %>% dplyr::select(samID, `RWR-NF` = RWRNF),
                    cluster_colors = cluster_colors_heatmap, 
                    output_path = paste0(home, "/Results/single_algorithm/", algorithm, "/Supplement"), 
                    title_add = "SNPs")

# Final affinity
pca_from_sim_matrix(sim_matrix = aff_final, algorithm = algorithm, 
                    clust_res = clust_annot_pheno %>% dplyr::select(samID, `RWR-NF` = RWRNF),
                    cluster_colors = cluster_colors_heatmap, 
                    output_path = paste0(home, "/Results/single_algorithm/", algorithm, "/Supplement"), 
                    title_add = "Fusion")

# Setup for barcharts ###
# Stage
scale_fill_stage = scale_fill_manual(values = c(`Stage I` = "#00C9FF", 
                                                `Stage II` = "#099CF5", 
                                                `Stage III` = "#097BF5", 
                                                `Stage IV` = "#0B5684", 
                                                `Unknown` = "grey40"))

# Lymph node status
scale_fill_lymph_node_status = scale_fill_manual(values = c(No = "grey75", 
                                                            Yes = "#4A0558", 
                                                            Unknown = "grey40"))

# ER status
scale_fill_ER_status = scale_fill_manual(values = c(Negative = "#C11D9C", 
                                                    Positive = "#0F1682", 
                                                    Unknown = "grey40"))

# PR status
scale_fill_PR_status = scale_fill_manual(values = c(Indeterminate = "aliceblue", 
                                                    Positive = "dodgerblue4", 
                                                    Negative = "#F0C6C3", 
                                                    Unknown = "grey40"))

# HER2 status
scale_fill_HER2_status = scale_fill_manual(values = c(Negative = "#0B9EF8", 
                                                      Positive = "#560DA7", 
                                                      Indeterminate = "mistyrose1", 
                                                      Equivocal = "hotpink4", 
                                                      Unknown = "grey40"))

# Vital status
scale_fill_vital_status = scale_fill_manual(values = c(Alive = "lightpink1", 
                                                       Dead = "black", 
                                                       Unknown = "grey40"))

# Ethnicity
scale_fill_ethnicity = scale_fill_manual(values = c(`Hispanic or latino` = "#E58606", 
                                                    `Not hispanic or latino` = "#24796C", 
                                                    Unknown = "grey40"))

# Race
scale_fill_race = scale_fill_manual(values = c(`American indian or alaska native` = "#E73F74", 
                                               Asian = "#3969AC", 
                                               `Black or african american` = "#666666", 
                                               White = "beige", 
                                               Unknown = "grey40"))

# Metastasis
scale_fill_metastasis = scale_fill_manual(values = c(Yes = "deeppink4", 
                                                     No = "cadetblue2", 
                                                     Unknown = "grey40"))

# Histology
scale_fill_histology = scale_fill_manual(values = c(`Infiltrating Carcinoma NOS` = "#88CCEE", 
                                                    `Infiltrating Ductal Carcinoma` = "#CC6677", 
                                                    `Infiltrating Lobular Carcinoma` = "#DDCC77", 
                                                    `Medullary Carcinoma` = "#117733", 
                                                    `Metaplastic Carcinoma` = "#332288", 
                                                    Mixed = "#AA4499", 
                                                    `Mucinous Carcinoma` = "#44AA99", 
                                                    Other = "#999933", 
                                                    Unknown = "grey40"))

# Menopausal status
scale_fill_menopausal_status = scale_fill_manual(values = c(Indeterminate = "mistyrose2", 
                                                            `Pre-menopausal` = "#FAA476", 
                                                            Perimenopausal = "#DC3977", 
                                                            `Post-menopausal` = "#7C1D6F", 
                                                            Unknown = "grey40"))

# Combine all scales into a list
barchart_scales = list(scale_fill_stage, scale_fill_lymph_node_status, scale_fill_ER_status, 
                       scale_fill_PR_status, scale_fill_HER2_status, scale_fill_vital_status, 
                       scale_fill_ethnicity, scale_fill_race, scale_fill_metastasis, 
                       scale_fill_histology, scale_fill_menopausal_status)

# Name the scales accordingly
names(barchart_scales) = c("Stage", "Lymph node status", "ER status", "PR status", "HER2 status", 
                           "Vital status", "Ethnicity", "Race", "Metastasis", "Histology", 
                           "Menopausal status")

# Chi-square tests ###
# Bias-corrected Cramer's V calculation using package rcompanion:
unbiased.cv.test = function(x, string, digits = 3) {
  CV = rcompanion::cramerV(x, bias.correct = TRUE)
  return(list(text = paste0("Bias-corrected Cramer's V / Phi for ", 
                            string, ": ", round(as.numeric(CV), digits)),
              value = round(as.numeric(CV), digits)))
}

clust_annot_pheno_nonas = clust_annot_pheno
for(i in 1:ncol(clust_annot_pheno_nonas)) {
  clust_annot_pheno_nonas[, i] = as.character(clust_annot_pheno_nonas[, i])
  nas = which(clust_annot_pheno_nonas[, i] == "Unknown")
  clust_annot_pheno_nonas[nas, i] = NA
  clust_annot_pheno_nonas[, i] = factor(clust_annot_pheno_nonas[, i])
}
rm(nas); gc()

voi = colnames(clust_annot_pheno_nonas)[1:11]
clust_annot_pheno_nonas = clust_annot_pheno_nonas %>% dplyr::rename(`RWR-NF` = RWRNF)
output = as.data.frame(matrix(NA, nrow = 0, ncol = 4))
for (v in 1:length(voi)){
  keepers = which(!is.na(clust_annot_pheno_nonas[, voi[v]]))
  test = suppressWarnings(chisq.test(table(clust_annot_pheno_nonas[keepers, algorithm], 
                                           clust_annot_pheno_nonas[keepers, voi[v]])))
  chifit_p = test$p.value
  chifit_xsq = test$statistic
  chifit_cv = suppressWarnings(unbiased.cv.test(table(clust_annot_pheno_nonas[keepers, algorithm], 
                                                      clust_annot_pheno_nonas[keepers, voi[v]]),
                                                string = voi[v],
                                                digits = 3)$value)
  comparison = paste0(voi[v], " vs ", algorithm, " cluster")
  output = rbind(output, c(comparison, chifit_p, chifit_xsq, chifit_cv))
  rm(test, comparison, chifit_p, chifit_xsq, chifit_cv, keepers)
}
colnames(output) = c("Comparison", "p-value", "Statistic", "Cramer's V")

rm(v); gc()
openxlsx::write.xlsx(output, 
                     paste0(home, 
                            "/Results/single_algorithm/", algorithm, "/Supplement/Chisq_tests.xlsx"),
                     overwrite = TRUE)

# Bar chart generation
RWRNF_barcharts = list()
plotdata_bar = clust_annot_pheno_nonas
plotdata_bar[[algorithm]] = factor(plotdata_bar[[algorithm]])
for (i in 1:length(voi)) {
  chifit = output
  loc = which(grepl(voi[i], chifit$Comparison))
  chifit = chifit[loc, ]
  RWRNF_barcharts[[i]] = create_annot_barchart(plotdata = plotdata_bar, fill = voi[i],
                                               chifit = chifit,
                                               na.action = "na.omit",
                                               algorithm = algorithm,
                                               barchart_ylim = 650,
                                               text_y = 600, rect_ymin = 500,
                                               rect_ymax = 620, x_annot = 1.5,
                                               v_gap = 35, rect_xmin = 1,
                                               rect_xmax = 2, 
                                               annot_text_size = 2.25,
                                               legend.text.size = 5,
                                               x.axis.text.size = 5) +
    barchart_scales[[voi[i]]]
  print(RWRNF_barcharts[[i]])
  ggsave(filename = paste0(algorithm, "_", voi[i], "_barchart.png"),
         path = paste0(home, 
                       "/Results/single_algorithm/", algorithm, "/Supplement"), 
         width = 2320, height = 2320, device = 'png', units = "px",
         dpi = 700)
  dev.off()
}
names(RWRNF_barcharts) = voi
rm(loc, chifit)

# Multiplot (PNG) - bar charts
library(ggpubr)
ggarrange(RWRNF_barcharts[[1]], RWRNF_barcharts[[2]], RWRNF_barcharts[[3]],
          RWRNF_barcharts[[4]], RWRNF_barcharts[[5]], RWRNF_barcharts[[6]],
          RWRNF_barcharts[[7]], RWRNF_barcharts[[8]], RWRNF_barcharts[[9]],
          RWRNF_barcharts[[10]], RWRNF_barcharts[[11]],
          ncol = 3, nrow = 4, labels = c("A", "B", "C", "D", "E", "F", "G", "H",
                                         "I", "J", "K"),
          font.label = list(size = 8, face = "bold", color ="black"))
ggsave(filename = paste0("Multiplot_", algorithm, "_barcharts.png"),
       path = paste0(home, 
                     "/Results/single_algorithm/", algorithm, "/Supplement"), 
       width = 7000, height = 8000, device = 'png', units = "px",
       dpi = 700)
dev.off()

# Just significant ones now
RWRNF_barcharts_sig = list()
plotdata_bar_sig = clust_annot_pheno_nonas %>% dplyr::select(`RWR-NF`, Race, Histology, 
                                                             `ER status`, `PR status`,
                                                             `Menopausal status`, Stage)
plotdata_bar_sig$`RWR-NF` = factor(plotdata_bar_sig$`RWR-NF`)
voi_sig = setdiff(colnames(plotdata_bar_sig), algorithm)
for (i in 1:length(voi_sig)) {
  chifit = output
  loc = which(grepl(voi_sig[i], chifit$Comparison))
  chifit = chifit[loc, ]
  RWRNF_barcharts_sig[[i]] = create_annot_barchart(plotdata = plotdata_bar_sig, fill = voi_sig[i],
                                                   chifit = chifit,
                                                   na.action = "na.omit",
                                                   algorithm = algorithm,
                                                   barchart_ylim = 650,
                                                   text_y = 600, rect_ymin = 500,
                                                   rect_ymax = 620, x_annot = 1.5,
                                                   v_gap = 35, rect_xmin = 1,
                                                   rect_xmax = 2, 
                                                   annot_text_size = 2.25,
                                                   legend.text.size = 5,
                                                   x.axis.text.size = 5) +
    barchart_scales[[voi_sig[i]]]
  print(RWRNF_barcharts_sig[[i]])
  ggsave(filename = paste0("sig_", algorithm, "_", voi_sig[i], "_barchart.png"),
         path = paste0(home, 
                       "/Results/single_algorithm/", algorithm, "/Supplement"), 
         width = 2320, height = 2320, device = 'png', units = "px",
         dpi = 700)
  dev.off()
}
names(RWRNF_barcharts_sig) = voi_sig
rm(loc, chifit)

# Multiplot (PNG) - bar charts
ggarrange(RWRNF_barcharts_sig[[1]], RWRNF_barcharts_sig[[2]], RWRNF_barcharts_sig[[3]],
          RWRNF_barcharts_sig[[4]], RWRNF_barcharts_sig[[5]], RWRNF_barcharts_sig[[6]],
          ncol = 2, nrow = 3, labels = c("A", "B", "C", "D", "E", "F"),
          font.label = list(size = 8, face = "bold", color ="black"))
ggsave(filename = paste0("sig_Multiplot_", algorithm, "_barcharts.png"),
       path = paste0(home, 
                     "/Results/single_algorithm/", algorithm, "/Supplement"), 
       width = 5500, height = 8500, device = 'png', units = "px",
       dpi = 700)
dev.off()

# Sunburst plot ###
library(plotly)
Pheno_sunburst_RWRNF = clust_annot_pheno
Pheno_sunburst_RWRNF$`ER status` = gsub("Unknown", "Unkn ER status", Pheno_sunburst_RWRNF$`ER status`)
Pheno_sunburst_RWRNF$`ER status` = gsub("Positive", "ER+", Pheno_sunburst_RWRNF$`ER status`)
Pheno_sunburst_RWRNF$`ER status` = gsub("Negative", "ER-", Pheno_sunburst_RWRNF$`ER status`)
Pheno_sunburst_RWRNF$`Menopausal status` = gsub("Indeterminate", "Indet", Pheno_sunburst_RWRNF$`Menopausal status`)
Pheno_sunburst_RWRNF$`Menopausal status` = gsub("Pre-menopausal", "Pre", Pheno_sunburst_RWRNF$`Menopausal status`)
Pheno_sunburst_RWRNF$`Menopausal status` = gsub("Perimenopausal", "Peri", Pheno_sunburst_RWRNF$`Menopausal status`)
Pheno_sunburst_RWRNF$`Menopausal status` = gsub("Post-menopausal", "Post", Pheno_sunburst_RWRNF$`Menopausal status`)
Pheno_sunburst_RWRNF$`Menopausal status` = gsub("Unknown", "Unkn Meno", Pheno_sunburst_RWRNF$`Menopausal status`)
Pheno_sunburst_RWRNF$Stage = gsub("Unknown", "Unkn Stage", Pheno_sunburst_RWRNF$Stage)
Pheno_sunburst_RWRNF = Pheno_sunburst_RWRNF %>%
  dplyr::select(RWRNF, `ER status`, `Menopausal status`, Stage) %>%
  group_by(RWRNF, `ER status`, `Menopausal status`, Stage) %>%
  summarise(Counts = n()) %>%
  as.data.frame()

sunburst_coloring_RWRNF = data.frame(stringsAsFactors = FALSE,
                                     colors = tolower(gplots::col2hex(c("#2EC4B6", "#E71D36", 
                                                                        "#C11D9C", "#0F1682",  "grey40",
                                                                        "mistyrose2", "#FAA476", "#DC3977", 
                                                                        "#7C1D6F", "grey40",
                                                                        "#00C9FF", "#099CF5", "#097BF5", 
                                                                        "#0B5684", "grey40"))),
                                     labels = c("RWR-NF1", "RWR-NF2",
                                                "ER-", "ER+", "Unkn ER status",
                                                "Indet", "Pre", "Peri", "Post", "Unkn Meno",
                                                "Stage I", "Stage II", "Stage III", "Stage IV",
                                                "Unkn Stage"))

sunburstDF_RWRNF = as.sunburstDF(Pheno_sunburst_RWRNF, value_column = "Counts", add_root = FALSE) %>%
  inner_join(sunburst_coloring_RWRNF, by = "labels")

pie_RWRNF = plot_ly() %>%
  add_trace(ids = sunburstDF_RWRNF$ids, labels= sunburstDF_RWRNF$labels, 
            parents = sunburstDF_RWRNF$parents, 
            values= sunburstDF_RWRNF$values, type='sunburst', branchvalues = 'total',
            insidetextorientation='radial', maxdepth = 5,
            marker = list(colors = sunburstDF_RWRNF$colors)) %>%
  layout(
    grid = list(columns =1, rows = 1),
    margin = list(l = 0, r = 0, b = 0, t = 0)
  )
pie_RWRNF
rm(Pheno_sunburst_RWRNF, sunburstDF_RWRNF, sunburst_coloring_RWRNF, pie_RWRNF); gc()

# Graphs ###
library(igraph)

aff_CNV_S = calculate_S(aff_CNV)
aff_rna_S = calculate_S(aff_rna)
aff_miRNA_S = calculate_S(aff_miRNA)
aff_Methyl_S = calculate_S(aff_Methyl)
aff_SNPs_S = calculate_S(aff_SNPs)
aff_final_S = calculate_S(aff_final)

list_aff_S = list(aff_CNV_S, aff_rna_S, aff_miRNA_S, aff_Methyl_S, 
                  aff_SNPs_S, aff_final_S)
names(list_aff_S) = c(paste0("CNV Original Affinity ", optN, "-NN Graph"),
                      paste0("RNAseq Original Affinity ", optN, "-NN Graph"),
                      paste0("miRNA Original Affinity ", optN, "-NN Graph"),
                      paste0("Methylation Original Affinity ", optN, "-NN Graph"),
                      paste0("SNPs Original Affinity ", optN, "-NN Graph"),
                      paste0("Final Fused Affinity ", optN, "-NN Graph"))

quantile_thresh = 0.75
for (i in seq_along(list_aff_S)) {
  g <- graph_from_adjacency_matrix(
    list_aff_S[[i]],
    mode     = "max",
    weighted = TRUE,
    diag     = FALSE)
  
  # drop the zero-weight edges
  g <- delete_edges(g, E(g)[E(g)$weight == 0])
  
  # Edge threshold for drawing
  thresh      <- quantile(E(g)$weight, quantile_thresh)          # 4th quartile
  keep_edge   <- E(g)$weight >= thresh               # logical mask
  
  ## edge-specific plotting attributes
  E(g)$plot_width  <- ifelse(keep_edge,
                             sqrt(E(g)$weight)*10,      # visible edges
                             0)                      # invisible edges
  E(g)$plot_color  <- ifelse(keep_edge,
                             "gray85",               # visible color
                             NA)                     # NA 
  
  nodes_data <- data.frame(name = V(g)$name) %>%
    inner_join(clust_annot_pheno %>% select(samID, `RWR-NF` = RWRNF),
               by = c("name" = "samID"))
  
  nodes_data[[algorithm]] <- as.factor(nodes_data[[algorithm]])
  V(g)$`RWR-NF` <- nodes_data[[algorithm]]
  
  V(g)$color <- fifelse(V(g)$`RWR-NF` == paste0(algorithm, "1"), "#2EC4B6", "#E71D36")
  
  png(paste0(home,
             "/Results/single_algorithm/", algorithm, "/Supplement/",
             names(list_aff_S)[i], ".png"),
      width = 6000, height = 6000, res = 700)
  
  par(mar = c(2, 2, 2, 5))
  
  plot(g,
       layout       = layout_with_fr(g),
       vertex.color = V(g)$color,
       vertex.size  = 4,
       vertex.label = NA,
       edge.width   = E(g)$plot_width,
       edge.color   = E(g)$plot_color,
       main         = "")
  
  title(main = names(list_aff_S)[i], cex.main = 1.7)
  
  legend("bottomright",
         title  = "Node Color Legend",
         legend = paste0(algorithm, 1:2),
         fill   = cluster_colors_heatmap,
         cex    = 0.7,
         box.lwd = 1)
  
  dev.off()
}
rm(g, nodes_data)

# Define the affinity matrices for each modality
# Function to determine the important modalities for each edge
determine_edge_support <- function(edge_weights) {
  sorted_weights <- sort(edge_weights, decreasing = TRUE)
  # If the highest weight is more than 10% greater than all others, it is supported by a single modality
  if (sorted_weights[1] > sorted_weights[2] * 1.1) {
    return(which(edge_weights == sorted_weights[1]))
  }
  # If the difference between the two highest weights is less than 10%, it is supported by those two modalities
  else if (sorted_weights[1] <= sorted_weights[2] * 1.1 && sorted_weights[2] > sorted_weights[3] * 1.1) {
    return(which(edge_weights >= sorted_weights[2]))
  }
  # If the difference between all weights is less than 10%, it is supported by all modalities
  else {
    return(which(edge_weights >= sorted_weights[1] * 0.9))
  }
}

# Prepare the graph object for the final affinity matrix
g <- graph_from_adjacency_matrix(aff_final_S, mode = "max", weighted = TRUE, diag = FALSE)
g <- delete_edges(g, E(g)[E(g)$weight == 0])

# Calculate the weights of each edge in individual networks
edge_weights_list <- list(aff_CNV_S, aff_rna_S, aff_miRNA_S, aff_Methyl_S, aff_SNPs_S)
edge_support_list <- lapply(E(g), function(e) {
  from <- ends(g, e)[1]
  to <- ends(g, e)[2]
  sapply(edge_weights_list, function(mat) mat[from, to])
})

# Assign color to each edge based on the modalities that support it
edge_colors <- c("#FF5733", "#33FF57", "#3357FF", "#FF33A1", "#F3FF33")  # Define unique colors for each modality
combinations_colors <- list()

for (i in 1:length(E(g))) {
  supported_modalities <- determine_edge_support(edge_support_list[[i]])
  if (length(supported_modalities) == 5) {
    combinations_colors[[i]] <- "black"  # Assign black color if all modalities support the edge
  } else if (length(supported_modalities) == 1) {
    combinations_colors[[i]] <- edge_colors[supported_modalities]
  } else {
    combination <- paste(supported_modalities, collapse = " + ")
    if (!is.null(combinations_colors[[combination]])) {
      combinations_colors[[i]] <- combinations_colors[[combination]]
    } else {
      new_color <- colorRampPalette(edge_colors[supported_modalities])(1)
      combinations_colors[[combination]] <- new_color
      combinations_colors[[i]] <- new_color
    }
  }
}

E(g)$color <- unlist(combinations_colors)

# Prepare the node data
nodes_data <- data.frame(name = V(g)$name) %>%
  inner_join(clust_annot_pheno %>% dplyr::select(samID, `RWR-NF` = RWRNF), by = c("name" = "samID"))
nodes_data[[algorithm]] <- as.factor(nodes_data[[algorithm]])
V(g)$`RWR-NF` <- nodes_data[[algorithm]]

# Set color based on RWRNF
V(g)$color <- fifelse(V(g)$`RWR-NF` == paste0(algorithm, "1"), "#2EC4B6", "#E71D36")

# Draw the graph for the final affinity matrix
png(paste0(home, "/Results/single_algorithm/", algorithm, "/Supplement/", "colored_Final_Fused_Affinity_", optN, "-NN_Graph.png"),
    width = 6000, height = 6000, res = 700)

par(mar = c(0, 0, 0, 17))  # Adjust right margin to make enough space for the legend

# Plot the graph
plot(g, vertex.color = V(g)$color,
     edge.width = E(g)$width * 0.1,  # Make the edges thinner
     vertex.size = 3,  # Make the nodes smaller
     vertex.label = NA, 
     edge.color = adjustcolor(E(g)$color, alpha.f = 0.35),  # Add transparency to edges
     layout = layout_with_fr(g, niter = 1000))

# Add a legend for the edge colors
legend(x = 1.2, y = 1,  # Manually adjust the position of the legend to the right of the plot
       xpd = TRUE,  # Allow legend to be drawn outside the plot area
       title="Edge Color Legend",
       legend=c("RNA-seq", "miRNA", "Methylation", "CNV", "SNPs", 
                "RNA-seq + miRNA", "RNA-seq + Methylation", "RNA-seq + CNV", "RNA-seq + SNPs", 
                "miRNA + Methylation", "miRNA + CNV", "miRNA + SNPs", 
                "Methylation + CNV", "Methylation + SNPs", 
                "CNV + SNPs", 
                "RNA-seq + miRNA + Methylation", "RNA-seq + miRNA + CNV", "RNA-seq + miRNA + SNPs", 
                "RNA-seq + Methylation + CNV", "RNA-seq + Methylation + SNPs", 
                "RNA-seq + CNV + SNPs", 
                "miRNA + Methylation + CNV", "miRNA + Methylation + SNPs", 
                "miRNA + CNV + SNPs", 
                "Methylation + CNV + SNPs", 
                "RNA-seq + miRNA + Methylation + CNV", "RNA-seq + miRNA + Methylation + SNPs", 
                "RNA-seq + miRNA + CNV + SNPs", "RNA-seq + Methylation + CNV + SNPs", 
                "miRNA + Methylation + CNV + SNPs", 
                "RNA-seq + miRNA + Methylation + CNV + SNPs"),
       fill=c("#33FF57", "#3357FF", "#FF5733", "#FF33A1", "#F3FF33", 
              "#abcdef", "#123456", "#654321", "#fedcba", 
              "#a1b2c3", "#b2c3d4", "#c3d4e5", 
              "#d4e5f6", "#e5f6a7", "#f6a7b8", 
              "#aabbcc", "#bbccdd", "#ccddee", 
              "#ddeeff", "#eeffaa", "#ffaabb", 
              "#aaffcc", "#bbffdd", "#ccffee", 
              "#ffccaa", 
              "#112233", "#223344", "#334455", 
              "#445566", "#556677", "#667788", 
              "black"),
       cex=0.7,
       box.lwd=1)

dev.off()

# Wrap up #####
# Hyperparameters and other parameters ###
hyperparameters = list(num_neighbors_min = min(num_neighbors_range),
                       num_neighbors_max = max(num_neighbors_range),
                       num_neighbors_step = neighbor_step,
                       sigma_min = min(sigma_range),
                       sigma_max = max(sigma_range),
                       gamma = RWRNF_gamma_fixed,
                       RWRNF_num_neighbors = RWRNF_num_neighbors_fixed,
                       RWRNF_alpha = RWRNF_alpha_fixed,
                       RWRNF_beta = RWRNF_beta_fixed,
                       sigma_step = sigma_step,
                       optimal_N = optN,
                       optimal_sigma = optSigma,
                       conclusion = conclusion, # if there is agreement, np_conclusion can also be used
                       max_iter = RWR_iteration_max,
                       feature_ranks_text = feature_ranks_text,
                       optk = optk
)

# Put all parameters in a list
params = list(algorithm = algorithm, data_source = data_source, data_types = data_types,
              evaluation_source = evaluation_source, title = title, subtitle = subtitle,
              description = description, in_a_nutshell = in_a_nutshell, optk_text = optk_text,
              citation = citation, NMI_to_MOVICS = NMI_to_MOVICS, ARI_to_MOVICS = ARI_to_MOVICS,
              hyperparameters = hyperparameters, ground_truth_k = ground_truth_k,
              transNEO_var2comp = transNEO_var2comp,
              sessionInfo = sessionInfo(), home = home)


# Render the R Markdown document with the parameters
rmarkdown::render(paste0(getwd(), "/Results/single_algorithm/", algorithm,
                         "/", algorithm, "_report.Rmd"), 
                  params = params, 
                  output_file = paste0(home, "/Results/single_algorithm/", 
                                       algorithm, "/", algorithm, "_report_",
                                       data_source, "_",
                                       data_types, "_eval_on_", evaluation_source,
                                       ".html"))

# Export session info as .txt
writeLines(capture.output(sessionInfo()), paste0("sessionInfo/",
                                                 algorithm, "_", data_source, "_",
                                                 data_types, "_eval_on_", evaluation_source,
                                                 "_sessionInfo.txt"))

# Save environment
save.image(paste0(home, "/Results/single_algorithm/", 
                  algorithm, "/", algorithm, "_", data_source, "_",
                  data_types, "_eval_on_", evaluation_source,
                  "_env.RData"))
