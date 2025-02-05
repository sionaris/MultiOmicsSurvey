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
algorithm = "ab-SNF"
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

# However abSNF prefers features in columns so we transpose the matrices.

# Extract the names of the modalities that will be used
modalities = unlist(strsplit(data_types, "-"))

# Replace with TRUE wherever features are in rows
features_in_rows = rep(TRUE, length(modalities))

# Run algorithm #####
library(abSNF)

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
COSMIC_BC_drivers = read.csv("Resources/COSMIC_CGC_Breast_somatic.csv")$Gene.Symbol %>%
  as.character()

# Setup ###
# In ab-SNF continuous features are firstly weighted by using associations statistics
# with an outcome of interest. Here, we want to discover patterns in an outcome-agnostic
# way, therefore, the features will be ranked by MAD (variance is equal for all)

# For binary modalities, e.g. SNPs, the authors suggest weighting based on whether the
# associated gene is a census driver for the cancer of interest or not. 
# They suggest to use weight = 1 for drivers and 0 for the rest.
# Here we will use 0.8 for drivers and 0.2 for the rest.

# Returns an vector with dimensions equal to the number of features of each modality
# of normalized variance-based weights for the columns.
mad_based_weights <- function(X) {
  v <- apply(X, 2, mad, na.rm = TRUE) # MAD of each feature (features in columns)
  v[!is.finite(v)] <- 0
  if (sum(v) == 0) {
    # all features are constant, fallback: uniform weighting or return zeros
    w <- rep(0, length(v))
  } else {
    w <- v / sum(v) # normalize so sum of all weights is 1
  }
  return(w)
}

# Weighted Euclidean distance, using precomputed per-feature weights
weighted_euclidean_dist <- function(X, w) {
  Xw <- sweep(X, 2, sqrt(w), `*`) # multiply each feature column by sqrt(w_k)
  dmat <- abSNF::dist2(as.matrix(Xw), as.matrix(Xw))
  return(dmat)
}

# For binary data, we weigh 0.8/0.2 drivers vs. non-drivers
weighted_hamming_fraction <- function(binary_mat, features, drivers,
                                      driver_weight, nondriver_weight) {
  w <- ifelse(features %in% drivers, driver_weight, nondriver_weight)
  
  # 2) We will compute for each pair (i, j):
  #       d(i, j) = [ sum_{k}  w_k * |SNP_mat[i,k] - SNP_mat[j,k]| ] / sum_{k} w_k
  #    This yields the *fraction* of mismatches, weighted by w_k.
  
  total_w <- sum(w)
  N <- nrow(binary_mat)
  dmat <- matrix(0, nrow = N, ncol = N)
  
  # The following is time consuming
  for (i in seq_len(N)) {
    for (j in seq_len(N)) {
      # Weighted sum of absolute differences across features
      mismatches_ij <- sum(w * abs(binary_mat[i, ] - binary_mat[j, ]))
      # Divide by total weight to get fraction of weighted mismatches
      dmat[i, j] <- mismatches_ij / total_w
    }
  }
  
  return(dmat)
}

# Modality types
continuous = c("RNAseq", "CNV", "Methylation", "miRNA")
categorical = c("SNPs")

# Calculate the pair-wise distance (Euclidean for continuous modalities)
input_dists = lapply(input[continuous], function(x) {
  weights = mad_based_weights(x)
  weighted_dists = weighted_euclidean_dist(x, weights)
  return(weighted_dists)
})

# Binary for SNPs (see ?dist for details)
input_dists[["SNPs"]] = weighted_hamming_fraction(binary_mat = input[["SNPs"]],
                                                  features = colnames(input[["SNPs"]]),
                                                  drivers = COSMIC_BC_drivers,
                                                  driver_weight = 0.8,
                                                  nondriver_weight = 0.2)
gc()

# Hyperparameter tuning ###
sigma_step = 0.1
# iter_step = 10
neighbor_step = 5
num_neighbors_range = seq(10, 50, neighbor_step) # number of neighbors, usually (10~30)
sigma_range = seq(0.3, 0.8, sigma_step) 	# hyperparameter, usually (0.3~0.8) -REFFERED as \mu in report text
# iterations = seq(10, 100, iter_step)      # Number of Iterations, usually (10~20)
n_iterations = 50     # Number of Iterations, usually (10~20)

# Affinity matrices ###
affinity_object = list()
for (nn in num_neighbors_range) {
  nn_list <- list()
  for (sigma in sigma_range) {
    sigma_list <- list()
    for (modality in modalities) {
      aff_mat <- abSNF::affinityMatrix(as.matrix(input_dists[[modality]]), K = nn, sigma = sigma)
      sigma_list[[modality]] <- list(
        num_neighbors = nn,
        regularization = sigma,
        affinity_matrix = aff_mat
      )
      rm(aff_mat)
    }
    nn_list[[paste0("sigma = ", sigma)]] <- sigma_list
    rm(sigma_list)
  }
  affinity_object[[paste0("NN = ", nn)]] <- nn_list
  rm(nn_list)
}
rm(nn, sigma)

# Fusions ###
# Try parallel
library(parallel)
library(foreach)
library(doParallel)

# Use 5 cores
cl = makeCluster(6) # or 3 (depending on resources)
registerDoParallel(cl)

Fusions = list()

# Use foreach to parallelize the computation (< 30 min)
Fusions <- foreach(nn = num_neighbors_range, .combine = 'c', .packages = 'abSNF') %:%
  foreach(sigma = sigma_range, .combine = 'c') %dopar% {
    sublist <- affinity_object[[paste0("NN = ", nn)]][[paste0("sigma = ", sigma)]]
    iter_matrices <- lapply(sublist, `[[`, "affinity_matrix")
    fusion_result <- abSNF::SNF(iter_matrices, K = nn, t = n_iterations)
    list(fusion_result)
  }

# Stop the cluster
stopCluster(cl)
gc()

# Restructure the Fusions list to match the desired output format
names(Fusions) <- unlist(lapply(num_neighbors_range, function(nn) {
  lapply(sigma_range, function(sigma) {
    paste0("NN = ", nn, ", sigma = ", sigma)
  })
}))

# Give appropriate colnames and rownames
column_names = colnames(affinity_object[["NN = 10"]][["sigma = 0.3"]][["RNAseq"]][["affinity_matrix"]])
row_names = rownames(affinity_object[["NN = 10"]][["sigma = 0.3"]][["RNAseq"]][["affinity_matrix"]])

for (i in 1:length(Fusions)) {
  colnames(Fusions[[i]]) = column_names
  rownames(Fusions[[i]]) = row_names
}

save.image(paste0(home, "/Results/single_algorithm/", 
                  algorithm, "/", algorithm, "_", data_source, "_",
                  data_types, "_eval_on_", evaluation_source,
                  "_env.RData"))

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
  optSigma = median(sigma_range) # 0.55
}

# We will therefore pick the median sigma, rounded down to 0.5, a value we actually used in affinities.
optSigma = 0.5

# The choice of number of neighbors affects the final matrix more than sigma
# according to both parametric and non-parametric tests

# Pearson and Frobenius histograms for different nn AND sigma = 0.5
Fusions_filt = Fusions[which(grepl("sigma = 0.5", names(Fusions)))]
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

# Based on Pearson similarities and Frobenius norms, the matrices for varying NN
# and sigma = 0.5, are very similar, so the choice of NN does not really matter in this case

# Spectral clustering for the estimated optimal k by SNFtool for different values of NN
RNGversion("4.2.2")
set.seed(123)
numc_df = data.frame(matrix(nrow = 0, ncol = 9))
colnames = c("Fusion", "K1", "K12", "K2", "K22", "eigengap_K1_score", "eigengap_K2_score",
             "rotation_K2_score", "rotation_K22_score")
for (i in 1:length(Fusions_filt)) {
  fusion = Fusions_filt[[i]]
  optks = as.data.frame(list(Fusion = names(Fusions_filt)[i], 
                             estimateNumberOfClustersGivenGraph_mod(fusion, NUMC=2:10))) # similar to the ab-SNF paper (but they start from >=6 clusters)
  numc_df = rbind(numc_df, optks)
}

# Inspection of numc_df shows that 2 is unequivocally the optimal number of clusters
optk = 2

### NOT RUN ###
# # We therefore proceed with picking the fused matrix with the highest separation
# # as per the heuristic below:
# 
# # We then choose the nn value for which the
# # fused similarity matrix has the "best" bimodal distribution of low and high values.
# # WE USE THIS APPROACH ONLY BECAUSE THE NUMBER OF CLUSTERS WE SEEK IS 2!
# 
# # We use combine two methodologies to do it:
# 
# # 1. Get the sum of variance and IQR for every matrix
# # 2. Get the sum of absolute skewness and kurtosis
# # 3. Find the nn matrix for which the sum of 1 and 2 is maximum
# 
# # Variance and IQR
# choose_matrix_contrasts <- function(similarity_matrices) {
#   contrast_values <- list()
#   
#   for (name in names(similarity_matrices)) {
#     similarity_values <- as.vector(similarity_matrices[[name]])
#     similarity_values <- similarity_values[similarity_values != 1] # remove self-similarities
#     
#     # Calculate measures of contrast
#     variance <- var(similarity_values)
#     iqr <- IQR(similarity_values)
#     contrast_metric <- variance + iqr
#     contrast_values[[name]] <- contrast_metric
#   }
#   
#   return(contrast_values)
# }
# 
# # Skewness and kurtosis
# library(e1071)
# 
# choose_matrix_skewness_kurtosis <- function(similarity_matrices) {
#   skewness_kurtosis_values <- list()
#   
#   for (name in names(similarity_matrices)) {
#     similarity_values <- as.vector(similarity_matrices[[name]])
#     similarity_values <- similarity_values[similarity_values != 1] # remove self-similarities
#     
#     # Calculate skewness and kurtosis
#     skewness_value <- skewness(similarity_values)
#     kurtosis_value <- kurtosis(similarity_values)
#     
#     # Calculate a combined metric: |skewness| + kurtosis
#     combined_metric <- abs(skewness_value) + kurtosis_value
#     skewness_kurtosis_values[[name]] <- combined_metric
#   }
#   
#   return(skewness_kurtosis_values)
# }
# 
# contrast_list <- choose_matrix_contrasts(Fusions_filt)
# skewness_kurtosis_list <- choose_matrix_skewness_kurtosis(Fusions_filt)
# 
# # Min-max normalization to [0,1]
# minmax_normalize_values <- function(values) {
#   min_value <- min(values)
#   max_value <- max(values)
#   
#   normalized_values <- (values - min_value) / (max_value - min_value)
#   return(normalized_values)
# }
# 
# contrast_values <- unlist(contrast_list)
# skewness_kurtosis_values <- unlist(skewness_kurtosis_list)
# 
# # Normalize the contrast values and skewness-kurtosis values
# normalized_contrast <- minmax_normalize_values(contrast_values)
# normalized_skewness_kurtosis <- minmax_normalize_values(skewness_kurtosis_values)
# 
# # Get the sum
# combined_scores <- normalized_contrast + normalized_skewness_kurtosis
# combined_scores_list <- setNames(as.list(combined_scores), names(contrast_list))
# print(combined_scores_list)
# 
# best_combined_matrix <- names(combined_scores_list)[which.max(combined_scores)]
# cat("Best similarity matrix based on combined normalized scores:", best_combined_matrix, "\n")

# We now run spectral clustering for the different values of NN and pick the one with
# the highest average silhouette index ###
RNGversion("4.2.2")
set.seed(123)

# Adapted source code to return eigenvectors as well
spectralClustering_eig <- function (affinity, K, type = 3) 
{
  library(cluster)
  d = rowSums(affinity)
  d[d == 0] = .Machine$double.eps
  D = diag(d)
  L = D - affinity
  if (type == 1) {
    NL = L
  }
  else if (type == 2) {
    Di = diag(1/d)
    NL = Di %*% L
  }
  else if (type == 3) {
    Di = diag(1/sqrt(d))
    NL = Di %*% L %*% Di
  }
  eig = eigen(NL)
  res = sort(abs(eig$values), index.return = TRUE)
  U = eig$vectors[, res$ix[1:K]]
  normalize <- function(x) x/sqrt(sum(x^2))
  if (type == 3) {
    U = t(apply(U, 1, normalize))
  }
  eigDiscrete = SNFtool:::.discretisation(U)
  eigDiscrete = eigDiscrete$discrete
  labels = apply(eigDiscrete, 1, which.max)
  U = as.data.frame(cbind(U, labels))
  U$Sample.ID = colnames(affinity)
  colnames(U)[(ncol(U)-1):ncol(U)] = c("Cluster", "Sample.ID")
  return(U)
}

NN_runs = vector("list", length(Fusions_filt))
names(NN_runs) = names(Fusions_filt)
for (i in 1:length(Fusions_filt)) {
  res = spectralClustering_eig(Fusions_filt[[i]], optk)
  col_index = ncol(res) - 2
  sil = silhouette(as.integer(res$Cluster),
                   dist = Rfast::Dist(res[, 1:col_index], method = "euclidean"))
  avg_width = summary(sil)$avg.width
  NN_runs[[names(Fusions_filt)[i]]][["Results"]] = res
  NN_runs[[names(Fusions_filt)[i]]][["Silhouette"]] = sil
  NN_runs[[names(Fusions_filt)[i]]][["Avg. sil. width"]] = avg_width
}
rm(res, col_index, sil, avg_width, i); gc()

# Optimal combination is for maximum avg. silhouette width
# NN = 10, sigma = 0.5
opt_index = which.max(lapply(NN_runs, function(x) x[["Avg. sil. width"]]))
optimal_abSNF = NN_runs[[opt_index]]
optN = as.numeric(substr(names(Fusions_filt)[opt_index], 6, 7))
final_affinity_matrix = Fusions[[paste0("NN = ", optN, ", sigma = ", optSigma)]]

# Concordance between final matrix and individual modality affinity matrices
conc_NMI = concordanceNetworkNMI(list(final_affinity_matrix, 
                                      affinity_object[[paste0("NN = ", optN)]][[paste0("sigma = ", optSigma)]]$RNAseq$affinity_matrix,
                                      affinity_object[[paste0("NN = ", optN)]][[paste0("sigma = ", optSigma)]]$CNV$affinity_matrix,
                                      affinity_object[[paste0("NN = ", optN)]][[paste0("sigma = ", optSigma)]]$Methylation$affinity_matrix,
                                      affinity_object[[paste0("NN = ", optN)]][[paste0("sigma = ", optSigma)]]$miRNA$affinity_matrix,
                                      affinity_object[[paste0("NN = ", optN)]][[paste0("sigma = ", optSigma)]]$SNPs$affinity_matrix),
                                 C = optk)
dimnames(conc_NMI) = list(c("Fusion", "RNAseq", "CNV", "Methylation", "miRNA", "SNPs"),
                          c("Fusion", "RNAseq", "CNV", "Methylation", "miRNA", "SNPs"))

abSNF_clusters = as.data.frame(list(Sample.ID = optimal_abSNF$Results$Sample.ID,
                                    Cluster = optimal_abSNF$Results$Cluster))
abSNF_clusters$Sample.ID = gsub("\\.", "-", abSNF_clusters$Sample.ID)
rownames(abSNF_clusters) = abSNF_clusters$Sample.ID


# Feature ranking
abSNF_feature_ranks = list()
binary_flags = c(TRUE, FALSE, FALSE, FALSE, FALSE)
for (i in 1:length(input)) {
  abSNF_feature_ranks[[i]] = rankFeaturesByNMI_parallely(data = list(input[[i]]), 
                                                         W = final_affinity_matrix,
                                                         ncores = 8,
                                                         binary = binary_flags[i])
  cat("Done with", names(input)[i], "\n")
}
names(abSNF_feature_ranks) = names(input)

rank_sum_nas = sapply(abSNF_feature_ranks, function(x) sum(is.na(x$NMI_ranks)))
names(rank_sum_nas) = names(abSNF_feature_ranks)

feature_ranks_text1 = paste0("Feature ranks could not be calculated for some features in the ",
                             paste(names(rank_sum_nas)[which(rank_sum_nas > 0)], collapse = ", "), 
                             ifelse(length(which(rank_sum_nas > 0)) > 1, " modalities (", " modality ("),
                             paste(rank_sum_nas[which(rank_sum_nas > 0)], collapse = ", "),
                             ").")
rm(rank_sum_nas); gc()

# Combine all modalities into a single data frame
for (i in 1:length(abSNF_feature_ranks)) {
  names(abSNF_feature_ranks[[i]]$NMI_scores) = colnames(input[[i]])
}

# Prepare data for plotting
# Combine all modalities into a single data frame
feature_data <- do.call(rbind, lapply(names(abSNF_feature_ranks), function(modality) {
  data.frame(
    Feature = paste(modality, seq_along(abSNF_feature_ranks[[modality]]$NMI_scores), sep = "_"),
    Modality = modality,
    NMI_Scores = abSNF_feature_ranks[[modality]]$NMI_scores
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
                                    cluster_df2 = abSNF_clusters,
                                    sample_col = "Sample.ID",
                                    clust_col = "Cluster",
                                    suffixes = c("_MOVICS_CS",
                                                 paste0("_", algorithm)))

NMI_to_MOVICS = calculate_nmi_index(cluster_df1 = ground_truth_labels,
                                    cluster_df2 = abSNF_clusters,
                                    sample_col = "Sample.ID",
                                    clust_col = "Cluster",
                                    suffixes = c("_MOVICS_CS",
                                                 paste0("_", algorithm)))

# Very low statistics when compared to the MOVICS. Results differ

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
  inner_join(abSNF_clusters, by = "Sample.ID") %>%
  tibble::column_to_rownames(var = "Sample.ID") %>%
  mutate(abSNF = paste0(algorithm, Cluster)) %>%
  dplyr::select(abSNF, everything()) %>%
  dplyr::select(-Cluster)
rm(scheme); gc()

# Silhouette
getSilhouette_ggplot(sil      = optimal_abSNF$Silhouette,
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

plot_object = list(clust.res = abSNF_clusters %>%
                     dplyr::rename(samID = Sample.ID, clust = Cluster))

# Export consensus clustering object
clust = as.data.frame(plot_object$clust.res)
colnames(clust) = c("Sample.ID", "Cluster")
clust$Cluster = paste0(algorithm, clust$Cluster)
openxlsx::write.xlsx(clust, paste0(home, "/Results/single_algorithm/", algorithm, "/", 
                                   algorithm, "_", data_source, "_",
                                   data_types, "_eval_on_", evaluation_source,
                                   "_clusterings.xlsx"))