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

# Load NEMO
# devtools::install_github('Shamir-Lab/NEMO/NEMO')
library(NEMO)
library(SNFtool) # required by NEMO

# Preamble
home = getwd()
algorithm = "NEMO"
alg_feature_pref = "rows" # Where does the algorithm expect the features to be
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

# However NEMO prefers features in columns so we transpose the matrices.

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

# Setup ###

# NEMO is using the SNFtool approach to construct an affinity matrix for each omic
# separately. However, it treats them differently downstream

# We can therefore modify the function nemo.affinity.graph() to allow for tuning
# the number of neighbors and sigma as we did in SNF. Additionally, we want to
# incorporate binary distance calculation instead of Euclidean for the SNPs kernels

# The individual omics affinity matrices should be identical to the SNF cases, but
# the final affinity matrix (Average Relative Similarity Matrix) - on which the 
# final clustering will be based - will be different to the Fusion matrix from SNF

# Hyperparameter tuning
sigma_step = 0.1
# iter_step = 10
neighbor_step = 5
num_neighbors_range = seq(10, 50, neighbor_step) # number of neighbors, usually (10~30)
sigma_range = seq(0.3, 0.8, sigma_step) 	# hyperparameter, usually (0.3~0.8)

# NEMO pipeline ###
# Parallel approach
library(parallel)
library(foreach)
library(doParallel)

# Set up a cluster using 6 cores
cl <- makeCluster(6)
registerDoParallel(cl)

# Define the list to store results
ARS_object <- list()

timestamp()
# Parallelized computation using foreach for both loops
ARS_object <- foreach(nn = num_neighbors_range, .combine = 'c', 
                      .packages = c("NEMO", "SNFtool")) %:%
  foreach(sigma = sigma_range, .combine = 'c') %dopar% {
    
    # Compute the affinity matrix using nemo.affinity.graph_mod function
    aff_mat <- nemo.affinity.graph_mod(input, k = nn, sigma = sigma, 
                                       binary_flags = c("Yes", "No", "No", "No", "No"),
                                       binary_distance = "binary")
    
    # Create the sublist to store in ARS_object
    result_list <- list(
      num_neighbors = nn,
      regularization = sigma,
      affinity_matrix = aff_mat
    )
    
    # Name the result using the key and return the result list
    named_result <- list(result_list)
    names(named_result) <- paste0("NN = ", nn, ", sigma = ", sigma)
    return(named_result)
  }

# Stop the cluster once processing is complete
stopCluster(cl)

timestamp() #~27min

# Compare similarity matrices similarly to what we did to SNF
# Check similarities for a given nn
nn_similarities = list()
for (nn in num_neighbors_range) {
  indices = grepl(paste0("NN = ", nn), names(ARS_object))
  matrices = ARS_object[indices]
  matrices = lapply(matrices, function(x) x[["affinity_matrix"]])
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
  indices = grepl(paste0("sigma = ", sigma), names(ARS_object))
  matrices = ARS_object[indices]
  matrices = lapply(matrices, function(x) x[["affinity_matrix"]])
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
all_similarities = compute_matrix_similarity(lapply(ARS_object, function(x) x[["affinity_matrix"]]))

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

# The choice of number of neighbors affects the final matrix more than sigma, which
# also plays a role. These results are verified by both parametric and non-parametric tests

# Pearson and Frobenius histograms for all afinity matrices
Pearson_hist_matrix = all_similarities$Pearson
Pearson_values <- Pearson_hist_matrix[lower.tri(Pearson_hist_matrix, diag = FALSE)]
mean_Pearson_value <- mean(Pearson_values)
median_Pearson_value <- median(Pearson_values)
sd_Pearson_value <- sd(Pearson_values)

Frobenius_hist_matrix = all_similarities$Frobenius
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
                                     sigma, " values"))), 
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
                     "/Results/single_algorithm/NEMO/Supplement"), 
       width = 2880, height = 1820, device = 'pdf', units = "px",
       dpi = 700)
dev.off()

# Plot histogram of Frobenius values
ggplot(data = data.frame(Frobenius_values), aes(x = Frobenius_values)) +
  geom_histogram(breaks = seq(0, 47, length.out = length(Frobenius_values)),
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
                                     sigma, " values"))), 
       x = "Affinity Matrix Frobenius Values", y = "Frequency") +
  scale_x_continuous(name = "Affinity Matrix Frobenius Values", limits = c(0, 47),
                     breaks = seq(0, 47, 5), expand = c(0, 0)) +
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
                     "/Results/single_algorithm/NEMO/Supplement"), 
       width = 2880, height = 1820, device = 'pdf', units = "px",
       dpi = 700)
dev.off()

# We then choose the nn and sigma values for which the
# fused similarity matrix has the "best" bimodal distribution of low and high values.

# We use combine two methodologies to do it:

# 1. Get the sum of variance and IQR for every matrix
# 2. Get the sum of absolute skewness and kurtosis
# 3. Find the nn matrix for which the sum of 1 and 2 is maximum

library(e1071)
contrast_list <- choose_matrix_contrasts(lapply(ARS_object, function(x) x[["affinity_matrix"]]))
skewness_kurtosis_list <- choose_matrix_skewness_kurtosis(lapply(ARS_object, function(x) x[["affinity_matrix"]]))
contrast_values <- unlist(contrast_list)
skewness_kurtosis_values <- unlist(skewness_kurtosis_list)

# Normalize the contrast values and skewness-kurtosis values
normalized_contrast <- minmax_normalize_values(contrast_values)
normalized_skewness_kurtosis <- minmax_normalize_values(skewness_kurtosis_values)

# Get the sum
combined_scores <- normalized_contrast + normalized_skewness_kurtosis
combined_scores_list <- setNames(as.list(combined_scores), names(contrast_list))
print(combined_scores_list)

best_combined_matrix <- names(combined_scores_list)[which.max(combined_scores)]
cat("Best similarity matrix based on combined normalized scores:", best_combined_matrix, "\n")

# We choose nn = 10, sigma = 0.3
optN = as.numeric(substr(best_combined_matrix, 6, 7))
optSigma = as.numeric(substr(best_combined_matrix, 18, 20))

conclusion2 = paste0("Best similarity matrix based on combined normalized scores is for $nn' = ", 
                     optN, "$ and $\\sigma = ", optSigma, "$. We therefore proceed with $nn' = ",
                     optN, "$ and $\\sigma = ", optSigma, "$.")
conclusion = paste(conclusion, conclusion2)

# Spectral clustering for k = ground_truth_k from MOVICS
RNGversion("4.2.2")
set.seed(123)

final_affinity_matrix = ARS_object[[paste0("NN = ", optN, ", sigma = ", optSigma)]]$affinity_matrix

# NEMO estimation of number of clusters:
num.clusters = nemo.num.clusters(final_affinity_matrix, NUMC = 2:10)

# 9 is the optimal number of clusters based on the NEMO process
group = spectralClustering(final_affinity_matrix, num.clusters)
names(group) = colnames(final_affinity_matrix)

NEMO_clusters = as.data.frame(list(Sample.ID = names(group),
                                  Cluster = group))
rownames(NEMO_clusters) = NEMO_clusters$Sample.ID

# Import resources
scheme = readRDS("Resources/scheme.rds")
annCol = scheme$annCol
annColors = scheme$annColors
cluster_colors = scheme$clust.colors
col.list = scheme$col.list
var2comp = scheme$var2comp %>%
  dplyr::select(-`Consensus Subtype`) %>%
  mutate(Sample.ID = rownames(.)) %>%
  inner_join(NEMO_clusters, by = "Sample.ID") %>%
  tibble::column_to_rownames(var = "Sample.ID") %>%
  mutate(NEMO = paste0(algorithm, Cluster)) %>%
  dplyr::select(NEMO, everything()) %>%
  dplyr::select(-Cluster)
rm(scheme); gc()
cluster_colors = c("#2EC4B6", "#E71D36", "#FF9F1C", "#BDD5EA", 
                   "#FFA5AB", "#011627", "#023E8A", "#9D4EDD", "#f09c6c")

# Main results ###
# Examine cluster similarity to MOVICS by measuring NMI and ARI indices #####
# (Jaccard may be misleading)

# Calculate ARI and NMI
library(mclust)
library(clue)

ARI_to_MOVICS = calculate_ari_index(cluster_df1 = ground_truth_labels,
                                    cluster_df2 = NEMO_clusters,
                                    sample_col = "Sample.ID",
                                    clust_col = "Cluster",
                                    suffixes = c("_MOVICS", "_NEMO"))

NMI_to_MOVICS = calculate_nmi_index(cluster_df1 = ground_truth_labels,
                                    cluster_df2 = NEMO_clusters,
                                    sample_col = "Sample.ID",
                                    clust_col = "Cluster",
                                    suffixes = c("_MOVICS", "_NEMO"))

# Very low statistics when compared to the MOVICS. Results differ

# MOVICS-like analysis #####
library(MOVICS)
library(ComplexHeatmap)

plotdata <- lapply(lapply(input, as.matrix), 
                   function(mat) mat[rowSums(mat != 0) > 0, ])
heatmap_plotdata = getStdiz(
  data = plotdata,
  halfwidth = c(NA, 3, 3, 3, 3), # No halfwidth for SNPs
  centerFlag = c(F, F, F, F, F),
  scaleFlag = c(F, F, F, F, F)
)

plot_object = list(clust.res = NEMO_clusters %>%
                     dplyr::rename(samID = Sample.ID, clust = Cluster))

# comprehensive heatmap (may take a while)
getMoHeatmap_single_algorithm(algorithm_name = algorithm,
                              data          = heatmap_plotdata,
                              row.title     = names(heatmap_plotdata),
                              is.binary     = c(T,F,F,F,F), 
                              legend.name   = c("SNPs",
                                                "Standardized RNAseq norm. counts",
                                                "Standardized CNV",
                                                "Standardized miRNA norm. counts",
                                                "Standardized Methylation M-values"
                              ),
                              clust.res     = plot_object$clust.res, # consensusMOIC-like results
                              clust.dend    = NULL, # show no dendrogram for samples
                              show.rownames = c(F,F,F,F,F), # specify for each omics data
                              show.colnames = FALSE, # show no sample names
                              show.row.dend = c(F,F,F,F,F), # show no dendrogram for features
                              annRow        = NULL, # no selected features
                              color         = col.list,
                              annCol        = annCol, # annotation for samples
                              annColors     = annColors, # annotation color
                              width         = 20, # width of each subheatmap
                              height        = 10, # height of each subheatmap
                              fig.path      = paste0(home, "/Results/single_algorithm/NEMO"),
                              fig.name      = paste0("default_", algorithm, "_Comprehensive_heatmap"))
dev.off()
gc()

# Clinical variables ###
# Remove unknown levels for statistical tests
var2comp_nonas = var2comp
for (i in 1:ncol(var2comp)) {
  nas = which(var2comp[, i] == "Unknown")
  var2comp_nonas[nas, i] = NA
}
rm(nas); gc()

# Statistical comparisons
clin_comp = compClinvar_single_algorithm(algorithm_name = algorithm,
                                         moic.res = plot_object,
                                         var2comp = var2comp_nonas,
                                         strata = algorithm,
                                         factorVars = c("vital_status", "race_list", "ethnicity",
                                                        "history_of_neoadjuvant_treatment",
                                                        "primary_lymph_node_presentation_assessment",
                                                        "histological_type", "menopause_status",
                                                        "breast_carcinoma_progesterone_receptor_status",
                                                        "breast_carcinoma_estrogen_receptor_status",
                                                        "lab_proc_her2_neu_immunohistochemistry_receptor_status",
                                                        "distant_metastasis_present_ind2",
                                                        "stage_event_pathologic_stage"),
                                         includeNA = FALSE,
                                         doWord = TRUE,
                                         tab.name = "Summary_of_clinical_variables",
                                         res.path = paste0(home, "/Results/single_algorithm/NEMO/"),
                                         output_pdf = TRUE,
                                         pdf_level_col_width = c("5em", "5em"),
                                         pdf_count_col_width = "5em",
                                         pdf_tab_font_size = 5)

# race_list, ER status, PR status, metastasis are sig

# Oncoprint ###
# Internally set the workspace argument in Fisher's exact test to 5e+09
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
                                      fig.path     = paste0(home, "/Results/single_algorithm/NEMO"),
                                      res.path     = paste0(home, "/Results/single_algorithm/NEMO"),
                                      simulate.p.value = TRUE)

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
                                                 fig.path = paste0(home, "/Results/single_algorithm/NEMO"),
                                                 width = 10,
                                                 notch     = TRUE)

# Agreement with other subtypes ###
subtype_agreement <- compAgree_single_algorithm(algorithm_name = algorithm,
                                                moic.res  = plot_object,
                                                subt2comp = annCol[, c("ER status", "PR status",
                                                                       "HER2 status", "Metastasis", "Stage")],
                                                doPlot    = TRUE,
                                                box.width = 0.2,
                                                fig.name  = "Classification_agreement",
                                                fig.path  = paste0(home, "/Results/single_algorithm/NEMO"),
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
              res.path = paste0(home, "/Results/single_algorithm/NEMO"),
              algorithm = algorithm)

# # Identify unique subtype biomarkers
# # 1. Up-regulated markers
dgea.marker.up <- runMarker_single_algorithm(algorithm_name = algorithm,
                                             moic.res = plot_object,
                                             dea.method    = "limma", # name of DEA method
                                             prefix        = "dgea_", # MUST be the same of argument in runDEA()
                                             dat.path      = paste0(home, "/Results/single_algorithm/NEMO"), # path of DEA files
                                             res.path      = paste0(home, "/Results/single_algorithm/NEMO"), # path to save marker files
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
                                             fig.path = paste0(home, "/Results/single_algorithm/NEMO"),
                                             width = 14,
                                             height = 12,
                                             fontsize_row = 0, # 3 default
                                             name = "normalized RNA-seq")
dev.off()

# # 2. Down-regulated markers
dgea.marker.down <- runMarker_single_algorithm(algorithm_name = algorithm,
                                               moic.res = plot_object,
                                               dea.method    = "limma", # name of DEA method
                                               prefix        = "dgea_", # MUST be the same of argument in runDEA()
                                               dat.path      = paste0(home, "/Results/single_algorithm/NEMO"), # path of DEA files
                                               res.path      = paste0(home, "/Results/single_algorithm/NEMO"), # path to save marker files
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
                                               fig.path = paste0(home, "/Results/single_algorithm/NEMO"),
                                               width = 14,
                                               height = 12,
                                               fontsize_row = 0, # 3 default
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
                  res.path = paste0(home, "/Results/single_algorithm/NEMO"),
                  algorithm = algorithm)

# # Identify unique subtype biomarkers
# # 1. Up-regulated markers
methyl.marker.up <- runMarker_single_algorithm(algorithm_name = algorithm,
                                               moic.res = plot_object,
                                               dea.method    = "limma", # name of DEA method
                                               prefix        = "dmea_", # MUST be the same of argument in runDEA()
                                               dat.path      = paste0(home, "/Results/single_algorithm/NEMO"), # path of DEA files
                                               res.path      = paste0(home, "/Results/single_algorithm/NEMO"), # path to save marker files
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
                                               fig.path = paste0(home, "/Results/single_algorithm/NEMO"),
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
                                                 dat.path      = paste0(home, "/Results/single_algorithm/NEMO"), # path of DEA files
                                                 res.path      = paste0(home, "/Results/single_algorithm/NEMO"), # path to save marker files
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
                                                 fig.path = paste0(home, "/Results/single_algorithm/NEMO"),
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
                    res.path = paste0(home, "/Results/single_algorithm/NEMO"),
                    algorithm = algorithm)

# # Identify unique subtype biomarkers
# # 1. Up-regulated markers
miRNA.marker.up <- runMarker_single_algorithm(algorithm_name = algorithm,
                                              moic.res = plot_object,
                                              dea.method    = "limma", # name of DEA method
                                              prefix        = "dmiRea_", # MUST be the same of argument in runDEA()
                                              dat.path      = paste0(home, "/Results/single_algorithm/NEMO"), # path of DEA files
                                              res.path      = paste0(home, "/Results/single_algorithm/NEMO"), # path to save marker files
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
                                              fig.path = paste0(home, "/Results/single_algorithm/NEMO"),
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
                                                dat.path      = paste0(home, "/Results/single_algorithm/NEMO"), # path of DEA files
                                                res.path      = paste0(home, "/Results/single_algorithm/NEMO"), # path to save marker files
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
                                                fig.path = paste0(home, "/Results/single_algorithm/NEMO"),
                                                width = 14,
                                                height = 12,
                                                fontsize_row = 0, # 3 default
                                                name = "normalized miRNA")
dev.off()

# DCNVA ###
dCNVea = runDEA_mod(dea.method = "limma", # we use normalized data as input
                    expr = plotdata$CNV,
                    moic.res = plot_object,
                    prefix = "dCNVea_",
                    sort.p = TRUE,
                    overwt = TRUE,
                    verbose = TRUE,
                    res.path = paste0(home, "/Results/single_algorithm/NEMO"),
                    algorithm = algorithm)

# # Identify unique subtype biomarkers
# # 1. Up-regulated markers
CNV.marker.up <- runMarker_single_algorithm(algorithm_name = algorithm,
                                            moic.res = plot_object,
                                            dea.method    = "limma", # name of DEA method
                                            prefix        = "dCNVea_", # MUST be the same of argument in runDEA()
                                            dat.path      = paste0(home, "/Results/single_algorithm/NEMO"), # path of DEA files
                                            res.path      = paste0(home, "/Results/single_algorithm/NEMO"), # path to save marker files
                                            p.cutoff      = 0.05, # p cutoff to identify significant DEGs
                                            p.adj.cutoff  = 0.05, # padj cutoff to identify significant DEGs
                                            dirct         = "up", # direction of dysregulation in expression
                                            n.marker      = 100, # number of biomarkers for each subtype
                                            doplot        = TRUE, # generate diagonal heatmap
                                            norm.expr     = plotdata$CNV, # use normalized expression as heatmap input
                                            annCol        = annCol, # sample annotation in heatmap
                                            annColors     = annColors, # colors for sample annotation
                                            show_rownames = TRUE, # show no rownames (biomarker name)
                                            centerFlag = F,
                                            scaleFlag = F,
                                            halfwidth = 3,
                                            fig.name      = "upregulated_CNV_biomarkers_heatmap",
                                            fig.path = paste0(home, "/Results/single_algorithm/NEMO"),
                                            width = 14,
                                            height = 12,
                                            fontsize_row = 0, # 3 default
                                            name = "normalized CNV")
dev.off()

# # 2. Down-regulated markers
CNV.marker.down <- runMarker_single_algorithm(algorithm_name = algorithm,
                                              moic.res = plot_object,
                                              dea.method    = "limma", # name of DEA method
                                              prefix        = "dCNVea_", # MUST be the same of argument in runDEA()
                                              dat.path      = paste0(home, "/Results/single_algorithm/NEMO"), # path of DEA files
                                              res.path      = paste0(home, "/Results/single_algorithm/NEMO"), # path to save marker files
                                              p.cutoff      = 0.05, # p cutoff to identify significant DEGs
                                              p.adj.cutoff  = 0.05, # padj cutoff to identify significant DEGs
                                              dirct         = "down", # direction of dysregulation in expression
                                              n.marker      = 100, # number of biomarkers for each subtype
                                              doplot        = TRUE, # generate diagonal heatmap
                                              norm.expr     = plotdata$CNV, # use normalized expression as heatmap input
                                              annCol        = annCol, # sample annotation in heatmap
                                              annColors     = annColors, # colors for sample annotation
                                              show_rownames = TRUE, # show no rownames (biomarker name)
                                              centerFlag = F,
                                              scaleFlag = F,
                                              halfwidth = 3,
                                              fig.name      = "downregulated_CNV_biomarkers_heatmap",
                                              fig.path = paste0(home, "/Results/single_algorithm/NEMO"),
                                              width = 14,
                                              height = 12,
                                              fontsize_row = 0, # 3 default
                                              name = "normalized CNV")
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
                                            dat.path      = paste0(home, "/Results/single_algorithm/NEMO"), # path of DEA files
                                            res.path      = paste0(home, "/Results/single_algorithm/NEMO"), # path to save marker files
                                            msigdb.path  = MSIGDB.FILE, # MUST be the ABSOLUTE path of msigdb file
                                            norm.expr    = plotdata$RNAseq, # use normalized expression to calculate enrichment score
                                            dirct        = "up", # direction of dysregulation in pathway
                                            n.path       = 10,
                                            p.cutoff     = 0.05, # p cutoff to identify significant pathways
                                            p.adj.cutoff = 0.1, # padj cutoff to identify significant pathways
                                            gsva.method  = "gsva", # method to calculate single sample enrichment score
                                            name         = "GSVA scores", # name for colorbar
                                            norm.method  = "mean", # normalization method to calculate subtype-specific enrichment score
                                            fig.name     = "upregulated_pathway_heatmap",
                                            nPerm = 10000,
                                            minGSSize = 5, # default: 10
                                            maxGSSize = 500,
                                            fig.path = paste0(home, "/Results/single_algorithm/NEMO"),
                                            width = 14, height = 18) # default 12

# GSEA down-regulated
RNGversion("4.2.2")
set.seed(123)
gsea.down <- runGSEA_mod_4.4_single_algorithm(algorithm_name = algorithm,
                                              moic.res     = plot_object,
                                              dea.method   = "limma", # name of DEA method
                                              prefix       = "dgea_", # MUST be the same of argument in runDEA()
                                              dat.path      = paste0(home, "/Results/single_algorithm/NEMO"), # path of DEA files
                                              res.path      = paste0(home, "/Results/single_algorithm/NEMO"), # path to save marker files
                                              msigdb.path  = MSIGDB.FILE, # MUST be the ABSOLUTE path of msigdb file
                                              norm.expr    = plotdata$RNAseq, # use normalized expression to calculate enrichment score
                                              dirct        = "down", # direction of dysregulation in pathway
                                              n.path       = 10,
                                              p.cutoff     = 0.05, # p cutoff to identify significant pathways
                                              p.adj.cutoff = 0.1, # padj cutoff to identify significant pathways
                                              gsva.method  = "gsva", # method to calculate single sample enrichment score
                                              name         = "GSVA scores", # name for colorbar
                                              norm.method  = "mean", # normalization method to calculate subtype-specific enrichment score
                                              fig.name     = "downregulated_pathway_heatmap",
                                              nPerm = 10000,
                                              minGSSize = 5, # default: 10
                                              maxGSSize = 500,
                                              fig.path = paste0(home, "/Results/single_algorithm/NEMO"),
                                              width = 14, height = 18) # default 12

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
                                            fig.path      = paste0(home, "/Results/single_algorithm/NEMO"),
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

# Fraction Genome Altered ###
readRDS("Resources/TCGA/fga_df.rds")

fga.NEMO <- compFGA_mod(moic.res     = plot_object,
                           segment      = fga_df,
                           iscopynumber = TRUE, 
                           test.method  = "nonparametric", # statistical testing method (Wilcoxon with asymptotic approximation. Consider Kruskall Wallis?)
                           fig.path     = paste0(home, "/Results/single_algorithm/NEMO"),
                           fig.name     = paste0("FGA_barplot_", algorithm),
                           prefix = algorithm,
                           width = 16,
                           ga_column = "ga", # genome altered column
                           clust.col = cluster_colors)

fga.NEMO.COSMIC <- compFGA_mod(moic.res     = plot_object,
                        segment      = fga_df,
                        iscopynumber = TRUE, 
                        test.method  = "nonparametric", # statistical testing method (Wilcoxon with asymptotic approximation. Consider Kruskall Wallis?)
                        fig.path     = paste0(home, "/Results/single_algorithm/NEMO"),
                        fig.name     = paste0("COSMIC_criteria_FGA_barplot_", algorithm),
                        prefix = algorithm,
                        width = 16,
                        ga_column = "COSMIC_ga", # genome altered column
                        clust.col = cluster_colors)
# Evaluation #####
# Run Nearest Template Prediction in transNEO cohort ###
# Load transNEO data
transNEO_mm_inputs = readRDS("Resources/transNEO/transNEO_multimodal_inputs.rds")
transcr = readRDS("Resources/transNEO/log2TPMplus1_transNEO.rds")

# get as many templates as possible
dgea.marker.up_full <- runMarker_single_algorithm_no_export(algorithm_name = algorithm,
                                             moic.res = plot_object,
                                             dea.method    = "limma", # name of DEA method
                                             prefix        = "dgea_", # MUST be the same of argument in runDEA()
                                             dat.path      = paste0(home, "/Results/single_algorithm/NEMO"), # path of DEA files
                                             p.cutoff      = 0.05, # p cutoff to identify significant DEGs
                                             p.adj.cutoff  = 0.05, # padj cutoff to identify significant DEGs
                                             dirct         = "up", # direction of dysregulation in expression
                                             n.marker      = 60000)

# 2. Down-regulated markers
dgea.marker.down_full <- runMarker_single_algorithm_no_export(algorithm_name = algorithm,
                                                              moic.res = plot_object,
                                                              dea.method    = "limma", # name of DEA method
                                                              prefix        = "dgea_", # MUST be the same of argument in runDEA()
                                                              dat.path      = paste0(home, "/Results/single_algorithm/NEMO"), # path of DEA files
                                                              p.cutoff      = 0.05, # p cutoff to identify significant DEGs
                                                              p.adj.cutoff  = 0.05, # padj cutoff to identify significant DEGs
                                                              dirct         = "down", # direction of dysregulation in expression
                                                              n.marker      = 60000)

# Filter templates if too long
number_of_ntp_markers = 1000

ntp_topn_up <- dgea.marker.up_full %>%
  dplyr::filter(probe %in% rownames(transcr)) %>%
  group_by(class) %>% 
  slice_max(order_by = -padj, n = number_of_ntp_markers, with_ties = FALSE) %>%
  ungroup()

ntp_topn_down <- dgea.marker.down_full %>%
  group_by(class) %>% 
  slice_max(order_by = -padj, n = number_of_ntp_markers, with_ties = FALSE) %>%
  ungroup()

# Up-regulated expression features
RNGversion("4.2.2")
timestamp()
transNEO_ntp_expr_up = runNTP_mod(
  expr = transcr,
  templates = ntp_topn_up,
  scaleFlag = TRUE,
  centerFlag = TRUE,
  nPerm = 10000,
  seed = 123,
  distance = "cosine", # default
  doPlot = TRUE,
  height = 8,
  width = 12,
  fig.path = paste0(home, "/Results/single_algorithm/NEMO"),
  fig.name = "ntp_expr_up_heatmap_transNEO")
timestamp() # 16 min

# down-regulated
RNGversion("4.2.2")
transNEO_ntp_expr_down = runNTP_mod(
  expr = transcr,
  templates = ntp_topn_down,
  scaleFlag = TRUE,
  centerFlag = TRUE,
  nPerm = 10000,
  seed = 123,
  distance = "cosine", # default
  doPlot = TRUE,
  height = 8,
  width = 12,
  fig.path = paste0(home, "/Results/single_algorithm/NEMO"),
  fig.name = "ntp_expr_down_heatmap_transNEO")


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
  inner_join(expr_conc %>% dplyr::select(Donor.ID = samID, NEMO = clust_up),
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

eval_moic_res_clinvar = transNEO_ntp_expr_up
eval_moic_res_clinvar$clust.res$clust = factor(gsub(pattern = algorithm, 
                                             x = eval_moic_res_clinvar$clust.res$clust,
                                             replacement = ""))
  
transNEO_clincomp = compClinvar_single_algorithm(algorithm_name = algorithm,
                                                 moic.res = eval_moic_res_clinvar,
                                                 var2comp = transNEO_var2comp,
                                                 strata = algorithm,
                                                 factorVars = c("ER.status", "HER2.status", "Grade.pre.NAT",
                                                                "NAT.regimen", 
                                                                "pCR.RD", "LN.status.at.diagnosis"),
                                                 includeNA = FALSE,
                                                 doWord = TRUE,
                                                 tab.name = "transNEO_Summary_of_clinical_variables",
                                                 res.path = paste0(home, "/Results/single_algorithm/NEMO"),
                                                 output_pdf = TRUE,
                                                 pdf_level_col_width = c("5em", "5em"),
                                                 pdf_count_col_width = "5em",
                                                 pdf_tab_font_size = 5)

# Run PAM ###
RNGversion("4.2.2.")
set.seed(123)
transNEO_pam = runPAM_single_algorithm(algorithm_name = algorithm,
                                       train.expr = plotdata$RNAseq,
                                       moic.res   = plot_object,
                                       test.expr  = transcr)

# Check consistency across methods

# Get predictions for TCGA (discovery cohort)
RNGversion("4.2.2.")
set.seed(123)
TCGA.ntp.pred = runNTP_mod(expr = plotdata$RNAseq[, plot_object$clust.res$samID],
                       templates = ntp_topn_up, distance = "cosine",
                       doPlot = F, nPerm = 10000)

TCGA.pam.pred = runPAM_single_algorithm(algorithm_name = algorithm,
                                        train.expr = plotdata$RNAseq[, plot_object$clust.res$samID],
                                        moic.res = plot_object,
                                        test.expr = plotdata$RNAseq[, plot_object$clust.res$samID])

# consensus TCGA vs NTP TCGA # FAILS
runKappa_single_algorithm(algorithm_name = algorithm,
                          subt1 = plot_object$clust.res$clust,
                          subt2 = gsub(algorithm, "", TCGA.ntp.pred$clust.res$clust),
                          subt1.lab = "NEMO",
                          subt2.lab = "NTP TCGA",
                          height = 8,
                          width = 8,
                          fig.path = paste0(home, "/Results/single_algorithm/NEMO"),
                          fig.name = paste0("kappa_", algorithm, "_vs_NTP_TCGA"))

# consensus TCGA vs PAM TCGA
runKappa_single_algorithm(algorithm_name = algorithm,
                          subt1 = plot_object$clust.res$clust,
                          subt2 = gsub(algorithm, "", TCGA.pam.pred$clust.res$clust),
                          subt1.lab = "NEMO",
                          subt2.lab = "PAM TCGA",
                          height = 8,
                          width = 8,
                          fig.path = paste0(home, "/Results/single_algorithm/NEMO"),
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
                          fig.path = paste0(home, "/Results/single_algorithm/NEMO"),
                          fig.name = "kappa_NTP_vs_PAM_transNEO")

# Export consensus clustering object
clust = as.data.frame(plot_object$clust.res)
colnames(clust) = c("Sample.ID", "Cluster")
clust$Cluster = paste0("NEMO", clust$Cluster)
openxlsx::write.xlsx(clust, paste0(home, "/Results/single_algorithm/NEMO/", 
                                   algorithm, "_", data_source, "_",
                                   data_types, "_eval_on_", evaluation_source,
                                   "_clusterings.xlsx"))

# Supplementary results #####

# Create subdirectory for supplementary plots
if (!dir.exists(paste0(home, "/Results/single_algorithm/NEMO/Supplement"))) {
  dir.create(paste0(home, "/Results/single_algorithm/NEMO/Supplement"))
}

# Setup for heatmaps
colors_heatmap = rev(colorRampPalette(viridisLite::magma(10))(255))
cluster_colors_heatmap = c("#2EC4B6", "#E71D36", 
                           "#FF9F1C", "#BDD5EA", "#FFA5AB", "#011627", "#023E8A", 
                           "#9D4EDD", "#f09c6c")
clust_annot_pheno = annCol %>% mutate(Sample.ID = rownames(.)) %>%
  inner_join(clust, by = "Sample.ID") %>%
  dplyr::rename(NEMO = Cluster, samID = "Sample.ID")
rownames(clust_annot_pheno) = clust_annot_pheno$samID
afh_colnames = colnames(annCol)

# Prepare final affinity matrix
aff_final = final_affinity_matrix

# Final affinity matrix
create_MO_heatmap(matrix = final_affinity_matrix, algorithm = algorithm, 
                  need.diag.zero = TRUE, 
                  clust_annot_pheno = clust_annot_pheno,
                  afh_colnames = afh_colnames, 
                  colors = colors_heatmap,
                  annColors = annColors,
                  heatmap_title = "Average Relative Similarity matrix heatmap",
                  cluster_colors = cluster_colors_heatmap,
                  legend_title = "Normalized affinity",
                  cluster_cols_flag = FALSE,
                  cluster_rows_flag = FALSE,
                  splits_flag = TRUE,
                  output_file_name = paste0(home, "/Results/single_algorithm/NEMO/Supplement/ARS_heatmap.png"))

# Final affinity matrix with clustered rows and columns
create_MO_heatmap(matrix = final_affinity_matrix, algorithm = algorithm, 
                  need.diag.zero = TRUE, 
                  clust_annot_pheno = clust_annot_pheno,
                  afh_colnames = afh_colnames, 
                  colors = colors_heatmap,
                  annColors = annColors,
                  heatmap_title = "Average Relative Similarity matrix heatmap",
                  cluster_colors = cluster_colors_heatmap,
                  legend_title = "Normalized affinity",
                  cluster_cols_flag = TRUE,
                  cluster_rows_flag = TRUE,
                  splits_flag = FALSE,
                  output_file_name = paste0(home, "/Results/single_algorithm/NEMO/Supplement/hclust_ARS_heatmap.png"))

# PCA from original matrices (features as rows) ###
# Cluster labels must be numbers (not paste0(algorithm#))
NEMO_clust_res = NEMO_clusters %>% dplyr::rename(samID = Sample.ID, NEMO = Cluster)

# RNA
pca_from_original_matrix(mydata = input$RNAseq, 
                         algorithm = "NEMO", 
                         clust_res = NEMO_clust_res,
                         cluster_colors = cluster_colors, 
                         output_path = paste0(home, "/Results/single_algorithm/NEMO/Supplement"),
                         title_add = "RNAseq")

# miRNA
pca_from_original_matrix(mydata = input$miRNA, 
                         algorithm = "NEMO", 
                         clust_res = NEMO_clust_res,
                         cluster_colors = cluster_colors, 
                         output_path = paste0(home, "/Results/single_algorithm/NEMO/Supplement"),
                         title_add = "miRNA")

# CNV
pca_from_original_matrix(mydata = input$CNV, 
                         algorithm = "NEMO", 
                         clust_res = NEMO_clust_res,
                         cluster_colors = cluster_colors, 
                         output_path = paste0(home, "/Results/single_algorithm/NEMO/Supplement"),
                         title_add = "CNV")

# Use multidimensional scaling for SNPs
# Features must be in rows
mds_from_original_matrix(matrix = input$SNPs, dist_method = "binary",
                         algorithm = "NEMO", 
                         clust_res = NEMO_clust_res,
                         cluster_colors = cluster_colors, 
                         output_path = paste0(home, "/Results/single_algorithm/NEMO/Supplement"),
                         title_add = "SNPs")

# Methylation
pca_from_original_matrix(mydata = input$Methylation, 
                         algorithm = "NEMO", 
                         clust_res = NEMO_clust_res,
                         cluster_colors = cluster_colors, 
                         output_path = paste0(home, "/Results/single_algorithm/NEMO/Supplement"),
                         title_add = "Methylation")

# Final affinity
pca_from_sim_matrix(sim_matrix = aff_final, algorithm = algorithm, 
                    clust_res = NEMO_clust_res,
                    cluster_colors = cluster_colors,
                    output_path = paste0(home, "/Results/single_algorithm/NEMO/Supplement"), 
                    title_add = "ARS value")

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
                            "/Results/single_algorithm/NEMO/Supplement/Chisq_tests.xlsx"),
                     overwrite = TRUE)

# Bar chart generation
NEMO_barcharts = list()
plotdata_bar = clust_annot_pheno_nonas
plotdata_bar[[algorithm]] = factor(plotdata_bar[[algorithm]])
for (i in 1:length(voi)) {
  chifit = output
  loc = which(grepl(voi[i], chifit$Comparison))
  chifit = chifit[loc, ]
  NEMO_barcharts[[i]] = create_annot_barchart(plotdata = plotdata_bar, fill = voi[i],
                                              na.action = "na.omit",
                                              chifit = chifit,
                                              algorithm = algorithm,
                                              barchart_ylim = 500,
                                              text_y = 420, rect_ymin = 320,
                                              rect_ymax = 450, x_annot = 5.5,
                                              v_gap = 35, rect_xmin = 4.5,
                                              rect_xmax = 6.5, 
                                              annot_text_size = 2.25,
                                              legend.text.size = 5,
                                              x.axis.text.size = 5) +
    barchart_scales[[voi[i]]]
  print(NEMO_barcharts[[i]])
  ggsave(filename = paste0(algorithm, "_", voi[i], "_barchart.png"),
         path = paste0(home, 
                       "/Results/single_algorithm/NEMO/Supplement"), 
         width = 5320, height = 2320, device = 'png', units = "px",
         dpi = 700)
  dev.off()
}
names(NEMO_barcharts) = voi
rm(loc, chifit)

# Multiplot (PNG) - bar charts
library(ggpubr)
ggarrange(NEMO_barcharts[[1]], NEMO_barcharts[[2]], NEMO_barcharts[[3]],
          NEMO_barcharts[[4]], NEMO_barcharts[[5]], NEMO_barcharts[[6]],
          NEMO_barcharts[[7]], NEMO_barcharts[[8]], NEMO_barcharts[[9]],
          NEMO_barcharts[[10]], NEMO_barcharts[[11]],
          ncol = 3, nrow = 4, labels = c("A", "B", "C", "D", "E", "F", "G", "H",
                                         "I", "J", "K"),
          font.label = list(size = 8, face = "bold", color ="black"))
ggsave(filename = paste0("Multiplot_", algorithm, "_barcharts.png"),
       path = paste0(home, 
                     "/Results/single_algorithm/NEMO/Supplement"), 
       width = 12000, height = 8000, device = 'png', units = "px",
       dpi = 700)
dev.off()

# Just significant ones now
NEMO_barcharts_sig = list()
plotdata_bar_sig = clust_annot_pheno_nonas %>% dplyr::select(NEMO, Race, Histology, 
                                                             `ER status`, `PR status`,
                                                             `HER2 status`)
plotdata_bar_sig$NEMO = factor(plotdata_bar_sig$NEMO)
voi_sig = setdiff(colnames(plotdata_bar_sig), "NEMO")
for (i in 1:length(voi_sig)) {
  chifit = output
  loc = which(grepl(voi_sig[i], chifit$Comparison))
  chifit = chifit[loc, ]
  NEMO_barcharts_sig[[i]] = create_annot_barchart(plotdata = plotdata_bar_sig, fill = voi_sig[i],
                                                  na.action = "na.omit",
                                                  chifit = chifit,
                                                  algorithm = algorithm,
                                                  barchart_ylim = 500,
                                                  text_y = 420, rect_ymin = 320,
                                                  rect_ymax = 450, x_annot = 5.5,
                                                  v_gap = 35, rect_xmin = 4.5,
                                                  rect_xmax = 6.5, 
                                                  annot_text_size = 2.25,
                                                  legend.text.size = 5,
                                                  x.axis.text.size = 5) +
    barchart_scales[[voi_sig[i]]]
  print(NEMO_barcharts_sig[[i]])
  ggsave(filename = paste0("sig_NEMO_", voi_sig[i], "_barchart.png"),
         path = paste0(home, 
                       "/Results/single_algorithm/NEMO/Supplement"), 
         width = 5320, height = 2320, device = 'png', units = "px",
         dpi = 700)
  dev.off()
}
names(NEMO_barcharts_sig) = voi_sig
rm(loc, chifit)

# Multiplot (PNG) - bar charts
ggarrange(NEMO_barcharts_sig[[1]], NEMO_barcharts_sig[[2]], NEMO_barcharts_sig[[3]],
          NEMO_barcharts_sig[[4]], NEMO_barcharts_sig[[5]],
          ncol = 2, nrow = 3, labels = c("A", "B", "C", "D", "E"),
          font.label = list(size = 8, face = "bold", color ="black"))
ggsave(filename = "sig_Multiplot_NEMO_barcharts.png",
       path = paste0(home, 
                     "/Results/single_algorithm/NEMO/Supplement"), 
       width = 8500, height = 5500, device = 'png', units = "px",
       dpi = 700)
dev.off()

# Sunburst plot ###
library(plotly)
Pheno_sunburst_NEMO = clust_annot_pheno
Pheno_sunburst_NEMO$`ER status` = gsub("Unknown", "Unkn ER status", Pheno_sunburst_NEMO$`ER status`)
Pheno_sunburst_NEMO$`ER status` = gsub("Positive", "ER+", Pheno_sunburst_NEMO$`ER status`)
Pheno_sunburst_NEMO$`ER status` = gsub("Negative", "ER-", Pheno_sunburst_NEMO$`ER status`)
Pheno_sunburst_NEMO$`HER2 status` = gsub("Unknown", "Unkn HER2 status", 
                                          Pheno_sunburst_NEMO$`HER2 status`)
Pheno_sunburst_NEMO$`HER2 status` = gsub("Positive", "HER2+", Pheno_sunburst_NEMO$`HER2 status`)
Pheno_sunburst_NEMO$`HER2 status` = gsub("Negative", "HER2-", Pheno_sunburst_NEMO$`HER2 status`)
Pheno_sunburst_NEMO = Pheno_sunburst_NEMO %>%
  dplyr::select(NEMO, `ER status`, `HER2 status`) %>%
  group_by(NEMO, `ER status`, `HER2 status`) %>%
  summarise(Counts = n()) %>%
  as.data.frame()

sunburst_coloring_NEMO = data.frame(stringsAsFactors = FALSE,
                                    colors = tolower(gplots::col2hex(c(cluster_colors, 
                                                                       "#C11D9C", "#0F1682",  "grey40",
                                                                       "#0B9EF8", "#560DA7", "mistyrose1", 
                                                                       "hotpink4", "grey40"))),
                                    labels = c(paste0("NEMO", c(1:9)),
                                               "ER-", "ER+", "Unkn ER status",
                                               "HER2-", "HER2+", "Indeterminate",
                                               "Equivocal", "Unkn HER2 status"))

sunburstDF_NEMO = as.sunburstDF(Pheno_sunburst_NEMO, value_column = "Counts", add_root = FALSE) %>%
  inner_join(sunburst_coloring_NEMO, by = "labels")

pie_NEMO = plot_ly() %>%
  add_trace(ids = sunburstDF_NEMO$ids, labels= sunburstDF_NEMO$labels, 
            parents = sunburstDF_NEMO$parents, 
            values= sunburstDF_NEMO$values, type='sunburst', branchvalues = 'total',
            insidetextorientation='radial', maxdepth = 5,
            marker = list(colors = sunburstDF_NEMO$colors)) %>%
  layout(
    grid = list(columns =1, rows = 1),
    margin = list(l = 0, r = 0, b = 0, t = 0)
  )
pie_NEMO
rm(Pheno_sunburst_NEMO, sunburstDF_NEMO, sunburst_coloring_NEMO, pie_NEMO); gc()

# Graphs ###
library(igraph)
list_aff = list(`ARS matrix graph` = final_affinity_matrix)

for (i in 1:length(list_aff)) {
  
  # Prepare the graph object
  g <- graph_from_adjacency_matrix(list_aff[[i]], 
                                   mode = "undirected", weighted = TRUE, diag = FALSE)
  g <- delete_edges(g, E(g)[weight == 0])
  E(g)$width <- sqrt(E(g)$weight) * 5  # Example transformation for visibility
  nodes_data <- data.frame(name = V(g)$name) %>%
    inner_join(clust_annot_pheno %>% dplyr::select(samID, NEMO),
               by = c("name" = "samID"))
  
  # Set NEMO as a factor for coloring
  nodes_data[[algorithm]] <- as.factor(nodes_data[[algorithm]])
  V(g)$NEMO <- nodes_data[[algorithm]] # modify `$NEMO` manually
  
  # Create a named vector that maps cluster names (paste0(algorithm, 1:9)) to colors
  cluster_color_map <- setNames(cluster_colors, paste0(algorithm, 1:9))
  
  # Assign colors to the vertices based on the NEMO cluster
  V(g)$color <- cluster_color_map[V(g)$NEMO]
  
  png(paste0(home, 
             "/Results/single_algorithm/NEMO/Supplement/",
             names(list_aff)[i], ".png"),
      width = 6000, height = 6000, res = 700)
  
  par(mar = c(2, 2, 2, 5))  # Adjust right margin to accommodate legend
  
  # Plot the graph with a layout that spreads nodes well
  plot(g, vertex.color = V(g)$color,
       edge.width = E(g)$width,
       vertex.size = 4, 
       vertex.label = NA, 
       edge.color = "gray85",
       layout = layout_with_fr(g),  # Use Fruchterman-Reingold layout
       main = "")
  
  # Add title with reduced size using title() function
  title(main = names(list_aff)[i], cex.main = 1.7)
  
  # Add a legend to the right of the plot
  legend("bottomright", 
         title="Node Color Legend",    
         legend=c(paste0(algorithm, "1"),
                  paste0(algorithm, "2")), 
         fill=cluster_colors_heatmap,  
         cex=0.7,      
         box.lwd=1)  
  
  dev.off() 
}
rm(g, nodes_data)

# Compare these NEMO results with the NEMO output from MOVICS ###
load("Results/MOVICS_baseline/MOVICS_TCGA_RNAseq-CNV-Methylation-miRNA-SNPs_eval_on_transNEO_moic.res.list.rda")
MOVICS_NEMO = moic.res.list$NEMO$clust.res

ARI_to_MOVICS_NEMO = calculate_ari_index(cluster_df1 = MOVICS_NEMO %>%
                                           dplyr::rename(Sample.ID = samID,
                                                         Cluster = clust),
                                         cluster_df2 = NEMO_clusters,
                                         sample_col = "Sample.ID",
                                         clust_col = "Cluster",
                                         suffixes = c("_MOVICS_NEMO", "_NEMO"))

NMI_to_MOVICS_NEMO = calculate_nmi_index(cluster_df1 = MOVICS_NEMO %>%
                                           dplyr::rename(Sample.ID = samID,
                                                         Cluster = clust),
                                         cluster_df2 = NEMO_clusters,
                                         sample_col = "Sample.ID",
                                         clust_col = "Cluster",
                                         suffixes = c("_MOVICS_NEMO", "_NEMO"))

# Wrap up #####
hyperparameters = list(num_neighbors_min = min(num_neighbors_range),
                       num_neighbors_max = max(num_neighbors_range),
                       num_neighbors_step = neighbor_step,
                       sigma_min = min(sigma_range),
                       sigma_max = max(sigma_range),
                       eigengap_k = num.clusters,
                       sigma_step = sigma_step,
                       optimal_N = optN,
                       optimal_sigma = optSigma,
                       conclusion = conclusion
)

# Put all parameters in a list
params = list(algorithm = algorithm, data_source = data_source, data_types = data_types,
              evaluation_source = evaluation_source, title = title, subtitle = subtitle,
              description = description, in_a_nutshell = in_a_nutshell, optk_text = optk_text,
              citation = citation, NMI_to_MOVICS = NMI_to_MOVICS, ARI_to_MOVICS = ARI_to_MOVICS,
              NMI_to_MOVICS_NEMO = NMI_to_MOVICS_NEMO, ARI_to_MOVICS_NEMO = ARI_to_MOVICS_NEMO,
              hyperparameters = hyperparameters, ground_truth_k = ground_truth_k,
              sessionInfo = sessionInfo(), home = home, transNEO_var2comp = transNEO_var2comp)

# Render the R Markdown document with the parameters
rmarkdown::render(paste0(getwd(), "/Results/single_algorithm/NEMO/NEMO_report.Rmd"), 
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
