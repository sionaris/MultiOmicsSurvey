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
algorithm = "MOFA"
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
optk_boolean = "FALSE" # either TRUE or FALSE. Answers whether the algorithm suggests an optimal k
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
library(MOFA2)

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

# Hyperaparameters ###
factor_candidates <- c(2:10, 15, 20, 25)
convergence_mode = "slow"
dropR2 = -1
iter = 20000
freqELBO = 5
startELBO = 1

# Import MOFA output objects #####
library(rhdf5)

# Example output object high-level inspection
h5ls("Python/MOFA/MOFA_output_factors9_iter20000_R2_-1.hdf5")

# Load results in list
mofa_results <- list()
for (nf in factor_candidates) {
  h5file <- sprintf("Python/MOFA/MOFA_output_factors%d_iter20000_R2_-1.hdf5", nf)
  Z_raw         <- h5read(h5file, "/expectations/Z")$group0
  r2_per_factor <- h5read(h5file, "/variance_explained/r2_per_factor")$group0
  r2_total      <- h5read(h5file, "/variance_explained/r2_total")$group0
  sample_names  <- h5read(h5file, "/samples")$group0
  view_names    <- h5read(h5file, "/views/views")
  
  # Add row/column names for clarity
  rownames(Z_raw) <- sample_names
  colnames(Z_raw) <- paste0("Factor", seq_len(nf))
  
  rownames(r2_per_factor) <- paste0("Factor", seq_len(nrow(r2_per_factor)))
  colnames(r2_per_factor) <- view_names
  
  names(r2_total) <- view_names
  
  # Store everything in a nested list object
  mofa_results[[paste0("num.factors = ", nf)]] <- list(
    n.factors         = nf,
    Z                 = Z_raw,
    r2_per_factor     = r2_per_factor,
    r2_total          = r2_total,
    sample_names      = sample_names,
    view_names        = view_names
  )
}
rm(h5file, Z_raw, r2_per_factor, r2_total, sample_names, view_names, nf); gc()

# Elbow plot
# Initialize a numeric vector to store total variance explained for each K
total_var_vec <- numeric(length(factor_candidates))

for (i in seq_along(factor_candidates)) {
  k_val <- factor_candidates[i]
  r2_tot <- mofa_results[[paste0("num.factors = ", k_val)]]$r2_total
  
  # Sum across all views
  total_var_vec[i] <- sum(r2_tot)
}
rm(i, k_val, r2_tot)

plot(factor_candidates, total_var_vec, type = "b",
     xlab = "Number of Factors",
     ylab = "Sum of total variance explained across all views",
     main = "Elbow Plot: MOFA Total Variance Explained")

# There is an elbow at 7
views <- mofa_results[["num.factors = 2"]]$view_names
num_views <- length(views)

# Make a matrix to store the total R2 for each view across all K
r2_by_view <- matrix(NA, nrow = length(factor_candidates), ncol = num_views,
                     dimnames = list(paste0("num.factors = ", factor_candidates), views))

for (i in seq_along(factor_candidates)) {
  k_val <- factor_candidates[i]
  r2_view_i <- mofa_results[[paste0("num.factors = ", k_val)]]$r2_total
  r2_by_view[i, ] <- r2_view_i
}
rm(i, k_val, r2_view_i)

# Per view plot (total R2)
matplot(factor_candidates, r2_by_view, type = "b", pch = 19,
        xlab = "Number of Factors",
        ylab = "Total R2 explained (per view)",
        main = "View-specific total variance explained by MOFA factors")
legend("bottomright", legend = colnames(r2_by_view), lty = 1:ncol(r2_by_view),
       pch = 19, col = 1:ncol(r2_by_view), cex = 0.8)

# Correlations
for (k_val in factor_candidates) {
  cat("\nExamining factor correlations for K=", k_val, ":\n")
  
  Z_mat <- mofa_results[[paste0("num.factors = ", k_val)]]$Z
  cor_mat <- cor(Z_mat)
  
  print(round(cor_mat, 2))
  # Optionally visualize a heatmap
  heatmap(cor_mat, main = paste0("Factor correlation: n.f. = ", k_val))
}
rm(k_val, Z_mat, cor_mat); gc()
