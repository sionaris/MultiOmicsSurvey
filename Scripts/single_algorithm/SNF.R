# Import data from gitignored "Resources/BRCA complete/" folder #####

# These datasets are pre-standardized
brca_cnv = read.csv("Resources/BRCA complete/BRCA_CNV.csv")
brca_met = read.csv("Resources/BRCA complete/BRCA_Methy.csv")
brca_exp = read.csv("Resources/BRCA complete/BRCA_mRNA.csv")
brca_mirna = read.csv("Resources/BRCA complete/BRCA_miRNA.csv")
data_object = list(CNV = brca_cnv, Methylation = brca_met,
            RNAseq = brca_exp, miRNA = brca_mirna)
rm(brca_cnv, brca_exp, brca_met, brca_mirna); gc()

# Setup environment variables for markdown #####

# Ensure reproducibility
RNGversion("4.2.2")
set.seed(123)

# Load custom helper functions
source("Scripts/automated_scripts/custom_functions.R")

# Preamble
home = getwd()
algorithm = "SNF"
alg_feature_pref = "cols" # Where does the algorithm expect the features to be
citation = fetch_citation(algorithm = algorithm)
data_source = "TCGA" # e.g. TCGA, TCGA-transNEO, transNEO-PARTNER
data_types = "RNAseq-CNV-Methylation-miRNA" # e.g. RNAseq, RNAseq-CNV-miRNA
evaluation_source = "transNEO" # e.g. PARTNER, transNEO-PARTNER 
title = paste0("Results from ", algorithm)
subtitle = paste0("<b>Train</b>: ", data_source, " ", data_types, " | <b>Evaluation</b>: ", evaluation_source)
in_a_nutshell = fetch_in_a_nutshell(algorithm = algorithm)
optk_boolean = "TRUE" # either TRUE or FALSE. Answers whether the algorithm suggests an optimal k
optk_text = ifelse(optk_boolean == TRUE,
                   "<u>suggests</u> an estimate of the optimal number of multi-omic clusters $k$",
                   "<u>does not suggest</u> an optimal number of multi-omic clusters $k$")

# Detailed description of the algorithm
description = paste(readLines("Resources/algorithm_descriptions/snf_description.Rmd"),
                    collapse = "\n") # File path to .Rmd file within Resources/algorithm_descriptions

# Preprocessing flags and code #####
library(stringr)
library(dplyr)

# Extract the names of the modalities that will be used
modalities = unlist(strsplit(data_types, "-"))

# The standardization boolean gets the names from the modalities vector
standardization_booleans = rep(NULL, length(modalities))

# Replace with TRUE wherever standardization is required or FALSE otherwise
standardization_booleans = rep(FALSE, length(modalities))

# The same logic follows for features in rows
features_in_rows = rep(NULL, length(modalities))

# Replace with TRUE wherever features are in rows. Each row will then be standardized
features_in_rows = rep(TRUE, length(modalities))

# If features are in columns, each column will be standardized

# The feature column vector is TRUE when features are not rownames, but a column
feature_column = rep(NULL, length(modalities))

# Replace with either a numeric value or a column name where applicable
# Leave NULL if features are in columns
# Use feature_column EVEN IF data are pre-standardized, as long as the data frame
# contains a feature column
feature_column = rep("X", length(modalities))

# Give names to the vectors
names(standardization_booleans) = names(features_in_rows) = names(feature_column) = modalities

# Perform standardization of features if required
input = data_object
for (i in 1:length(modalities)) {
  if (standardization_booleans[modalities[i]] == TRUE) {
    if (features_in_rows[i] == TRUE && !is.null(feature_column[i])) {
      
      # Check if numerical or character string was used for the feature column indicator
      feature_column_indicator_type = ifelse(is.numeric(feature_column[i]),
                                             "numeric", "character")
      
      # Extract features
      features = data_object[[i]][, feature_column[i]]
      
      # Remove feature column based on indicator type
      if (feature_column_indicator_type == "numeric") {
        z_data = data_object[[i]][, -feature_column[i]]
      } else {
        z_data = data_object[[i]] %>% dplyr::select(-feature_column[i])
      }
      
      # Convert to matrix and standardize
      cols = colnames(z_data) # Extract colnames
      z_data = as.matrix(z_data)
      z_data = t(apply(z_data, 1, scale))
      rownames(z_data) = features # Set/ensure rownames
      colnames(z_data) = cols # Set/ensure colnames
      
      # Save the transformed matrix in the input object
      input[[modalities[i]]] = z_data
      rm(feature_column_indicator_type, features, z_data, cols)
    } else if (features_in_rows[i] == TRUE && is.null(feature_column[i])) {
      
      # We assume the features are the rownames
      features = rownames(data_object[[i]])
      cols = colnames(data_object[[i]])
      z_data = as.matrix(data_object[[i]])
      z_data = t(apply(z_data, 1, scale))
      rownames(z_data) = features # Set/ensure rownames
      colnames(z_data) = cols # Set/ensure colnames
      
      # Save the transformed matrix in the input object
      input[[modalities[i]]] = z_data
      rm(features, z_data, cols)
    } else if (features_in_rows[i] == FALSE) {
      
      # We assume that all columns are numeric and are to be standardized
      rows = rownames(data_object[[i]])
      features = colnames(data_object[[i]])
      z_data = as.matrix(data_object[[i]])
      z_data = apply(z_data, 1, scale)
      rownames(z_data) = rows # Set/ensure rownames
      colnames(z_data) = features # Set/ensure colnames
      
      # Save the transformed matrix in the input object
      input[[modalities[i]]] = z_data
      rm(features, z_data, rows)
    } 
  } else if (standardization_booleans[modalities[i]] == FALSE && !is.null(feature_column[i])) {
    # Only used for pre-standardized data which have the feature names in one column
    
    # Check if numerical or character string was used for the feature column indicator
    feature_column_indicator_type = ifelse(is.numeric(feature_column[i]),
                                           "numeric", "character")
    
    # Extract the rownames from the designated column
    features = data_object[[i]][, feature_column[i]]
    # Remove feature column based on indicator type
    if (feature_column_indicator_type == "numeric") {
      z_data = data_object[[i]][, -feature_column[i]]
    } else {
      z_data = data_object[[i]] %>% dplyr::select(-feature_column[i])
    }
    
    rownames(z_data) = features # Set/ensure colnames
    
    # Save the matrix in the input object
    input[[modalities[i]]] = z_data
    rm(features, z_data)
  }
}

rm(i); gc()

# Run algorithm #####
library(SNFtool)

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

# Keep samples that have measurements in all modalities
if (alg_feature_pref == "rows") {
  # Sample names are in the columns
  overlap = Reduce(intersect, lapply(input, colnames))
  
  # Filter inputs
  input = lapply(input, function(x) {
    x = x[, overlap]
  })
} else {
  # Sample names are in the rows
  overlap = Reduce(intersect, lapply(input, rownames))
  
  # Filter inputs
  input = lapply(input, function(x) {
    x = x[overlap, ]
  })
}

# Download clinical data for the TCGA samples of interest
library(TCGAbiolinks)
tcga_samples = gsub("\\.", "-", overlap)
tcga_samples = str_sub(tcga_samples, 1, 12)

# Create a query to retrieve clinical data for the specified samples
query <- GDCquery(
  project = "TCGA-BRCA",  # Replace with the appropriate TCGA project ID
  data.category = "Clinical",
  data.type = "Clinical Supplement",
  barcode = tcga_samples
)

# Execute the query
GDCdownload(query)

# Prepare the clinical data
clinical_data = GDCprepare_clinic(query, clinical.info = "patient")

# Setup ###
# Hyperparameter tuning
sigma_step = 0.1
# iter_step = 10
neighbor_step = 5
num_neighbors_range = seq(10, 30, neighbor_step) # number of neighbors, usually (10~30)
sigma_range = seq(0.3, 0.8, sigma_step) 	# hyperparameter, usually (0.3~0.8) -REFFERED as \mu in report text
# iterations = seq(10, 100, iter_step)      # Number of Iterations, usually (10~20)
n_iterations = 20     # Number of Iterations, usually (10~20)

# Here, the simulation data (Data1, Data2) has two data types. They are complementary to each other. And two data types have the same number of points. The first half data belongs to the first cluster; the rest belongs to the second cluster.
# rcb_label = multimodal_inputs$Labels$RCB.category # 1: pCR, 0: RCB-I/-II/-III
# names(rcb_label) = multimodal_inputs$Labels$Donor.ID
# rcb_label = rcb_label[overlap]

# Calculate the pair-wise distance
inputs_dists = lapply(input, function(x) {
  x = as.matrix(x)
  x = dist2(x, x)
})

# Affinity matrices ###
affinity_object = list()
for (nn in num_neighbors_range) {
  nn_list <- list()
  for (sigma in sigma_range) {
    sigma_list <- list()
    for (modality in modalities) {
      aff_mat <- affinityMatrix(as.matrix(inputs_dists[[modality]]), K = nn, sigma = sigma)
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
Fusions = list()

# # Try parallel
# library(parallel)
# 
# # Use 5 cores
# cl = makeCluster(5)
# 
# # Define the function to perform the fusion for a given combination of nn and sigma
# process_fusion <- function(params_parallel) {
#   nn <- params_parallel$nn
#   sigma <- params_parallel$sigma
#   sublist <- affinity_object[[paste0("NN = ", nn)]][[paste0("sigma = ", sigma)]]
#   iter_matrices <- lapply(sublist, `[[`, "affinity_matrix")
#   result <- list()
#   
#   for (iter in iterations) {
#     run_name <- paste0("NN = ", nn, ", sigma = ", sigma, ", n_iter = ", iter)
#     result[[run_name]] <- SNF(iter_matrices, K = nn, t = iter, parallel = FALSE)
#   }
#   # Print update after completing all iterations for the current (nn, sigma)
#   cat("Completed NN =", nn, ", sigma =", sigma, "\n")
#   
#   return(result)
# }
# 
# # Create a list of parameter combinations
# params_parallel_combinations <- expand.grid(nn = num_neighbors_range, sigma = sigma_range)
# params_parallel_list <- split(params_parallel_combinations, 
#                               seq(nrow(params_parallel_combinations)))
# 
# # Export the necessary variables and functions to the cluster
# clusterExport(cl, varlist = c("affinity_object", "iterations", "SNF"))
# 
# # Use parLapply to parallelize the process_fusion function
# results <- parLapply(cl, params_parallel_list, process_fusion)
# 
# # Stop the cluster
# stopCluster(cl)
# 
# # Combine the results into a single list
# Fusions <- do.call(c, results)

# Try parallel
library(parallel)
library(foreach)
library(doParallel)

# Use 5 cores
cl = makeCluster(6) # or 3 (depending on resources)
registerDoParallel(cl)

Fusions = list()

# Use foreach to parallelize the computation (< 30 min)
Fusions <- foreach(nn = num_neighbors_range, .combine = 'c', .packages = 'SNFtool') %:%
  foreach(sigma = sigma_range, .combine = 'c') %dopar% {
    sublist <- affinity_object[[paste0("NN = ", nn)]][[paste0("sigma = ", sigma)]]
    iter_matrices <- lapply(sublist, `[[`, "affinity_matrix")
    fusion_result <- SNF(iter_matrices, K = nn, t = n_iterations, parallel = FALSE)
    list(fusion_result)
  }

# Stop the cluster
stopCluster(cl)

# Restructure the Fusions list to match the desired output format
names(Fusions) <- unlist(lapply(num_neighbors_range, function(nn) {
  lapply(sigma_range, function(sigma) {
    paste0("NN = ", nn, ", sigma = ", sigma)
  })
}))

# # Create loop that will create lists of affinity matrices and run fusions
# for (nn in num_neighbors_range) {
#   for (sigma in sigma_range) {
#   sublist = affinity_object[[paste0("NN = ", nn)]][[paste0("sigma = ", sigma)]]
#   iter_matrices = lapply(sublist, `[[`, "affinity_matrix")
#   Fusions[[paste0("NN = ", nn)]][[paste0("sigma = ", sigma)]] = SNF(iter_matrices, 
#                                                                     K = nn, 
#                                                                     t = n_iterations,
#                                                                     parallel = FALSE)
#   rm(sublist)
#   }
# }
# rm(iter_matrices); gc()

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

# Overall tests
# Reshape data for ANOVA
nn_reshape <- reshape_SNF_Pearson_for_anova(nn_similarities)
sigma_reshape <- reshape_SNF_Pearson_for_anova(sigma_similarities)

# Perform ANOVA for nn
anova_nn <- aov(Value ~ Factor, data = nn_reshape)
summary(anova_nn)

# Perform ANOVA for sigma
anova_sigma <- aov(Value ~ Factor, data = sigma_reshape)
summary(anova_sigma)

# Conclusion
if (summary(anova_nn)[[1]][["Pr(>F)"]][1] < 0.05) {
  conclusion1 = "Overall, the choice of sigma significantly affects the results for a given nn."
  cat(conclusion1)
} else {
  conclusion1 = "Overall, the choice of sigma does not significantly affect the results for a given nn."
  cat(conclusion1)
}

if (summary(anova_sigma)[[1]][["Pr(>F)"]][1] < 0.05) {
  conclusion2 = "Overall, the choice of nn significantly affects the results for a given sigma."
  cat(conclusion2)
} else {
  conclusion2 = "Overall, the choice of nn does not significantly affect the results for a given sigma."
  cat(conclusion2)
}

# Determine overall effect
if (summary(anova_nn)[[1]][["Pr(>F)"]][1] < 0.05 &&
    summary(anova_sigma)[[1]][["Pr(>F)"]][1] < 0.05) {
  nn_mean_diff <- max(nn_summary$mean) - min(nn_summary$mean)
  sigma_mean_diff <- max(sigma_summary$mean) - min(sigma_summary$mean)
  
  # Means
  if (nn_mean_diff > sigma_mean_diff) {
    conclusion3 = paste0("The nn effect is stronger than the sigma effect based on mean",
                         " Pearson similarities (", nn_mean_diff, " vs. ", sigma_mean_diff, ").")
    cat(conclusion3)
  } else if (nn_mean_diff < sigma_mean_diff) {
    conclusion3 = paste0("The sigma effect is stronger than the nn effect based on mean",
                         " Pearson similarities (", sigma_mean_diff, " vs. ", nn_mean_diff, ").")
    cat(conclusion3)
  }
  
  # Standard deviations
  nn_sd_diff <- max(nn_summary$sd) - min(nn_summary$sd)
  sigma_sd_diff <- max(sigma_summary$sd) - min(sigma_summary$sd)
  
  if (nn_sd_diff > sigma_sd_diff) {
    conclusion4 = paste0("The nn effect is stronger than the sigma effect based on", 
                         " the standard deviation of Pearson similarities (",
                         nn_sd_diff, " vs. ", sigma_sd_diff, ").")
    cat(conclusion4)
  } else if (nn_sd_diff < sigma_sd_diff) {
    conclusion4 = paste0("The sigma effect is stronger than the nn effect based on",  
                         " the standard deviation of Pearson similarities (",
                         sigma_sd_diff, " vs. ", nn_sd_diff, ").")
    cat(conclusion4)
  }
} else {
  conclusion3 = "The nn effect and sigma effect are practically equal based on mean Pearson similarities"
  cat(conclusion3)
}

if (exists("conclusion4")) {
  conclusion = paste0(conclusion1, conclusion2, conclusion3, conclusion4, collapse = " ")
  rm(conclusion1, conclusion2, conclusion3, conclusion4)
  sig_status = TRUE
} else {
  conclusion = paste0(conclusion1, conclusion2, conclusion3, collapse = " ")
  rm(conclusion1, conclusion2, conclusion3)
  sig_status = FALSE
}

# If no significant differences are shown between/across hyperparameters then pick median values
if (sig_status) {
  # Code to pick best hyperparameters
} else {
  optN = 20 # median(num_neighbors_range)
  optSigma = 0.5 # ~median(sigma_range)
}

# Monte Carlo Consensus Clustering using M3C and optimal hyperparameters matrix
RNGversion("4.2.2")
set.seed(123)
final_affinity_matrix = Fusions[[paste0("NN = ", optN, ", sigma = ", optSigma)]]
cc_pheno = clinical_data %>%
  dplyr::select(ID = bcr_patient_barcode, histological_type, menopause_status,
                breast_carcinoma_progesterone_receptor_status,
                breast_carcinoma_estrogen_receptor_status,
                her2_immunohistochemistry_level_result)
cc_in = t(input$RNAseq)
newcols = gsub("\\.", "-", str_sub(colnames(cc_in), 1, 12))
colnames(cc_in) = newcols
cc_pheno = cc_pheno[cc_pheno$ID %in% newcols, ] %>%
  distinct(ID, .keep_all = TRUE)
cc = M3C_mod(cc_in, des = cc_pheno, iters = 100, repsref = 250, 
    repsreal = 250, seed = 123, fsize = 18, lthick = 2, dotsize = 1.25,
    my_affinity_matrix = final_affinity_matrix, clusteralg = "spectral")

## With this unified graph W of size n x n, you can do either spectral clustering or Kernel NMF.
group = spectralClustering(Fusion, 2) 	# spectral clustering with K = 2

## you can evaluate the goodness of the obtained clustering results by 
# calculating Normalized mutual information (NMI): if NMI is close to 1, it 
# indicates that the obtained clustering is very close to the "true" cluster information;
# if NMI is close to 0, it indicates the obtained clustering is not similar to the "true"
# cluster information.

displayClusters(Fusion, group)
SNFNMI = calNMI(group, rcb_label)

## you can also find the concordance between each individual network and the fused network
ConcordanceMatrix = concordanceNetworkNMI(list(Fusion,
                                               SNF_affinity_matrices[[1]],
                                               SNF_affinity_matrices[[2]],
                                               SNF_affinity_matrices[[3]],
                                               SNF_affinity_matrices[[4]]), 2)

# Export main results #####


# Evaluation #####


# Wrap up #####
hyperparameters = list(num_neighbors_min = min(num_neighbors_range),
                       num_neighbors_max = max(num_neighbors_range),
                       num_neighbors_step = neighbor_step,
                       sigma_min = min(sigma_range),
                       sigma_max = max(sigma_range),
                       sigma_step = sigma_step,
                       optimal_N = optN,
                       optimal_sigma = optSigma,
                       conclusion = conclusion
                       )

# Put all parameters in a list
params = list(algorithm = algorithm, data_source = data_source, data_types = data_types,
              evaluation_source = evaluation_source, title = title, subtitle = subtitle,
              description = description, in_a_nutshell = in_a_nutshell, optk_text = optk_text,
              hyperparams_text = generate_hyperparams_text(algorithm = algorithm,
                                                           hyperparameters = hyperparameters),
              criterion = criterion)

# Create algorithm directory if it doesn't exist
if (!dir.exists(paste0(home, "/Results/single_algorithm/", 
                       algorithm))) {
  dir.create(paste0(home, "/Results/single_algorithm/", 
                    algorithm))
}

# Render the R Markdown document with the parameters
rmarkdown::render(paste0(getwd(), "/Scripts/automated_scripts/single_algorithm_results_report.Rmd"), 
                  params = params, 
                  output_file = paste0(home, "/Results/single_algorithm/", 
                                       algorithm, "/", algorithm, "_report_",
                                       data_source, "_",
                                       data_types, "_eval_on_", evaluation_source,
                                       ".html"))

# Save environment
save.image(paste0(home, "/Results/single_algorithm/", 
                  algorithm, "/", algorithm, "_", data_source, "_",
                  data_types, "_eval_on_", evaluation_source,
                  "_env.RData"))

# Session info
# Capture the output of sessionInfo() to a variable
session_info <- capture.output(sessionInfo())

# Write the captured output to a .txt file
writeLines(session_info, paste0(home, "/Results/single_algorithm/", 
                                algorithm, "/", algorithm, "_", data_source, "_",
                                data_types, "_eval_on_", evaluation_source,
                                "_session_info.txt"))
