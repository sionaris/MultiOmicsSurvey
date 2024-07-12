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
evaluation_source = "PARTNER" # e.g. PARTNER, transNEO-PARTNER 
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
rm(rogue_indices); gc()

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

# Create a query to retrieve clinical data for the specified samples
query <- GDCquery(
  project = "TCGA-BRCA",  # Replace with the appropriate TCGA project ID
  data.category = "Clinical",
  barcode = tcga_samples
)

# Execute the query
GDCdownload(query)

# Prepare the clinical data
clinical_data = GDCprepare_clinic(query, clinical.info = "patient")

# Setup ###
# Hyperparameter tuning
sigma_step = 0.1
iter_step = 5
num_neighbors_range = 10:30		    # number of neighbors, usually (10~30)
sigma_range = seq(0.3, 0.8, sigma_step) 	# hyperparameter, usually (0.3~0.8)
iterations = seq(10, 100, iter_step)      # Number of Iterations, usually (10~20)

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
for (modality in modalities) {
  affinity_mod = list()
  for (nn in num_neighbors_range) {
    for (sigma in sigma_range) {
      aff_mat = affinityMatrix(as.matrix(inputs_dists[[modality]]), K = nn, sigma = sigma)
      run_name = paste0("num_neighbors = ", nn, ", sigma = ", sigma)
      affinity_mod[[run_name]] = list(num_neighbors = nn,
                                      regularization = sigma,
                                      affinity_matrix = aff_mat)
      rm(run_name, aff_mat)
    }
  }
  affinity_object[[modality]] = affinity_mod
  rm(affinity_mod)
}
names(affinity_object) = modalities

# Fusions ###
Fusions = list()

# Extract total number of fusions to be performed
num_fusions = as.numeric(length(sigma_range) * length(num_neighbors_range))

# Create loop that will create lists of affinity matrices and run fusions
for (i in 1:num_fusions) {
  iter_matrices = lapply(lapply(affinity_object, `[[`, i), `[[`, "affinity_matrix")
  for (nn in num_neighbors_range) {
    for (iter in iterations) {
      run_name = paste0("num_neighbors = ", nn, ", n_iter = ", iter)
      Fusions[[run_name]] = SNF(iter_matrices, K = nn, t = iter, parallel = FALSE)
      rm(run_name)
    }
  }
}
rm(iter_matrices); gc()

Fusion = SNF(SNF_affinity_matrices, K = K, t = t)

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
                       sigma_min = min(sigma_range),
                       sigma_max = max(sigma_range),
                       sigma_step = sigma_step,
                       iter_min = min(iterations),
                       iter_max = max(iterations),
                       iter_step = iter_step,
                       optimal_N = optN,
                       optimal_sigma = optSigma,
                       optimal_iter = optIter)

# Put all parameters in a list
params = list(algorithm = algorithm, data_source = data_source, data_types = data_types,
              evaluation_source = evaluation_source, title = title, subtitle = subtitle,
              description = description, in_a_nutshell = in_a_nutshell, optk_text = optk_text,
              hyperparams_text = generate_hyperparams_text(algorithm = algorithm,
                                                           hyperparameters = hyperparameters))

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
