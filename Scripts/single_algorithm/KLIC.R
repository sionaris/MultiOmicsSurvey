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
algorithm = "KLIC"
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

# However KLIC prefers features in columns so we transpose the matrices.

# Extract the names of the modalities that will be used
modalities = unlist(strsplit(data_types, "-"))

# Replace with TRUE wherever features are in rows
features_in_rows = rep(TRUE, length(modalities))

# Run algorithm #####
library(coca)
library(klic)

# IMPORTANT: klic requires Rmosek. To install it in Windows (trickier as always):
# Firstly download and install Rtools (if you don't have it already)
# Download the Rmosek distr from here: https://www.mosek.com/downloads/
# Unzip it in C:/Users/username so that there is a
# C:/Users/username/mosek directory
# Add the following to your .Rprofile:

# # Mosek
# mosek_path = "C:/Users/username/mosek/10.2/tools/platform/win64x86/bin"
# 
# if (!mosek_path %in% paths) {
#   new_path <- c(mosek_path, paths)
#   Sys.setenv(PATH = paste(new_path, collapse = .Platform$path.sep))
# }

# rm(new_path, mosek_path)
# gc()

# Exit and then re-open RStudio

# Run the following:
# source("C:/Users/username/mosek/10.2/tools/platform/win64x86/rmosek/builder.R")
# attachbuilder()
# install.rmosek()

# If a license is required, that can be obtained quickly for academics from:
# https://www.mosek.com/products/academic-licenses/
# once obtained, it can be placed in C:/Users/username/mosek/

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

# Export input object for HPC run
saveRDS(input, "Resources/KLIC_input.rds")

# Modality types
continuous = c("RNAseq", "CNV", "Methylation", "miRNA")
categorical = c("SNPs")

# Calculate the pair-wise distance (Euclidean for continuous modalities)
input_dists = lapply(input[continuous], function(x) {
  x = as.matrix(x)
  x = SNFtool::dist2(x, x)
})

# Binary for SNPs (see ?dist for details)
input_dists[["SNPs"]] = as.matrix(dist(as.matrix(input$SNPs),
                                       as.matrix(input$SNPs),
                                       method = "binary"))
gc()

saveRDS(input_dists, "Resources/KLIC_distances.rds")

# Hyperparameters and tuning ###
maxK = 10
B = 250
pItem = 0.8

# Build consensus matrices with COCA
allCM = list()
nSamples <- nrow(input[[1]])
nDatasets <- length(input)

timestamp()
for (k in 2:maxK) {
  # 1) Initialize CM for this k
  CM <- array(NA, dim = c(nSamples, nSamples, nDatasets))
  
  # 2) Loop over datasets
  for (i in 1:nDatasets) {
    temp_mat <- coca_cc_mod(
      dist         = input_dists[[i]],
      clMethod     = "hclust",
      B            = B,
      K            = k,
      pItem        = pItem,
      hclustMethod = "average"
    )
    
    # 3) Shift eigenvalues to ensure PSD
    CM[, , i] <- klic::spectrumShift(temp_mat, verbose = FALSE)
  }
  allCM[[paste0("k = ", k)]] = CM
}
timestamp() # ~3 min
gc()

# Determine best global k using lmmkmeans
library(cluster)

kRange    <- 2:10 # local
globalK   <- 2:10 # global

# Create all combinations of per-dataset k
# For M=5, this is {2..10}^5 = 59,049 combos
kCombos <- expand.grid(
  k1 = kRange,
  k2 = kRange,
  k3 = kRange,
  k4 = kRange,
  k5 = kRange
)

# Set up local kernel k-means parameters
km_parameters <- list()
km_parameters$iteration_count <- 100  # max iterations

bestSil        <- -Inf
bestCombo      <- NULL  # will store the (k1, k2, ..., k5)
bestGlobalK    <- NULL
bestClustering <- NULL

# Main loop over combos
for (comboIdx in seq_len(nrow(kCombos))) {
  kvals <- as.numeric(kCombos[comboIdx, ])  # e.g. c(k1, k2, k3, k4, k5)
  
  # Build a 3D array of dimension (nSamples x nSamples x nDatasets)
  # by picking the co-clustering matrix for dataset i from allCM[[k_i - 1]]
  CMcombo <- array(0, dim = c(nSamples, nSamples, nDatasets))
  for (i in seq_len(nDatasets)) {
    # k_i is kvals[i]
    # allCM[[1]] => k=2, allCM[[2]] => k=3, so the index in allCM is (k_i - 1).
    CMcombo[, , i] <- allCM[[kvals[i] - 1]][, , i]
  }
  
  # Now try each global k in 2..10
  for (gk in globalK) {
    km_parameters$cluster_count <- gk
    
    # Run local kernel k-means
    res <- klic::lmkkmeans(CMcombo, km_parameters)
    
    # Weighted kernel: WKM = sum_{dataset j} (Theta[, j] %*% t(Theta[, j])) * CMcombo[,, j]
    WKM <- matrix(0, nrow = nSamples, ncol = nSamples)
    for (j in seq_len(nDatasets)) {
      WKM <- WKM + (res$Theta[, j] %*% t(res$Theta[, j])) * CMcombo[,, j]
    }
    
    # Compute silhouette => we must convert WKM (similarities) to a dist/dissimilarity
    #  Simple approach: dissimilarity = 1 - similarity (since WKM in [0,1]).
    diss_mat <- 1 - WKM
    dd <- as.dist(diss_mat)
    
    sil_obj <- cluster::silhouette(res$clustering, dd)
    avgSil  <- summary(sil_obj)$avg.width
    
    # Check if this is the best silhouette so far
    if (avgSil > bestSil) {
      bestSil        <- avgSil
      bestCombo      <- kvals      # (k1, k2, k3, k4, k5)
      bestGlobalK    <- gk
      bestClustering <- res$clustering
    }
  }
  cat("Done for ", comboIdx)
}

# Inspect results
bestSil
#> e.g. 0.55 ...
bestCombo
#> e.g. c(2, 2, 4, 10, 9)
bestGlobalK
#> e.g. 4
#> 
#> 


# Parallel version
library(foreach)
library(doParallel)
library(doRNG)       # for reproducible parallel loops

# Suppose you have:
# 1. 'allCM' = list of length 9 => co-clustering mats for k=2..10
# 2. 'kCombos' = expand.grid(...) => all combos of per-dataset k
# 3. 'nSamples' = number of samples (e.g., 625)
# 4. 'nDatasets' = 5 (for 5 omics)
# 5. 'globalK' = 2:10

# Prepare a cluster of 9 cores (match 'globalK' length or as you wish)
nCores <- 9
cl <- makeCluster(nCores)
registerDoParallel(cl)

# Optionally fix a seed for reproducibility in parallel:
registerDoRNG(seed = 123)

# Local kernel k-means parameters
km_parameters <- list()
km_parameters$iteration_count <- 100

bestSil        <- -Inf
bestCombo      <- NULL  # will store the best (k1, k2, ..., k5)
bestGlobalK    <- NULL
bestClustering <- NULL

# Main loop over combos
for (comboIdx in seq_len(nrow(kCombos))) {
  kvals <- as.numeric(kCombos[comboIdx, ])  # e.g. c(k1, k2, k3, k4, k5)
  
  # 1) Build a 3D array (nSamples x nSamples x nDatasets)
  #    by picking the co-clustering matrices from allCM
  #    allCM[[1]] => k=2, allCM[[2]] => k=3, ...
  CMcombo <- array(NA, dim = c(nSamples, nSamples, nDatasets))
  for (i in seq_len(nDatasets)) {
    # k_i in {2..10}, index in allCM is (k_i - 2 + 1) => (k_i - 1)
    CMcombo[, , i] <- allCM[[kvals[i] - 1]][, , i]
  }
  
  # 2) Parallelize the loop over globalK
  #    We gather results in a data.frame or matrix
  #    using .combine='rbind' so we can pick the best silhouette
  results_gk <- foreach(gk = globalK, 
                        .combine = rbind,        # row-bind results
                        .packages = c("klic","cluster")) %dopar% {
                          local_params <- km_parameters
                          local_params$cluster_count <- gk
                          
                          # 2a) Run local kernel k-means
                          res <- klic::lmkkmeans(CMcombo, local_params)
                          
                          # 2b) Weighted kernel: sum_{j} (Theta[, j] %*% t(Theta[, j])) * CMcombo[,, j]
                          WKM <- matrix(0, nrow = nSamples, ncol = nSamples)
                          for (j in seq_len(nDatasets)) {
                            WKM <- WKM + (res$Theta[, j] %*% t(res$Theta[, j])) * CMcombo[,, j]
                          }
                          
                          # 2c) Silhouette => transform WKM => dissimilarity
                          #     If WKM is 0..1, do 1 - WKM
                          diss_mat <- 1 - WKM
                          dd       <- as.dist(diss_mat)
                          
                          sil_obj <- silhouette(res$clustering, dd)
                          avgSil  <- summary(sil_obj)$avg.width
                          
                          # Return info: [gk, silhouette, plus if needed the cluster labels, etc.]
                          # We only do numeric columns in .combine='rbind'; store cluster as NA or skip
                          c(global_k = gk, avgSil = avgSil)
                        }
  
  # 3) among the 9 globalK's tested, find which yields the best silhouette
  bestIndex   <- which.max(results_gk[,"avgSil"])
  gk_best     <- results_gk[bestIndex, "global_k"]
  localBestSil <- results_gk[bestIndex, "avgSil"]
  
  # If that best is better than our overall best, re-run to store final labels
  if (localBestSil > bestSil) {
    # Re-run local kernel k-means once more *in serial* or keep it from above
    # to store the final clustering. We'll do a quick re-run:
    km_parameters$cluster_count <- gk_best
    res <- klic::lmkkmeans(CMcombo, km_parameters)
    
    # Weighted kernel and silhouette re-check if you want
    # ...
    
    bestSil        <- localBestSil
    bestCombo      <- kvals
    bestGlobalK    <- gk_best
    bestClustering <- res$clustering
  }
  
  cat("Done comboIdx =", comboIdx, "=> Best local silhouette:", localBestSil, "\n")
}

# Stop cluster
stopCluster(cl); gc()

# Print final results
cat("Overall best silhouette =", bestSil, "\n")
cat("Best combo of dataset-level k's =", bestCombo, "\n")
cat("Best global k =", bestGlobalK, "\n")
