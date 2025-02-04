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
maxK = 5 # max individual modality k
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
timestamp() # ~1 min
gc()

# Determine best global k using lmmkmeans
library(cluster)
library(foreach)
library(doParallel)
library(doRNG)       # for reproducible parallel loops

kRange    <- 2:5 # range of individual omic clusters
globalK   <- 2:10 # range of global clusters
nSamples  <- dim(allCM[[1]])[1] # 625
nDatasets <- dim(allCM[[1]])[3] # 5

kCombos <- expand.grid(
  k1 = kRange,
  k2 = kRange,
  k3 = kRange,
  k4 = kRange,
  k5 = kRange
)
cat("Total combos of (k1..k5):", nrow(kCombos), "\n") # 1024

# Set local kernel k-means parameters
km_parameters <- list()
km_parameters$iteration_count <- 100

bestSil        <- -Inf
bestCombo      <- NULL
bestGlobalK    <- NULL
bestClustering <- NULL

# Loop over each per-dataset k combination
timestamp()
for (comboIdx in seq_len(nrow(kCombos))) {
  kvals <- as.numeric(kCombos[comboIdx, ])
  
  # Build the 3D array (nSamples x nSamples x nDatasets) 
  # by extracting from 'allCM' for each dataset i
  CMcombo <- array(0, dim = c(nSamples, nSamples, nDatasets))
  for (i in seq_len(nDatasets)) {
    # If k_i = 2 => index in allCM is (2 - 2 + 1) = 1, if k_i = 3 => 2, etc.
    idxInAllCM <- kvals[i] - 2 + 1
    CMcombo[, , i] <- allCM[[idxInAllCM]][, , i]
  }
  
  cl <- makeCluster(9)
  registerDoParallel(cl)
  
  # Ensure reproducibility across workers:
  registerDoRNG(seed = 123)
  
  # Parallelize the inner loop over global K
  results_gk <- foreach(gk = globalK,
                        .combine = rbind,
                        .packages = c("klic","cluster")) %dopar% {
                          # Copy the parameters locally
                          local_params <- km_parameters
                          local_params$cluster_count <- gk
                          
                          # Run local kernel k-means
                          res <- klic::lmkkmeans(CMcombo, local_params)
                          
                          # Weighted kernel
                          WKM <- matrix(0, nrow = nSamples, ncol = nSamples)
                          for (j in seq_len(nDatasets)) {
                            WKM <- WKM + (res$Theta[, j] %*% t(res$Theta[, j])) * CMcombo[,, j]
                          }
                          
                          # Convert to dissimilarity => (1 - similarity)
                          diss_mat <- 1 - WKM
                          dd <- as.dist(diss_mat)
                          
                          # Silhouette
                          sil_obj <- silhouette(res$clustering, dd)
                          avgSil  <- summary(sil_obj)$avg.width
                          
                          # Return numeric row: (globalK, silhouette)
                          c(gkVal = gk, silhouette = avgSil)
                        }
  
  stopCluster(cl)
  
  # Among the 9 tested globalK's, find which yields best silhouette
  bestLocalIdx <- which.max(results_gk[, "silhouette"])
  localBestSil <- results_gk[bestLocalIdx, "silhouette"]
  bestLocalGk  <- results_gk[bestLocalIdx, "gkVal"]
  
  # Update global best if improved
  if (localBestSil > bestSil) {
    # Re-run to store cluster labels
    km_parameters$cluster_count <- bestLocalGk
    res <- klic::lmkkmeans(CMcombo, km_parameters)
    
    bestSil        <- localBestSil
    bestCombo      <- kvals
    bestGlobalK    <- bestLocalGk
    bestClustering <- res$clustering
  }
  
  cat(sprintf("Done combo %d/%d, local best: K=%d, silhouette=%.4f\n",
              comboIdx, nrow(kCombos), bestLocalGk, localBestSil))
}
timestamp()
# Stop the cluster
stopCluster(cl)

# Print or save final results
cat("======= FINAL RESULTS ========\n")
cat("Overall best silhouette =", bestSil, "\n")
cat("Best combo of dataset-level k's =", paste(bestCombo, collapse=", "), "\n")
cat("Best global k =", bestGlobalK, "\n")

# Export
final_KLIC = list(
  bestSil        = bestSil,
  bestCombo      = bestCombo,
  bestGlobalK    = bestGlobalK,
  bestClustering = bestClustering
)
saveRDS(
  final_KLIC,
  file = "Results/single_algorithm/KLIC/KLIC_finalResults.rds"
)

# Optimal k and cluster df
optk = bestGlobalK
KLIC_clusters = as.data.frame(list(Sample.ID = rownames(input$SNPs), Cluster = bestClustering))

# Main results ###
# Examine cluster similarity to MOVICS by measuring NMI and ARI indices #####
# (Jaccard may be misleading)

# Calculate ARI and NMI
library(mclust)
library(clue)

ARI_to_MOVICS = calculate_ari_index(cluster_df1 = ground_truth_labels,
                                    cluster_df2 = KLIC_clusters,
                                    sample_col = "Sample.ID",
                                    clust_col = "Cluster",
                                    suffixes = c("_MOVICS_CS",
                                                 paste0("_", algorithm)))

NMI_to_MOVICS = calculate_nmi_index(cluster_df1 = ground_truth_labels,
                                    cluster_df2 = KLIC_clusters,
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
  inner_join(KLIC_clusters, by = "Sample.ID") %>%
  tibble::column_to_rownames(var = "Sample.ID") %>%
  mutate(KLIC = paste0(algorithm, Cluster)) %>%
  dplyr::select(KLIC, everything()) %>%
  dplyr::select(-Cluster)
rm(scheme); gc()

# Compute a similarity matrix based on the methodology
final_kvals = bestCombo
final_CM = array(0, dim = c(nSamples, nSamples, nDatasets))
for (i in seq_len(nDatasets)) {
  # If k_i = 2 => index in allCM is (2 - 2 + 1) = 1, if k_i = 3 => 2, etc.
  idxInAllCM <- final_kvals[i] - 2 + 1
  final_CM[, , i] <- allCM[[idxInAllCM]][, , i]
}

RNGversion("4.2.2")
set.seed(123)

local_params = list(iteration_count = 100, cluster_count = 2)
final_res = klic::lmkkmeans(final_CM, local_params)
final_WKM <- matrix(0, nrow = nSamples, ncol = nSamples)
for (j in seq_len(nDatasets)) {
  final_WKM <- final_WKM + (res$Theta[, j] %*% t(res$Theta[, j])) * final_CM[,, j]
}
dimnames(final_WKM) = list(rownames(input$SNPs), rownames(input$SNPs))

# Silhouette
sil = compute_silhouette(cluster_df = KLIC_clusters %>% dplyr::rename(samID = Sample.ID),
                         similarity_matrix = final_WKM,
                         normalize_matrix = TRUE)

getSilhouette_ggplot(sil      = sil,
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

plot_object = list(clust.res = KLIC_clusters %>%
                     dplyr::rename(samID = Sample.ID, clust = Cluster))

# Export consensus clustering object
clust = as.data.frame(plot_object$clust.res)
colnames(clust) = c("Sample.ID", "Cluster")
clust$Cluster = paste0(algorithm, clust$Cluster)
openxlsx::write.xlsx(clust, paste0(home, "/Results/single_algorithm/", algorithm, "/", 
                                   algorithm, "_", data_source, "_",
                                   data_types, "_eval_on_", evaluation_source,
                                   "_clusterings.xlsx"))
