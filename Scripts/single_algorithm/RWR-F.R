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
algorithm = "SNF"
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

# However SNF prefers features in columns so we transpose the matrices.

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
source("Scripts/single_algorithm/RWR-F source/RWR-F source.R")
library(SNFtool)
library(clValid)       # Provides the 'dunn' index
library(kernlab)       # For spectral clustering (specc)
library(foreach)
library(doParallel)

# Parallel set-up (the RWR functions contain parallelization within them)
cl <- makeCluster(8)  
registerDoParallel(cl) 

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


timestamp()
fusion_results <- list()
for (num_neighbors in num_neighbors_range) {
  for (sigma in sigma_range) {
    run_name = paste0("NN = ", num_neighbors, " sigma = ", sigma)
    # Construct affinity matrices for SNF
    affinityL <- list()
    for (modnm in names(distL)) {
      aff <- affinityMatrix(distL[[modnm]], K = num_neighbors, sigma = sigma)
      affinityL[[modnm]] <- aff
    }
    
    # RWRF fusion (fixed hyperparams)
    fused_rwrf <- RWR_fusion(
      sim_list      = affinityL,
      iteration_max = RWR_iteration_max,
      gama          = RWRF_gamma_fixed
    )
    
    # RWRNF fusion (fixed hyperparams)
    fused_rwrnf <- RWR_fusion_neighbor(
      sim_list      = affinityL,
      iteration_max = RWR_iteration_max,
      gama          = RWRNF_gamma_fixed,
      neighbor_num  = RWRNF_num_neighbors_fixed,
      alpha         = RWRNF_alpha_fixed,
      beta          = RWRNF_beta_fixed
    )
    
    # Store
    fusion_results[[run_name]] <- list(
      aff_num_neighbors  = num_neighbors,
      aff_sigma          = sigma,
      affinity           = affinityL,
      fused_rwrf         = fused_rwrf,
      fused_rwrnf        = fused_rwrnf
    )
  }
  cat("Done for combo ", run_name, ".")
}
timestamp()
rm(run_name, affinityL, fused_rwrf, fused_rwrnf)
gc()

## Next four slots: continuous data
for (k in 2:5) {
  # Option 1: If you use SNFtool's dist2 function:
  #   distL[[k]] <- dist2(input[[k]], input[[k]])
  #
  # Option 2: If you rely on base R's dist for Euclidean distance:
  tmp <- as.matrix(dist(input[[k]], method = "euclidean"))
  distL[[k]] <- tmp
}

## Convert distances to similarities
affinityL <- list()
for (k in 1:5) {
  # In the original script, SNFtool uses:
  #   affinityMatrix(Dist, K=20, alpha=0.5)
  # If the distance matrix is raw, you can do:
  aff <- affinityMatrix(distL[[k]], K_neighbors, alpha)
  affinityL[[k]] <- aff
}

##################################################
##  Compute fusion via RWR, RWR_neighbor, and SNF
##################################################
similarity_fusion_list <- list()

# 1) RWR_fusion
temp <- RWR_fusion(sim_list = affinityL
                   # Possibly pass your RWR hyperparameters, e.g.:
                   # restart_prob = rwr_restart_prob,
                   # max_iter     = rwr_max_iter
)
similarity_fusion_list <- c(similarity_fusion_list, list(temp))
rm(temp)

# 2) RWR_fusion_neighbor
temp <- RWR_fusion_neighbor(sim_list = affinityL
                            # Possibly pass neighbor-based parameters
)
similarity_fusion_list <- c(similarity_fusion_list, list(temp))
rm(temp)

# 3) SNF
time1 <- Sys.time()
temp  <- SNF(affinityL, K_neighbors, T_snf)  # K=20, T=20 by default in your script
time2 <- Sys.time()
print(time2 - time1)
similarity_fusion_list <- c(similarity_fusion_list, list(temp))
rm(temp)

# 4) Optionally include a simple concatenation or any other method
#    E.g., "Concatenation" from your original code
#    If you'd like it, build it similarly:
data_concat <- do.call(cbind, input)
dist_concat <- as.matrix(dist(data_concat, method="euclidean"))
# Rescale to [0,1] if needed:
dist_concat <- (dist_concat - min(dist_concat)) / (max(dist_concat) - min(dist_concat))
aff_concat  <- 1 - dist_concat
similarity_fusion_list <- c(similarity_fusion_list, list(aff_concat))

# 5) You could optionally add the raw affinity of each data type again
for (k in 1:5) {
  similarity_fusion_list <- c(similarity_fusion_list, list(affinityL[[k]]))
}

# If needed, remove diagonal self-similarities
for (i in seq_along(similarity_fusion_list)) {
  diag(similarity_fusion_list[[i]]) <- 0
}

######################################################
##  Clustering on each fused / single-similarity matrix
######################################################
set.seed(7)
cluster_list <- list()

for (k in k_range) {
  message("Clustering for K = ", k)
  tmp_result <- list()
  
  # Example of calling spectral clustering from kernlab or SNFtool
  # Indices:
  #   1 => RWR_fusion
  #   2 => RWR_fusion_neighbor
  #   3 => SNF
  #   4 => Concatenation
  #   5..9 => single affinity matrices
  #
  # Modify indices as needed, since you now have more or fewer slots in similarity_fusion_list
  tmp_result[[1]] <- specc(similarity_fusion_list[[1]], centers = k)@.Data
  tmp_result[[2]] <- specc(similarity_fusion_list[[2]], centers = k)@.Data
  tmp_result[[3]] <- spectralClustering(similarity_fusion_list[[3]], k)
  tmp_result[[4]] <- spectralClustering(similarity_fusion_list[[4]], k)
  for (i in 5:length(similarity_fusion_list)) {
    tmp_result[[i]] <- specc(similarity_fusion_list[[i]], centers = k)@.Data
  }
  
  cluster_list[[k]] <- tmp_result
}

######################################################
##   Survival analysis & cluster validation (example)
######################################################
# Prepare your survival_data accordingly:
# survival_data[, 1] = survival_data[, 1] / 30  # Example: convert days to months
# survival_data[, 3] = cluster labels placeholder

# result_matrix can be expanded to hold the p-value and Dunn index
num_fusion_methods <- length(similarity_fusion_list)
result_matrix <- matrix(0, 
                        nrow = num_fusion_methods, 
                        ncol = 2 * length(k_range))
colnames_vec <- c()
for (k in k_range) {
  colnames_vec <- c(colnames_vec, paste0("p-value_", k), paste0("Dunn_", k))
}
colnames(result_matrix) <- colnames_vec
rownames(result_matrix) <- paste0("Method_", 1:num_fusion_methods)

# Evaluate each method
for (k in k_range) {
  for (i in seq_along(similarity_fusion_list)) {
    # Assign cluster labels to survival_data
    survival_data[, "Cluster"] <- cluster_list[[k]][[i]]
    
    # Fit survival
    y   <- Surv(time = as.numeric(survival_data[,1]),
                event = as.numeric(survival_data[,2]))
    sdf <- survdiff(y ~ as.character(survival_data[,"Cluster"]))
    
    # p-value
    p_val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
    result_matrix[i, paste0("p-value_", k)] <- p_val
    
    # Dunn index (requires a distance matrix, so pick appropriately)
    # If i is your fused method, pick the corresponding distance
    # For example, if i = 1 or 2 => use something akin to -log(similarity),
    # if i = 3 => SNF, etc.
    # You may have to track your distances in a separate list as in the original code.
    
    # Example:
    #   result_matrix[i, paste0("Dunn_", k)] <- dunn(distX, cluster_list[[k]][[i]])
    
    rm(y, sdf)
  }
}

# Write output
write.table(result_matrix,
            file = paste0("Resources/RWR-F_result_matrix_", subpath, ".txt"),
            sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)

################################
##  Stop the parallel cluster ##
################################
stopCluster(cl)
