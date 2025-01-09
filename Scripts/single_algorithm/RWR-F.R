# Source RWR-F code
source("Scripts/single_algorithm/RWR-F source/RWR-F source.R")
library(SNFtool)
library(clValid)       # Provides the 'dunn' index
library(kernlab)       # For spectral clustering (specc)
library(foreach)
library(doParallel)

##########################
##   Parallelization    ##
##########################
cl <- makeCluster(8)  
registerDoParallel(cl) 

##########################
##   Set parameters     ##
##########################
subpath <- "MY_PROJECT"       # Adjust as needed
CLUSTER_NUM_LIST <- 2:10

# (Hyper)parameters for computing affinities:
K_neighbors <- 20       # e.g. number of neighbors for the affinityMatrix
alpha       <- 0.5      # e.g. parameter for the affinityMatrix
T_snf       <- 20       # e.g. number of iterations in SNF

# If your RWR_fusion functions allow any parameters like restart probability,
# alpha, or number of iterations, set them here. (Example placeholders)
rwr_restart_prob <- 0.7
rwr_max_iter     <- 100

# Example survival data & placeholders
# survival_data <- read.table("clinicalData.txt", header=TRUE, sep="\t", ...)
# Suppose survival_data has columns: (time_in_days, event, cluster_placeholder)

###################################
##  Prepare your 5 input matrices
###################################
# input is a list of 5 matrices:
#   input[[1]] is binary
#   input[[2]] to input[[5]] are continuous
#
# Many scripts require rows = samples and columns = features.
# If your data have features in rows, you might transpose:
#   input[[k]] <- t(input[[k]])
# for each k.

for (k in seq_along(input)) {
  # Transpose if necessary; depends on how your RWR/SNF scripts are implemented
  # input[[k]] <- t(input[[k]])
  #
  # The key is that each matrix should have the same number of samples (rows).
}

##############################################
##  Construct distance & affinity matrices  ##
##############################################
distL <- list()

## First slot: binary data
distL[[1]] <- as.matrix(dist(input[[1]], method = "binary"))

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

for (CLUSTER_NUM in CLUSTER_NUM_LIST) {
  message("Clustering for K = ", CLUSTER_NUM)
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
  tmp_result[[1]] <- specc(similarity_fusion_list[[1]], centers = CLUSTER_NUM)@.Data
  tmp_result[[2]] <- specc(similarity_fusion_list[[2]], centers = CLUSTER_NUM)@.Data
  tmp_result[[3]] <- spectralClustering(similarity_fusion_list[[3]], CLUSTER_NUM)
  tmp_result[[4]] <- spectralClustering(similarity_fusion_list[[4]], CLUSTER_NUM)
  for (i in 5:length(similarity_fusion_list)) {
    tmp_result[[i]] <- specc(similarity_fusion_list[[i]], centers = CLUSTER_NUM)@.Data
  }
  
  cluster_list[[CLUSTER_NUM]] <- tmp_result
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
                        ncol = 2 * length(CLUSTER_NUM_LIST))
colnames_vec <- c()
for (k in CLUSTER_NUM_LIST) {
  colnames_vec <- c(colnames_vec, paste0("p-value_", k), paste0("Dunn_", k))
}
colnames(result_matrix) <- colnames_vec
rownames(result_matrix) <- paste0("Method_", 1:num_fusion_methods)

# Evaluate each method
for (CLUSTER_NUM in CLUSTER_NUM_LIST) {
  for (i in seq_along(similarity_fusion_list)) {
    # Assign cluster labels to survival_data
    survival_data[, "Cluster"] <- cluster_list[[CLUSTER_NUM]][[i]]
    
    # Fit survival
    y   <- Surv(time = as.numeric(survival_data[,1]),
                event = as.numeric(survival_data[,2]))
    sdf <- survdiff(y ~ as.character(survival_data[,"Cluster"]))
    
    # p-value
    p_val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
    result_matrix[i, paste0("p-value_", CLUSTER_NUM)] <- p_val
    
    # Dunn index (requires a distance matrix, so pick appropriately)
    # If i is your fused method, pick the corresponding distance
    # For example, if i = 1 or 2 => use something akin to -log(similarity),
    # if i = 3 => SNF, etc.
    # You may have to track your distances in a separate list as in the original code.
    
    # Example:
    #   result_matrix[i, paste0("Dunn_", CLUSTER_NUM)] <- dunn(distX, cluster_list[[CLUSTER_NUM]][[i]])
    
    rm(y, sdf)
  }
}

# Write output
write.table(result_matrix,
            file = paste0("result_matrix_", subpath, ".txt"),
            sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)

################################
##  Stop the parallel cluster ##
################################
stopCluster(cl)
