
# Load libraries
library(SNFtool)
library(foreach)
library(doParallel)

# Source RWR-F
source("RWR-F_source.R")

# Load input data
distL <- readRDS("RWR-F_input_dists.rds")

# Define hyperparameters
num_neighbors <- 50
sigma_val     <- 0.5

# Prepare the cluster for parallelization (up to 20 cores)
num_cores <- 20
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Construct affinity matrices
affinityL <- list()
for (modnm in names(distL)) {
  aff <- affinityMatrix(distL[[modnm]], K = num_neighbors, sigma = sigma_val)
  affinityL[[modnm]] <- aff
}

# RWRF
fused_rwrf <- RWR_fusion(
  sim_list      = affinityL,
  iteration_max = 1000,
  gama          = 0.7
)

# RWRNF
fused_rwrnf <- RWR_fusion_neighbor(
  sim_list      = affinityL,
  iteration_max = 1000,
  gama          = 0.7,
  neighbor_num  = 10,
  alpha         = 0.9,
  beta          = 0.9
)

stopCluster(cl)

# Create a fused object to save
fusion <- list(
  aff_num_neighbors  = num_neighbors,
  aff_sigma          = sigma_val,
  affinity           = affinityL,
  fused_rwrf         = fused_rwrf,
  fused_rwrnf        = fused_rwrnf
)

# Save the results
saveRDS(fusion, paste0("NN_", num_neighbors, "_sigma_", sigma_val, "_fusion.rds"))

# Also record sessionInfo
writeLines(capture.output(sessionInfo()), paste0("NN_", num_neighbors, "_sigma_", sigma_val, "_sessionInfo.txt"))

