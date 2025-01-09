# Load LRAcluster code
source(paste0(getwd(), "/full_LRAcluster_source.R"))

# Read in the data input
input <- readRDS(paste0(getwd(), "/LRAcluster_input.rds"))

# Set the types vector accordingly: first one binary, the rest Gaussian
types <- c("binary", "gaussian", "gaussian", "gaussian", "gaussian")

# Dimensions to try
dimensions_to_try <- c(2:10)

# Run LRAcluster for each dimension and store the potential value
library(doParallel)
library(foreach)

# Number of cores
num_cores <- 9
cl <- makeCluster(num_cores)
registerDoParallel(cl)

timestamp()
# Run LRAcluster in parallel using foreach
results_list <- foreach(dim_val = dimensions_to_try) %dopar% {
  LRAcluster(data = input, 
             types = types, 
             dimension = dim_val, 
             names = names(input))
}

stopCluster(cl)
timestamp()

# Name the results by the dimension for clarity
names(results_list) <- paste0("dim_", dimensions_to_try)
saveRDS(results_list, paste0(getwd(), "/LRAcluster_results_list.rds"))
