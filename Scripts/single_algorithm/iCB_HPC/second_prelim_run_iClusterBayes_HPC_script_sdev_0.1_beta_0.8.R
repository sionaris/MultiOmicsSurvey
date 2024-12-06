
# Load necessary packages
library(iClusterPlus)
library(parallel)

# Read in the data input
input <- readRDS(paste0(getwd(), "/iCB_input.rds"))

# Define hyperparameters
sdev <- 0.1
beta_var_scale <- 0.8

# Generate an informative suffix for output files
suffix <- paste0("sdev_", sdev, "_beta_", beta_var_scale)

# Define the range of K
K_values <- 1:9

# Set up the number of cores for parallelization
num_cores <- 9  # Adjust based on available cores

# Define data types for each dataset
data_types <- c("binomial", "gaussian", "gaussian", "gaussian", "gaussian")

# Fix other hyperparameters
n_burnin <- 1000
n_draw <- 1200
prior_gamma <- rep(0.1, length(data_types))
thin <- 1
pp_cutoff <- 0.5

# Run tune.iClusterBayes in parallel
tune_results <- tune.iClusterBayes(
  cpus = num_cores,
  dt1 = input$SNPs,
  dt2 = input$RNAseq,
  dt3 = input$CNV,
  dt4 = input$miRNA,
  dt5 = input$Methylation,
  type = data_types,
  K = K_values,
  n.burnin = n_burnin,
  n.draw = n_draw,
  prior.gamma = prior_gamma,
  sdev = sdev,
  beta.var.scale = beta_var_scale,
  thin = thin,
  pp.cutoff = pp_cutoff
)

# Save the results with an informative filename
saveRDS(tune_results, paste0(getwd(), "/tune_results_", suffix, ".rds"))

