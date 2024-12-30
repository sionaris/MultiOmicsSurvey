library(Spectrum)

# Load modified Spectrum functions
# (updates: binary distances and parallelization)
source("Spectrum_functions.R")

# Import input
readRDS("Spectrum_input.rds")

# Setup
parameters = list(
  method = 2, # multimodality gap method (Gaussian/non-Gaussian clusters)
  diffusion = TRUE, # whether to perform graph diffusion
  kerneltype = 'density', # adaptive density-aware kernel from the paper
  maxk = 10, # maximum number of expected clusters
  # Higher values for NN and NN2 will prefer global structures
  NN = 3, # the number of nearest neighbors to use sigma parameters (default=3)
  NN2 = 7, # the number of nearest neighbors to use for the common nearest neighbors
  frac = 2, # optk search param, fraction to find the last substantial drop (multimodality gap method param)
  thresh = 7, # optk search param, how many points ahead to keep searching (multimodality gap method param)
  tunekernel = TRUE, # whether to tune the kernel, only applies for method 2 (default=FALSE)
  clusteralg = 'km', # or GMM
  diffusion_iters = 5, # default 4
  KNNs_p = 10, # number of KNNs when making KNN graph (default=10, suggested=10-20)
  fontsize = 18, # plot parameter
  dotsize = 1.25 # plot parameter
)

Spectrum_results = Spectrum_bin_and_par(input,
                                           method = parameters$method,
                                           diffusion = parameters$diffusion,
                                           kerneltype = parameters$kerneltype,
                                           maxk = parameters$maxk,
                                           NN = parameters$NN,
                                           NN2 = parameters$NN2,
                                           frac = parameters$frac,
                                           thresh = parameters$thresh,
                                           tunekernel = parameters$tunekernel,
                                           clusteralg = parameters$clusteralg,
                                           diffusion_iters = parameters$diffusion_iters,
                                           KNNs_p = parameters$KNNs_p,
                                           fontsize = parameters$fontsize,
                                           dotsize = parameters$dotsize,
                                           distances = c("binary",
                                                         rep("euclidean", 4)),
                                           cores = 5)

# Export results object
saveRDS(Spectrum_results, "Spectrum_results.rds")
