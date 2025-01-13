
suppressPackageStartupMessages(library(mixOmics))

# Load input data
input = readRDS("RGCCA_input.rds")

RNGversion("4.2.2")
set.seed(123)

# Hyperparameter setup ###
penalty = 0.9
design = 1 - diag(length(input))
ncomp = 10
scheme = "horst"
keepX = NULL
scale = FALSE
max.iter = 1000 # default
init = "svd.single"
near.zero.var = FALSE

# Run RGCCA
rgcca = wrapper.rgcca(
  input,
  penalty = penalty,
  design = design,
  ncomp = ncomp,
  keepX = keepX,
  scheme = scheme,
  scale = scale,
  init = init,
  tol = .Machine$double.eps, # default
  max.iter = max.iter,
  near.zero.var = near.zero.var,
  all.outputs = TRUE
)

# Save the results
saveRDS(rgcca, "rgcca_results.rds")

# Also record sessionInfo
writeLines(capture.output(sessionInfo()), 
           "rgcca_sessionInfo.txt")

