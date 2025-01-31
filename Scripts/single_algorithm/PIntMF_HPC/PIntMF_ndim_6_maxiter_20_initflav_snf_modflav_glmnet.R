
suppressPackageStartupMessages(library(PintMF))

# Load input data
input = readRDS("PIntNMF_input.rds")

RNGversion("4.2.2")
set.seed(123)

# Hyperparameter setup ###
p = 6
max.it = 20
init_flavor = snf
flavor_mod = glmnet

# Run PIntNMF
pintmf = SolveInt(
  Y=input, 
  p=p, 
  max.it=max.it, 
  verbose=FALSE, 
  init_flavor=init_flavor, 
  flavor_mod=flavor_mod
)

# Save the results
saveRDS(pintmf, "PIntMF_ndim_6_maxiter_20_initflav_snf_modflav_glmnet_results.rds")

# Also record sessionInfo
writeLines(capture.output(sessionInfo()), 
           "PIntMF_ndim_6_maxiter_20_initflav_snf_modflav_glmnet_sessionInfo.txt")

