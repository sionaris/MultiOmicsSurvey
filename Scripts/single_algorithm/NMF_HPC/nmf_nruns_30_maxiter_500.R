#
# Scripts/single_algorithm/NMF_HPC/nmf_nruns_30_maxiter_500.R
library(foreach)
library(doParallel)
library(mclust)
library(MASS)

source('nmf_integrative_functions.R') # Ensure this file has nmf.opt.k.integrative defined

# Load data
dat <- readRDS('NMF_input.rds')

res <- nmf.opt.k.integrative(
  dat = dat,
  is.binary = c(TRUE, FALSE, FALSE, FALSE, FALSE),
  n.runs = 30,
  n.fold = 10,
  k.range = 2:10,
  maxiter = 500,
  allowParallel = TRUE,
  n.cores = 10,
  make.plot = FALSE,
  result = TRUE
)

saveRDS(res, file = 'nmf_nruns_30_maxiter_500_results.rds')
