#
# Scripts/single_algorithm/NMF_HPC/nmf_nruns_75_maxiter_1000.R
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
  n.runs = 75,
  n.fold = 10,
  k.range = 2:10,
  maxiter = 1000,
  allowParallel = TRUE,
  n.cores = 10,
  make.plot = FALSE,
  result = TRUE
)

saveRDS(res, file = 'nmf_nruns_75_maxiter_1000_results.rds')
