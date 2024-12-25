
# Scripts/single_algorithm/NMF_HPC/nmf_nruns_200_maxiter_200_lr_1e-04.R
library(foreach)
library(doParallel)
library(mclust)
library(MASS)

source('nmf_integrative_functions.R')

dat <- readRDS('NMF_input.rds')

res <- nmf.opt.k.integrative(
  dat = dat,
  is.binary = c(TRUE, FALSE, FALSE, FALSE, FALSE),
  n.runs = 200,
  n.fold = 10,
  k.range = 2:10,
  maxiter = 200,
  lr = 1e-04,
  allowParallel = TRUE,
  n.cores = 10,
  make.plot = FALSE,
  result = TRUE
)

saveRDS(res, file = 'nmf_nruns_200_maxiter_200_lr_1e-04_results.rds')
