
# Scripts/single_algorithm/NMF_HPC/nmf_nruns_75_maxiter_1000_lr_0.001.R
library(foreach)
library(doParallel)
library(mclust)
library(MASS)

source('nmf_integrative_functions.R')

dat <- readRDS('NMF_input.rds')

res <- nmf.opt.k.integrative(
  dat = dat,
  is.binary = c(TRUE, FALSE, FALSE, FALSE, FALSE),
  n.runs = 75,
  n.fold = 10,
  k.range = 2:10,
  maxiter = 1000,
  lr = 0.001,
  allowParallel = TRUE,
  n.cores = 10,
  make.plot = FALSE,
  result = TRUE
)

saveRDS(res, file = 'nmf_nruns_75_maxiter_1000_lr_0.001_results.rds')
