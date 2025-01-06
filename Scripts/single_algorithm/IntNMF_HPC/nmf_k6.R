
# HPC R script to run IntNMF parallely for k=6 (single k-range).

set.seed(123)

suppressPackageStartupMessages({
  library(parallel)
  library(doParallel)
  library(foreach)
  library(MASS)
  library(mclust)
  library(IntNMF)
  source('parallel_IntNMF.R')
})

dat <- readRDS('IntNMF_input.rds')

res <- nmf.opt.k_parallel(
  dat       = dat,
  n.runs    = 50,
  n.fold    = 10,
  k.range   = c(6, 6),
  maxiter   = 1000,
  st.count  = 10,
  num_cores = 50,
  init_boolean = TRUE,
  n.init = 30,
  result    = TRUE,
  make.plot = FALSE,
  progress  = TRUE
)

saveRDS(res, file = 'nmf_k6_results.rds')
writeLines(capture.output(sessionInfo()), 'nmf_k6_sessionInfo.txt')
cat('IntNMF parallely finished for k=6.\n')
