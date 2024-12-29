library(foreach)
library(doParallel)
library(parallel)
library(mclust)
library(MASS)
library(IntNMF)
library(doRNG)

dat <- readRDS("IntNMF_input.rds")
source("parallel_IntNMF.R")

res <- nmf.opt.k_parallel(dat = dat,
                 n.runs = 5,
                 n.fold = 5,
                 k.range = 9:10,
                 result = TRUE,
                 make.plot = TRUE,
                 progress = FALSE,
                 maxiter = 50,
                 num_cores = 5)

saveRDS(res, file = '", job_name, "_results.rds')