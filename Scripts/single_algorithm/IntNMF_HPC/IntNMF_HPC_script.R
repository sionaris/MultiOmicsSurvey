library(foreach)
library(doParallel)
library(mclust)
library(MASS)
library(IntNMF)

dat <- readRDS("NMF_input.rds")

res <- nmf.opt.k(dat = dat,
                 n.runs = 300,
                 n.fold = 10,
                 k.range = 2:10,
                 result = TRUE,
                 make.plot = TRUE,
                 progress = FALSE,
                 maxiter = 2000)

saveRDS(res, file = '", job_name, "_results.rds')