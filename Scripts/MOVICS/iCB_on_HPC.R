.libPaths(c("/home/as3582/rds/hpc-work/Rlibs/win-library/4.4",
            "/home/as3582/rds/hpc-work/Rlibs/R-4.4.0/library"))

library(MOVICS)

load("/home/as3582/rds/hpc-work/MultiOmicsSurvey/movics_env.RData")

iClusterBayes.res = getMOIC(data        = input,
                            N.clust     = 3,
                            methodslist = "iClusterBayes",
                            type        = c("gaussian",
                                            "gaussian",
                                            "gaussian",
                                            "gaussian",
                                            "binomial"),
                            n.burnin    = 18000,
                            n.draw      = 12000,
                            prior.gamma = c(0.5, 0.5, 0.5, 0.5, 0.5),
                            sdev        = 0.05,
                            thin        = 3)

saveRDS(iClusterBayes.res, "/home/as3582/rds/hpc-work/MultiOmicsSurvey/iCB_full.rds")