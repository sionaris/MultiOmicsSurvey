suppressPackageStartupMessages({
  library(coca)
  library(klic)
  library(foreach)
  library(doParallel)
})

# Import input
input = readRDS("KLIC_input.rds")

# Hyperparameters and tuning ###
maxK = 10
B = 100

# Build consensus matrices with COCA
allCM = vector("list", maxK - 1)

# Distance types
distances = c("binary", rep("euclidean", 4))

# Set up parallel loop
cl = makeCluster(9) # 9 cores
registerDoParallel(cl)

nSamples <- nrow(input[[1]])
nDatasets <- length(input)

allCM <- vector("list", maxK - 1)  # Will store results for K=2..maxK

res <- foreach(k = 2:maxK, 
               .packages = c("coca", "klic")) %dopar% {
                 # 1) Initialize CM for this k
                 CM <- array(NA, dim = c(nSamples, nSamples, nDatasets))
                 
                 # 2) Loop over datasets
                 for (i in seq_len(nDatasets)) {
                   temp_mat <- coca::consensusCluster(
                     data         = input[[i]], 
                     dist         = distances[i],
                     clMethod     = "hclust",
                     B            = B,
                     K            = k,
                     hclustMethod = "average"
                   )
                   
                   # 3) Shift eigenvalues to ensure PSD
                   CM[, , i] <- klic::spectrumShift(temp_mat, verbose = FALSE)
                 }
                 
                 # Return the CM array for this 'k'
                 CM
               }

# Export objects
saveRDS(res, "KLIC_CM_res.rds")
writeLines(capture.output(sessionInfo()), "KLIC_CM_sessionInfo.txt")