suppressPackageStartupMessages({
  library(mixOmics)
  library(mixKernel)
  library(foreach)
  library(doParallel)
})

# Load input
input = readRDS("mixKernel_input.rds")
kernel_matrices = list()

# Set seed if necessary
RNGversion("4.2.2")
set.seed(123)

# Initialize a cluster of 5 cores
cl <- makeCluster(5)
registerDoParallel(cl)

timestamp()

kernel_matrices <- foreach(j = seq_along(input),
                           .packages = c("mixKernel", "mixOmics")) %dopar% {
                             # Determine the method
                             if (names(input)[j] == "SNPs") {
                               method <- "jaccard"
                             } else {
                               method <- "euclidean"
                             }
                             
                             # Compute kernel using mixKernel::compute.kernel
                             compute.kernel(
                               input[[j]],
                               kernel.func = "abundance",
                               method = method,
                               test.pos.semidef = FALSE
                             )
                           }

timestamp()

# Stop the cluster
stopCluster(cl)

saveRDS(kernel_matrices, "kernel_matrices.rds")
writeLines(capture.output(sessionInfo()), "kernel_matrices_sessionInfo.txt")