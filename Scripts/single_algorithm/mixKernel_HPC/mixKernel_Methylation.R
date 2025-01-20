suppressPackageStartupMessages({
  library(mixOmics)
  library(mixKernel)
})

# Load input
input = readRDS("mixKernel_input.rds")[["Methylation"]]

# Set seed if necessary
RNGversion("4.2.2")
set.seed(123)

Methylation_kernel = compute.kernel(input,
                              kernel.func = "abundance",
                              method = "euclidean",
                              test.pos.semidef = TRUE)

saveRDS(Methylation_kernel, "Methylation_kernel.rds")
writeLines(capture.output(sessionInfo()), "Methylation_kernel_sessionInfo.txt")