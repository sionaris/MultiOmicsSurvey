suppressPackageStartupMessages({
  library(mixOmics)
  library(mixKernel)
})

# Load input
input = readRDS("mixKernel_input.rds")[["SNPs"]]

# Set seed if necessary
RNGversion("4.2.2")
set.seed(123)

SNPs_kernel = compute.kernel(input,
                             kernel.func = "abundance",
                             method = "jaccard",
                             test.pos.semidef = TRUE)

saveRDS(SNPs_kernel, "SNPs_kernel.rds")
writeLines(capture.output(sessionInfo()), "SNPs_kernel_sessionInfo.txt")