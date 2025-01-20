suppressPackageStartupMessages({
  library(mixOmics)
  library(mixKernel)
})

# Load input
input = readRDS("mixKernel_input.rds")[["miRNA"]]

# Set seed if necessary
RNGversion("4.2.2")
set.seed(123)

miRNA_kernel = compute.kernel(input,
                             kernel.func = "abundance",
                             method = "euclidean",
                             test.pos.semidef = TRUE)

saveRDS(miRNA_kernel, "miRNA_kernel.rds")
writeLines(capture.output(sessionInfo()), "miRNA_kernel_sessionInfo.txt")