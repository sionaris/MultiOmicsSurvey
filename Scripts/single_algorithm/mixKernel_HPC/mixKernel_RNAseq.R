suppressPackageStartupMessages({
  library(mixOmics)
  library(mixKernel)
})

# Load input
input = readRDS("mixKernel_input.rds")[["RNAseq"]]

# Set seed if necessary
RNGversion("4.2.2")
set.seed(123)

RNAseq_kernel = compute.kernel(input,
                            kernel.func = "abundance",
                            method = "euclidean",
                            test.pos.semidef = TRUE)

saveRDS(RNAseq_kernel, "RNAseq_kernel.rds")
writeLines(capture.output(sessionInfo()), "RNAseq_kernel_sessionInfo.txt")