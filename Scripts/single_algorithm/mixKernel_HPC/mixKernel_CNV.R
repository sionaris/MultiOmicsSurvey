suppressPackageStartupMessages({
  library(mixOmics)
  library(mixKernel)
})

# Load input
input = readRDS("mixKernel_input.rds")[["CNV"]]

# Set seed if necessary
RNGversion("4.2.2")
set.seed(123)

CNV_kernel = compute.kernel(input,
                                    kernel.func = "abundance",
                                    method = "euclidean",
                                    test.pos.semidef = FALSE)

saveRDS(CNV_kernel, "CNV_kernel.rds")
writeLines(capture.output(sessionInfo()), "CNV_kernel_sessionInfo.txt")