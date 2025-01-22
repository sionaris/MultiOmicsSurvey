library(FactoMineR)

# Set seed if necessary
RNGversion("4.2.2")
set.seed(123)

# Load input
mfa_input = readRDS("MFA_input.rds")

# Run for ncp = 10
mfa = MFA(mfa_input,
          group = c(12926, 56717, 35858, 1568, 36761), # dimensionalities of the datasets
          type = c("n", rep("c", 4)),
          excl = NULL,
          ncp = 10,
          name.group = c("SNPs", "RNAseq", "CNV", "miRNA", "Methylation"),
          graph = FALSE,
          axes = c(1,2))

# Export
saveRDS(mfa, "MFA_results.rds")
writeLines(capture.output(sessionInfo()), "MFA_sessionInfo.txt")

