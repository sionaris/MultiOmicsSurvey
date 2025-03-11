
library(FactoMineR)

# Set seed if necessary
RNGversion("4.2.2")
set.seed(123)

# Load input
mfa_input = readRDS("MFA_input.rds")

# Define number of components
ncp = 25

# Run for ncp = 25
snp_cols = sum(grep("SNPs_", colnames(mfa_input)))
rna_cols = sum(grep("RNAseq_", colnames(mfa_input)))
cnv_cols = sum(grep("CNV_", colnames(mfa_input)))
mirna_cols = sum(grep("miRNA_", colnames(mfa_input)))
methyl_cols = sum(grep("Methylation_", colnames(mfa_input)))

t1 = Sys.time()
mfa = MFA(mfa_input,
          group = c(snp_cols, rna_cols, cnv_cols, mirna_cols, methyl_cols), # dimensionalities of the datasets
          type = c("n", rep("c", 4)),
          excl = NULL,
          ncp = ncp,
          name.group = c("SNPs", "RNAseq", "CNV", "miRNA", "Methylation"),
          graph = FALSE,
          axes = c(1,2))
dt = Sys.time() - t1

# Export
saveRDS(mfa, paste0("MFA_ncp_", ncp, "_results.rds"))
writeLines(capture.output(sessionInfo()), paste0("MFA_ncp_", ncp, "_sessionInfo.txt"))

