
library(FactoMineR)

# Set seed if necessary
RNGversion("4.2.2")
set.seed(123)

# Load input
mfa_input = readRDS("MFA_input.rds")

# Define number of components
ncp = 20

# Run for ncp = 20
snp_cols = sum(grepl("SNPs_", colnames(mfa_input)))
rna_cols = sum(grepl("RNAseq_", colnames(mfa_input)))
cnv_cols = sum(grepl("CNV_", colnames(mfa_input)))
mirna_cols = sum(grepl("miRNA_", colnames(mfa_input)))
methyl_cols = sum(grepl("Methylation_", colnames(mfa_input)))

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
results = list(MFA = mfa, dt = dt)
rm(mfa); gc()

# Export
saveRDS(results, paste0("MFA_ncp_", ncp, "_results.rds"))
writeLines(capture.output(sessionInfo()), paste0("MFA_ncp_", ncp, "_sessionInfo.txt"))

