# This script is used to determine the optimal choice of pseudocount for transNEO
# TPM data, using the approach that we applied on the TCGA data

# Import transNEO TPM, clinical and raw count data from the Resources/transNEO folder
library(dplyr)
library(DESeq2)
transNEO_mm_inputs = readRDS("Resources/transNEO/transNEO_multimodal_inputs.rds")
tpm_transNEO = readRDS("Resources/transNEO/TPM.rds")
raw_counts_transNEO = readRDS("Resources/transNEO/raw_counts.rds")

# Convert raw counts to matrix
rnames = raw_counts_transNEO$Hugo
count_matrix = as.matrix(raw_counts_transNEO[, 1:153])
rownames(count_matrix) = rnames

# https://www.biorxiv.org/content/10.1101/404962v1.full -> Process description

# Get clinical data
clinical_data = transNEO_mm_inputs$`Full pheno` %>%
  dplyr::rename(Sample.ID = Donor.ID)

# We need a grouping for the breast cancer samples. We are going to go for 
# pCR vs RD: response to neoadjuvant treatment in transNEO
pCR_samples = clinical_data$Sample.ID[clinical_data$pCR.RD == "pCR"]
RD_samples = clinical_data$Sample.ID[clinical_data$pCR.RD == "RD"]

# Get size factors similarly to TCGA
# Mock sample conditions
sample_conditions <- data.frame(
  row.names = colnames(count_matrix),
  condition = factor(rep("condition", ncol(count_matrix)))
)

# Create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = count_matrix, 
                              colData = sample_conditions, 
                              design = ~ 1) 

# Estimate size factors
dds <- estimateSizeFactors(dds)

# Save size factors for later if needed
size_factors_transNEO = dds$sizeFactor
names(size_factors_transNEO) = dds@colData@rownames

# Get the normalized counts
normalized_counts <- counts(dds, normalized = TRUE)

# Assign proper row and column names to the normalized counts
rownames(normalized_counts) = rownames(count_matrix)
colnames(normalized_counts) = colnames(count_matrix)

# Estimate averages
avg_pCR_size_factor = mean(size_factors_transNEO[pCR_samples], na.rm = TRUE)
avg_RD_size_factor = mean(size_factors_transNEO[RD_samples], na.rm = TRUE)
sf = c(avg_pCR_size_factor, avg_RD_size_factor)

# Lun et al. suggest pseudocount = max{1, r|1/smin - 1/smax|}, where r = 1 (suggestion)
pseudocount_RNA = max(c(1, abs(1/min(sf) - 1/max(sf)))) # 1

# We proceed with the log2(norm.counts + 1) transformation
log2TPMplus1_transNEO = log2(tpm_transNEO + pseudocount_RNA)
saveRDS(log2TPMplus1_transNEO, "Resources/transNEO/log2TPMplus1_transNEO.rds")
