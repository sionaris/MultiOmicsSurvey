# Import data from gitignored "Resources/BRCA complete/" folder #####

# This RDS object was produced using the Scripts/MOVICS/MOVICS_baseline.R script
input = readRDS("Resources/TCGA/mm_input.rds")

# Setup environment variables for markdown #####

# Ensure reproducibility
RNGversion("4.2.2")
set.seed(123)

# Load custom helper functions
source("Scripts/automated_scripts/custom_functions.R")
source("Scripts/automated_scripts/modified_MOVICS_functions.R")

# Preamble
home = getwd()
algorithm = "MFA"
alg_feature_pref = "cols" # Where does the algorithm expect the features to be
citation = fetch_citation(algorithm = algorithm)
data_source = "TCGA" # e.g. TCGA, TCGA-transNEO, transNEO-PARTNER
data_types = "RNAseq-CNV-Methylation-miRNA-SNPs" # e.g. RNAseq, RNAseq-CNV-miRNA
evaluation_source = "transNEO" # e.g. PARTNER, transNEO-PARTNER 
title = paste0("Results from ", algorithm)
subtitle = paste0("<b>Train</b>: ", data_source, " ", data_types, 
                  " | <b>Evaluation</b>: ", evaluation_source)
in_a_nutshell = fetch_in_a_nutshell(algorithm = algorithm)
ground_truth_labels = openxlsx::read.xlsx("Results/MOVICS_baseline/MOVICS_TCGA_RNAseq-CNV-Methylation-miRNA-SNPs_eval_on_transNEO_clusterings.xlsx")
ground_truth_k = 2 # optk from MOVICS
optk_boolean = "FALSE" # either TRUE or FALSE. Answers whether the algorithm suggests an optimal k
optk_text = ifelse(optk_boolean == TRUE,
                   "<u>suggests</u> an estimate of the optimal number of multi-omic clusters $k$",
                   "<u>does not suggest</u> an optimal number of multi-omic clusters $k$")

# Detailed description of the algorithm
description = paste(readLines(paste0("Resources/algorithm_descriptions/", algorithm,
                                     "_description.Rmd")),
                    collapse = "\n") # File path to .Rmd file within Resources/algorithm_descriptions

# Create algorithm directory if it doesn't exist
if (!dir.exists(paste0(home, "/Results/single_algorithm/", 
                       algorithm))) {
  dir.create(paste0(home, "/Results/single_algorithm/", 
                    algorithm))
}

# Preprocessing flags and code #####
library(stringr)
library(dplyr)

# All preprocessing for this input has already been performed using the 
# Scripts/MOVICS/MOVICS_baseline.R script

# Extract the names of the modalities that will be used
modalities = unlist(strsplit(data_types, "-"))

# Replace with TRUE wherever features are in rows
features_in_rows = rep(TRUE, length(modalities))

# Run algorithm #####
library(FactoMineR)

# Check what kind of arrangement the algorithm requires as input 
# (i.e. features in rows or columns)?
if (alg_feature_pref == "rows") {
  rogue_indices = which(features_in_rows == FALSE)
  if (length(rogue_indices >= 1)) {
    for (index in rogue_indices) {
      cols = colnames(input[[index]])
      rows = rownames(input[[index]])
      input[[index]] = t(input[[index]])
      rownames(input[[index]]) = cols # transpose names
      colnames(input[[index]]) = rows # transpose names
      rm(rows, cols)
    }
  }
} else if (alg_feature_pref == "cols") {
  rogue_indices = which(features_in_rows == TRUE)
  if (length(rogue_indices >= 1)) {
    for (index in rogue_indices) {
      cols = colnames(input[[index]])
      rows = rownames(input[[index]])
      input[[index]] = t(input[[index]])
      rownames(input[[index]]) = cols # transpose names
      colnames(input[[index]]) = rows # transpose names
      rm(rows, cols)
    }
  }
}
rm(rogue_indices, index); gc()

# Import clinical data for the TCGA samples of interest
clinical_data = openxlsx::read.xlsx("Resources/TCGA/clinical_data.xlsx")

# Keep the top 33% features for each continuous dataset based on MAD
mfa_input = input

# for continuous datasets
mfa_input[c("RNAseq", "CNV", "Methylation", "miRNA")] = lapply(mfa_input[c("RNAseq", "CNV", "Methylation", "miRNA")], function(x) {
  colMAD <- apply(x, 2, mad, na.rm = TRUE)
  ordered_cols <- order(colMAD, decreasing = TRUE)
  n_top <- ceiling(0.33 * ncol(x))
  x[, ordered_cols[1:n_top]]
})

# Keep cancer drivers and the top 33%  most mutated genes of the rest for SNPs
COSMIC_BC_drivers = read.csv("Resources/COSMIC_CGC_Breast_somatic.csv")$Gene.Symbol %>%
  as.character()

snp_mat <- input[["SNPs"]]
mutation_counts <- colSums(snp_mat, na.rm = TRUE)
non_drivers <- setdiff(colnames(snp_mat), COSMIC_BC_drivers)
ordered_non_drivers <- non_drivers[order(mutation_counts[non_drivers], decreasing = TRUE)]
n_top <- ceiling(0.33 * length(non_drivers))
top_non_drivers <- ordered_non_drivers[1:n_top]
genes_to_keep <- unique(c(COSMIC_BC_drivers, top_non_drivers))
mfa_input[["SNPs"]] <- snp_mat[, intersect(colnames(snp_mat), genes_to_keep)]

# Export input for HPC
for (i in 1:length(mfa_input)) {
  colnames(mfa_input[[i]]) = paste0(names(mfa_input)[i], "_", colnames(mfa_input[[i]]))
}
snps_no = as.numeric(length(intersect(colnames(snp_mat), genes_to_keep)))

mfa_input = do.call(cbind, mfa_input)
snp_cols = sum(grepl("SNPs_", colnames(mfa_input)))
rna_cols = sum(grepl("RNAseq_", colnames(mfa_input)))
cnv_cols = sum(grepl("CNV_", colnames(mfa_input)))
mirna_cols = sum(grepl("miRNA_", colnames(mfa_input)))
methyl_cols = sum(grepl("Methylation_", colnames(mfa_input)))

mfa_input = as.data.frame(mfa_input)
mfa_input[, 1:snps_no] <- lapply(mfa_input[, 1:snps_no], as.factor)

saveRDS(mfa_input, "Resources/MFA_input.rds")
rm(i, snps_no, snp_mat, mutation_counts, non_drivers, ordered_non_drivers,
   n_top, top_non_drivers, genes_to_keep); gc()

# Setup ###
# Hyperparameter tuning - example with ncp = 10

# t1 = Sys.time()
# mfa = MFA(mfa_input,
#           group = c(as.numeric(snp_cols), 
#                     as.numeric(rna_cols), 
#                     as.numeric(cnv_cols), 
#                     as.numeric(mirna_cols), 
#                     as.numeric(methyl_cols)), # dimensionalities of the datasets
#           type = c("n", rep("c", 4)),
#           excl = NULL,
#           ncp = 10,
#           name.group = c("SNPs", "RNAseq", "CNV", "miRNA", "Methylation"),
#           graph = FALSE,
#           axes = c(1,2))
# 
# dt = Sys.time() - t1
# print(dt)