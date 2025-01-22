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
data_source = "TCGA" # e.g. TCGA, TCGA-transNEO, transNEO-PARTNER
data_types = "RNAseq-CNV-Methylation-miRNA-SNPs" # e.g. RNAseq, RNAseq-CNV-miRNA
evaluation_source = "transNEO" # e.g. PARTNER, transNEO-PARTNER 
title = paste0("Results from comparisons between algorithms")
subtitle = paste0("<b>Train</b>: ", data_source, " ", data_types, 
                  " | <b>Evaluation</b>: ", evaluation_source)
ground_truth_labels = openxlsx::read.xlsx("Results/MOVICS_baseline/MOVICS_TCGA_RNAseq-CNV-Methylation-miRNA-SNPs_eval_on_transNEO_clusterings.xlsx")
ground_truth_k = 2 # optk from MOVICS

# Detailed description of the algorithm
description = paste(readLines(paste0("Resources/Comparisons_description.Rmd")),
                    collapse = "\n") # File path to .Rmd file within Resources/algorithm_descriptions

# Create algorithm directory if it doesn't exist
if (!dir.exists(paste0(home, "/Results/Comparisons/", 
                       algorithm))) {
  dir.create(paste0(home, "/Results/Comparisons/", 
                    algorithm))
}

# Imports #####
library(dplyr)
library(openxlsx)


# Create list of clustering outputs
clusterings = list()
for (algorithm in algorithms) {
  dir = paste0("Results/single_algorithm/", algorithm, "/")
  files = list.files(dir)
  clusterfile = files[grepl("clusterings.xlsx", files)]
  clusterings[[algorithm]] = read.xlsx(paste0(dir, clusterfile))
}

rm(algorith, dir, files, clusterfile); gc()

