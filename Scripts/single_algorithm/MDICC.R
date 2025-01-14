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
algorithm = "MDICC"
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
optk_boolean = "TRUE" # either TRUE or FALSE. Answers whether the algorithm suggests an optimal k
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

# However MDICC prefers features in columns so we transpose the matrices.

# Extract the names of the modalities that will be used
modalities = unlist(strsplit(data_types, "-"))

# Replace with TRUE wherever features are in rows
features_in_rows = rep(TRUE, length(modalities))

# Run algorithm #####
# In order to use MDICC, you must follow guidance from:
# https://github.com/yushanqiu/MDICC

# Download MDICC-main.zip from GitHub and place into a directory of preference
# I choose .libPaths()[1]

# MDICC dir
MDICC_dir = paste0(.libPaths()[1], '/MDICC')

# To avoid known overwhelming warnings from MDICClabel() you can set some environment
# variables as below:

# Avoid “Could not find the number of physical cores”
# scikit-learn’s joblib backend (loky) tries to detect the number of physical cores.
# On some Windows systems, it cannot find certain OS files, so it defaults to 
# logical cores and warns the user.
Sys.setenv(LOKY_MAX_CPU_COUNT = 8)  # or any integer >= 1

# KMeans + MKL Memory Leak Warning
# KMeans on Windows with MKL can leak memory if it parallelizes using more 
# threads than chunks of data. scikit-learn suggests either limiting the number 
# of threads or chunks so you avoid the leak.
Sys.setenv(OMP_NUM_THREADS = 3)

# Load MDICC source code
library(Rcpp)
library(parallel)
library(Matrix)
source(paste0(MDICC_dir, '/NetworkFusion.R'))

# Dynamically load the .dll file from the same directory
dyn.load(paste0(MDICC_dir, '/projsplx_R.dll'))

# You also need Anaconda and reticulate for MDICC
library(reticulate)
# use_virtualenv("base")
Sys.setenv(RETICULATE_PYTHON="C:/ProgramData/anaconda3/python.exe")
use_python("C:/ProgramData/anaconda3/python.exe")
py_config()
py_available()

# > py_config()
# python:         C:/ProgramData/anaconda3/python.exe
# libpython:      C:/ProgramData/anaconda3/python312.dll
# pythonhome:     C:/ProgramData/anaconda3
# version:        3.12.4 | packaged by Anaconda, Inc. | (main, Jun 18 2024, 15:03:56) [MSC v.1929 64 bit (AMD64)]
# Architecture:   64bit
# numpy:          C:/ProgramData/anaconda3/Lib/site-packages/numpy
# numpy_version:  1.26.4
# 
# NOTE: Python version was forced by RETICULATE_PYTHON
# > py_available()
# [1] TRUE

# Source MDICC Python files
source_python(paste0(MDICC_dir, "/LocalAffinityMatrix.py"))
source_python(paste0(MDICC_dir, "/score.py"))
source_python(paste0(MDICC_dir, "/label.py"))

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

# Distance calculations ###
# Modality types
continuous = c("RNAseq", "CNV", "Methylation", "miRNA")
categorical = c("SNPs")

# Calculate the pair-wise distance (Euclidean for continuous modalities)
input_dists = lapply(input[continuous], function(x) {
  x = as.matrix(x)
  x = SNFtool::dist2(x, x)
})

# Binary for SNPs (see ?dist for details)
input_dists[["SNPs"]] = as.matrix(dist(as.matrix(input$SNPs),
                                       as.matrix(input$SNPs),
                                       method = "binary"))
gc()

# Hyperparameter setup
# Some hyperparameters are fixed, based on authors' recommendations here:
# https://academic.oup.com/bib/article/23/3/bbac132/6569541
aff_matrix_neighbors = 18 # author suggestion in the paper
k2 = 41:44 # suggested by the authors
k3 = 2:10 # possible numbers of clusters
cc = 2:10 # subspace dimensions


# Construct affinity matrices
aff_input = list()
for (modality in modalities) {
  aff_input[[modality]] <- testaff(as.matrix(input_dists[[modality]]),
                                   aff_matrix_neighbors)
}
gc()

# # Set up parallel background
# library(doParallel)
# library(foreach)
# 
# nCores = 9 # set depending on your machine's capabilities. Here: one core for each k
# cl <- makeCluster(nCores)
# registerDoParallel(cl)
# 
# # We want the following functions to be accessible on each worker:
# func_vec <- c("dominate.set", "transition.fields", "dn", 
#               "eig1", "L2_distance_1", "umkl", "Hbeta", "MDICC",
#               "MDICClabel", "MDICCscore")
# 
# # MDICC runs
# timestamp()
# mdicc_results_parallel <- foreach(
#   cc = 2:10, 
#   .combine = 'rbind', 
#   .export = func_vec,  
#   .packages = c("Matrix", "reticulate")  
# ) %:%
#   foreach(
#     k2_val = k2, 
#     .combine = 'rbind',
#     .export = func_vec,
#     .packages = c("Matrix", "reticulate")
#   ) %:%
#   foreach(
#     k3_val = 2:10,
#     .combine = 'rbind',
#     .export = func_vec,
#     .packages = c("Matrix", "reticulate")
#   ) %dopar% {
#     # Running MDICC
#     S <- MDICC(aff_input, c = cc, k = k2_val)
#     label_vec <- MDICClabel(S, k3_val)
#     
#     data.frame(
#       c              = cc,
#       k2             = k2_val,
#       k3             = k3_val,
#       cluster_labels = I(list(label_vec)),
#       stringsAsFactors = FALSE
#     )
#   }
# 
# stopCluster(cl)
# timestamp()

# ---- Serial version (no parallel, nested for loops) ----

mdicc_results_serial <- data.frame(
  c              = numeric(),
  k2             = numeric(),
  k3             = numeric(),
  cluster_labels = I(list()),
  stringsAsFactors = FALSE
)

S_matrices = list()

timestamp()
for (cc in 2:10) {
  cat("---------------", "\n")
  cat("Starting cc =", cc, "\n")
  for (k2_val in k2) {
    cat("---------------", "\n")
    cat("Starting k2 =", k2_val, "\n")
    for (k3_val in 2:10) {
      S <- suppressMessages(MDICC(aff_input, c = cc, k = k2_val))
      S <- as.matrix(S)
      label_vec <- MDICClabel(S, k3_val)

      # Combine current result into a data.frame row
      new_row <- data.frame(
        c              = cc,
        k2             = k2_val,
        k3             = k3_val,
        cluster_labels = I(list(label_vec)),
        stringsAsFactors = FALSE
      )

      # Append to the main result object
      mdicc_results_serial <- rbind(mdicc_results_serial, new_row)
      S_matrices[[paste0("cc = ", cc, ", k2 = ", k2_val, ", k3 = ", k3_val)]] = S
      cat("Done for k3 =", k3_val, "\n")
    }
    cat("---------------", "\n")
    cat("Done for k2 =", k2_val, "\n")
  }
  cat("---------------", "\n")
  cat("Done for cc =", cc, "\n")
}
timestamp()

# Save environment
save.image(paste0(home, "/Results/single_algorithm/", 
                  algorithm, "/", algorithm, "_", data_source, "_",
                  data_types, "_eval_on_", evaluation_source,
                  "_env.RData"))

