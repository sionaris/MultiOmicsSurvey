# This script generates .sh and .R scripts for each combination of NMF parameters.

# Define parameter grids
is.binary <- c(TRUE, FALSE, FALSE, FALSE, FALSE)
# Create a string that represents is.binary correctly in R code:
is_binary_str <- paste0("c(", paste(is.binary, collapse = ", "), ")")

n.runs.values <- c(50, 75, 100, 150, 200, 300)    
maxiter.values <- c(100, 200, 500, 1000, 2000) 
k.range.values <- list(2:10)
lr.values <- c(1e-4, 5e-4, 1e-3, 5e-3, 1e-2)
n.fold <- 10
allowParallel <- TRUE
n.cores <- 10

# We'll assume input_data.rds has a list of matrices for dat
input_dat_path <- "NMF_input.rds"

# Directories for generated scripts
script_dir <- "Scripts/single_algorithm/NMF_HPC/"
if (!dir.exists(script_dir)) {
  dir.create(script_dir, recursive = TRUE)
}

counter <- 0
for (lr in lr.values) {
  for (nr in n.runs.values) {
    for (mi in maxiter.values) {
      for (kr in k.range.values) {
        counter <- counter + 1
        job_name <- paste0("nmf_nruns_", nr, "_maxiter_", mi, "_lr_", lr)
        
        r_script_basename <- paste0(job_name, ".R")
        sh_script_basename <- paste0(job_name, ".sh")
        r_script_name <- paste0(script_dir, r_script_basename)
        sh_script_name <- paste0(script_dir, sh_script_basename)
        
        cat("
# ", r_script_name, "
library(foreach)
library(doParallel)
library(mclust)
library(MASS)

source('nmf_integrative_functions.R')

dat <- readRDS('", input_dat_path, "')

res <- nmf.opt.k.integrative(
  dat = dat,
  is.binary = ", is_binary_str, ",
  n.runs = ", nr, ",
  n.fold = ", n.fold, ",
  k.range = ", deparse(kr), ",
  maxiter = ", mi, ",
  lr = ", lr, ",
  allowParallel = ", allowParallel, ",
  n.cores = ", n.cores, ",
  make.plot = FALSE,
  result = TRUE
)

saveRDS(res, file = '", job_name, "_results.rds')
", file = r_script_name, sep="")

#
# Write the SLURM (.sh) script
#
cat("#!/bin/bash
#SBATCH --partition=icelake-himem
#SBATCH --nodes=1
#SBATCH --job-name=", job_name, "
#SBATCH --output=logs/", job_name, "_%j.out
#SBATCH --error=logs/", job_name, "_%j.err
#SBATCH --time=11:59:59
#SBATCH --cpus-per-task=", n.cores, "
#SBATCH --mem=50G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=as3582@cam.ac.uk

# Load R module
module load r-4.0.2-gcc-5.4.0-xyx46xb

# Set working directory
cd $HOME/MO_survey/IntNMF

# Set R_LIBS_USER to install packages in your home directory
export R_LIBS_USER=$HOME/MO_survey/R/library

# 1) Install required packages (except MASS) if not already installed
Rscript -e \"packages <- c('parallel','doParallel','foreach','mclust'); install.packages(setdiff(packages, installed.packages()[,'Package']), repos='https://cran.r-project.org')\"

# 2) Check if MASS 7.3-60.0.1 is installed; if not, install it
Rscript -e \"if (!requireNamespace('MASS', quietly=TRUE) || packageVersion('MASS') != '7.3.60.0.1') {
  if(!requireNamespace('remotes', quietly=TRUE)) {
    install.packages('remotes', repos='https://cran.r-project.org')
  }
  remotes::install_version('MASS', version='7.3-60.0.1', repos='https://cran.r-project.org')
}\"

# 3) Run R script
Rscript ", r_script_basename, "
", sep = "", file = sh_script_name)
      }
    }
  }
}

message(counter, " job scripts generated.")
