# This script generates .sh and .R scripts for each combination of NMF parameters.

# Define parameter grids
is.binary <- c(TRUE, FALSE, FALSE, FALSE, FALSE)
# Create a string that represents is.binary correctly in R code:
is_binary_str <- paste0("c(", paste(is.binary, collapse = ", "), ")")

n.runs.values <- c(30, 50, 75, 100)    
maxiter.values <- c(50, 100, 200, 500, 1000) 
k.range.values <- list(2:10)
n.fold <- 10
allowParallel <- TRUE
n.cores <- 10

# We'll assume input_data.rds has a list of matrices for dat
input_dat_path <- "NMF_input.rds"

# Directories for generated scripts
script_dir <- "Scripts/single_algorithm/NMF_HPC/"
if(!dir.exists(script_dir)) {
  dir.create(script_dir, recursive = TRUE)
}

counter <- 0
for (nr in n.runs.values) {
  for (mi in maxiter.values) {
    for (kr in k.range.values) {
      counter <- counter + 1
      job_name <- paste0("nmf_nruns_", nr, "_maxiter_", mi)
      
      # Filenames for scripts
      r_script_basename <- paste0(job_name, ".R")
      sh_script_basename <- paste0(job_name, ".sh")
      r_script_name <- paste0(script_dir, r_script_basename)
      sh_script_name <- paste0(script_dir, sh_script_basename)
      
      # Write the R script
      # Use cat with careful formatting for readability:
      cat(
        "#
# ", r_script_name, "
library(foreach)
library(doParallel)
library(mclust)
library(MASS)

source('nmf_integrative_functions.R') # Ensure this file has nmf.opt.k.integrative defined

# Load data
dat <- readRDS('", input_dat_path, "')

res <- nmf.opt.k.integrative(
  dat = dat,
  is.binary = ", is_binary_str, ",
  n.runs = ", nr, ",
  n.fold = ", n.fold, ",
  k.range = ", deparse(kr), ",
  maxiter = ", mi, ",
  allowParallel = ", allowParallel, ",
  n.cores = ", n.cores, ",
  make.plot = FALSE,
  result = TRUE
)

saveRDS(res, file = '", job_name, "_results.rds')
", sep="", file = r_script_name)

# Write the SLURM (.sh) script
# Note the quoting around Rscript -e command
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
cd $HOME/MO_survey

# Set R_LIBS_USER to install packages in your home directory
export R_LIBS_USER=$HOME/MO_survey/R/library

# Install required packages if not already installed
Rscript -e \"packages <- c('parallel','doParallel','foreach','mclust','MASS'); install.packages(setdiff(packages, installed.packages()[,'Package']), repos='https://cran.r-project.org')\"

# Run R script
Rscript ", r_script_basename, "
", sep="", file = sh_script_name)

    }
  }
}

message(counter, " job scripts generated.")
