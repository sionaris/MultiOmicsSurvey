# Define hyperparameter values
sdev_values <- c(0.01, 0.025, 0.05)
beta_var_scale_values <- c(0.3, 0.5, 0.8)

# Template for R script
r_script_template <- '
# Load necessary packages
library(iClusterPlus)
library(parallel)

# Read in the data input
input <- readRDS(paste0(getwd(), "/iCB_input.rds"))

# Define hyperparameters
sdev <- <sdev_value>
beta_var_scale <- <beta_var_scale_value>

# Generate an informative suffix for output files
suffix <- paste0("sdev_", sdev, "_beta_", beta_var_scale)

# Define the range of K
K_values <- 1:9

# Set up the number of cores for parallelization
num_cores <- 9  # Adjust based on available cores

# Define data types for each dataset
data_types <- c("binomial", "gaussian", "gaussian", "gaussian", "gaussian")

# Fix other hyperparameters
n_burnin <- 1000
n_draw <- 1200
prior_gamma <- rep(0.1, length(data_types))
thin <- 1
pp_cutoff <- 0.5

# Run tune.iClusterBayes in parallel
tune_results <- tune.iClusterBayes(
  cpus = num_cores,
  dt1 = input$SNPs,
  dt2 = input$RNAseq,
  dt3 = input$CNV,
  dt4 = input$miRNA,
  dt5 = input$Methylation,
  type = data_types,
  K = K_values,
  n.burnin = n_burnin,
  n.draw = n_draw,
  prior.gamma = prior_gamma,
  sdev = sdev,
  beta.var.scale = beta_var_scale,
  thin = thin,
  pp.cutoff = pp_cutoff
)

# Save the results with an informative filename
saveRDS(tune_results, paste0(getwd(), "/tune_results_", suffix, ".rds"))
'

# Template for SLURM script
slurm_script_template <- '#!/bin/bash
#SBATCH --partition=icelake-himem
#SBATCH --nodes=1
#SBATCH --job-name=iCB_sdev_<sdev>_beta_<beta>
#SBATCH --output=logs/iCB_sdev_<sdev>_beta_<beta>_%j.out
#SBATCH --error=logs/iCB_sdev_<sdev>_beta_<beta>_%j.err
#SBATCH --time=11:59:59
#SBATCH --cpus-per-task=9
#SBATCH --mem=100G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=as3582@cam.ac.uk

# Load R module
module load r-4.0.2-gcc-5.4.0-xyx46xb

# Set working directory
cd $HOME/MO_survey

# Set R_LIBS_USER to install packages in your home directory
export R_LIBS_USER=$HOME/MO_survey/R/library

# Install required packages if not already installed
Rscript -e \'packages <- c("parallel", "doParallel"); install.packages(setdiff(packages, installed.packages()[,"Package"]), repos="https://cran.r-project.org")\'

# Install iClusterPlus from Bioconductor
Rscript -e \'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos="https://cran.r-project.org"); BiocManager::install("iClusterPlus", lib=Sys.getenv("R_LIBS_USER"))\'

# Run R script
Rscript iClusterBayes_HPC_script_sdev_<sdev>_beta_<beta>.R
'

# Loop over hyperparameter combinations
for (sdev in sdev_values) {
  for (beta_var_scale in beta_var_scale_values) {
    # Skip the combination already running (sdev = 0.05, beta_var_scale = 0.5)
    if (sdev == 0.05 && beta_var_scale == 0.5) {
      next
    }
    
    # Create R script content
    r_script <- r_script_template
    r_script <- gsub("<sdev_value>", sdev, r_script)
    r_script <- gsub("<beta_var_scale_value>", beta_var_scale, r_script)
    
    # Write R script to file
    r_script_filename <- paste0("iClusterBayes_HPC_script_sdev_", sdev, "_beta_", beta_var_scale, ".R")
    writeLines(r_script, con = r_script_filename)
    
    # Create SLURM script content
    slurm_script <- slurm_script_template
    slurm_script <- gsub("<sdev>", sdev, slurm_script)
    slurm_script <- gsub("<beta>", beta_var_scale, slurm_script)
    
    # Write SLURM script to file
    slurm_script_filename <- paste0("iClusterBayes_sdev_", sdev, "_beta_", beta_var_scale, ".slurm")
    writeLines(slurm_script, con = slurm_script_filename)
  }
}
