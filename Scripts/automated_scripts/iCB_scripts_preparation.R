# Define hyperparameter values to be exhaustively tested
sdev_values <- c(0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.05)
beta_var_scale_values <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.8, 1.0)

# Fixed hyperparameters based on the paper
thin <- 1
pp_cutoff <- 0.5
# Assuming we have 5 data sets as per your previous scripts
prior_gamma <- rep(0.1, 5)
n_burnin <- 3000
n_draw <- 4000

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

# Fixed hyperparameters
thin <- <thin_value>
pp_cutoff <- <pp_cutoff_value>
n_burnin <- <n_burnin_value>
n_draw <- <n_draw_value>
prior_gamma <- c(<prior_gamma_values>)

# Generate an informative suffix for output files
suffix <- paste0("sdev_", sdev, "_beta_", beta_var_scale)

# Define the range of K
K_values <- 1:9  # Adjust if needed based on your data

# Set up the number of cores for parallelization
num_cores <- 9  # Adjust based on available cores and HPC capacity

# Define data types for each dataset (adjust if needed)
data_types <- c("binomial", "gaussian", "gaussian", "gaussian", "gaussian")

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
cd $HOME/MO_survey/iCB_HPC

# Set R_LIBS_USER to install packages in your home directory
export R_LIBS_USER=$HOME/MO_survey/R/library

# Install required packages if not already installed
Rscript -e \'packages <- c("parallel", "doParallel"); install.packages(setdiff(packages, installed.packages()[,"Package"]), repos="https://cran.r-project.org")\'

# Install iClusterPlus from Bioconductor
Rscript -e \'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos="https://cran.r-project.org"); BiocManager::install("iClusterPlus", lib=Sys.getenv("R_LIBS_USER"))\'

# Run R script
Rscript iClusterBayes_HPC_script_sdev_<sdev>_beta_<beta>.R
'

# Create directories if not exist
if(!dir.exists("Scripts/single_algorithm/iCB_HPC")) {
  dir.create("Scripts/single_algorithm/iCB_HPC", recursive = TRUE)
}
if(!dir.exists("logs")) {
  dir.create("logs")
}

# Convert prior_gamma to a string
pg_str <- paste(prior_gamma, collapse=",")

# Generate R and SLURM scripts for all combinations
for (sdev in sdev_values) {
  for (beta_var_scale in beta_var_scale_values) {
    # Create R script content
    r_script <- r_script_template
    r_script <- gsub("<sdev_value>", sdev, r_script)
    r_script <- gsub("<beta_var_scale_value>", beta_var_scale, r_script)
    r_script <- gsub("<thin_value>", thin, r_script)
    r_script <- gsub("<pp_cutoff_value>", pp_cutoff, r_script)
    r_script <- gsub("<n_burnin_value>", n_burnin, r_script)
    r_script <- gsub("<n_draw_value>", n_draw, r_script)
    r_script <- gsub("<prior_gamma_values>", pg_str, r_script)
    
    r_script_filename <- paste0("Scripts/single_algorithm/iCB_HPC/iClusterBayes_HPC_script_sdev_",
                                sdev, "_beta_", beta_var_scale, ".R")
    writeLines(r_script, con = r_script_filename)
    
    # Create SLURM script content
    slurm_script <- slurm_script_template
    slurm_script <- gsub("<sdev>", sdev, slurm_script)
    slurm_script <- gsub("<beta>", beta_var_scale, slurm_script)
    
    slurm_script_filename <- paste0("Scripts/single_algorithm/iCB_HPC/iClusterBayes_sdev_",
                                    sdev, "_beta_", beta_var_scale, ".sh")
    writeLines(slurm_script, con = slurm_script_filename)
  }
}
