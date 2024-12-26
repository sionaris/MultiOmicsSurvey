#!/bin/bash
#SBATCH --partition=icelake-himem
#SBATCH --nodes=1
#SBATCH --job-name=iCB_sdev_0.025_beta_2.5
#SBATCH --output=logs/iCB_sdev_0.025_beta_2.5_%j.out
#SBATCH --error=logs/iCB_sdev_0.025_beta_2.5_%j.err
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
Rscript -e 'packages <- c("parallel", "doParallel"); install.packages(setdiff(packages, installed.packages()[,"Package"]), repos="https://cran.r-project.org")'

# Install iClusterPlus from Bioconductor
Rscript -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos="https://cran.r-project.org"); BiocManager::install("iClusterPlus", lib=Sys.getenv("R_LIBS_USER"))'

# Run R script
Rscript iClusterBayes_HPC_script_sdev_0.025_beta_2.5.R

