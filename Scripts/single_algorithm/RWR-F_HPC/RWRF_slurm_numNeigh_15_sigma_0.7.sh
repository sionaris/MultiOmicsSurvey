#!/bin/bash
#SBATCH --partition=icelake-himem
#SBATCH --nodes=1
#SBATCH --job-name=RWRF_15_0.7
#SBATCH --output=logs/RWRF_15_0.7.out
#SBATCH --error=logs/RWRF_15_0.7.err
#SBATCH --time=11:59:59
#SBATCH --cpus-per-task=20
#SBATCH --mem=100G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=as3582@cam.ac.uk

# Load R module
module purge
module load rhel8/cclake/base
module load r/4.4.0/gcc/3z5n2tph

# Set R_LIBS_USER to install packages in your home directory
export R_LIBS_USER=$HOME/MO_survey/R/library

# Create the library directory if it does not exist
mkdir -p "$R_LIBS_USER"

# Verify .libPaths()
echo "Current R library paths:"
Rscript -e "print(.libPaths())"

# Optional: install any missing packages
Rscript -e "packages <- c('SNFtool', 'foreach', 'doParallel'); to_install <- setdiff(packages, installed.packages()[,'Package']); if(length(to_install)>0) install.packages(to_install, repos='https://cran.r-project.org')"

# Run the R script
Rscript RWRF_script_numNeigh_15_sigma_0.7.R

