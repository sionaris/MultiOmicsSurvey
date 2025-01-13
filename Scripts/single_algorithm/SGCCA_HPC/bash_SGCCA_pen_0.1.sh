#!/bin/bash
#SBATCH --partition=icelake-himem
#SBATCH --nodes=1
#SBATCH --job-name=SGCCA_pen_0.1
#SBATCH --output=logs/SGCCA_pen_0.1.out
#SBATCH --error=logs/SGCCA_pen_0.1.err
#SBATCH --time=11:59:59
#SBATCH --cpus-per-task=5
#SBATCH --mem=200G
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
Rscript -e "packages <- c('mixOmics'); to_install <- setdiff(packages, installed.packages()[,'Package']); if(length(to_install)>0) BiocManager::install(to_install)"

# Run the R script
Rscript SGCCA_pen_0.1.R

