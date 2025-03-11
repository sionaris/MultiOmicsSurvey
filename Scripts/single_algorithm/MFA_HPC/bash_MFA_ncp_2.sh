#!/bin/bash
#SBATCH --partition=icelake-himem
#SBATCH --nodes=1
#SBATCH --job-name=MFA_ncp_2
#SBATCH --output=logs/MFA_ncp_2.out
#SBATCH --error=logs/MFA_ncp_2.err
#SBATCH --time=11:59:59
#SBATCH --cpus-per-task=20
#SBATCH --mem=300G
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
Rscript -e "packages <- c('FactoMineR'); to_install <- setdiff(packages, installed.packages()[,'Package']); if(length(to_install)>0) install.packages(to_install, repos='https://cloud.r-project.org')"

# Run R script
Rscript MFA_ncp_2.R

