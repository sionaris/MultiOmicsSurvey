#!/bin/bash
#SBATCH --partition=icelake-himem
#SBATCH --nodes=1
#SBATCH --job-name=KLIC_hlcust_input_down_KLIC2
#SBATCH --output=logs/KLIC_hlcust_input_down_KLIC2.out
#SBATCH --error=logs/KLIC_hlcust_input_down_KLIC2.err
#SBATCH --time=11:59:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=20G
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
Rscript -e "packages <- c('pathfindR', 'fastcluster'); to_install <- setdiff(packages, installed.packages()[,'Package']); if(length(to_install)>0) install.packages(to_install)"

# Run script
Rscript KLIC_hclust_input_down_KLIC2.R