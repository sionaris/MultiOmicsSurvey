#!/bin/bash
#SBATCH --partition=icelake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=LRAcluster
#SBATCH --output=logs/LRAcluster_%j.out
#SBATCH --error=logs/LRAcluster_%j.err
#SBATCH --time=11:59:59
#SBATCH --cpus-per-task=10
#SBATCH --mem=90G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=as3582@cam.ac.uk

# Set working directory
cd $HOME/MO_survey/LRAcluster

# Load R module
module load r-4.0.2-gcc-5.4.0-xyx46xb

# Set R_LIBS_USER to install packages in your home directory
export R_LIBS_USER=$HOME/MO_survey/R/library

# Install required packages if not already installed
Rscript -e "packages <- c('parallel','doParallel','foreach'); install.packages(setdiff(packages, installed.packages()[,'Package']), repos='https://cran.r-project.org')"

# Run R script
Rscript LRAcluster_HPC.R
