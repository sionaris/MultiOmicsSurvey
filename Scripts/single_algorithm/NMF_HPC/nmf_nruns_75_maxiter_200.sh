#!/bin/bash
#SBATCH --partition=icelake-himem
#SBATCH --nodes=1
#SBATCH --job-name=nmf_nruns_75_maxiter_200
#SBATCH --output=logs/nmf_nruns_75_maxiter_200_%j.out
#SBATCH --error=logs/nmf_nruns_75_maxiter_200_%j.err
#SBATCH --time=11:59:59
#SBATCH --cpus-per-task=10
#SBATCH --mem=50G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=as3582@cam.ac.uk

# Load R module
module load r-4.0.2-gcc-5.4.0-xyx46xb

# Set working directory
cd $HOME/MO_survey/IntNMF

# Set R_LIBS_USER to install packages in your home directory
export R_LIBS_USER=$HOME/MO_survey/R/library

# Install required packages if not already installed
Rscript -e "packages <- c('parallel','doParallel','foreach','mclust','MASS'); install.packages(setdiff(packages, installed.packages()[,'Package']), repos='https://cran.r-project.org')"

# Run R script
Rscript nmf_nruns_75_maxiter_200.R
