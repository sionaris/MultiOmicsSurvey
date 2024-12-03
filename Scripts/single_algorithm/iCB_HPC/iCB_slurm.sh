#!/bin/bash
#SBATCH --partition=icelake
#SBATCH --nodes=1                # Use 1 nodes
#SBATCH --job-name=iClusterBayes
#SBATCH --output=iClusterBayes_%j.out
#SBATCH --error=iClusterBayes_%j.err
#SBATCH --time=23:59:59             # Max runtime
#SBATCH --cpus-per-task=9          # Number of cores
#SBATCH --mem=100G                  # Memory required
#SBATCH --mail-type=ALL             # Get email notifications
#SBATCH --mail-user=as3582@cam.ac.uk  # Your email address

# Load R module (adjust module name based on your HPC)
module load R/4.1.0

# Set R_LIBS_USER to install packages in your home directory
export R_LIBS_USER=$HOME/MO_survey/R/library

# Install required packages if not already installed
Rscript -e 'packages <- c("iClusterPlus", "parallel", "doParallel"); \
             install.packages(setdiff(packages, installed.packages()[,"Package"]), repos="http://cran.r-project.org")'

# Run R script
Rscript iClusterBayes_HPC_script.R
