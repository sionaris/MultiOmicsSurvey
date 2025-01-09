#!/bin/bash
#SBATCH --partition=icelake-himem
#SBATCH --nodes=1
#SBATCH --job-name=Spectrum_HPC_run
#SBATCH --output=logs/Spectrum_HPC_run.out
#SBATCH --error=logs/Spectrum_HPC_run.err
#SBATCH --time=11:59:59
#SBATCH --cpus-per-task=5
#SBATCH --mem=100G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=as3582@cam.ac.uk

# ----------------------------------------
# 1. Set Up R Library Path
# ----------------------------------------
module purge
module load rhel8/cclake/base
module load r/4.4.0/gcc/3z5n2tph

# Set R_LIBS_USER to install packages in your home directory
export R_LIBS_USER=$HOME/MO_survey/R/library

# Create the library directory if it doesn't exist
mkdir -p "$R_LIBS_USER"

# Verify .libPaths()
echo "Current R library paths:"
Rscript -e "print(.libPaths())"

# Run the script

echo "Running Spectrum_HPC_script.R..."
Rscript Spectrum_HPC_script.R

echo "Spectrum analysis completed successfully."
