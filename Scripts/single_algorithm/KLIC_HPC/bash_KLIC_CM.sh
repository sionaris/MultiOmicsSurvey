#!/bin/bash
#SBATCH --partition=icelake-himem
#SBATCH --nodes=1
#SBATCH --job-name=klic_cm_run
#SBATCH --output=logs/klic_cm_run.out
#SBATCH --error=logs/klic_cm_run.err
#SBATCH --time=11:59:59
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=9
#SBATCH --mem=100G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=as3582@cam.ac.uk

module purge
module load rhel8/cclake/base
module load r/4.4.0/gcc/3z5n2tph

# For reproducibility, set library path in your home directory
export R_LIBS_USER=$HOME/MO_survey/R/library
mkdir -p $R_LIBS_USER

echo '--- Starting job: klic_installation  ---'
echo 'Working directory:' `pwd`
echo 'R version:'
R --version

# Run script
Rscript KLIC_CM.R