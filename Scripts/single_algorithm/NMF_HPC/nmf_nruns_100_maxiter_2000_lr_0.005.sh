#!/bin/bash
#SBATCH --partition=icelake-himem
#SBATCH --nodes=1
#SBATCH --job-name=nmf_nruns_100_maxiter_2000_lr_0.005
#SBATCH --output=logs/nmf_nruns_100_maxiter_2000_lr_0.005_%j.out
#SBATCH --error=logs/nmf_nruns_100_maxiter_2000_lr_0.005_%j.err
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

# 1) Install required packages (except MASS) if not already installed
Rscript -e "packages <- c('parallel','doParallel','foreach','mclust'); install.packages(setdiff(packages, installed.packages()[,'Package']), repos='https://cran.r-project.org')"

# 2) Check if MASS 7.3-60.0.1 is installed; if not, install it
Rscript -e "if (!requireNamespace('MASS', quietly=TRUE) || packageVersion('MASS') != '7.3.60.0.1') {
  if(!requireNamespace('remotes', quietly=TRUE)) {
    install.packages('remotes', repos='https://cran.r-project.org')
  }
  remotes::install_version('MASS', version='7.3-60.0.1', repos='https://cran.r-project.org')
}"

# 3) Run R script
Rscript nmf_nruns_100_maxiter_2000_lr_0.005.R
