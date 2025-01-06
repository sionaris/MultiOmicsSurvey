#!/bin/bash
#SBATCH --partition=icelake-himem
#SBATCH --nodes=1
#SBATCH --job-name=nmf_k9
#SBATCH --output=logs/nmf_k9_%j.out
#SBATCH --error=logs/nmf_k9_%j.err
#SBATCH --time=11:59:59
#SBATCH --cpus-per-task=50
#SBATCH --mem=500G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=as3582@cam.ac.uk

module purge
module load rhel8/cclake/base
module load r/4.4.0/gcc/3z5n2tph

# For reproducibility, set library path in your home directory
export R_LIBS_USER=$HOME/MO_survey/R/library
mkdir -p $R_LIBS_USER

echo '--- Starting job: nmf_k9 (k=9) ---'
echo 'Working directory:' `pwd`
echo 'R version:'
R --version

# Verify minimal needed packages are present, otherwise stop:
Rscript -e "reqs <- c('IntNMF','mclust','MASS','doParallel','doRNG'); \
missing_pkgs <- setdiff(reqs, rownames(installed.packages())); \
if(length(missing_pkgs) > 0) {
  stop(paste0('ERROR: The following package(s) are missing: ', paste(missing_pkgs, collapse=', ')))
} else {
  cat('Packages OK\n')
}"

# Now run the R script
Rscript nmf_k9.R

echo '--- Job done: nmf_k9 (k=9) ---'
