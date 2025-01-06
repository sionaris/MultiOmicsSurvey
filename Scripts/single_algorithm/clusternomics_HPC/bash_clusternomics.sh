#!/bin/bash
#SBATCH --partition=icelake-himem
#SBATCH --nodes=1
#SBATCH --job-name=clusternomics
#SBATCH --output=logs/clusternomics_%j.out
#SBATCH --error=logs/clusternomics_%j.err
#SBATCH --time=11:59:59
#SBATCH --cpus-per-task=30
#SBATCH --mem=300G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=as3582@cam.ac.uk

module purge
module load rhel8/cclake/base
module load r/4.4.0/gcc/3z5n2tph

echo '--- Starting job: clusternomics ---'
echo 'Working directory:' `pwd`
echo 'R version:'
R --version

# Navigate to HOME directory
cd $HOME/MO_survey/clusternomics

# For reproducibility, set library path in your home directory
export R_LIBS_USER=$HOME/MO_survey/R/library
mkdir -p $R_LIBS_USER

# Verify .libPaths()
Rscript -e "print(.libPaths())"

# Install parallel and clusternomics if not already installed
Rscript -e "
packages <- c('parallel', 'devtools')
# Exclude packages that are already installed
packages_to_install <- setdiff(packages, rownames(installed.packages(lib.loc='$R_LIBS_USER')))
if(length(packages_to_install) > 0){
  install.packages(packages_to_install, repos='https://cran.r-project.org', dependencies=TRUE, lib='$R_LIBS_USER')
}
"

Rscript -e "
packages <- c('clusternomics')
# Exclude packages that are already installed
packages_to_install <- setdiff(packages, rownames(installed.packages(lib.loc='$R_LIBS_USER')))
if(length(packages_to_install) > 0){
  devtools::install_github('evelinag/clusternomics')
}
"

# Verify installations and load necessary libraries
Rscript -e "
packages <- c('parallel', 'devtools', 'clusternomics')
for(pkg in packages){
  if(!require(pkg, character.only=TRUE, lib.loc='$R_LIBS_USER')){
    stop(paste('Package', pkg, 'failed to install.'))
  }
}
"

# Run script
Rscript run_clusternomics.R