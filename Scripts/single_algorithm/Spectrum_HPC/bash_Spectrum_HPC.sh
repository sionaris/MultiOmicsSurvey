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

# Load the correct R module
module load r/4.4.0/gcc/3z5n2tph  # Ensure this is the full R module

# Verify R version
Rscript -e "cat('R version:', R.version.string, '\n')"

# Set working directory
cd $HOME/MO_survey/IntNMF

# Set R_LIBS_USER to install packages in your home directory
export R_LIBS_USER=$HOME/MO_survey/R/library

# Create the library directory if it doesn't exist
mkdir -p $R_LIBS_USER

# Verify .libPaths()
Rscript -e "print(.libPaths())"

# 1) Install BiocManager and Bioconductor packages
Rscript -e "
if (!requireNamespace('BiocManager', quietly=TRUE)) {
  install.packages('BiocManager', repos='https://cran.r-project.org', lib='$R_LIBS_USER')
}
BiocManager::install(version='3.20', update=FALSE, ask=FALSE, lib='$R_LIBS_USER')  # Bioc 3.20 for R 4.4
BiocManager::install('Biobase', update=FALSE, ask=FALSE, lib='$R_LIBS_USER')
"

# 2) Install required CRAN packages, including NMF and its dependencies
Rscript -e "
packages <- c('parallel', 'foreach', 'Spectrum')
# Exclude packages that are already installed
packages_to_install <- setdiff(packages, rownames(installed.packages(lib.loc='$R_LIBS_USER')))
if(length(packages_to_install) > 0){
  install.packages(packages_to_install, repos='https://cran.r-project.org', dependencies=TRUE, lib='$R_LIBS_USER')
}
"

# 3) Check and install specific version of MASS if required
Rscript -e "
if (!requireNamespace('MASS', quietly=TRUE) || packageVersion('MASS') != '7.3-60.0.1') {
  if(!requireNamespace('remotes', quietly=TRUE)) {
    install.packages('remotes', repos='https://cran.r-project.org', lib='$R_LIBS_USER')
  }
  remotes::install_version('MASS', version='7.3-60.0.1', repos='https://cran.r-project.org', lib='$R_LIBS_USER')
}
"

# 4) Verify installations and load necessary libraries
Rscript -e "
packages <- c('parallel', 'foreach', 'Spectrum', 'MASS')
for(pkg in packages){
  if(!require(pkg, character.only=TRUE, lib.loc='$R_LIBS_USER')){
    stop(paste('Package', pkg, 'failed to install.'))
  }
}
"

# 5) Run R script
Rscript Spectrum_HPC_script.R