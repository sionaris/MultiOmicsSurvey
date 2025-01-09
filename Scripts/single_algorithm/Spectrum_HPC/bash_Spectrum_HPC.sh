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

module load r/4.4.0/gcc/3z5n2tph

# Set R_LIBS_USER to install packages in your home directory
export R_LIBS_USER=$HOME/MO_survey/R/library

# Create the library directory if it doesn't exist
mkdir -p "$R_LIBS_USER"

# Verify .libPaths()
echo "Current R library paths:"
Rscript -e "print(.libPaths())"

# ----------------------------------------
# 2. Install Specific Version of MASS Package
# ----------------------------------------

echo "Checking and installing MASS package if necessary..."
Rscript -e "
if (!requireNamespace('MASS', quietly=TRUE) || packageVersion('MASS') != '7.3-60.0.1') {
  if(!requireNamespace('remotes', quietly=TRUE)) {
    install.packages('remotes', repos='https://cran.r-project.org', lib='$R_LIBS_USER')
  }
  remotes::install_version('MASS', version='7.3-60.0.1', repos='https://cran.r-project.org', lib='$R_LIBS_USER')
}
"

echo "Installing BiocManager and Biobase (if not already installed)..."
Rscript -e "
if (!requireNamespace('BiocManager', quietly=TRUE)) {
  install.packages('BiocManager', repos='https://cran.r-project.org', lib='$R_LIBS_USER')
}
BiocManager::install(version='3.20', update=FALSE, ask=FALSE, lib='$R_LIBS_USER')  # Bioc 3.20 for R 4.4
BiocManager::install('Biobase', update=FALSE, ask=FALSE, lib='$R_LIBS_USER')
"

# ----------------------------------------
# 3. Install Required CRAN Packages
# ----------------------------------------

echo "Installing required CRAN packages (parallel, foreach, ClusterR, Rfast, Spectrum)..."
Rscript -e "
required_packages <- c('parallel', 'foreach', 'ggplot2', 'ClusterR', 'Rfast', 'Spectrum')
installed_packages <- rownames(installed.packages(lib.loc='$R_LIBS_USER'))
packages_to_install <- setdiff(required_packages, installed_packages)

if(length(packages_to_install) > 0){
  install.packages(packages_to_install, repos='https://cran.r-project.org', dependencies=TRUE, lib='$R_LIBS_USER')
}
"

# ----------------------------------------
# 4. Verify Installations and Load Necessary Libraries
# ----------------------------------------

echo "Verifying installations and loading necessary libraries..."
Rscript -e "
required_packages <- c('parallel', 'foreach', 'ClusterR', 'Rfast', 'Spectrum', 'MASS')
for(pkg in required_packages){
  if(!require(pkg, character.only=TRUE, lib.loc='$R_LIBS_USER')){
    stop(paste('Package', pkg, 'failed to install or load.'))
  } else {
    cat(paste('Package', pkg, 'loaded successfully.\n'))
  }
}
"

# ----------------------------------------
# 5. Run the Spectrum R Script
# ----------------------------------------

echo "Running Spectrum_HPC_script.R..."
Rscript Spectrum_HPC_script.R

echo "Spectrum analysis completed successfully."

