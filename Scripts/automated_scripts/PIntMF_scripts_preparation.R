# Parameters
dimensions = c(2:10)
max.iters = c(20)
init_flavors = c("snf")
flavor_mods = c("glmnet")

r_script_template <- '
suppressPackageStartupMessages(library(PintMF))

# Load input data
input = readRDS("PIntNMF_input.rds")

RNGversion("4.2.2")
set.seed(123)

# Hyperparameter setup ###
p = <ndim>
max.it = <max.iter>
init_flavor = <flav_init>
flavor_mod = <flav_mod>

# Run PIntNMF
pintmf = SolveInt(
  Y=input, 
  p=p, 
  max.it=max.it, 
  verbose=FALSE, 
  init_flavor=init_flavor, 
  flavor_mod=flavor_mod
)

# Save the results
saveRDS(pintmf, "PIntMF_ndim_<ndim>_maxiter_<max.iter>_initflav_<flav_init>_modflav_<flav_mod>_results.rds")

# Also record sessionInfo
writeLines(capture.output(sessionInfo()), 
           "PIntMF_ndim_<ndim>_maxiter_<max.iter>_initflav_<flav_init>_modflav_<flav_mod>_sessionInfo.txt")
'

# SLURM
slurm_script_template <- '#!/bin/bash
#SBATCH --partition=icelake-himem
#SBATCH --nodes=1
#SBATCH --job-name=PIntMF_ndim_<ndim>_maxiter_<max.iter>_initflav_<flav_init>_modflav_<flav_mod>
#SBATCH --output=logs/PIntMF_ndim_<ndim>_maxiter_<max.iter>_initflav_<flav_init>_modflav_<flav_mod>.out
#SBATCH --error=logs/PIntMF_ndim_<ndim>_maxiter_<max.iter>_initflav_<flav_init>_modflav_<flav_mod>.err
#SBATCH --time=11:59:59
#SBATCH --cpus-per-task=10
#SBATCH --mem=100G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=as3582@cam.ac.uk

# Load R module
module purge
module load rhel8/cclake/base
module load r/4.4.0/gcc/3z5n2tph

# Set R_LIBS_USER to install packages in your home directory
export R_LIBS_USER=$HOME/MO_survey/R/library

# Create the library directory if it does not exist
mkdir -p "$R_LIBS_USER"

# Verify .libPaths()
echo "Current R library paths:"
Rscript -e "print(.libPaths())"

# Optional: install any missing packages
Rscript -e "packages <- c(\'PintMF\'); to_install <- setdiff(packages, installed.packages()[,\'Package\']); if(length(to_install)>0) devtools::install_github(\'mpierrejean/pintmf\')"

# Run the R script
Rscript PIntMF_ndim_<ndim>_maxiter_<max.iter>_initflav_<flav_init>_modflav_<flav_mod>.R
'
# Make sure directories exist
if(!dir.exists("Scripts/single_algorithm/PIntMF_HPC")) {
  dir.create("Scripts/single_algorithm/PIntMF_HPC", recursive = TRUE)
}

for (ndim in dimensions) {
  for (max.iter in max.iters) {
    for (flav_init in init_flavors) {
      for(flav_mod in flavor_mods) {
        # 1) Create the R script
        r_script <- r_script_template
        r_script <- gsub("<ndim>", ndim, r_script)
        r_script <- gsub("<max.iter>", max.iter, r_script)
        r_script <- gsub("<flav_init>", flav_init, r_script)
        r_script <- gsub("<flav_mod>", flav_mod, r_script)
        
        r_script_filename <- paste0("Scripts/single_algorithm/PIntMF_HPC/PIntMF_ndim_",
                                    ndim, "_maxiter_", max.iter, "_initflav_",
                                    flav_init, "_modflav_", flav_mod, ".R")
        writeLines(r_script, con = r_script_filename)
        
        # 2) Create the SLURM script
        slurm_script <- slurm_script_template
        slurm_script <- gsub("<ndim>", ndim, slurm_script)
        slurm_script <- gsub("<max.iter>", max.iter, slurm_script)
        slurm_script <- gsub("<flav_init>", flav_init, slurm_script)
        slurm_script <- gsub("<flav_mod>", flav_mod, slurm_script)
        
        slurm_script_filename <- paste0("Scripts/single_algorithm/PIntMF_HPC/bash_PIntMF_ndim_",
                                        ndim, "_maxiter_", max.iter, "_initflav_",
                                        flav_init, "_modflav_", flav_mod, ".sh")
        writeLines(slurm_script, con = slurm_script_filename)
      }
    }
  }
}
