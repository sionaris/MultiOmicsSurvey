###########################
## 1) Penalties for SGCCA #
###########################
penalties = seq(0.1, 0.9, 0.1)

###############################
## 2) Template for R script  ##
###############################
r_script_template <- '
suppressPackageStartupMessages(library(mixOmics))

# Load input data
input = readRDS("SGCCA_input.rds")

RNGversion("4.2.2")
set.seed(123)

# Hyperparameter setup ###
penalty = <penalty_value>
design = 1 - diag(length(input))
ncomp = 10
scheme = "horst"
keepX = NULL
scale = FALSE
max.iter = 1000 # default
init = "svd.single"
near.zero.var = FALSE

# Run SGCCA
sgcca = wrapper.sgcca(
  input,
  penalty = penalty,
  design = design,
  ncomp = ncomp,
  keepX = keepX,
  scheme = scheme,
  scale = scale,
  init = init,
  tol = .Machine$double.eps, # default
  max.iter = max.iter,
  near.zero.var = near.zero.var,
  all.outputs = TRUE
)

# Save the results
saveRDS(sgcca, "sgcca_pen_<penalty_value>_results.rds")

# Also record sessionInfo
writeLines(capture.output(sessionInfo()), 
           "sgcca_pen_<penalty_value>_sessionInfo.txt")
'

###################################
## 3) Template for SLURM script  ##
###################################
slurm_script_template <- '#!/bin/bash
#SBATCH --partition=icelake-himem
#SBATCH --nodes=1
#SBATCH --job-name=SGCCA_pen_<penalty_value>
#SBATCH --output=logs/SGCCA_pen_<penalty_value>.out
#SBATCH --error=logs/SGCCA_pen_<penalty_value>.err
#SBATCH --time=11:59:59
#SBATCH --cpus-per-task=5
#SBATCH --mem=200G
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
Rscript -e "packages <- c(\'mixOmics\'); to_install <- setdiff(packages, installed.packages()[,\'Package\']); if(length(to_install)>0) BiocManager::install(to_install)"

# Run the R script
Rscript SGCCA_pen_<penalty_value>.R
'

# Make sure directories exist
if(!dir.exists("Scripts/single_algorithm/SGCCA_HPC")) {
  dir.create("Scripts/single_algorithm/SGCCA_HPC", recursive = TRUE)
}

for (pen in penalties) {
  # 1) Create the R script
  r_script <- r_script_template
  r_script <- gsub("<penalty_value>", pen, r_script)
  
  r_script_filename <- paste0("Scripts/single_algorithm/SGCCA_HPC/SGCCA_pen_",
                              pen, ".R")
  writeLines(r_script, con = r_script_filename)
  
  # 2) Create the SLURM script
  slurm_script <- slurm_script_template
  slurm_script <- gsub("<penalty_value>", pen, slurm_script)
  
  slurm_script_filename <- paste0("Scripts/single_algorithm/SGCCA_HPC/bash_SGCCA_pen_",
                                  pen, ".sh")
  writeLines(slurm_script, con = slurm_script_filename)
}
