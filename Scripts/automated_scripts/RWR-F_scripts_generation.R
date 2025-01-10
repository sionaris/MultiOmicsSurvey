##########################
## 1) Define the ranges  ##
##########################
num_neighbors_range = seq(10, 50, 5)
sigma_range = seq(0.3, 0.8, 0.1)

# RWR-F hyperparameters (fixed)
RWR_iteration_max     <- 1000
RWRF_gamma_fixed      <- 0.7
RWRNF_gamma_fixed     <- 0.7
RWRNF_num_neighbors_fixed <- 10
RWRNF_alpha_fixed     <- 0.9
RWRNF_beta_fixed      <- 0.9

###############################
## 2) Template for R script  ##
###############################
r_script_template <- '
# Load libraries
library(SNFtool)
library(foreach)
library(doParallel)

# Source RWR-F
source("RWR-F_source.R")

# Load input data
distL <- readRDS("RWR-F_input_dists.rds")

# Define hyperparameters
num_neighbors <- <neighbors_value>
sigma_val     <- <sigma_value>

# Prepare the cluster for parallelization (up to 20 cores)
num_cores <- 20
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Construct affinity matrices
affinityL <- list()
for (modnm in names(distL)) {
  aff <- affinityMatrix(distL[[modnm]], K = num_neighbors, sigma = sigma_val)
  affinityL[[modnm]] <- aff
}

# RWRF
fused_rwrf <- RWR_fusion(
  sim_list      = affinityL,
  iteration_max = <rwr_iter_max>,
  gama          = <rwrf_gamma>
)

# RWRNF
fused_rwrnf <- RWR_fusion_neighbor(
  sim_list      = affinityL,
  iteration_max = <rwr_iter_max>,
  gama          = <rwrnf_gamma>,
  neighbor_num  = <rwrnf_neighbors>,
  alpha         = <rwrnf_alpha>,
  beta          = <rwrnf_beta>
)

stopCluster(cl)

# Create a fused object to save
fusion <- list(
  aff_num_neighbors  = num_neighbors,
  aff_sigma          = sigma_val,
  affinity           = affinityL,
  fused_rwrf         = fused_rwrf,
  fused_rwrnf        = fused_rwrnf
)

# Save the results
saveRDS(fusion, paste0("NN_", num_neighbors, "_sigma_", sigma_val, "_fusion.rds"))

# Also record sessionInfo
writeLines(capture.output(sessionInfo()), paste0("NN_", num_neighbors, "_sigma_", sigma_val, "_sessionInfo.txt"))
'

###################################
## 3) Template for SLURM script  ##
###################################
slurm_script_template <- '#!/bin/bash
#SBATCH --partition=icelake-himem
#SBATCH --nodes=1
#SBATCH --job-name=RWRF_<neighbors>_<sigma>
#SBATCH --output=logs/RWRF_<neighbors>_<sigma>.out
#SBATCH --error=logs/RWRF_<neighbors>_<sigma>.err
#SBATCH --time=11:59:59
#SBATCH --cpus-per-task=20
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
Rscript -e "packages <- c(\'SNFtool\', \'foreach\', \'doParallel\'); to_install <- setdiff(packages, installed.packages()[,\'Package\']); if(length(to_install)>0) install.packages(to_install, repos=\'https://cran.r-project.org\')"

# Run the R script
Rscript RWRF_script_numNeigh_<neighbors>_sigma_<sigma>.R
'

# Make sure directories exist
if(!dir.exists("Scripts/single_algorithm/RWR-F_HPC")) {
  dir.create("Scripts/single_algorithm/RWR-F_HPC", recursive = TRUE)
}

for (nn in num_neighbors_range) {
  for (sig in sigma_range) {
    # 1) Create the R script
    r_script <- r_script_template
    r_script <- gsub("<neighbors_value>", nn, r_script)
    r_script <- gsub("<sigma_value>", sig, r_script)
    r_script <- gsub("<rwr_iter_max>", RWR_iteration_max, r_script)
    r_script <- gsub("<rwrf_gamma>", RWRF_gamma_fixed, r_script)
    r_script <- gsub("<rwrnf_gamma>", RWRNF_gamma_fixed, r_script)
    r_script <- gsub("<rwrnf_neighbors>", RWRNF_num_neighbors_fixed, r_script)
    r_script <- gsub("<rwrnf_alpha>", RWRNF_alpha_fixed, r_script)
    r_script <- gsub("<rwrnf_beta>", RWRNF_beta_fixed, r_script)
    
    r_script_filename <- paste0("Scripts/single_algorithm/RWR-F_HPC/RWRF_script_numNeigh_",
                                nn, "_sigma_", sig, ".R")
    writeLines(r_script, con = r_script_filename)
    
    # 2) Create the SLURM script
    slurm_script <- slurm_script_template
    slurm_script <- gsub("<neighbors>", nn, slurm_script)
    slurm_script <- gsub("<sigma>", sig, slurm_script)
    
    slurm_script_filename <- paste0("Scripts/single_algorithm/RWR-F_HPC/RWRF_slurm_numNeigh_",
                                    nn, "_sigma_", sig, ".sh")
    writeLines(slurm_script, con = slurm_script_filename)
  }
}
