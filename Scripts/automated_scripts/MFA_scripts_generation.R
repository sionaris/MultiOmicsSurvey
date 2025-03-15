# Define hyperparameter values to be exhaustively tested
ncp_values = c(2:10, 15, 20, 25)

# Template for R script
r_script_template <- '
library(FactoMineR)

# Set seed if necessary
RNGversion("4.2.2")
set.seed(123)

# Load input
mfa_input = readRDS("MFA_input.rds")

# Define number of components
ncp = <ncp_value>

# Run for ncp = <ncp_value>
snp_cols = sum(grepl("SNPs_", colnames(mfa_input)))
rna_cols = sum(grepl("RNAseq_", colnames(mfa_input)))
cnv_cols = sum(grepl("CNV_", colnames(mfa_input)))
mirna_cols = sum(grepl("miRNA_", colnames(mfa_input)))
methyl_cols = sum(grepl("Methylation_", colnames(mfa_input)))

t1 = Sys.time()
mfa = MFA(mfa_input,
          group = c(snp_cols, rna_cols, cnv_cols, mirna_cols, methyl_cols), # dimensionalities of the datasets
          type = c("n", rep("c", 4)),
          excl = NULL,
          ncp = ncp,
          name.group = c("SNPs", "RNAseq", "CNV", "miRNA", "Methylation"),
          graph = FALSE,
          axes = c(1,2))
dt = Sys.time() - t1
results = list(MFA = mfa, dt = dt)
rm(mfa); gc()

# Export
saveRDS(results, paste0("MFA_ncp_", ncp, "_results.rds"))
writeLines(capture.output(sessionInfo()), paste0("MFA_ncp_", ncp, "_sessionInfo.txt"))
'

# Template for SLURM script
slurm_script_template <- '#!/bin/bash
#SBATCH --partition=icelake-himem
#SBATCH --nodes=1
#SBATCH --job-name=MFA_ncp_<ncp_value>
#SBATCH --output=logs/MFA_ncp_<ncp_value>.out
#SBATCH --error=logs/MFA_ncp_<ncp_value>.err
#SBATCH --time=11:59:59
#SBATCH --cpus-per-task=20
#SBATCH --mem=300G
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
Rscript -e "packages <- c(\'FactoMineR\'); to_install <- setdiff(packages, installed.packages()[,\'Package\']); if(length(to_install)>0) install.packages(to_install, repos=\'https://cloud.r-project.org\')"

# Run R script
Rscript MFA_ncp_<ncp_value>.R
'

# Create directories if not exist
if(!dir.exists("Scripts/single_algorithm/MFA_HPC")) {
  dir.create("Scripts/single_algorithm/MFA_HPC", recursive = TRUE)
}

for (ncp_value in ncp_values) {
  # Create R script content
  r_script <- r_script_template
  r_script <- gsub("<ncp_value>", ncp_value, r_script)
  
  r_script_filename <- paste0("Scripts/single_algorithm/MFA_HPC/MFA_ncp_",
                              ncp_value, ".R")
  writeLines(r_script, con = r_script_filename)
  
  # Create SLURM script content
  slurm_script <- slurm_script_template
  slurm_script <- gsub("<ncp_value>", ncp_value, slurm_script)
  slurm_script_filename <- paste0("Scripts/single_algorithm/MFA_HPC/bash_MFA_ncp_",
                                  ncp_value, ".sh")
  writeLines(slurm_script, con = slurm_script_filename)
}