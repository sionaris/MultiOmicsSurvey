#!/usr/bin/env bash
set -euo pipefail

# --- modules (adjust if your cluster uses different names) ---
module purge
module load rhel8/cclake/base
module load r/4.4.0/gcc/3z5n2tph

# --- user R library ---
export R_LIBS_USER="${R_LIBS_USER:-$HOME/MO_survey/R/library}"
mkdir -p "$R_LIBS_USER"

echo "Using R_LIBS_USER=$R_LIBS_USER"
Rscript -e 'print(R.version.string); print(.libPaths())'

# --- install devtools (if missing), then install wMKL from GitHub ---
Rscript -e 'if (!requireNamespace("devtools", quietly=TRUE)) install.packages("devtools", repos="https://cran.r-project.org")'
Rscript -e 'devtools::install_github("Shamir-Lab/NEMO/NEMO")'

echo "Done."

