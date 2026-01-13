#!/usr/bin/env bash

module load rhel8/cclake/base
module load r/4.4.0/gcc/3z5n2tph

export R_LIBS_USER="$PWD/R/library"
mkdir -p "$R_LIBS_USER"

Rscript -e 'install.packages(c("data.table","jsonlite","openxlsx","mclust","cluster","Rfast"), repos="https://cran.r-project.org", lib=Sys.getenv("R_LIBS_USER"))'

# M3C (CRAN on many setups; if not, install from GitHub/other as you normally do)
Rscript -e 'install.packages("M3C", repos="https://cran.r-project.org", lib=Sys.getenv("R_LIBS_USER"))'

# MOFA2 (Bioc)
Rscript -e 'if(!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cran.r-project.org"); BiocManager::install("MOFA2", ask=FALSE, update=FALSE)'

