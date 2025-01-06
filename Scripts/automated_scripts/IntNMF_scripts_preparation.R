###############################################################################
## File: generate_nmf_parallel_k_scripts.R
##
## Purpose:
##   Generate HPC job scripts (one for each k in 2..10) that each run
##   "IntNMF parallely" for a single k-value in parallel.
##   We fix: n.runs=50, n.fold=10, maxiter=1000, num_cores=50, seed=123
##
###############################################################################

# Where to store the generated scripts:
script_dir <- "Scripts/single_algorithm/IntNMF_HPC/"
if (!dir.exists(script_dir)) {
  dir.create(script_dir, recursive = TRUE)
}

# Slurm HPC specs:
partition_name <- "icelake-himem"
time_limit     <- "11:59:59"
cpus_per_task  <- 50
memory_req     <- "500G"
email_user     <- "as3582@cam.ac.uk"

# NMF parallel function parameters:
n_runs   <- 50
n_fold   <- 10
max_iter <- 1000
num_cores <- 50
seed_val  <- 123  # for reproducibility
st_count  <- 10   # matches default logic from your function

# Range of k-values (each job handles exactly one k)
k_values <- 2:10

# Path to input data (adjust as needed)
input_data_path <- "IntNMF_input.rds"

# We'll create job files named "nmf_k2.R", "nmf_k2.sh", etc.
job_prefix <- "nmf_k"

counter <- 0
for (kval in k_values) {
  
  counter <- counter + 1
  job_name <- paste0(job_prefix, kval)   # e.g. "nmf_k2"
  
  # Filenames
  r_script_basename <- paste0(job_name, ".R")
  sh_script_basename <- paste0(job_name, ".sh")
  r_script_path <- file.path(script_dir, r_script_basename)
  sh_script_path <- file.path(script_dir, sh_script_basename)
  
  ############################################################################
  # 1) Generate the R script that calls "IntNMF parallely" for a single k
  ############################################################################
  cat(sprintf("
# HPC R script to run IntNMF parallely for k=%d (single k-range).

set.seed(%d)

suppressPackageStartupMessages({
  library(parallel)
  library(doParallel)
  library(foreach)
  library(MASS)
  library(mclust)
  library(IntNMF)
  source('parallel_IntNMF.R')
})

dat <- readRDS('%s')

res <- nmf.opt.k_parallel(
  dat       = dat,
  n.runs    = %d,
  n.fold    = %d,
  k.range   = c(%d, %d),
  maxiter   = %d,
  st.count  = %d,
  num_cores = %d,
  init_boolean = TRUE,
  n.init = 30,
  result    = TRUE,
  make.plot = FALSE,
  progress  = TRUE
)

saveRDS(res, file = '%s_results.rds')
writeLines(capture.output(sessionInfo()), '%s_sessionInfo.txt')
cat('IntNMF parallely finished for k=%d.\\n')
",
      kval,               # 1) %d (k)
      seed_val,           # 2) %d (seed)
      input_data_path,    # 3) %s
      n_runs,             # 4) %d
      n_fold,             # 5) %d
      kval,               # 6) %d
      kval,               # 7) %d
      max_iter,           # 8) %d
      st_count,           # 9) %d
      num_cores,          # 10) %d
      job_name,           # 11) %s -> '%s_results.rds'
      job_name,           # 12) %s -> '%s_sessionInfo.txt'
      kval                # 13) %d
  ),
  file = r_script_path, sep=""
  )

############################################################################
# 2) Generate the .sh script (Slurm)
############################################################################
cat(sprintf("#!/bin/bash
#SBATCH --partition=%s
#SBATCH --nodes=1
#SBATCH --job-name=%s
#SBATCH --output=logs/%s_%%j.out
#SBATCH --error=logs/%s_%%j.err
#SBATCH --time=%s
#SBATCH --cpus-per-task=%d
#SBATCH --mem=%s
#SBATCH --mail-type=ALL
#SBATCH --mail-user=%s

module purge
module load rhel8/cclake/base
module load r/4.4.0/gcc/3z5n2tph

# For reproducibility, set library path in your home directory
export R_LIBS_USER=$HOME/MO_survey/R/library
mkdir -p $R_LIBS_USER

echo '--- Starting job: %s (k=%d) ---'
echo 'Working directory:' `pwd`
echo 'R version:'
R --version

# Verify minimal needed packages are present, otherwise stop:
Rscript -e \"reqs <- c('IntNMF','mclust','MASS','doParallel','doRNG'); \\
missing_pkgs <- setdiff(reqs, rownames(installed.packages())); \\
if(length(missing_pkgs) > 0) {
  stop(paste0('ERROR: The following package(s) are missing: ', paste(missing_pkgs, collapse=', ')))
} else {
  cat('Packages OK\\n')
}\"

# Now run the R script
Rscript %s

echo '--- Job done: %s (k=%d) ---'
",
            partition_name,
            job_name,
            job_name,
            job_name,
            time_limit,
            cpus_per_task,
            memory_req,
            email_user,
            job_name, kval,
            r_script_basename,
            job_name, kval
), file = sh_script_path, sep="")
}

message(counter, ' job scripts generated in: ', script_dir)
