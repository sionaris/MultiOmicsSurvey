#!/bin/bash
#SBATCH --partition=icelake-himem
#SBATCH -A SIMIDJIEVSKI-SL3-CPU
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --job-name=MSNE_feat_merge
#SBATCH --output=logs/MSNE_feat_merge_%j.out
#SBATCH --error=logs/MSNE_feat_merge_%j.err

set -euo pipefail

export SLURM_TRES_PER_TASK="cpu=${SLURM_CPUS_PER_TASK}"

: "${ARRAY_JOB_ID:?Must pass ARRAY_JOB_ID via --export=ALL,ARRAY_JOB_ID=<id>}"

METHOD="MSNE"
RUN_DIR="${SLURM_SUBMIT_DIR}/Results/Performance/Feature_perturbations/${METHOD}/array_${ARRAY_JOB_ID}"
mkdir -p "${RUN_DIR}/out" "${RUN_DIR}/meta" "${SLURM_SUBMIT_DIR}/logs"

export BENCH_RUN_DIR="$(readlink -f "${RUN_DIR}")"

module purge
module load rhel8/cclake/base
module load r/4.4.0/gcc/3z5n2tph

export R_LIBS_USER="${SLURM_SUBMIT_DIR}/R/library"
mkdir -p "$R_LIBS_USER"

srun --hint=nomultithread --cpu-bind=cores \
  Rscript "${SLURM_SUBMIT_DIR}/Jobs/MSNE/MSNE_feature_experiments.R" merge

