#!/usr/bin/env bash
#SBATCH -p icelake-himem
#SBATCH -A SIMIDJIEVSKI-SL3-CPU
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --job-name=LRAcluster_feat_merge
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail

METHOD="LRAcluster"
SCRIPT="Jobs/LRAcluster/LRAcluster_feature_experiments.R"
RUN_ROOT="${RUN_ROOT:-Results/Performance/Feature_perturbations}"

module purge
module load rhel8/cclake/base
module load r/4.4.0/gcc/3z5n2tph

export R_LIBS_USER="$HOME/MO_survey/R/library"
mkdir -p "$R_LIBS_USER" logs

: "${ARRAY_JOB_ID:?Set ARRAY_JOB_ID to the SLURM_ARRAY_JOB_ID of the feature array job}"

RUN_DIR="${RUN_ROOT}/${METHOD}/array_${ARRAY_JOB_ID}"
OUT_DIR="${RUN_DIR}/out"
mkdir -p "${OUT_DIR}"

export BENCH_RUN_DIR="${RUN_DIR}"
export BENCH_OUT_DIR="${OUT_DIR}"
export LRA_MODE="merge"

cd "${SLURM_SUBMIT_DIR}"
Rscript "${SCRIPT}"

