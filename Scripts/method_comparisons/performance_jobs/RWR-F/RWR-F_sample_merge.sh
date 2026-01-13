#!/usr/bin/env bash
#SBATCH -p icelake-himem
#SBATCH -A SIMIDJIEVSKI-SL3-CPU
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --job-name=RWR-F_samp_merge
#SBATCH --output=logs/RWR-F_samp_merge_%x_%j.out
#SBATCH --error=logs/RWR-F_samp_merge_%x_%j.err

set -euo pipefail

: "${ARRAY_JOB_ID:?Must pass ARRAY_JOB_ID via --export=ALL,ARRAY_JOB_ID=<id>}"

ROOT="${SLURM_SUBMIT_DIR}"

module purge
module load rhel8/cclake/base
module load r/4.4.0/gcc/3z5n2tph

export R_LIBS_USER="${ROOT}/R/library"
mkdir -p "${R_LIBS_USER}" "${ROOT}/logs"

RUN_DIR="${ROOT}/Results/Performance/Sample_perturbations/RWR-F/array_${ARRAY_JOB_ID}"
export BENCH_RUN_DIR="$(readlink -f "${RUN_DIR}")"

cd "${ROOT}"

srun --hint=nomultithread --cpu-bind=cores --export=ALL \
  Rscript "Jobs/RWR-F/RWR-F_merge.R" --mode sample

