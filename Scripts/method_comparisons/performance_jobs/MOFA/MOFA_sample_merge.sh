#!/bin/bash
#SBATCH --partition=icelake-himem
#SBATCH -A SIMIDJIEVSKI-SL3-CPU
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:30:00
#SBATCH --job-name=MOFA_samp_merge
#SBATCH --output=logs/MOFA_samp_merge_%j.out
#SBATCH --error=logs/MOFA_samp_merge_%j.err

set -euo pipefail
export SLURM_TRES_PER_TASK="cpu:${SLURM_CPUS_PER_TASK}"

: "${ARRAY_JOB_ID:?Must pass ARRAY_JOB_ID via --export=ALL,ARRAY_JOB_ID=<id>}"

METHOD="MOFA"
RUN_DIR="${SLURM_SUBMIT_DIR}/Results/Performance/Sample_perturbations/${METHOD}/array_${ARRAY_JOB_ID}"
OUT_DIR="${RUN_DIR}/out"
META_DIR="${RUN_DIR}/meta"
mkdir -p "${OUT_DIR}" "${META_DIR}" "${SLURM_SUBMIT_DIR}/logs"

module purge
module load rhel8/cclake/base
module load r/4.4.0/gcc/3z5n2tph

export R_LIBS_USER="${SLURM_SUBMIT_DIR}/R/library"
mkdir -p "${R_LIBS_USER}"

srun --hint=nomultithread --cpu-bind=cores --export=ALL \
  Rscript --vanilla "${SLURM_SUBMIT_DIR}/Jobs/MOFA/MOFA_merge_perf.R" \
    --run_dir "${RUN_DIR}" \
    --mode sample \
    --out_file "${OUT_DIR}/MOFA_sample_perturbations_performance.tsv" \
  > "${OUT_DIR}/merge.stdout.txt" 2> "${OUT_DIR}/merge.stderr.txt"

echo "MOFA sample merge done: ${OUT_DIR}/MOFA_sample_perturbations_performance.tsv"

