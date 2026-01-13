#!/usr/bin/env bash
#SBATCH -p icelake-himem
#SBATCH -A SIMIDJIEVSKI-SL3-CPU
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=15000M
#SBATCH --time=11:59:00
#SBATCH --job-name=RWR-NF_samp
#SBATCH --output=logs/RWR-NF_samp_%A_%a.out
#SBATCH --error=logs/RWR-NF_samp_%A_%a.err
#SBATCH --exclusive
#SBATCH --array=1-50%5

set -euo pipefail

export SRUN_CPUS_PER_TASK="${SLURM_CPUS_PER_TASK:-1}"
export SLURM_TRES_PER_TASK="cpu:${SLURM_CPUS_PER_TASK:-1}"

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

ROOT="${SLURM_SUBMIT_DIR}"

module purge
module load rhel8/cclake/base
module load r/4.4.0/gcc/3z5n2tph

export R_LIBS_USER="${ROOT}/R/library"
mkdir -p "${R_LIBS_USER}" "${ROOT}/logs"

METHOD_DIR="RWR-NF"
RUN_ROOT="${RUN_ROOT:-Results/Performance/Sample_perturbations}"
[[ "${RUN_ROOT}" = /* ]] || RUN_ROOT="${ROOT}/${RUN_ROOT}"

RUN_DIR="${RUN_ROOT}/${METHOD_DIR}/array_${SLURM_ARRAY_JOB_ID}"
TASK_DIR="${RUN_DIR}/task_${SLURM_ARRAY_TASK_ID}"
OUT_DIR="${TASK_DIR}/out"
META_DIR="${TASK_DIR}/meta"
mkdir -p "${OUT_DIR}" "${META_DIR}"

OUT_DIR="$(readlink -f "${OUT_DIR}")"
META_DIR="$(readlink -f "${META_DIR}")"
RUN_DIR="$(readlink -f "${RUN_DIR}")"

export BENCH_RUN_DIR="${RUN_DIR}"
export BENCH_OUT_DIR="${OUT_DIR}"
export BENCH_META_DIR="${META_DIR}"

cd "${ROOT}"

set +e
/usr/bin/time -v -o "${META_DIR}/usrbin_time_v.txt" \
  srun --hint=nomultithread --cpu-bind=cores --export=ALL \
  Rscript "Jobs/RWR-NF/RWR-NF_sample_experiments.R" \
    --task_index "${SLURM_ARRAY_TASK_ID}"
rc=$?
set -e

echo "${rc}" > "${META_DIR}/exit_code.txt"
exit "${rc}"

