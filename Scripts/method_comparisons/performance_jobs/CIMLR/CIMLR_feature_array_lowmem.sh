#!/usr/bin/env bash
#SBATCH -p icelake-himem
#SBATCH -A SIMIDJIEVSKI-SL3-CPU
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=300G
#SBATCH --time=06:00:00
#SBATCH --job-name=CIMLR_feat_lo
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --exclusive
#SBATCH --array=1-2%1

set -euo pipefail

# ---- Slurm env normalisation (prevents srun fatal mismatch) ----
if [[ -n "${SLURM_TRES_PER_TASK:-}" ]]; then
  _cpu="$(echo "${SLURM_TRES_PER_TASK}" | sed -n 's/.*cpu[:=]\([0-9]\+\).*/\1/p')"
  if [[ -n "${_cpu}" ]]; then
    export SLURM_CPUS_PER_TASK="${_cpu}"
  fi
fi
export SLURM_TRES_PER_TASK="cpu:${SLURM_CPUS_PER_TASK}"
# --------------------------------------------------------------

# Reduce hidden BLAS threading (stability + avoids oversubscription)
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

ROOT="${SLURM_SUBMIT_DIR}"

module purge
module load rhel8/cclake/base
module load r/4.4.0/gcc/3z5n2tph

export R_LIBS_USER="${ROOT}/R/library"
mkdir -p "${R_LIBS_USER}" "${ROOT}/logs"

RUN_ROOT="${RUN_ROOT:-Results/Performance/Feature_perturbations}"
[[ "${RUN_ROOT}" = /* ]] || RUN_ROOT="${ROOT}/${RUN_ROOT}"

# Shared run id so multiple arrays can contribute to the same folder
RUN_ID="${RUN_ID:-${SLURM_ARRAY_JOB_ID}}"

RUN_DIR="${RUN_ROOT}/CIMLR/array_${RUN_ID}"
TASK_DIR="${RUN_DIR}/task_${SLURM_ARRAY_TASK_ID}"
OUT_DIR="${TASK_DIR}/out"
META_DIR="${TASK_DIR}/meta"
mkdir -p "${OUT_DIR}" "${META_DIR}"

export BENCH_RUN_DIR="${RUN_DIR}"
export BENCH_OUT_DIR="${OUT_DIR}"
export BENCH_META_DIR="${META_DIR}"
export CIMLR_MODE="worker"

cd "${ROOT}"

TASK_META_DIR="${BENCH_RUN_DIR}/meta/task_${SLURM_ARRAY_TASK_ID}"
mkdir -p "${TASK_META_DIR}"

set +e
/usr/bin/time -v -o "${TASK_META_DIR}/timev.txt" \
  srun --hint=nomultithread --cpu-bind=cores --export=ALL \
  Rscript --vanilla "Jobs/CIMLR/CIMLR_feature_experiments.R"
rc=$?
set -e

echo "${rc}" > "${TASK_META_DIR}/exit_code.txt"
exit "${rc}"

