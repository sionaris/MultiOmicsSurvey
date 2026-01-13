#!/usr/bin/env bash
#SBATCH -p icelake-himem
#SBATCH -A SIMIDJIEVSKI-SL3-CPU
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:45:00
#SBATCH --job-name=CIMLR_samp_merge
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

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

ROOT="${SLURM_SUBMIT_DIR}"

module purge
module load rhel8/cclake/base
module load r/4.4.0/gcc/3z5n2tph

export R_LIBS_USER="${ROOT}/R/library"
mkdir -p "${R_LIBS_USER}" "${ROOT}/logs"

: "${ARRAY_JOB_ID:?Set ARRAY_JOB_ID to the SLURM_ARRAY_JOB_ID of the sample array job}"

RUN_ROOT="${RUN_ROOT:-Results/Performance/Sample_perturbations}"
[[ "${RUN_ROOT}" = /* ]] || RUN_ROOT="${ROOT}/${RUN_ROOT}"

RUN_DIR="${RUN_ROOT}/CIMLR/array_${ARRAY_JOB_ID}"
OUT_DIR="${RUN_DIR}/out"
mkdir -p "${OUT_DIR}"

export BENCH_RUN_DIR="${RUN_DIR}"
export BENCH_OUT_DIR="${OUT_DIR}"
export CIMLR_MODE="merge"

cd "${ROOT}"

META_DIR="${RUN_DIR}/meta/merge_${SLURM_JOB_ID}"
mkdir -p "${META_DIR}"

set +e
/usr/bin/time -v -o "${META_DIR}/timev.txt" \
  srun --hint=nomultithread --cpu-bind=cores --export=ALL \
  Rscript --vanilla "Jobs/CIMLR/CIMLR_sample_experiments.R"
rc=$?
set -e

echo "${rc}" > "${META_DIR}/exit_code.txt"
exit "${rc}"

