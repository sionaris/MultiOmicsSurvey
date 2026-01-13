#!/usr/bin/env bash
#SBATCH -p icelake-himem
#SBATCH -A SIMIDJIEVSKI-SL3-CPU
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=15000M
#SBATCH --time=12:00:00
#SBATCH --job-name=LRAcluster_samp
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --exclusive
#SBATCH --array=1-50%5

set -euo pipefail

METHOD="LRAcluster"
SCRIPT="Jobs/LRAcluster/LRAcluster_sample_experiments.R"
RUN_ROOT="${RUN_ROOT:-Results/Performance/Sample_perturbations}"

module purge
module load rhel8/cclake/base
module load r/4.4.0/gcc/3z5n2tph

export R_LIBS_USER="$HOME/MO_survey/R/library"
mkdir -p "$R_LIBS_USER" logs

# avoid SLURM env mismatch
if [[ -n "${SLURM_TRES_PER_TASK:-}" ]]; then
  CPUS_FROM_TRES="$(echo "${SLURM_TRES_PER_TASK}" | sed -n 's/.*cpu:\([0-9]\+\).*/\1/p')"
  if [[ -n "${CPUS_FROM_TRES}" ]]; then
    export SLURM_CPUS_PER_TASK="${CPUS_FROM_TRES}"
  fi
fi

# fairness/stability: keep BLAS threads at 1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export PYTHONHASHSEED=0

# keep M3C at 1 core (not benchmarked, reduces memory spikes)
export M3C_CORES=1

JOBID="${SLURM_ARRAY_JOB_ID}"
TASKID="${SLURM_ARRAY_TASK_ID}"

RUN_DIR="${RUN_ROOT}/${METHOD}/array_${JOBID}"
TASK_DIR="${RUN_DIR}/task_${TASKID}"
META_DIR="${TASK_DIR}/meta"
OUT_DIR="${TASK_DIR}/out"
mkdir -p "${META_DIR}" "${OUT_DIR}"

export BENCH_RUN_DIR="${RUN_DIR}"
export BENCH_OUT_DIR="${OUT_DIR}"
export LRA_MODE="worker"

cd "${SLURM_SUBMIT_DIR}"

TIME_FILE="${META_DIR}/usrbin_time_v.txt"
set +e
/usr/bin/time -v -o "${TIME_FILE}" \
  srun --hint=nomultithread --cpu-bind=cores --export=ALL \
  Rscript "${SCRIPT}"
RC=$?
set -e

echo "${RC}" > "${META_DIR}/exit_code.txt"
exit "${RC}"

