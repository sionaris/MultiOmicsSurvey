#!/bin/bash
#SBATCH --job-name=iCB_feat
#SBATCH --output=logs/iCB_feat_%a.out
#SBATCH --error=logs/iCB_feat_%a.err
#SBATCH --account=SIMIDJIEVSKI-SL3-CPU
#SBATCH --partition=icelake-himem
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=11:59:59
#SBATCH --array=1-5%1

set -euo pipefail

# SLURM bugfix for mismatched TRES/cpus-per-task
export SLURM_TRES_PER_TASK="cpu=${SLURM_CPUS_PER_TASK}"

module purge
module load rhel8/cclake/base
module load r/4.4.0/gcc/3z5n2tph

# user library (avoid re-installing into repo)
export R_LIBS_USER="$HOME/MO_survey/R/library"
mkdir -p "$R_LIBS_USER"

# benchmarking fairness
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

METHOD="iClusterBayes"
RUN_ROOT="Results/Performance/Feature_perturbations/${METHOD}"
RUN_DIR="${RUN_ROOT}/array_${SLURM_ARRAY_JOB_ID}"
TASK_DIR="${RUN_DIR}/task_${SLURM_ARRAY_TASK_ID}"
OUT_DIR="${TASK_DIR}/out"
META_DIR="${TASK_DIR}/meta"

mkdir -p "${OUT_DIR}" "${META_DIR}" logs

env | sort > "${META_DIR}/env.txt"
module list 2>&1 > "${META_DIR}/modules.txt"
{
  echo "SLURM_ARRAY_JOB_ID=${SLURM_ARRAY_JOB_ID}"
  echo "SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"
  echo "SLURM_JOB_ID=${SLURM_JOB_ID}"
  echo "SLURM_JOB_NODELIST=${SLURM_JOB_NODELIST}"
} > "${META_DIR}/slurm_ids.txt"

export ICB_MODE="worker"
export ICB_CENTILE_INDEX="${SLURM_ARRAY_TASK_ID}"
export ICB_CPUS="1"
export BENCH_OUT_DIR="${OUT_DIR}"
export BENCH_META_DIR="${META_DIR}"
export BENCH_RUN_DIR="${RUN_DIR}"

echo "[INFO] Running ${METHOD} feature worker centile index ${ICB_CENTILE_INDEX}"

# capture time even if R fails (and keep logs separate)
TIME_FILE="${META_DIR}/time.txt"
set +e
/usr/bin/time -v -o "${TIME_FILE}" \
  srun --hint=nomultithread --cpu-bind=cores \
  Rscript Jobs/${METHOD}/${METHOD}_feature_experiments.R
RC=$?
set -e
exit "${RC}"

