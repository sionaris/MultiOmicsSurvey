#!/bin/bash
#SBATCH --job-name=iCB_feat_merge
#SBATCH --output=logs/iCB_feat_merge.out
#SBATCH --error=logs/iCB_feat_merge.err
#SBATCH --account=SIMIDJIEVSKI-SL3-CPU
#SBATCH --partition=icelake-himem
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00

set -euo pipefail

export SLURM_TRES_PER_TASK="cpu=${SLURM_CPUS_PER_TASK}"

module purge
module load rhel8/cclake/base
module load r/4.4.0/gcc/3z5n2tph

export R_LIBS_USER="$HOME/MO_survey/R/library"
mkdir -p "$R_LIBS_USER"

METHOD="iClusterBayes"
: "${ICB_ARRAY_JOB_ID:?You must export ICB_ARRAY_JOB_ID to this merge job}"

RUN_ROOT="Results/Performance/Feature_perturbations/${METHOD}"
RUN_DIR="${RUN_ROOT}/array_${ICB_ARRAY_JOB_ID}"
OUT_DIR="${RUN_DIR}/out"
META_DIR="${RUN_DIR}/meta"

mkdir -p "${OUT_DIR}" "${META_DIR}" logs

env | sort > "${META_DIR}/env_merge.txt"
module list 2>&1 > "${META_DIR}/modules_merge.txt"

export ICB_MODE="merge"
export BENCH_OUT_DIR="${OUT_DIR}"
export BENCH_META_DIR="${META_DIR}"
export BENCH_RUN_DIR="${RUN_DIR}"

TIME_FILE="${META_DIR}/time_merge.txt"
set +e
/usr/bin/time -v -o "${TIME_FILE}" \
  srun --hint=nomultithread --cpu-bind=cores \
  Rscript Jobs/${METHOD}/${METHOD}_feature_experiments.R
RC=$?
set -e
exit "${RC}"

