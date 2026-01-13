#!/bin/bash
#SBATCH --partition=icelake-himem
#SBATCH -A SIMIDJIEVSKI-SL3-CPU
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=200G
#SBATCH --time=11:59:59
#SBATCH --array=1-50%5
#SBATCH --job-name=MSNE_samp
#SBATCH --output=logs/MSNE_samp_%A_%a.out
#SBATCH --error=logs/MSNE_samp_%A_%a.err

set -euo pipefail

export SRUN_CPUS_PER_TASK="${SLURM_CPUS_PER_TASK}"
export SLURM_TRES_PER_TASK="cpu:${SLURM_CPUS_PER_TASK}"

METHOD="MSNE"
RUN_ROOT="${SLURM_SUBMIT_DIR}/Results/Performance/Sample_perturbations/${METHOD}"
RUN_DIR="${RUN_ROOT}/array_${SLURM_ARRAY_JOB_ID}"
TASK_DIR="${RUN_DIR}/task_${SLURM_ARRAY_TASK_ID}"
OUT_DIR="${TASK_DIR}/out"
META_DIR="${TASK_DIR}/meta"

mkdir -p "${OUT_DIR}" "${META_DIR}" "${SLURM_SUBMIT_DIR}/logs"

OUT_DIR="$(readlink -f "${OUT_DIR}")"
META_DIR="$(readlink -f "${META_DIR}")"
RUN_DIR="$(readlink -f "${RUN_DIR}")"

export BENCH_RUN_DIR="${RUN_DIR}"
export BENCH_OUT_DIR="${OUT_DIR}"
export BENCH_META_DIR="${META_DIR}"

module purge
module load rhel8/cclake/base
module load r/4.4.0/gcc/3z5n2tph

export R_LIBS_USER="${SLURM_SUBMIT_DIR}/R/library"
mkdir -p "$R_LIBS_USER"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

export SEED=123
export PYTHONHASHSEED=123

echo "[$(date)] ${METHOD} sample task ${SLURM_ARRAY_TASK_ID} starting"

srun --hint=nomultithread --cpu-bind=cores \
  Rscript "${SLURM_SUBMIT_DIR}/Jobs/MSNE/MSNE_sample_experiments.R" prepare \
  > "${OUT_DIR}/prepare.stdout.txt" 2> "${OUT_DIR}/prepare.stderr.txt"

CONDA="${HOME}/miniconda/bin/conda"
MSNE_ENV_NAME="${MSNE_ENV_NAME:-MSNE}"

MSNE_ENV_PREFIX="$("${CONDA}" env list | awk -v e="${MSNE_ENV_NAME}" '$1==e {print $NF; exit}')"
if [[ -z "${MSNE_ENV_PREFIX}" || ! -x "${MSNE_ENV_PREFIX}/bin/python" ]]; then
  echo "ERROR: conda env '${MSNE_ENV_NAME}' not found (or python missing)." >&2
  "${CONDA}" env list >&2
  exit 2
fi
MSNE_PY="${MSNE_ENV_PREFIX}/bin/python"

MSNE_SRC_DIR="${SLURM_SUBMIT_DIR}/MSNE/MSNE"
export PYTHONPATH="${MSNE_SRC_DIR}:${PYTHONPATH:-}"

"${MSNE_PY}" -c 'import sys, numpy as np; print(sys.executable); print(np.__version__)' \
  > "${META_DIR}/python_env_check.txt" 2>&1

cd "${MSNE_SRC_DIR}"

/usr/bin/time -v -o "${META_DIR}/timev_msne.txt" \
  srun --hint=nomultithread --cpu-bind=cores --export=ALL \
  "${MSNE_PY}" "${SLURM_SUBMIT_DIR}/Jobs/MSNE/MSNE_run.py" \
    --dists_dir "${OUT_DIR}/MO_Dists" \
    --out_dir "${OUT_DIR}" \
    --seed 123 \
    --n_clusters 5 \
    --k 50 \
    --workers 1 \
    --walk_length 40 \
    --num_walks 200 \
    --embed_size 50 \
    --window_size 5 \
  > "${OUT_DIR}/msne.stdout.txt" 2> "${OUT_DIR}/msne.stderr.txt"

cd "${SLURM_SUBMIT_DIR}"

srun --hint=nomultithread --cpu-bind=cores \
  Rscript "${SLURM_SUBMIT_DIR}/Jobs/MSNE/MSNE_sample_experiments.R" post \
  > "${OUT_DIR}/post.stdout.txt" 2> "${OUT_DIR}/post.stderr.txt"

echo "[$(date)] ${METHOD} sample task ${SLURM_ARRAY_TASK_ID} done"

