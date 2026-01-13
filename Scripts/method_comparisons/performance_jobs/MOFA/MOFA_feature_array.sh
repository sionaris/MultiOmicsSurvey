#!/bin/bash
#SBATCH --partition=ampere
#SBATCH -A SIMIDJIEVSKI-SL3-GPU
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --mem=64G
#SBATCH --time=11:59:59
#SBATCH --array=1-5%1
#SBATCH --job-name=MOFA_feat
#SBATCH --output=logs/MOFA_feat_%A_%a.out
#SBATCH --error=logs/MOFA_feat_%A_%a.err

set -euo pipefail

export SRUN_CPUS_PER_TASK="${SLURM_CPUS_PER_TASK}"
export SLURM_TRES_PER_TASK="cpu:${SLURM_CPUS_PER_TASK}"

METHOD="MOFA"
RUN_ROOT="${SLURM_SUBMIT_DIR}/Results/Performance/Feature_perturbations/${METHOD}"
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
module load rhel8/default-amp
module load cuda/11.4
module load rhel8/cclake/base
module load r/4.4.0/gcc/3z5n2tph

export R_LIBS_USER="${SLURM_SUBMIT_DIR}/R/library"
mkdir -p "${R_LIBS_USER}"

source "${HOME}/miniconda3/etc/profile.d/conda.sh"
conda activate mofa_env_gpu

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

CENTILES=(0.10 0.20 0.50 0.75 0.90)
MOFA_CENTILE="${CENTILES[$((SLURM_ARRAY_TASK_ID-1))]}"
PCT="$(awk -v c="${MOFA_CENTILE}" 'BEGIN{printf("%d", c*100+0.5)}')"

MOFA_IN_DIR="${OUT_DIR}/MOFA_input"
mkdir -p "${MOFA_IN_DIR}"

PREP_META="${META_DIR}/MOFA_prepare_meta_${PCT}pct.tsv"
SAMPLE_IDS="${META_DIR}/MOFA_sample_ids_${PCT}pct.tsv"

MOFA_MODEL="${OUT_DIR}/MOFA_model_factors7_iter20000_R2_-1.hdf5"
MOFA_METRICS_JSON="${OUT_DIR}/MOFA_train_metrics.json"

echo "[$(date)] ${METHOD} feature task ${SLURM_ARRAY_TASK_ID} starting (centile=${MOFA_CENTILE})"

# 1) R: export per-view CSVs for MOFA
srun --hint=nomultithread --cpu-bind=cores --export=ALL \
  Rscript --vanilla "${SLURM_SUBMIT_DIR}/Jobs/MOFA/MOFA_export_inputs_feature.R" \
    --centile "${MOFA_CENTILE}" \
    --out_dir "${MOFA_IN_DIR}" \
    --prep_meta "${PREP_META}" \
    --sample_ids "${SAMPLE_IDS}" \
  > "${OUT_DIR}/prepare.stdout.txt" 2> "${OUT_DIR}/prepare.stderr.txt"

# Sanity check: expect 5 CSVs
n_csv="$(ls -1 "${MOFA_IN_DIR}"/*.csv 2>/dev/null | wc -l | tr -d ' ')"
if [[ "${n_csv}" -ne 5 ]]; then
  echo "ERROR: Expected 5 MOFA input CSVs in ${MOFA_IN_DIR}, found ${n_csv}" >&2
  ls -lah "${MOFA_IN_DIR}" >&2 || true
  exit 3
fi

# 2) Python GPU: train MOFA (BENCHMARKED: build+run measured inside; CPU RSS via time -v)
TIMEV="${META_DIR}/timev_mofa.txt"
/usr/bin/time -v -o "${TIMEV}" \
  srun --hint=nomultithread --cpu-bind=cores --export=ALL \
  python "${SLURM_SUBMIT_DIR}/Jobs/MOFA/run_MOFA_train_gpu.py" \
    --input_dir "${MOFA_IN_DIR}" \
    --out_hdf5 "${MOFA_MODEL}" \
    --metrics_json "${MOFA_METRICS_JSON}" \
    --seed 123 \
  > "${OUT_DIR}/mofa.stdout.txt" 2> "${OUT_DIR}/mofa.stderr.txt"

# 3) R: postprocess -> clusters + perf row
OUT_CLUSTERS="${OUT_DIR}/mofa_clusters.tsv.gz"
OUT_PERF_ROW="${OUT_DIR}/mofa_perf_row.tsv"

srun --hint=nomultithread --cpu-bind=cores --export=ALL \
  Rscript --vanilla "${SLURM_SUBMIT_DIR}/Jobs/MOFA/MOFA_postprocess.R" \
    --mode feature \
    --prep_meta "${PREP_META}" \
    --sample_ids "${SAMPLE_IDS}" \
    --mofa_hdf5 "${MOFA_MODEL}" \
    --metrics_json "${MOFA_METRICS_JSON}" \
    --timev_file "${TIMEV}" \
    --out_clusters "${OUT_CLUSTERS}" \
    --out_perf_row "${OUT_PERF_ROW}" \
  > "${OUT_DIR}/post.stdout.txt" 2> "${OUT_DIR}/post.stderr.txt"

echo "[$(date)] ${METHOD} feature task ${SLURM_ARRAY_TASK_ID} done"

