#!/bin/bash
#SBATCH --partition=icelake-himem
#SBATCH -A SIMIDJIEVSKI-SL3-CPU
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=200G
#SBATCH --time=11:59:59
#SBATCH --array=1-5%2
#SBATCH --job-name=MONET_feat
#SBATCH --output=logs/MONET_feat_%A_%a.out
#SBATCH --error=logs/MONET_feat_%A_%a.err

set -euo pipefail

export SRUN_CPUS_PER_TASK="${SLURM_CPUS_PER_TASK}"
export SLURM_TRES_PER_TASK="cpu:${SLURM_CPUS_PER_TASK}"

METHOD="MONET"
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
module load rhel8/cclake/base
module load r/4.4.0/gcc/3z5n2tph

export R_LIBS_USER="${SLURM_SUBMIT_DIR}/R/library"
mkdir -p "$R_LIBS_USER"

CONDA="${HOME}/miniconda/bin/conda"
MONET_ENV_NAME="${MONET_ENV_NAME:-MONET}"

MONET_ENV_PREFIX="$("${CONDA}" env list | awk -v e="${MONET_ENV_NAME}" '$1==e {print $NF; exit}')"
if [[ -z "${MONET_ENV_PREFIX}" || ! -x "${MONET_ENV_PREFIX}/bin/python" ]]; then
  echo "ERROR: conda env '${MONET_ENV_NAME}' not found (or python missing)." >&2
  "${CONDA}" env list >&2
  exit 2
fi
MONET_PY="${MONET_ENV_PREFIX}/bin/python"

# Feature centiles by array index
CENTILES=(0.10 0.20 0.50 0.75 0.90)
MONET_CENTILE="${CENTILES[$((SLURM_ARRAY_TASK_ID-1))]}"
PCT="$(awk -v c="${MONET_CENTILE}" 'BEGIN{printf("%d", c*100+0.5)}')"

CORR_DIR="${OUT_DIR}/MO_Corrs"
PREP_META="${META_DIR}/MONET_prepare_meta_${PCT}pct.tsv"

MONET_OUT_DIR="${OUT_DIR}/monet_results"
mkdir -p "${MONET_OUT_DIR}"

echo "[$(date)] ${METHOD} feature task ${SLURM_ARRAY_TASK_ID} starting (centile=${MONET_CENTILE}, env=${MONET_ENV_NAME})"

# 0) Preflight: verify python + networkx on the compute node
"${MONET_PY}" -c 'import sys; import networkx as nx; print(sys.executable); print(nx.__version__)' \
  > "${META_DIR}/python_env_check.txt" 2>&1

# 1) R: make correlations (creates MO_Corrs/*.csv)
srun --hint=nomultithread --cpu-bind=cores --export=ALL \
  Rscript --vanilla "${SLURM_SUBMIT_DIR}/Jobs/MONET/MONET_make_correlations_feature.R" \
    --centile "${MONET_CENTILE}" \
    --corr_dir "${CORR_DIR}" \
    --prep_meta "${PREP_META}" \
  > "${OUT_DIR}/prepare.stdout.txt" 2> "${OUT_DIR}/prepare.stderr.txt"

# Sanity check: expect 5 CSVs
n_csv="$(ls -1 "${CORR_DIR}"/*.csv 2>/dev/null | wc -l | tr -d ' ')"
if [[ "${n_csv}" -ne 5 ]]; then
  echo "ERROR: Expected 5 correlation CSVs in ${CORR_DIR}, found ${n_csv}" >&2
  ls -lah "${CORR_DIR}" >&2 || true
  exit 3
fi

# 2) Python: MONET mainloop (BENCHMARKED)
/usr/bin/time -v -o "${META_DIR}/timev_monet.txt" \
  srun --hint=nomultithread --cpu-bind=cores --export=ALL \
  "${MONET_PY}" "${SLURM_SUBMIT_DIR}/Jobs/MONET/run_MONET_mainloop.py" \
    --corr_dir "${CORR_DIR}" \
    --out_dir "${MONET_OUT_DIR}" \
    --seed 123 \
  > "${OUT_DIR}/monet.stdout.txt" 2> "${OUT_DIR}/monet.stderr.txt"

# 3) R: postprocess into benchmark-friendly outputs
CLUST_TSV="${MONET_OUT_DIR}/monet_clustering.tsv"
AVG_ADJ_GZ="${MONET_OUT_DIR}/monet_avg_adjacency.tsv.gz"
METRICS_JSON="${MONET_OUT_DIR}/monet_metrics.json"

OUT_CLUSTERS="${OUT_DIR}/monet_clusters.tsv.gz"
OUT_PERF_ROW="${OUT_DIR}/monet_perf_row.tsv"

srun --hint=nomultithread --cpu-bind=cores --export=ALL \
  Rscript --vanilla "${SLURM_SUBMIT_DIR}/Jobs/MONET/MONET_postprocess.R" \
    --mode feature \
    --prep_meta "${PREP_META}" \
    --metrics_json "${METRICS_JSON}" \
    --timev_file "${META_DIR}/timev_monet.txt" \
    --clustering_tsv "${CLUST_TSV}" \
    --avg_adj_tsv_gz "${AVG_ADJ_GZ}" \
    --out_clusters "${OUT_CLUSTERS}" \
    --out_perf_row "${OUT_PERF_ROW}" \
  > "${OUT_DIR}/post.stdout.txt" 2> "${OUT_DIR}/post.stderr.txt"

echo "[$(date)] ${METHOD} feature task ${SLURM_ARRAY_TASK_ID} done"

