#!/bin/bash
#SBATCH --partition=icelake-himem
#SBATCH -A SIMIDJIEVSKI-SL3-CPU
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --job-name=MONET_samp_merge
#SBATCH --output=logs/MONET_samp_merge_%j.out
#SBATCH --error=logs/MONET_samp_merge_%j.err

set -euo pipefail
shopt -s nullglob

: "${ARRAY_JOB_ID:?Must pass ARRAY_JOB_ID via --export=ALL,ARRAY_JOB_ID=<id>}"

METHOD="MONET"
RUN_DIR="${SLURM_SUBMIT_DIR}/Results/Performance/Sample_perturbations/${METHOD}/array_${ARRAY_JOB_ID}"
OUT_MERGE="${RUN_DIR}/out"
META_MERGE="${RUN_DIR}/meta"

mkdir -p "${OUT_MERGE}" "${META_MERGE}" "${SLURM_SUBMIT_DIR}/logs"

echo "[$(date)] ${METHOD} sample merge starting for array_${ARRAY_JOB_ID}"
echo "RUN_DIR=${RUN_DIR}"

perf_files=( "${RUN_DIR}"/task_*/out/monet_perf_row.tsv )
if [[ ${#perf_files[@]} -eq 0 ]]; then
  echo "ERROR: No per-task monet_perf_row.tsv files found under ${RUN_DIR}/task_*/out" >&2
  exit 2
fi

# Merge perf rows (header once)
merged_perf="${OUT_MERGE}/monet_perf_rows.tsv"
awk 'FNR==1 && NR!=1 {next} {print}' "${perf_files[@]}" > "${merged_perf}"
echo "Wrote merged perf rows: ${merged_perf}"

# Copy clusters + prep meta + metrics for convenience
clust_dir="${OUT_MERGE}/monet_clusters_by_task"
prep_dir="${OUT_MERGE}/monet_prepare_meta_by_task"
metrics_dir="${OUT_MERGE}/monet_metrics_by_task"
mkdir -p "${clust_dir}" "${prep_dir}" "${metrics_dir}"

for tdir in "${RUN_DIR}"/task_*; do
  task="$(basename "${tdir}")"

  cfile="${tdir}/out/monet_clusters.tsv.gz"
  if [[ -f "${cfile}" ]]; then
    cp -f "${cfile}" "${clust_dir}/${task}_monet_clusters.tsv.gz"
  fi

  mfile="${tdir}/meta/MONET_prepare_meta_task_$(echo "${task}" | sed 's/^task_//').tsv"
  if [[ -f "${mfile}" ]]; then
    cp -f "${mfile}" "${prep_dir}/${task}_$(basename "${mfile}")"
  fi

  jfile="${tdir}/out/monet_results/monet_metrics.json"
  if [[ -f "${jfile}" ]]; then
    cp -f "${jfile}" "${metrics_dir}/${task}_monet_metrics.json"
  fi
done

echo "[$(date)] ${METHOD} sample merge done"

