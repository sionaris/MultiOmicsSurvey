#!/usr/bin/env bash
#SBATCH -p icelake-himem
#SBATCH -A SIMIDJIEVSKI-SL3-CPU
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=60000M
#SBATCH --time=12:00:00
#SBATCH --job-name=LRAcluster_feat
#SBATCH --output=logs/LRAcluster_feat.out
#SBATCH --error=logs/LRAcluster_feat.err
#SBATCH --exclusive

set -euo pipefail

METHOD="${METHOD:-LRAcluster}"
SCRIPT="${SCRIPT:-Jobs/LRAcluster/LRAcluster_feature_experiments.R}"
RUN_ROOT="${RUN_ROOT:-Results/Performance/Feature_perturbations}"

module purge
module load rhel8/cclake/base
module load r/4.4.0/gcc/3z5n2tph

export R_LIBS_USER="$HOME/MO_survey/R/library"
mkdir -p "$R_LIBS_USER" logs

echo "Current R library paths:"
Rscript -e "print(.libPaths())"

# Avoid SLURM env mismatches (prefer TRES if present)
if [[ -n "${SLURM_TRES_PER_TASK:-}" ]]; then
  CPUS_FROM_TRES="$(echo "${SLURM_TRES_PER_TASK}" | sed -n 's/.*cpu:\([0-9]\+\).*/\1/p')"
  if [[ -n "${CPUS_FROM_TRES}" ]]; then
    export SLURM_CPUS_PER_TASK="${CPUS_FROM_TRES}"
  fi
fi

THREADS="${SLURM_CPUS_PER_TASK:-1}"
export OMP_NUM_THREADS="${THREADS}"
export MKL_NUM_THREADS="${THREADS}"
export OPENBLAS_NUM_THREADS="${THREADS}"
export NUMEXPR_NUM_THREADS="${THREADS}"
export PYTHONHASHSEED=0

TS="$(date -u +%Y%m%dT%H%M%SZ)"
JOBID="${SLURM_JOB_ID:-nojobid}"
RUN_DIR="${RUN_ROOT}/${METHOD}/job_${JOBID}_${TS}"
META_DIR="${RUN_DIR}/meta"
OUT_DIR="${RUN_DIR}/out"
mkdir -p "${META_DIR}" "${OUT_DIR}"

export BENCH_RUN_DIR="${RUN_DIR}"
export BENCH_OUT_DIR="${OUT_DIR}"

cat > "${META_DIR}/ids.tsv" <<EOF
job_id  ${JOBID}
timestamp_utc   ${TS}
method  ${METHOD}
partition       ${SLURM_JOB_PARTITION:-}
nodelist        ${SLURM_JOB_NODELIST:-}
cpus_per_task   ${SLURM_CPUS_PER_TASK:-}
mem_per_cpu     ${SLURM_MEM_PER_CPU:-}
mem_per_node    ${SLURM_MEM_PER_NODE:-}
EOF

echo "${SCRIPT}" > "${META_DIR}/script.txt"

{
  date -u
  echo "hostname: $(hostname)"
  echo "submit_dir: ${SLURM_SUBMIT_DIR:-}"
  echo "pwd_before_cd: $(pwd)"
  echo
  scontrol show job "${JOBID}" || true
  echo
  lscpu || true
  echo
  free -h || true
  echo
  env | egrep '^(BENCH_RUN_DIR|BENCH_OUT_DIR|R_LIBS_USER|OMP_NUM_THREADS|MKL_NUM_THREADS|OPENBLAS_NUM_THREADS|NUMEXPR_NUM_THREADS|PYTHONHASHSEED)=' || true
  echo
  module -t list 2>&1 || true
} > "${META_DIR}/environment.txt"

{
  command -v R && R --version || true
} > "${META_DIR}/versions.txt"

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

