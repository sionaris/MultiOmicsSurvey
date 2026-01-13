#!/usr/bin/env bash
#SBATCH -p icelake-himem
#SBATCH -A SIMIDJIEVSKI-SL3-CPU
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=15000M
#SBATCH --time=12:00:00
#SBATCH --job-name=LRAcluster_feat
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --exclusive
#SBATCH --array=1-5%1

set -euo pipefail

METHOD="LRAcluster"
SCRIPT="Jobs/LRAcluster/LRAcluster_feature_experiments.R"
RUN_ROOT="${RUN_ROOT:-Results/Performance/Feature_perturbations}"

module purge
module load rhel8/cclake/base
module load r/4.4.0/gcc/3z5n2tph

export R_LIBS_USER="$HOME/MO_survey/R/library"
mkdir -p "$R_LIBS_USER" logs

# Normalize SLURM_CPUS_PER_TASK if SLURM_TRES_PER_TASK is set (prevents mismatch errors)
if [[ -n "${SLURM_TRES_PER_TASK:-}" ]]; then
  CPUS_FROM_TRES="$(echo "${SLURM_TRES_PER_TASK}" | sed -n 's/.*cpu:\([0-9]\+\).*/\1/p')"
  if [[ -n "${CPUS_FROM_TRES}" ]]; then
    export SLURM_CPUS_PER_TASK="${CPUS_FROM_TRES}"
  fi
fi

THREADS="${SLURM_CPUS_PER_TASK:-1}"

# Fair benchmarking: keep BLAS threads at 1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# But allow M3C to use available cores (not benchmarked)
export M3C_CORES="${THREADS}"

TS="$(date -u +%Y%m%dT%H%M%SZ)"
JOBID="${SLURM_ARRAY_JOB_ID:-nojobid}"
TASKID="${SLURM_ARRAY_TASK_ID:-notask}"

RUN_DIR="${RUN_ROOT}/${METHOD}/array_${JOBID}"
TASK_DIR="${RUN_DIR}/task_${TASKID}"
META_DIR="${TASK_DIR}/meta"
OUT_DIR="${TASK_DIR}/out"
mkdir -p "${META_DIR}" "${OUT_DIR}"

export BENCH_RUN_DIR="${RUN_DIR}"
export BENCH_OUT_DIR="${OUT_DIR}"
export LRA_MODE="worker"

cat > "${META_DIR}/ids.tsv" <<EOF
job_id  ${JOBID}
task_id ${TASKID}
timestamp_utc   ${TS}
method  ${METHOD}
partition       ${SLURM_JOB_PARTITION:-}
nodelist        ${SLURM_JOB_NODELIST:-}
cpus_per_task   ${THREADS}
mem_per_cpu     ${SLURM_MEM_PER_CPU:-}
mem_per_node    ${SLURM_MEM_PER_NODE:-}
EOF

{
  date -u
  echo "hostname: $(hostname)"
  echo "submit_dir: ${SLURM_SUBMIT_DIR:-}"
  echo "pwd_before_cd: $(pwd)"
  echo
  module -t list 2>&1 || true
  echo
  env | egrep '^(BENCH_RUN_DIR|BENCH_OUT_DIR|R_LIBS_USER|OMP_NUM_THREADS|MKL_NUM_THREADS|OPENBLAS_NUM_THREADS|NUMEXPR_NUM_THREADS|M3C_CORES)=' || true
} > "${META_DIR}/environment.txt"

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

