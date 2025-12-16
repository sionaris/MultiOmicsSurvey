#!/usr/bin/env bash
#SBATCH -p icelake-himem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=6500M
#SBATCH --time=06:00:00
#SBATCH --job-name=bench
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --exclusive

set -euo pipefail

###############################################################################
# USER-SET VARIABLES (set via sbatch --export=ALL,KEY=VALUE,... or edit here)
#
# Required:
#   METHOD      e.g. "MethodA"
#   PROPERTY    e.g. "robust_feat" | "robust_samp" | "seed_stab" | "scalability_n" | ...
#   POINT       e.g. "feat_10pct" or "samp_70pct_rep03" or "seed_7"
#   CMD         full command to run (quoted string), e.g.
#              CMD="python run_method.py --config cfg.yaml --seed 1 --out outdir"
#
# Optional:
#   RUN_ROOT    base output directory (default: results)
###############################################################################
METHOD="${METHOD:-UNSET_METHOD}"
PROPERTY="${PROPERTY:-UNSET_PROPERTY}"
POINT="${POINT:-UNSET_POINT}"
CMD="${CMD:-}"
RUN_ROOT="${RUN_ROOT:-results}"

if [[ -z "${CMD}" ]]; then
  echo "ERROR: CMD is empty. Provide CMD via: sbatch --export=ALL,METHOD=...,PROPERTY=...,POINT=...,CMD=\"...\" $0"
  exit 2
fi

###############################################################################
# CONSISTENT THREADING CONTROLS
###############################################################################
THREADS="${SLURM_CPUS_PER_TASK:-1}"

export OMP_NUM_THREADS="${THREADS}"
export MKL_NUM_THREADS="${THREADS}"
export OPENBLAS_NUM_THREADS="${THREADS}"
export NUMEXPR_NUM_THREADS="${THREADS}"
export PYTHONHASHSEED=0

###############################################################################
# OUTPUT STRUCTURE (one directory per job)
###############################################################################
TS="$(date -u +%Y%m%dT%H%M%SZ)"
JOBID="${SLURM_JOB_ID:-nojobid}"

RUN_DIR="${RUN_ROOT}/${METHOD}/${PROPERTY}/${POINT}/job_${JOBID}_${TS}"
META_DIR="${RUN_DIR}/meta"
mkdir -p "${META_DIR}" logs

###############################################################################
# RECORD IDENTIFIERS FOR LATER sacct HARVESTING
###############################################################################
cat > "${META_DIR}/ids.tsv" <<EOF
job_id	${JOBID}
timestamp_utc	${TS}
method	${METHOD}
property	${PROPERTY}
point	${POINT}
partition	${SLURM_JOB_PARTITION:-}
nodelist	${SLURM_JOB_NODELIST:-}
cpus_per_task	${SLURM_CPUS_PER_TASK:-}
mem_per_cpu	${SLURM_MEM_PER_CPU:-}
mem_per_node	${SLURM_MEM_PER_NODE:-}
EOF

echo "${CMD}" > "${META_DIR}/cmd.txt"

###############################################################################
# RECORD ENVIRONMENT / HARDWARE (for cross-job comparability)
###############################################################################
{
  echo "=== BASIC ==="
  date -u
  echo "hostname: $(hostname)"
  echo "pwd: $(pwd)"
  echo "submit_dir: ${SLURM_SUBMIT_DIR:-}"
  echo

  echo "=== SLURM JOB INFO (scontrol) ==="
  scontrol show job "${JOBID}" || true
  echo

  echo "=== NODE INFO ==="
  echo "SLURM_JOB_NODELIST=${SLURM_JOB_NODELIST:-}"
  lscpu || true
  echo
  free -h || true
  echo

  echo "=== THREAD ENV ==="
  env | egrep '^(OMP_NUM_THREADS|MKL_NUM_THREADS|OPENBLAS_NUM_THREADS|NUMEXPR_NUM_THREADS|PYTHONHASHSEED)=' || true
  echo

  echo "=== MODULES ==="
  module -t list 2>&1 || true
  echo

  echo "=== ULIMIT ==="
  ulimit -a || true
  echo
} > "${META_DIR}/environment.txt"

{
  echo "=== python ==="
  command -v python && python --version || true
  command -v python && python -c "import sys; print(sys.executable)" || true
  echo
  echo "=== R ==="
  command -v R && R --version || true
} > "${META_DIR}/versions.txt"

###############################################################################
# RUN WITH srun + CORE BINDING
###############################################################################
echo "RUN_DIR=${RUN_DIR}"
echo "METHOD=${METHOD} PROPERTY=${PROPERTY} POINT=${POINT}"
echo "THREADS=${THREADS}"
echo "CMD=${CMD}"

TIME_FILE="${META_DIR}/usrbin_time_v.txt"

set +e
/usr/bin/time -v -o "${TIME_FILE}" \
  srun --hint=nomultithread --cpu-bind=cores --export=ALL \
  bash -lc "${CMD}"
RC=$?
set -e

echo "${RC}" > "${META_DIR}/exit_code.txt"
exit "${RC}"
