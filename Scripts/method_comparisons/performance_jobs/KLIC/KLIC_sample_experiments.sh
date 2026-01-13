#!/usr/bin/env bash
#SBATCH -p icelake-himem
#SBATCH -A SIMIDJIEVSKI-SL3-CPU
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=15000M
#SBATCH --time=12:00:00
#SBATCH --job-name=KLIC_samp
#SBATCH --output=logs/KLIC_samp.out
#SBATCH --error=logs/KLIC_samp.err
#SBATCH --exclusive

set -euo pipefail

METHOD="${METHOD:-KLIC}"
SCRIPT="${SCRIPT:-Jobs/KLIC/KLIC_sample_experiments.R}"
RUN_ROOT="${RUN_ROOT:-Results/Performance/Sample_perturbations}"

module purge
module load rhel8/cclake/base
module load r/4.4.0/gcc/3z5n2tph

export R_LIBS_USER="${R_LIBS_USER:-$HOME/MO_survey/R/library}"
mkdir -p "$R_LIBS_USER" logs

# SLURM mismatch fix
if [[ -n "${SLURM_TRES_PER_TASK:-}" ]]; then
  CPUS_FROM_TRES="$(echo "${SLURM_TRES_PER_TASK}" | sed -n 's/.*cpu:\([0-9]\+\).*/\1/p')"
  if [[ -n "${CPUS_FROM_TRES}" ]]; then
    export SLURM_CPUS_PER_TASK="${CPUS_FROM_TRES}"
  fi
fi

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export PYTHONHASHSEED=0

export MOSEK_PLATFORM_DIR="$HOME/mosek/10.2/tools/platform/linux64x86"
export MOSEK_BINDIR="$MOSEK_PLATFORM_DIR/bin"
export PATH="$MOSEK_BINDIR:$PATH"

LD_ADD=""
[[ -d "$MOSEK_PLATFORM_DIR/bin" ]] && LD_ADD="$MOSEK_PLATFORM_DIR/bin"
[[ -d "$MOSEK_PLATFORM_DIR/lib" ]] && LD_ADD="${LD_ADD:+$LD_ADD:}$MOSEK_PLATFORM_DIR/lib"
export LD_LIBRARY_PATH="${LD_ADD}${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

export MOSEKLM_LICENSE_FILE="$SLURM_SUBMIT_DIR/mosek.lic"
export MOSEK_LICENSE_FILE="$MOSEKLM_LICENSE_FILE"

TS="$(date -u +%Y%m%dT%H%M%SZ)"
JOBID="${SLURM_JOB_ID:-nojobid}"
RUN_DIR="${RUN_ROOT}/${METHOD}/job_${JOBID}_${TS}"
META_DIR="${RUN_DIR}/meta"
OUT_DIR="${RUN_DIR}/out"
mkdir -p "${META_DIR}" "${OUT_DIR}"

export BENCH_RUN_DIR="${RUN_DIR}"
export BENCH_OUT_DIR="${OUT_DIR}"

test -d "$MOSEK_PLATFORM_DIR" || { echo "Missing MOSEK_PLATFORM_DIR=$MOSEK_PLATFORM_DIR"; exit 2; }
test -x "$MOSEK_BINDIR/mosek" || { echo "mosek binary not found in $MOSEK_BINDIR"; exit 2; }
test -f "$MOSEKLM_LICENSE_FILE" || { echo "Missing license: $MOSEKLM_LICENSE_FILE"; exit 2; }

if [[ -x "$MOSEK_BINDIR/msktestlic" ]]; then
  set +e
  "$MOSEK_BINDIR/msktestlic" > "${META_DIR}/msktestlic.txt" 2>&1
  RC_LIC=$?
  set -e
  if [[ $RC_LIC -ne 0 ]]; then
    echo "msktestlic FAILED (rc=$RC_LIC). See ${META_DIR}/msktestlic.txt"
    cat "${META_DIR}/msktestlic.txt" >&2
    exit 3
  fi
else
  echo "msktestlic not found at $MOSEK_BINDIR/msktestlic (continuing without it)" > "${META_DIR}/msktestlic.txt"
fi

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

{
  date -u
  echo "hostname: $(hostname)"
  echo "submit_dir: ${SLURM_SUBMIT_DIR:-}"
  echo
  env | egrep '^(MOSEK_PLATFORM_DIR|MOSEK_BINDIR|MOSEKLM_LICENSE_FILE|MOSEK_LICENSE_FILE|LD_LIBRARY_PATH|PATH|OMP_NUM_THREADS|MKL_NUM_THREADS|OPENBLAS_NUM_THREADS|NUMEXPR_NUM_THREADS|PYTHONHASHSEED)=' || true
  echo
  module -t list 2>&1 || true
} > "${META_DIR}/environment.txt"

{
  command -v R && R --version || true
  command -v Rscript && Rscript --version || true
  command -v mosek && mosek -v || true
} > "${META_DIR}/versions.txt"

cd "${SLURM_SUBMIT_DIR}"

TIME_FILE="${META_DIR}/usrbin_time_v.txt"
set +e
/usr/bin/time -v -o "${TIME_FILE}" \
  srun --hint=nomultithread --cpu-bind=cores \
  Rscript "${SCRIPT}"
RC=$?
set -e

echo "${RC}" > "${META_DIR}/exit_code.txt"
exit "${RC}"

