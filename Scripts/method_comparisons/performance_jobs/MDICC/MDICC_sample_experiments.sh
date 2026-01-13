#!/usr/bin/env bash
#SBATCH -p icelake-himem
#SBATCH -A SIMIDJIEVSKI-SL3-CPU
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20000M
#SBATCH --time=12:00:00
#SBATCH --job-name=MDICC_samp
#SBATCH --output=logs/MDICC_samp.out
#SBATCH --error=logs/MDICC_samp.err
#SBATCH --exclusive

set -euo pipefail

METHOD="${METHOD:-MDICC}"
RUN_ROOT="Results/Performance/Sample_perturbations"

module purge
module load rhel8/cclake/base
module load r/4.4.0/gcc/3z5n2tph

export R_LIBS_USER="${R_LIBS_USER:-$HOME/MO_survey/R/library}"
mkdir -p "$R_LIBS_USER" logs

echo "Current R library paths:"
Rscript -e "print(.libPaths())"

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

# --- conda / python for reticulate ---
unset RETICULATE_PYTHON || true
unset RETICULATE_PYTHON_FALLBACK || true

CONDA_SH=""
for CAND in \
  "$HOME/miniconda/etc/profile.d/conda.sh" \
  "$HOME/miniconda3/etc/profile.d/conda.sh" \
  "$HOME/anaconda3/etc/profile.d/conda.sh" \
  "$HOME/mambaforge/etc/profile.d/conda.sh"
do
  if [[ -f "$CAND" ]]; then
    CONDA_SH="$CAND"
    break
  fi
done

if [[ -z "$CONDA_SH" ]]; then
  echo "ERROR: Could not find conda.sh under $HOME (tried miniconda/miniconda3/anaconda3/mambaforge)."
  exit 2
fi

# shellcheck disable=SC1090
source "$CONDA_SH"

conda activate mdicc_py
if [[ "${CONDA_DEFAULT_ENV:-}" != "mdicc_py" ]]; then
  echo "ERROR: conda activate mdicc_py failed (CONDA_DEFAULT_ENV=${CONDA_DEFAULT_ENV:-})."
  exit 2
fi

PY="${CONDA_PREFIX}/bin/python"
if [[ ! -x "$PY" ]]; then
  echo "ERROR: expected python at $PY but it is not executable."
  exit 2
fi

export PATH="${CONDA_PREFIX}/bin:${PATH}"
export RETICULATE_PYTHON="$PY"
export MDICC_PYTHON="$PY"

echo "RETICULATE_PYTHON=$RETICULATE_PYTHON"
"$PY" --version

"$PY" - <<'PY'
import numpy, pandas, scipy, sklearn
print("Python imports OK:", numpy.__version__, pandas.__version__, scipy.__version__, sklearn.__version__)
PY

TS="$(date -u +%Y%m%dT%H%M%SZ)"
JOBID="${SLURM_JOB_ID:-nojobid}"
RUN_DIR="${RUN_ROOT}/${METHOD}/job_${JOBID}_${TS}"
META_DIR="${RUN_DIR}/meta"
OUT_DIR="${RUN_DIR}/out"
mkdir -p "${META_DIR}" "${OUT_DIR}"

export BENCH_RUN_DIR="${RUN_DIR}"
export BENCH_META_DIR="${META_DIR}"
export BENCH_OUT_DIR="${OUT_DIR}"

{
  echo "UTC_TS=${TS}"
  echo "JOBID=${JOBID}"
  echo "METHOD=${METHOD}"
  echo "HOSTNAME=$(hostname)"
  echo "PWD=$(pwd)"
  echo "RETICULATE_PYTHON=${RETICULATE_PYTHON}"
} > "${META_DIR}/run_meta.txt"

env | sort > "${META_DIR}/env.txt"
module list 2>&1 | cat > "${META_DIR}/modules.txt" || true
Rscript -e "writeLines(capture.output(sessionInfo()), '${META_DIR}/sessionInfo.txt')" || true

echo "Running sample perturbations for ${METHOD}..."
/usr/bin/time -v srun --hint=nomultithread --cpu-bind=cores \
  Rscript Jobs/MDICC/MDICC_sample_experiments.R \
  > "${OUT_DIR}/MDICC_sample_stdout.txt" 2> "${OUT_DIR}/MDICC_sample_stderr.txt"

echo "Done. Output in ${RUN_DIR}"

