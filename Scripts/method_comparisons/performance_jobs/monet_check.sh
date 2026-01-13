# 1) Get an interactive shell on a compute node
srun --pty -p icelake-himem -A SIMIDJIEVSKI-SL3-CPU \
  --ntasks=1 --cpus-per-task=1 --mem=2G --time=00:10:00 bash -l

# 2) Inside the compute node shell, define MONET_PY exactly like your sbatch script does
CONDA="${HOME}/miniconda/bin/conda"
MONET_ENV_NAME="MONET"   # or: export MONET_ENV_NAME=MONET
MONET_ENV_PREFIX="$("${CONDA}" env list | awk -v e="${MONET_ENV_NAME}" '$1==e {print $NF; exit}')"
MONET_PY="${MONET_ENV_PREFIX}/bin/python"

# 3) Force Python to import monet from /home/as3582/MONET (your patched source tree)
export PYTHONPATH="/home/as3582/MONET:${PYTHONPATH:-}"

# 4) Run the check
"${MONET_PY}" - <<'PY'
import monet, networkx as nx
import monet.monet as mm
print("monet imported from:", monet.__file__)
print("monet.monet imported from:", mm.__file__)
print("networkx version:", nx.__version__)
print("has from_numpy_array:", hasattr(nx, "from_numpy_array"))
print("has from_numpy_matrix:", hasattr(nx, "from_numpy_matrix"))
PY

