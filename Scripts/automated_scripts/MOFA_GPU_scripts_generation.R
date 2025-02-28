# generate_mofa_scripts.R

# Define the vector of latent factor values
factor_values <- c(2:10, 15, 20, 25)

# Template for the Python script (with placeholders %FACTOR% and %OUTFILE%)
py_template <- '
import glob
import os
import pandas as pd
import numpy as np
from mofapy2.run.entry_point import entry_point
import cupy as cp

print(cp.cuda.runtime.runtimeGetVersion())
print(cp.show_config())

# ------------------------------------------------------------------------------
# 1. Import CSV files, each file is treated as one \'view\' in MOFA
# ------------------------------------------------------------------------------
input_folder = "Omics_input"
csv_files = glob.glob(os.path.join(input_folder, "*.csv"))

# Read each CSV: first column is feature names and subsequent columns are samples
data_matrices = []
view_names = []
for file_path in csv_files:
    view_name = os.path.splitext(os.path.basename(file_path))[0] # give appropriate names to views
    view_names.append(view_name)

    # Read the CSV file; first column = row labels (features), rest of columns = sample measurements
    df = pd.read_csv(file_path, header=0, index_col=0)

    # Transpose so that rows = samples, columns = features
    mat = df.T.values  # shape: [number_of_samples, number_of_features]
    data_matrices.append(mat)

# Wrap in a nested list to signal that we have M views and 1 group: data_mat[m][g]
# For one group (g=0), each "view" is data_mat[m][0]
data_mat_nested = [[m] for m in data_matrices]

# ------------------------------------------------------------------------------
# 2. Initialize the MOFA entry point and set data options
# ------------------------------------------------------------------------------
ent = entry_point()
ent.set_data_options(
    scale_views=False,
    scale_groups=False,
    center_groups=False
)

# ------------------------------------------------------------------------------
# 3. Add data to the model
# ------------------------------------------------------------------------------
# Provide one likelihood for each view
likelihoods = ["gaussian","gaussian","bernoulli","gaussian","gaussian"]

ent.set_data_matrix(
    data_mat_nested,
    likelihoods=likelihoods,
    views_names=view_names
)

# ------------------------------------------------------------------------------
# 4. Set model options
# ------------------------------------------------------------------------------
ent.set_model_options(
    factors=%FACTOR%,
    spikeslab_weights=False,
    spikeslab_factors=False,
    ard_weights=True,
    ard_factors=False
)

# ------------------------------------------------------------------------------
# 5. Set training options
# ------------------------------------------------------------------------------
ent.set_train_options(
    convergence_mode="slow",
    dropR2=-1,
    iter=20000,
    gpu_mode=True,
    freqELBO=5,
    startELBO=1,
    seed=123
)

# ------------------------------------------------------------------------------
# 6. Build and train the model
# ------------------------------------------------------------------------------
ent.build()
ent.run()

# ------------------------------------------------------------------------------
# 7. Save the trained model output
# ------------------------------------------------------------------------------
ent.save(outfile="%OUTFILE%")

print("MOFA training complete. Results saved to \\"%OUTFILE%\\"")
'

# Template for the Bash script (with placeholders %FACTOR% and %FACTOROUT%)
sh_template <- '#!/bin/bash
#SBATCH --nodes=1
#SBATCH --gres=gpu:1
#SBATCH --ntasks=1
#SBATCH -p ampere
#SBATCH -A SIMIDJIEVSKI-SL3-GPU
#SBATCH --job-name=MOFA_ncp%FACTOR%
#SBATCH --output=logs/MOFA_ncp%FACTOR%.out
#SBATCH --error=logs/MOFA_ncp%FACTOR%.err
#SBATCH --time=11:59:59
#SBATCH --mail-type=ALL
#SBATCH --mail-user=as3582@cam.ac.uk

#! Number of nodes and tasks per node allocated by SLURM (do not change):
numnodes=$SLURM_JOB_NUM_NODES
numtasks=$SLURM_NTASKS
mpi_tasks_per_node=$(echo "$SLURM_TASKS_PER_NODE" | sed -e  \'s/^\\([0-9][0-9]*\\).*$/\\1/\')

. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
module purge                               # Removes all modules still loaded
module load rhel8/default-amp              # REQUIRED - loads the basic environment
module load cuda/11.4

# Insert additional module load commands after this line if needed:

#! Full path to application executable:
application="python MOFA_ncp%FACTOR%.py"

#! Run options for the application:
options=""

#! Work directory (i.e. where the job will run):
workdir="$SLURM_SUBMIT_DIR"  # The value of SLURM_SUBMIT_DIR sets workdir to the directory
                             # in which sbatch is run.

#! Are you using OpenMP (NB this is unrelated to OpenMPI)? If so increase this
#! safe value to no more than 128:
export OMP_NUM_THREADS=1

#! Number of MPI tasks to be started by the application per node and in total (do not change):
np=$[${numnodes}*${mpi_tasks_per_node}]

#! Choose this for a pure shared-memory OpenMP parallel program on a single node:
#! (OMP_NUM_THREADS threads will be created):
CMD="$application $options"

#! Choose this for a MPI code using OpenMPI:
#CMD="mpirun -npernode $mpi_tasks_per_node -np $np $application $options"


###############################################################
### You should not have to change anything below this line ####
###############################################################

cd $workdir
echo -e "Changed directory to `pwd`.\n"

JOBID=$SLURM_JOB_ID

echo -e "JobID: $JOBID\n======"
echo "Time: `date`"
echo "Running on master node: `hostname`"
echo "Current directory: `pwd`"

if [ "$SLURM_JOB_NODELIST" ]; then
        #! Create a machine file:
        export NODEFILE=`generate_pbs_nodefile`
        cat $NODEFILE | uniq > machine.file.$JOBID
        echo -e "\\nNodes allocated:\\n================"
        echo `cat machine.file.$JOBID | sed -e \'s/\\..*$//g\'`
fi

echo -e "\\nnumtasks=$numtasks, numnodes=$numnodes, mpi_tasks_per_node=$mpi_tasks_per_node (OMP_NUM_THREADS=$OMP_NUM_THREADS)"

echo -e "\\nExecuting command:\\n==================\\n$CMD\\n"

source /home/as3582/miniconda3/etc/profile.d/conda.sh
conda activate mofa_env_gpu
mkdir -p logs
eval $CMD
'

# Create the output directory if desired (optional)
out_dir <- "Scripts/single_algorithm/MOFA_GPU"
if (!dir.exists(out_dir)) {
  dir.create(out_dir)
}

for (f in factor_values) {
  # Create Python script content
  py_content <- py_template
  py_content <- gsub("%FACTOR%", f, py_content)
  py_content <- gsub("%OUTFILE%", paste0("MOFA_output_factors", f, "_iter20000_R2_-1.hdf5"), py_content)
  
  # Create Bash script content
  sh_content <- sh_template
  sh_content <- gsub("%FACTOR%", f, sh_content)
  
  # Write the files
  py_filename <- file.path(out_dir, paste0("MOFA_ncp", f, ".py"))
  sh_filename <- file.path(out_dir, paste0("bash_MOFA_ncp", f, ".sh"))
  
  cat(py_content, file = py_filename)
  cat(sh_content, file = sh_filename)
  
  message("Generated scripts for factor = ", f)
}
