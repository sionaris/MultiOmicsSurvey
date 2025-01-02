#!/bin/bash
#SBATCH --job-name=MONET_HPC_run            # Job name
#SBATCH --output=logs/MONET_HPC_run_%j.out   # Standard output log
#SBATCH --error=logs/MONET_HPC_run_%j.err    # Standard error log
#SBATCH --ntasks=1                           # Number of tasks (processes)
#SBATCH --cpus-per-task=4                    # Number of CPU cores per task
#SBATCH --mem=32G                            # Total memory
#SBATCH --time=11:59:59                       # Time limit hrs:min:sec
#SBATCH --partition=icelake-himem             # Partition name
#SBATCH --mail-type=ALL                       # Send email on all events
#SBATCH --mail-user=as3582@cam.ac.uk           # Email address

# Activate conda environment
source ~/miniconda3/bin/activate
conda activate MONET 

# Navigate to the project directory
cd ~/MO_survey/MONET/MONET  # Replace with your actual project path 

# Export environment variables for reproducibility
export SEED=123
export PYTHONHASHSEED=123

# Create logs directory if it doesn't exist
mkdir -p logs

# Run the Python script
echo "Starting MONET run with SEED=$SEED and PYTHONHASHSEED=$PYTHONHASHSEED"
python run_MONET.py
echo "MONET run completed."

