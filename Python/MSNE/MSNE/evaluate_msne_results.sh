#!/bin/bash
#SBATCH --partition=icelake-himem
#SBATCH --nodes=1
#SBATCH --job-name=MSNE_eval
#SBATCH --output=logs/MSNE_eval.out
#SBATCH --error=logs/MSNE_eval.err
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=as3582@cam.ac.uk

cd /home/as3582/MO_survey/MSNE/MSNE
poetry env activate
mkdir -p logs

export SEED=123
export PYTHONHASHSEED=123

echo "Starting evaluation of MSNE runs..."

poetry run python evaluate_msne_results.py

echo "Evaluation completed!"
