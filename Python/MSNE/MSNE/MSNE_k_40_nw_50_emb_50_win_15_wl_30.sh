#!/bin/bash
#SBATCH --partition=icelake-himem
#SBATCH --nodes=1
#SBATCH --job-name=MSNE_k_40_nw_50_emb_50_win_15_wl_30
#SBATCH --output=logs/MSNE_k_40_nw_50_emb_50_win_15_wl_30.out
#SBATCH --error=logs/MSNE_k_40_nw_50_emb_50_win_15_wl_30.err
#SBATCH --time=11:59:59
#SBATCH --cpus-per-task=10
#SBATCH --mem=200G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=as3582@cam.ac.uk

cd /home/as3582/MO_survey/MSNE/MSNE
poetry env activate
mkdir -p logs

export SEED=123
export PYTHONHASHSEED=123

echo "Starting MSNE_k_40_nw_50_emb_50_win_15_wl_30 with SEED=$SEED"

poetry run python run_MSNE_k_40_nw_50_emb_50_win_15_wl_30.py

echo "MSNE_k_40_nw_50_emb_50_win_15_wl_30 run completed."
