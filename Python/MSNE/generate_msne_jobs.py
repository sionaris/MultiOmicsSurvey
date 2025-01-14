#!/usr/bin/env python3

import itertools
import os

def generate_msne_jobs():
    # Parameter ranges
    k_list = [10, 15, 20, 25, 30, 35, 40, 45, 50]
    num_walks_list = [20, 50, 100, 150, 200]
    embed_size_list = [50, 100, 200]
    window_size_list = [5, 10, 15]
    walk_length_list = [20, 30, 40]

    workers = 10  # always

    # We'll keep a base directory as the working dir
    base_dir = os.getcwd()

    # For reproducibility
    seed = 123

    # Create a "logs" dir if needed
    os.makedirs(os.path.join(base_dir, "logs"), exist_ok=True)

    # Generate all combos
    for k in k_list:
        for nw in num_walks_list:
            for emb in embed_size_list:
                for win in window_size_list:
                    for wl in walk_length_list:

                        job_name = f"MSNE_k_{k}_nw_{nw}_emb_{emb}_win_{win}_wl_{wl}"
                        sh_file = f"{job_name}.sh"
                        py_file = f"run_{job_name}.py"

                        # -------------------------
                        # (1) Create the .py script
                        # -------------------------
                        py_content = f"""#!/usr/bin/env python3
import os
import numpy as np
import random
import pandas as pd

# We'll assume your code/ folder is in the current working directory
# and contains embedding.py with MSNE defined inside it.
from code.embedding import MSNE

def set_seeds(seed):
    np.random.seed(seed)
    random.seed(seed)

if __name__ == "__main__":
    set_seeds({seed})

    # Prepare a list to hold each distance matrix
    views = []

    # Load all CSVs from /MO_Dists/ folder
    # Each CSV has sample IDs in the first column, and column headers
    # are sample IDs as well.
    mo_dists_dir = os.path.join(os.getcwd(), "MO_Dists")
    for file_name in os.listdir(mo_dists_dir):
        if file_name.endswith(".csv"):
            path = os.path.join(mo_dists_dir, file_name)
            # read_csv with index_col=0 => first column are sample IDs
            df_temp = pd.read_csv(path, index_col=0)
            # We now have a (samples x samples) matrix
            views.append(df_temp)

    # We call MSNE with input_type="symmetric_relationships"
    # so we skip the internal distance computation
    result = MSNE(
        views=views,
        n_clusters=5,
        k={k},
        workers={workers},
        walk_length={wl},
        num_walks={nw},
        embed_size={emb},
        window_size={win},
        input_type="symmetric_relationships"
    )

    # Extract embeddings & clusters from the MSNE result
    embeddings = result["embeddings"]
    clusters = result["group"]

    # Save them to CSV
    embeddings.to_csv("output_{job_name}_embeddings.csv")
    clusters.to_csv("output_{job_name}_clusters.csv")
"""

                        # ---------------------------
                        # (2) Create the SLURM .sh script
                        # ---------------------------
                        sh_content = f"""#!/bin/bash
#SBATCH --partition=icelake-himem
#SBATCH --nodes=1
#SBATCH --job-name={job_name}
#SBATCH --output=logs/{job_name}.out
#SBATCH --error=logs/{job_name}.err
#SBATCH --time=11:59:59
#SBATCH --cpus-per-task={workers}
#SBATCH --mem=200G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=as3582@cam.ac.uk

# Activate your conda environment
source ~/miniconda3/bin/activate
conda activate MSNE

# Navigate to your main working directory (the MSNE folder)
cd {base_dir}
mkdir -p logs

export SEED={seed}
export PYTHONHASHSEED={seed}

echo "Starting {job_name} with SEED=$SEED"

python {py_file}

echo "{job_name} run completed."
"""

                        # Write out the .py file
                        with open(os.path.join(base_dir, py_file), "w") as f_py:
                            f_py.write(py_content)
                        os.chmod(os.path.join(base_dir, py_file), 0o755)

                        # Write out the .sh file
                        with open(os.path.join(base_dir, sh_file), "w") as f_sh:
                            f_sh.write(sh_content)
                        os.chmod(os.path.join(base_dir, sh_file), 0o755)

    print("Job scripts generated successfully!")

if __name__ == "__main__":
    generate_msne_jobs()
