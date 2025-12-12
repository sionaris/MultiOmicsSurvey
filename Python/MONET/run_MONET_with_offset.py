#!/usr/bin/env python3
import os
import numpy as np
import pandas as pd
import random
import glob
import time
import pickle

from monet import monet  # Adjust the import based on your project structure
from monet.monet import Monet


def set_seeds(seed):
    """Set random seeds for reproducibility across NumPy and Python's random module.

    Ensures deterministic behavior for stochastic operations in MONET clustering
    by setting identical seeds for both NumPy's and Python's random number generators.

    Args:
        seed (int): The random seed value to set for reproducibility.

    Returns:
        None

    Example:
        >>> set_seeds(123)
        # All subsequent random operations will be reproducible
    """
    np.random.seed(seed)
    random.seed(seed)


def load_correlation_matrices(directory):
    """Load pre-computed correlation matrices from CSV files for multi-omics integration.

    Reads all CSV files from the specified directory, expecting exactly 5 files
    corresponding to different omics data types (e.g., RNAseq, CNV, Methylation,
    miRNA, SNPs). Each CSV should be a square correlation matrix with sample IDs
    as both row and column indices.

    Args:
        directory (str): Path to the directory containing the CSV correlation
            matrix files. Each file represents one omic view.

    Returns:
        dict: A dictionary mapping omic names (derived from filenames without
            extension) to pandas DataFrames containing the correlation matrices.

    Raises:
        ValueError: If the directory does not contain exactly 5 CSV files.

    Example:
        >>> data = load_correlation_matrices("MO_Corrs_offset")
        >>> print(data.keys())
        dict_keys(['RNAseq', 'CNV', 'Methylation', 'miRNA', 'SNPs'])
    """
    data = {}
    csv_files = glob.glob(os.path.join(directory, "*.csv"))
    if len(csv_files) != 5:
        raise ValueError(f"Expected 5 CSV files in {directory}, found {len(csv_files)}.")
    
    for file_path in csv_files:
        omic_name = os.path.splitext(os.path.basename(file_path))[0]
        df = pd.read_csv(file_path, index_col=0)
        data[omic_name] = df
    return data

def main():
    # ---------------------------
    # 1. Seed for reproducibility
    # ---------------------------
    seed = 123
    set_seeds(seed)
    print(f"Random seed set to: {seed}")

    # --------------------------------
    # 2. Load your pre-computed dists
    # --------------------------------
    corr_dir = os.path.join(os.getcwd(), "MO_Corrs_offset")
    if not os.path.isdir(corr_dir):
        raise FileNotFoundError(f"Correlation {corr_dir} does not exist.")

    print("Loading correlation matrices...")
    data = load_correlation_matrices(corr_dir)
    print(f"Loaded {len(data)} correlation matrices.")

    # ---------------------------
    # 3. Initialize and run MONET
    # ---------------------------
    monet_instance = Monet()

    # Example parameters for ~625 samples; tweak as needed
    iters = 10000
    num_of_seeds = 100
    num_of_samples_in_seed = 10
    min_mod_size = 10
    max_pats_per_action = 10
    percentile_shift = None
    percentile_remove_edge = 80

    print("Starting MONET main loop...")
    start_time = time.time()
    # main_loop returns: glob_var, total_time, super_g, iteration_times, total_weight
    print("------ Using MONET parameters ------")
    print("iters =", iters)
    print("num_of_seeds =", num_of_seeds)
    print("num_of_samples_in_seed =", num_of_samples_in_seed)
    print("min_mod_size =", min_mod_size)
    print("max_pats_per_action =", max_pats_per_action)
    print("percentile_shift =", percentile_shift)
    print("percentile_remove_edge =", percentile_remove_edge)
    print("------------------------------------")

    glob_var, total_time, super_g, iteration_times, total_weight = monet_instance.main_loop(
        data=data,
        is_input_raw=False,
        init_modules=None,
        iters=iters,
        num_of_seeds=num_of_seeds,
        num_of_samples_in_seed=num_of_samples_in_seed,
        min_mod_size=min_mod_size,
        max_pats_per_action=max_pats_per_action,
        percentile_shift=percentile_shift,
        percentile_remove_edge=percentile_remove_edge
    )
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"MONET main loop completed in {elapsed_time:.2f} seconds.")

    # -------------------------------------------------------
    # 4. Build a numeric clustering from the final modules
    # -------------------------------------------------------
    all_modules = glob_var.modules   # final dictionary {module_name -> module_object}
    numeric_clustering = {}
    for i, (mod_name, module) in enumerate(all_modules.items(), start=1):
        for pat in module.patients.keys():
            numeric_clustering[pat] = i

    # -------------------------------------------------------
    # 5. Build a dictionary with EVERYTHING we care about
    # -------------------------------------------------------
    monet_ret_dict = {
        "glob_var": glob_var,               # The entire MONET state object
        "super_g": super_g,                 # If you need it
        "total_time": total_time,           # The total time from MONET's perspective
        "iteration_times": iteration_times, # Per-iteration runtimes
        "total_weight": total_weight,       # The final aggregated score
        "clustering": numeric_clustering,   # Sample -> integer cluster ID
        "all_modules": all_modules         # The final modules dictionary
        # Optional: you can include any other fields you want
        # e.g. the raw data or 'data' if you want to re-run partial steps
    }

    # -------------------------------------------------------
    # 6. Save the comprehensive object to a pickle
    # -------------------------------------------------------
    output_dir = os.path.join(os.getcwd(), "monet_results")
    os.makedirs(output_dir, exist_ok=True)
    results_path = os.path.join(output_dir, "monet_full_output_with_offset.pkl")
    with open(results_path, 'wb') as f:
        pickle.dump(monet_ret_dict, f)

    print(f"Results saved to {results_path}")

if __name__ == "__main__":
    main()
