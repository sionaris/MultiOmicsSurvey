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
    """
    Set the random seeds for reproducibility.
    """
    np.random.seed(seed)
    random.seed(seed)

def load_distance_matrices(directory):
    """
    Load all CSV distance matrices from the specified directory.
    
    :param directory: Path to the directory containing CSV files.
    :return: Dictionary mapping omic names to pandas DataFrames.
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
    # Set the seed for reproducibility
    seed = 123
    set_seeds(seed)
    print(f"Random seed set to: {seed}")

    # Define the path to the distance matrices
    distance_dir = os.path.join(os.getcwd(), "MO_Dists")
    if not os.path.isdir(distance_dir):
        raise FileNotFoundError(f"Directory {distance_dir} does not exist.")
    
    # Load the distance matrices
    print("Loading distance matrices...")
    data = load_distance_matrices(distance_dir)
    print(f"Loaded {len(data)} distance matrices.")

    # Initialize MONET
    monet = Monet()
    
    # Define parameters tailored for 625 samples
    iters = 5000
    num_of_seeds = 100  # A reasonable number of initial seeds; MONET will adjust as needed
    num_of_samples_in_seed = 25
    min_mod_size = 20
    max_pats_per_action = 10
    percentile_shift = None
    percentile_remove_edge = 80

    # Run the main loop
    print("Starting MONET main loop...")
    start_time = time.time()
    glob_var, total_time, _, iteration_times, total_weight = monet.main_loop(
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

    # Save the results
    output_dir = os.path.join(os.getcwd(), "monet_results")
    os.makedirs(output_dir, exist_ok=True)
    results_path = os.path.join(output_dir, "monet_output.pkl")
    
    with open(results_path, 'wb') as f:
        pickle.dump({
            'glob_var': glob_var,
            'total_time': total_time,
            'iteration_times': iteration_times,
            'total_weight': total_weight
        }, f)
    print(f"Results saved to {results_path}")

if __name__ == "__main__":
    main()
