#!/usr/bin/env python3
import argparse
import glob
import json
import os
import random
import time

import networkx as nx
import numpy as np
import pandas as pd

from monet.monet import Monet


def set_seeds(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)


def _load_one_csv(file_path: str) -> pd.DataFrame:
    df = pd.read_csv(file_path, index_col=0)
    df = df.apply(pd.to_numeric, errors="coerce")
    if df.shape[0] != df.shape[1]:
        raise ValueError(f"{file_path}: correlation matrix is not square: {df.shape}")
    # Ensure columns match index (same sample set)
    if set(df.index) != set(df.columns):
        raise ValueError(f"{file_path}: index/columns sample sets differ")
    # Reorder columns to match index order
    df = df.loc[df.index, df.index]
    return df


def load_correlation_matrices(directory: str) -> dict[str, pd.DataFrame]:
    csv_files = sorted(glob.glob(os.path.join(directory, "*.csv")))
    if len(csv_files) != 5:
        raise ValueError(f"Expected 5 CSV files in {directory}, found {len(csv_files)}.")

    data: dict[str, pd.DataFrame] = {}
    for file_path in csv_files:
        omic_name = os.path.splitext(os.path.basename(file_path))[0]
        data[omic_name] = _load_one_csv(file_path)
    return data


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--corr_dir", required=True, help="Directory with 5 correlation CSVs (omics).")
    ap.add_argument("--out_dir", required=True, help="Output directory for MONET exports.")
    ap.add_argument("--seed", type=int, default=123)

    # Fixed MONET params (your original/default)
    ap.add_argument("--iters", type=int, default=10000)
    ap.add_argument("--num_of_seeds", type=int, default=100)
    ap.add_argument("--num_of_samples_in_seed", type=int, default=10)
    ap.add_argument("--min_mod_size", type=int, default=10)
    ap.add_argument("--max_pats_per_action", type=int, default=10)
    ap.add_argument("--percentile_shift", default="None")  # keep None unless explicitly set
    ap.add_argument("--percentile_remove_edge", type=int, default=80)

    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    set_seeds(args.seed)

    if not os.path.isdir(args.corr_dir):
        raise FileNotFoundError(f"corr_dir does not exist: {args.corr_dir}")

    data = load_correlation_matrices(args.corr_dir)

    # Harmonise sample set across omics (intersection) + consistent order
    common = set.intersection(*[set(df.index) for df in data.values()])
    if len(common) == 0:
        raise ValueError("No common sample IDs across correlation matrices.")
    sample_ids = sorted(common)
    data = {k: df.loc[sample_ids, sample_ids] for k, df in data.items()}
    n = len(sample_ids)

    monet_instance = Monet()

    percentile_shift = None if str(args.percentile_shift).lower() == "none" else float(args.percentile_shift)

    t0 = time.time()
    glob_var, total_time, super_g, iteration_times, total_weight = monet_instance.main_loop(
        data=data,
        is_input_raw=False,
        init_modules=None,
        iters=args.iters,
        num_of_seeds=args.num_of_seeds,
        num_of_samples_in_seed=args.num_of_samples_in_seed,
        min_mod_size=args.min_mod_size,
        max_pats_per_action=args.max_pats_per_action,
        percentile_shift=percentile_shift,
        percentile_remove_edge=args.percentile_remove_edge,
    )
    t1 = time.time()
    main_loop_elapsed = float(t1 - t0)

    # Numeric clustering
    all_modules = getattr(glob_var, "modules")
    numeric_clustering: dict[str, int] = {}
    for i, (_mod_name, module) in enumerate(all_modules.items(), start=1):
        for pat in module.patients.keys():
            numeric_clustering[str(pat)] = int(i)

    # Assign unclustered samples (if any) to a separate cluster at the end
    assigned = set(numeric_clustering.keys())
    max_cluster = max(numeric_clustering.values()) if numeric_clustering else 0
    for sid in sample_ids:
        if str(sid) not in assigned:
            numeric_clustering[str(sid)] = int(max_cluster + 1)

    # Average adjacency across omics graphs (for silhouette in R)
    omics = getattr(glob_var, "omics")
    omic_names = list(omics.keys())

    avg_adj = np.zeros((n, n), dtype=float)
    for om in omic_names:
        G = omics[om].graph
        A = nx.to_numpy_array(G, nodelist=sample_ids, weight="weight", dtype=float)
        np.fill_diagonal(A, 1.0)
        avg_adj += A
    avg_adj /= float(len(omic_names))

    # Save clustering
    clust_df = pd.DataFrame(
        {"Sample.ID": [str(s) for s in sample_ids], "Cluster": [int(numeric_clustering[str(s)]) for s in sample_ids]}
    )
    clust_path = os.path.join(args.out_dir, "monet_clustering.tsv")
    clust_df.to_csv(clust_path, sep="\t", index=False)

    # Save avg adjacency as gzipped TSV
    avg_df = pd.DataFrame(avg_adj, index=[str(s) for s in sample_ids], columns=[str(s) for s in sample_ids])
    avg_path = os.path.join(args.out_dir, "monet_avg_adjacency.tsv.gz")
    avg_df.to_csv(avg_path, sep="\t", compression="gzip")

    metrics = {
        "seed": args.seed,
        "params": {
            "iters": args.iters,
            "num_of_seeds": args.num_of_seeds,
            "num_of_samples_in_seed": args.num_of_samples_in_seed,
            "min_mod_size": args.min_mod_size,
            "max_pats_per_action": args.max_pats_per_action,
            "percentile_shift": percentile_shift,
            "percentile_remove_edge": args.percentile_remove_edge,
        },
        "n_samples": n,
        "n_modules": len(all_modules),
        "main_loop_elapsed_seconds": main_loop_elapsed,
        "monet_total_time": float(total_time) if total_time is not None else None,
        "total_weight": float(total_weight) if total_weight is not None else None,
        "iteration_times_len": len(iteration_times) if iteration_times is not None else None,
        "outputs": {"clustering_tsv": clust_path, "avg_adjacency_tsv_gz": avg_path},
    }

    metrics_path = os.path.join(args.out_dir, "monet_metrics.json")
    with open(metrics_path, "w", encoding="utf-8") as f:
        json.dump(metrics, f, indent=2)

    print(f"[MONET] main_loop_elapsed_seconds={main_loop_elapsed:.4f}")
    print(f"[MONET] wrote: {clust_path}")
    print(f"[MONET] wrote: {avg_path}")
    print(f"[MONET] wrote: {metrics_path}")


if __name__ == "__main__":
    main()
