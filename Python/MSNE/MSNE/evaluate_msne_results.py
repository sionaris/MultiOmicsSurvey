#!/usr/bin/env python3
"""MSNE Hyperparameter Search Results Evaluator.

This script evaluates the results of an MSNE hyperparameter grid search
by computing silhouette scores for each clustering result and identifying
the best-performing hyperparameter combination.
"""

import os
import glob
import pandas as pd
import numpy as np
from sklearn.metrics import silhouette_score


def main():
    """Evaluate MSNE clustering results and identify optimal hyperparameters.

    Scans the current directory for MSNE output files (embeddings and cluster
    assignments), computes silhouette scores for each hyperparameter combination,
    and generates a comprehensive evaluation report.

    The function:
        1. Discovers all embedding/cluster file pairs matching 'output_*_embeddings.csv'
        2. Computes silhouette scores using sklearn.metrics.silhouette_score
        3. Ranks all hyperparameter combinations by silhouette score
        4. Saves detailed results to 'all_silhouette_scores.csv'
        5. Generates 'best_clustering_report.txt' with the optimal configuration

    File Naming Convention:
        - Embeddings: output_MSNE_k_{k}_nw_{nw}_emb_{emb}_win_{win}_wl_{wl}_embeddings.csv
        - Clusters: output_MSNE_k_{k}_nw_{nw}_emb_{emb}_win_{win}_wl_{wl}_clusters.csv

    Returns:
        None. Results are written to files in the current directory.

    Side Effects:
        - Creates 'all_silhouette_scores.csv' with columns:
            [job_name, silhouette, n_samples, n_clusters, embeddings_file, clusters_file]
        - Creates 'best_clustering_report.txt' with the optimal hyperparameter summary
        - Prints summary to stdout

    Notes:
        - Skips results with fewer than 2 clusters (silhouette undefined)
        - Silhouette score ranges from -1 (poor) to 1 (excellent)
        - Higher silhouette indicates better-defined cluster separation

    Example:
        >>> main()
        Best clustering job: MSNE_k_20_nw_100_emb_100_win_10_wl_30
        Silhouette: 0.4523
        N samples: 625 N clusters: 5
    """
    # Look for "output_*.csv" files in the current directory.
    # Each job produces:
    #   output_<job_name>_embeddings.csv
    #   output_<job_name>_clusters.csv
    # We'll match them by job_name and compute silhouette scores.

    # Gather all embeddings
    emb_files = sorted(glob.glob("output_*_embeddings.csv"))

    # We'll store a dict of results
    results = []  # list of dicts: {job_name, silhouette, n_samples, n_clusters, embeddings_file, clusters_file}

    for emb_file in emb_files:
        # Derive job name by removing suffix
        # e.g. "output_MSNE_k_50_nw_50_emb_50_win_5_wl_20_embeddings.csv" -> "output_MSNE_k_50_nw_50_emb_50_win_5_wl_20"
        base_name = emb_file.rsplit("_embeddings.csv", 1)[0]
        # That means the cluster file is base_name + "_clusters.csv"
        clust_file = base_name + "_clusters.csv"

        if not os.path.isfile(clust_file):
            continue  # skip if we can't find matching cluster file

        # Load embeddings
        emb_df = pd.read_csv(emb_file, index_col=0)
        # Load clusters
        clust_df = pd.read_csv(clust_file, index_col=0)

        # Both must have the same number of rows
        if emb_df.shape[0] != clust_df.shape[0]:
            # Possibly handle mismatch
            continue

        labels = clust_df.iloc[:, 0].values  # cluster labels in the first column
        if len(set(labels)) < 2:
            # silhouette requires at least 2 clusters
            continue

        # Compute silhouette
        # emb_df is shape (n_samples, n_features)
        try:
            score = silhouette_score(emb_df.values, labels)
        except ValueError:
            # Possibly occurs if there's only one sample or some weird condition
            continue

        # Record
        job_name = base_name.replace("output_", "")  # or just base_name, up to you
        n_clusters = len(set(labels))
        results.append({
            "job_name": job_name,
            "silhouette": score,
            "n_samples": emb_df.shape[0],
            "n_clusters": n_clusters,
            "embeddings_file": emb_file,
            "clusters_file": clust_file
        })

    # If no results found
    if len(results) == 0:
        print("No valid results found!")
        return

    # Convert to DataFrame
    df_results = pd.DataFrame(results)
    df_results.sort_values(by="silhouette", ascending=False, inplace=True)

    # Best result is the top row
    best = df_results.iloc[0]

    # Print a short summary
    print("Best clustering job:", best["job_name"])
    print("Silhouette:", best["silhouette"])
    print("N samples:", best["n_samples"], "N clusters:", best["n_clusters"])

    # Save a more detailed report
    df_results.to_csv("all_silhouette_scores.csv", index=False)

    # Also produce a best-clustering text summary
    with open("best_clustering_report.txt", "w") as f:
        f.write(f"Best Clustering:\n")
        f.write(f"Job Name: {best['job_name']}\n")
        f.write(f"Silhouette Score: {best['silhouette']:.4f}\n")
        f.write(f"N Samples: {best['n_samples']}  N Clusters: {best['n_clusters']}\n")
        f.write(f"Embeddings File: {best['embeddings_file']}\n")
        f.write(f"Clusters File: {best['clusters_file']}\n")

if __name__ == "__main__":
    main()
