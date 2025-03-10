#!/usr/bin/env python3
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
    set_seeds(123)

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
        k=10,
        workers=10,
        walk_length=30,
        num_walks=50,
        embed_size=200,
        window_size=10,
        input_type="symmetric_relationships"
    )

    # Extract embeddings & clusters from the MSNE result
    embeddings = result["embeddings"]
    clusters = result["group"]

    # Save them to CSV
    embeddings.to_csv("output_MSNE_k_10_nw_50_emb_200_win_10_wl_30_embeddings.csv")
    clusters.to_csv("output_MSNE_k_10_nw_50_emb_200_win_10_wl_30_clusters.csv")
