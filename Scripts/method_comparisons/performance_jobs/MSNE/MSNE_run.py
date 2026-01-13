#!/usr/bin/env python3
import argparse
import glob
import os

# --- Fairness: hard-cap implicit threading (BLAS/OpenMP) unless the environment already specifies it
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")

import random
import time

import numpy as np
import pandas as pd

# MSNE repo code expects to run from repo root (where "code/" lives)
from code.embedding import MSNE


def set_seeds(seed: int) -> None:
    np.random.seed(seed)
    random.seed(seed)


def load_views_from_csv(dists_dir: str) -> list[pd.DataFrame]:
    csvs = sorted(glob.glob(os.path.join(dists_dir, "*.csv")))
    if len(csvs) != 5:
        raise ValueError(f"Expected 5 distance CSVs in {dists_dir}, found {len(csvs)}")
    views = []
    for fp in csvs:
        df = pd.read_csv(fp, index_col=0)
        views.append(df)
    return views


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dists_dir", required=True)
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--seed", type=int, default=123)

    ap.add_argument("--n_clusters", type=int, default=5)
    ap.add_argument("--k", type=int, default=50)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--walk_length", type=int, default=40)
    ap.add_argument("--num_walks", type=int, default=200)
    ap.add_argument("--embed_size", type=int, default=50)
    ap.add_argument("--window_size", type=int, default=5)

    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    set_seeds(args.seed)

    views = load_views_from_csv(args.dists_dir)

    t0 = time.time()
    res = MSNE(
        views=views,
        n_clusters=args.n_clusters,
        k=args.k,
        workers=args.workers,
        walk_length=args.walk_length,
        num_walks=args.num_walks,
        embed_size=args.embed_size,
        window_size=args.window_size,
        input_type="symmetric_relationships",
    )
    t1 = time.time()

    # Expected keys from your prior runs:
    embeddings = res["embeddings"]  # DataFrame (n_samples x embed_size)
    clusters = res["group"]         # Series / DF with labels

    # Normalise cluster output -> single-column CSV with index as Sample.ID
    if isinstance(clusters, pd.DataFrame):
        if clusters.shape[1] >= 1:
            lab = clusters.iloc[:, 0]
        else:
            raise ValueError("Clusters dataframe has no columns.")
        lab = pd.Series(lab.values, index=clusters.index, name="Cluster_pred")
    else:
        lab = pd.Series(clusters, name="Cluster_pred")

    emb_path = os.path.join(args.out_dir, "msne_embeddings.csv")
    clu_path = os.path.join(args.out_dir, "msne_clusters.csv")
    meta_path = os.path.join(args.out_dir, "msne_python_runtime_seconds.txt")

    embeddings.to_csv(emb_path)
    lab.to_csv(clu_path, header=True)
    with open(meta_path, "w") as f:
        f.write(f"{t1 - t0:.6f}\n")

    print("Wrote:", emb_path)
    print("Wrote:", clu_path)
    print(f"MSNE call wall time (python-internal): {t1 - t0:.2f}s")


if __name__ == "__main__":
    main()
