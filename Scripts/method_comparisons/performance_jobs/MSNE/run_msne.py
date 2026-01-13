#!/usr/bin/env python3
import argparse
import os

# --- Fairness: hard-cap implicit threading (BLAS/OpenMP) unless the environment already specifies it
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")

import random
import numpy as np
import pandas as pd

from code.embedding import MSNE

BEST = dict(
    n_clusters=5,
    k=50,
    walk_length=40,
    num_walks=200,
    embed_size=50,
    window_size=5,
)

def set_seeds(seed: int) -> None:
    os.environ["PYTHONHASHSEED"] = str(seed)
    random.seed(seed)
    np.random.seed(seed)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mo_dists", required=True, help="Folder containing *.csv distance matrices")
    ap.add_argument("--out_prefix", required=True, help="Prefix for output files")
    ap.add_argument("--seed", type=int, default=123)
    ap.add_argument("--workers", type=int, default=1)
    args = ap.parse_args()

    set_seeds(args.seed)

    views = []
    for fn in sorted(os.listdir(args.mo_dists)):
        if fn.endswith(".csv"):
            df = pd.read_csv(os.path.join(args.mo_dists, fn), index_col=0)
            views.append(df)

    res = MSNE(
        views=views,
        input_type="symmetric_relationships",
        workers=args.workers,
        seed=args.seed,
        **BEST
    )

    emb = res["embeddings"]
    clu = res["group"]

    emb.to_csv(args.out_prefix + "_embeddings.csv")
    clu.to_csv(args.out_prefix + "_clusters.csv")

if __name__ == "__main__":
    main()
