#!/usr/bin/env python3
import argparse
import json
import os
import threading
import time
from typing import Optional

import numpy as np
import pandas as pd

from mofapy2.run.entry_point import entry_point


VIEW_ORDER = ["RNAseq", "CNV", "SNPs", "miRNA", "Methylation"]
LIKELIHOODS = ["gaussian", "gaussian", "bernoulli", "gaussian", "gaussian"]


def _get_first_visible_gpu_index() -> int:
    cvd = os.environ.get("CUDA_VISIBLE_DEVICES", "")
    if not cvd:
        return 0
    first = cvd.split(",")[0].strip()
    if first == "":
        return 0
    try:
        return int(first)
    except ValueError:
        return 0


class GPUMemMonitor:
    def __init__(self, poll_s: float = 0.05) -> None:
        self.poll_s = poll_s
        self._stop = threading.Event()
        self.peak_mib: Optional[float] = None
        self._t: Optional[threading.Thread] = None
        self._use_nvml = False
        self._nvml = None
        self._handle = None

    def start(self) -> None:
        try:
            import pynvml  # type: ignore

            pynvml.nvmlInit()
            gpu_index = _get_first_visible_gpu_index()
            self._handle = pynvml.nvmlDeviceGetHandleByIndex(gpu_index)
            self._nvml = pynvml
            self._use_nvml = True
        except Exception:
            self._use_nvml = False

        self._t = threading.Thread(target=self._loop, daemon=True)
        self._t.start()

    def stop(self) -> None:
        self._stop.set()
        if self._t is not None:
            self._t.join(timeout=2.0)
        if self._use_nvml and self._nvml is not None:
            try:
                self._nvml.nvmlShutdown()
            except Exception:
                pass

    def _loop(self) -> None:
        peak = 0.0
        while not self._stop.is_set():
            try:
                if self._use_nvml and self._nvml is not None and self._handle is not None:
                    info = self._nvml.nvmlDeviceGetMemoryInfo(self._handle)
                    used_mib = float(info.used) / (1024.0**2)
                    if used_mib > peak:
                        peak = used_mib
                # If NVML isn't available, we just don't record GPU peak.
            except Exception:
                pass
            time.sleep(self.poll_s)
        self.peak_mib = peak if peak > 0 else None


def load_views(input_dir: str) -> tuple[list[np.ndarray], list[str], list[str]]:
    mats: list[np.ndarray] = []
    view_names: list[str] = []
    feature_names: list[str] = []

    sample_ids: Optional[list[str]] = None

    for view in VIEW_ORDER:
        fp = os.path.join(input_dir, f"{view}.csv")
        if not os.path.exists(fp):
            raise FileNotFoundError(f"Missing MOFA view file: {fp}")

        df = pd.read_csv(fp, header=0, index_col=0)
        df = df.apply(pd.to_numeric, errors="coerce").fillna(0.0)

        # Enforce consistent sample set and order across views (columns are samples)
        if sample_ids is None:
            sample_ids = list(map(str, df.columns.tolist()))
        else:
            df = df.loc[:, [c for c in sample_ids if c in df.columns]]

        # For bernoulli view, enforce {0,1}
        if view == "SNPs":
            x = df.to_numpy(dtype=float)
            x = np.where(np.isfinite(x), x, 0.0)
            x = (x != 0).astype(float)
            df = pd.DataFrame(x, index=df.index, columns=df.columns)

        # MOFA wants samples x features
        df_sf = df.T
        mats.append(df_sf.to_numpy(dtype=float))
        view_names.append(view)
        feature_names.append(str(df_sf.shape[1]))

    if sample_ids is None or len(sample_ids) == 0:
        raise ValueError("No samples found in MOFA inputs.")

    return mats, view_names, sample_ids


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input_dir", required=True)
    ap.add_argument("--out_hdf5", required=True)
    ap.add_argument("--metrics_json", required=True)
    ap.add_argument("--seed", type=int, default=123)
    args = ap.parse_args()

    os.makedirs(os.path.dirname(args.out_hdf5), exist_ok=True)
    os.makedirs(os.path.dirname(args.metrics_json), exist_ok=True)

    mats, views, sample_ids = load_views(args.input_dir)
    data_nested = [[m] for m in mats]  # M views, 1 group

    ent = entry_point()
    ent.set_data_options(scale_views=False, scale_groups=False, center_groups=False)

    # Set data
    ent.set_data_matrix(
        data_nested,
        likelihoods=LIKELIHOODS,
        views_names=views,
    )

    # Fixed model options (as requested)
    ent.set_model_options(
        factors=7,
        spikeslab_weights=False,
        spikeslab_factors=False,
        ard_weights=True,
        ard_factors=False,
    )

    # Fixed training options (as requested)
    ent.set_train_options(
        convergence_mode="slow",
        dropR2=-1,
        iter=20000,
        gpu_mode=True,
        freqELBO=5,
        startELBO=1,
        seed=args.seed,
    )

    mon = GPUMemMonitor(poll_s=0.05)
    mon.start()

    t0 = time.perf_counter()
    tb0 = time.perf_counter()
    ent.build()
    tb1 = time.perf_counter()

    tr0 = time.perf_counter()
    ent.run()
    tr1 = time.perf_counter()
    t1 = time.perf_counter()

    mon.stop()

    ent.save(outfile=args.out_hdf5)

    metrics = {
        "seed": args.seed,
        "params": {
            "factors": 7,
            "iter": 20000,
            "dropR2": -1,
            "convergence_mode": "slow",
            "gpu_mode": True,
            "freqELBO": 5,
            "startELBO": 1,
            "likelihoods": LIKELIHOODS,
            "views": views,
        },
        "n_samples": int(len(sample_ids)),
        "build_seconds": float(tb1 - tb0),
        "run_seconds": float(tr1 - tr0),
        "build_plus_run_seconds": float(t1 - t0),
        "peak_gpu_mem_mib": float(mon.peak_mib) if mon.peak_mib is not None else None,
        "outputs": {"hdf5": args.out_hdf5},
    }

    with open(args.metrics_json, "w", encoding="utf-8") as f:
        json.dump(metrics, f, indent=2)

    print(f"[MOFA] build_seconds={metrics['build_seconds']:.4f}")
    print(f"[MOFA] run_seconds={metrics['run_seconds']:.4f}")
    print(f"[MOFA] peak_gpu_mem_mib={metrics['peak_gpu_mem_mib']}")
    print(f"[MOFA] wrote: {args.out_hdf5}")
    print(f"[MOFA] wrote: {args.metrics_json}")


if __name__ == "__main__":
    main()

