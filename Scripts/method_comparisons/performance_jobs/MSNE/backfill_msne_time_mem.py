#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import glob
import os
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple


ELAPSED_RE = re.compile(
    r"^Elapsed \(wall clock\) time \(h:mm:ss or m:ss\):\s*([0-9:.]+)\s*$"
)
MAXRSS_RE = re.compile(
    r"^Maximum resident set size \(kbytes\):\s*([0-9]+)\s*$"
)


def to_seconds(hms: str) -> Optional[float]:
    """
    Accepts m:ss.xx or h:mm:ss.xx (as emitted by GNU time -v) and returns seconds.
    """
    hms = hms.strip()
    parts = hms.split(":")
    try:
        if len(parts) == 2:
            m = float(parts[0])
            s = float(parts[1])
            return m * 60.0 + s
        if len(parts) == 3:
            h = float(parts[0])
            m = float(parts[1])
            s = float(parts[2])
            return h * 3600.0 + m * 60.0 + s
    except ValueError:
        return None
    return None


def parse_timev(path: Path) -> Tuple[Optional[float], Optional[float]]:
    """
    Returns (elapsed_seconds, peak_mib) from a GNU time -v output file.
    Strips leading whitespace because your files are tab-indented.
    """
    if not path.exists():
        return None, None

    elapsed_s: Optional[float] = None
    maxrss_kb: Optional[int] = None

    for raw in path.read_text(errors="replace").splitlines():
        line = raw.lstrip()  # <-- critical for your tab-indented files
        m1 = ELAPSED_RE.match(line)
        if m1:
            elapsed_s = to_seconds(m1.group(1))
            continue
        m2 = MAXRSS_RE.match(line)
        if m2:
            try:
                maxrss_kb = int(m2.group(1))
            except ValueError:
                maxrss_kb = None

    peak_mib = (maxrss_kb / 1024.0) if maxrss_kb is not None else None
    return elapsed_s, peak_mib


def read_tsv_single_row(path: Path) -> Tuple[List[str], Dict[str, str]]:
    with path.open(newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        rows = list(r)
        if len(rows) != 1:
            raise ValueError(f"Expected exactly 1 row in {path}, found {len(rows)}")
        return r.fieldnames or [], rows[0]


def write_tsv_single_row(path: Path, fieldnames: List[str], row: Dict[str, str]) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        w.writeheader()
        # Ensure all keys exist
        out = {k: row.get(k, "NA") for k in fieldnames}
        w.writerow(out)
    tmp.replace(path)


def is_missing(v: Optional[str]) -> bool:
    if v is None:
        return True
    v = v.strip()
    return v == "" or v.upper() == "NA"


def backfill_task(task_dir: Path, dry_run: bool, make_bak: bool) -> int:
    """
    Backfills all MSNE_perf_row*.tsv in task_dir/out using task_dir/meta/timev_msne.txt.
    Returns number of row files updated.
    """
    meta = task_dir / "meta" / "timev_msne.txt"
    out_dir = task_dir / "out"
    if not out_dir.is_dir():
        return 0

    elapsed_s, peak_mib = parse_timev(meta)
    if elapsed_s is None and peak_mib is None:
        return 0

    row_files = sorted(out_dir.glob("MSNE_perf_row*.tsv"))
    updated = 0

    for rf in row_files:
        fieldnames, row = read_tsv_single_row(rf)

        changed = False
        if "MSNE_Elapsed_Seconds" in row and not is_missing(str(elapsed_s)) and is_missing(row.get("MSNE_Elapsed_Seconds")):
            row["MSNE_Elapsed_Seconds"] = f"{elapsed_s:.6f}"
            changed = True
        if "MSNE_PeakRAM_MiB" in row and not is_missing(str(peak_mib)) and is_missing(row.get("MSNE_PeakRAM_MiB")):
            row["MSNE_PeakRAM_MiB"] = f"{peak_mib:.6f}"
            changed = True

        if changed:
            updated += 1
            if not dry_run:
                if make_bak:
                    bak = rf.with_suffix(rf.suffix + ".bak")
                    if not bak.exists():
                        bak.write_bytes(rf.read_bytes())
                write_tsv_single_row(rf, fieldnames, row)

    return updated


def union_fieldnames(dicts: List[Dict[str, str]]) -> List[str]:
    seen = set()
    out: List[str] = []
    for d in dicts:
        for k in d.keys():
            if k not in seen:
                seen.add(k)
                out.append(k)
    return out


def read_any_rows(paths: List[Path]) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []
    for p in paths:
        fns, row = read_tsv_single_row(p)
        rows.append({k: row.get(k, "NA") for k in fns})
    return rows


def sort_feature_rows(rows: List[Dict[str, str]]) -> List[Dict[str, str]]:
    def keyfun(r: Dict[str, str]):
        v = r.get("Feature_Percent", "NA")
        try:
            return int(float(v))
        except Exception:
            return 10**9
    return sorted(rows, key=keyfun)


def sort_sample_rows(rows: List[Dict[str, str]]) -> List[Dict[str, str]]:
    def keyfun(r: Dict[str, str]):
        a = r.get("Sample_Fraction", "NA")
        b = r.get("Replicate", "NA")
        try:
            aa = int(float(a))
        except Exception:
            aa = 10**9
        try:
            bb = int(float(b))
        except Exception:
            bb = 10**9
        return (aa, bb)
    return sorted(rows, key=keyfun)


def write_table(path: Path, rows: List[Dict[str, str]]) -> None:
    if not rows:
        return
    fieldnames = union_fieldnames(rows)
    tmp = path.with_suffix(path.suffix + ".tmp")
    path.parent.mkdir(parents=True, exist_ok=True)
    with tmp.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "NA") if not is_missing(r.get(k)) else "NA" for k in fieldnames})
    tmp.replace(path)


def discover_array_dirs(root: Path) -> List[Path]:
    if root.name.startswith("array_") and root.is_dir():
        return [root]
    return sorted([Path(p) for p in glob.glob(str(root / "array_*")) if Path(p).is_dir()])


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("path", help="An array_XXXX directory, or a directory containing array_* subdirs.")
    ap.add_argument("--dry-run", action="store_true", help="Scan and report, but do not modify files.")
    ap.add_argument("--no-bak", action="store_true", help="Do not create .bak backups of modified perf_row files.")
    args = ap.parse_args()

    root = Path(args.path).resolve()
    array_dirs = discover_array_dirs(root)
    if not array_dirs:
        raise SystemExit(f"No array_* directories found under: {root}")

    total_rows_updated = 0

    for arr in array_dirs:
        # Patch all task dirs with a meta/timev_msne.txt
        task_dirs = sorted([p.parent.parent for p in arr.glob("task_*/meta/timev_msne.txt")])
        updated_here = 0
        for td in task_dirs:
            updated_here += backfill_task(td, dry_run=args.dry_run, make_bak=not args.no_bak)
        total_rows_updated += updated_here

        # Rebuild aggregates exactly like your merge steps do (just from perf_row TSVs)
        out_dir = arr / "out"

        feat_row_paths = sorted(
            [Path(p) for p in glob.glob(str(arr / "**" / "MSNE_perf_row_*pct.tsv"), recursive=True)]
        )
        feat_row_paths = [p for p in feat_row_paths if "samples" not in p.name]  # avoid sample rows
        if feat_row_paths:
            feat_rows = sort_feature_rows(read_any_rows(feat_row_paths))
            write_table(out_dir / "MSNE_feature_perturbations_performance.tsv", feat_rows)

        samp_row_paths = sorted(
            [Path(p) for p in glob.glob(str(arr / "**" / "MSNE_perf_row_samples_*.tsv"), recursive=True)]
        )
        if samp_row_paths:
            samp_rows = sort_sample_rows(read_any_rows(samp_row_paths))
            write_table(out_dir / "MSNE_sample_perturbations_performance.tsv", samp_rows)

        print(f"[{arr.name}] updated perf_row files: {updated_here} (dry_run={args.dry_run})")

    print(f"TOTAL updated perf_row files: {total_rows_updated} (dry_run={args.dry_run})")


if __name__ == "__main__":
    main()
