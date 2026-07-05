#!/usr/bin/env python3
"""
gxcnv1_predict.py  —  Wisecondor-core CNV prediction with gxcnv1 enhancements.

Core algorithm: Wisecondor test (Python2 wisecondor.py, same WC reference.npz).
Enhancements added on top:
  • MAD-based robust z-score (mad_z column)  — more outlier-resistant than WC z-score
  • MAPD noise metric (Median Absolute Pairwise Difference)
  • Per-bin log2(ratio) preserved in output
  • TSV output schema compatible with gxcnv2 (bins / segments / calls / qcmetrics)

Strategy:
  1. sample.npz already produced by GXCNV1_CONVERT (wisecondor.py convert)
  2. out.npz  already produced by wisecondor.py test (called from GXCNV1_PREDICT shell)
  3. Load out.npz to get per-bin z-scores (results_z) and ratios (results_r)
  4. Re-annotate with MAD z-score and MAPD
  5. Write gxcnv2-compatible TSV outputs

Wisecondor out.npz layout (Python2-generated, loaded with allow_pickle=True):
  binsize        int   — bin size in bp
  results_z      list  — per-chromosome per-bin z-score arrays (chr1..chr22)
  results_r      list  — per-chromosome per-bin (ratio - 1) arrays (chr1..chr22)
  results_calls  list  — [[chrom_1based, start_bin, end_bin, zscore, ratio-1], ...]
  threshold_z    float — z-score threshold determined during test

Reference format: unchanged Wisecondor NPZ — same files used by RUN_WC.

Usage:
  gxcnv1_predict.py  sample.npz  out.npz  -o PREFIX  [--zscore FLOAT]
"""

import argparse
import os
import sys
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)

import logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [gxcnv1] %(levelname)s %(message)s",
)
logger = logging.getLogger("gxcnv1")

# Wisecondor analyses chr1-chr22 by default (no chrX/Y)
CHROM_NAMES = [f"chr{i}" for i in range(1, 23)]
CHROM_ORDER = {c: i for i, c in enumerate(CHROM_NAMES)}
CHROM_ORDER["chrX"] = 22
CHROM_ORDER["chrY"] = 23


# ── Load Wisecondor out.npz ────────────────────────────────────────────────────

def load_wisecondor_result(out_npz_path: str) -> dict | None:
    """
    Load wisecondor.py test output NPZ (Python2-generated).

    Returns dict with keys: binsize, results_z, results_r, results_calls, threshold_z
    Returns None on failure.
    """
    try:
        data = np.load(out_npz_path, allow_pickle=True, encoding="bytes")
    except Exception as e:
        logger.error("Error parsing out.npz: %s", e)
        logger.error("Failed to load Wisecondor results — aborting")
        return None

    try:
        binsize     = int(data["binsize"])
        results_z   = list(data["results_z"])    # list of per-chrom z-score arrays
        results_r   = list(data["results_r"])    # list of per-chrom (ratio - 1) arrays
        threshold_z = float(data["threshold_z"])

        calls_raw = data["results_calls"]
        if hasattr(calls_raw, "tolist"):
            calls_raw = calls_raw.tolist()
        results_calls = list(calls_raw) if calls_raw is not None else []

    except KeyError as e:
        logger.error("out.npz missing expected key: %s", e)
        return None
    except Exception as e:
        logger.error("Error parsing out.npz: %s", e)
        return None

    logger.info(
        "Loaded out.npz: binsize=%d  threshold_z=%.2f  n_chrom=%d  n_calls=%d",
        binsize, threshold_z, len(results_z), len(results_calls),
    )
    return {
        "binsize":       binsize,
        "results_z":     results_z,
        "results_r":     results_r,
        "results_calls": results_calls,
        "threshold_z":   threshold_z,
    }


# ── Build per-bin DataFrame ────────────────────────────────────────────────────

def build_bins_df(wc: dict) -> pd.DataFrame:
    """
    Convert Wisecondor per-chromosome arrays to a flat per-bin DataFrame.

    results_z[i]  = per-bin z-score array for chromosome i+1 (chr1=0 .. chr22=21)
    results_r[i]  = per-bin (ratio - 1) array for chromosome i+1

    log2_ratio = log2(results_r + 1)
      0     → diploid
      +0.585 → trisomy-like (3 copies)
      -1.0  → monosomy-like (1 copy)

    Bins masked in the Wisecondor test are set to 0 by inflateArrayMulti;
    they appear as log2_ratio=0 / z_score=0 (treated as diploid in plots).
    """
    binsize   = wc["binsize"]
    results_z = wc["results_z"]
    results_r = wc["results_r"]

    rows = []
    n_chrom = min(len(results_z), len(results_r), len(CHROM_NAMES))
    for chrom_idx in range(n_chrom):
        chrom = CHROM_NAMES[chrom_idx]
        z_arr = np.asarray(results_z[chrom_idx], dtype=float)
        r_arr = np.asarray(results_r[chrom_idx], dtype=float)   # ratio - 1
        ratio = r_arr + 1.0
        log2r = np.where(ratio > 0, np.log2(ratio), float("nan"))

        for bin_idx in range(len(z_arr)):
            rows.append({
                "chrom":      chrom,
                "start":      bin_idx * binsize,
                "end":        (bin_idx + 1) * binsize,
                "log2_ratio": float(log2r[bin_idx]),
                "z_score":    float(z_arr[bin_idx]),
            })

    if not rows:
        return pd.DataFrame(columns=["chrom", "start", "end", "log2_ratio", "z_score"])

    df = pd.DataFrame(rows)
    df["_ci"] = df["chrom"].map(CHROM_ORDER).fillna(99).astype(int)
    df = df.sort_values(["_ci", "start"]).drop(columns="_ci").reset_index(drop=True)
    return df


# ── Enhancements (same as gxcnv2) ─────────────────────────────────────────────

def add_mad_z(df: pd.DataFrame, col: str = "log2_ratio") -> pd.DataFrame:
    """Add MAD-based robust z-score (mad_z)."""
    lr = df[col].dropna().values
    if len(lr) < 2:
        df = df.copy()
        df["mad_z"] = float("nan")
        return df
    median = np.median(lr)
    mad    = np.median(np.abs(lr - median))
    df = df.copy()
    if mad == 0:
        df["mad_z"] = 0.0
    else:
        df["mad_z"] = (df[col] - median) / (1.4826 * mad)
    return df


def compute_mapd(df: pd.DataFrame, col: str = "log2_ratio") -> float:
    """MAPD — Median Absolute Pairwise Difference. Values >0.35 indicate a noisy sample."""
    vals = df[col].dropna().values
    if len(vals) < 2:
        return float("nan")
    return float(np.median(np.abs(np.diff(vals))))


# ── Build calls DataFrame ──────────────────────────────────────────────────────

def build_calls_df(wc: dict, df_bins: pd.DataFrame) -> pd.DataFrame:
    """
    Convert Wisecondor results_calls to a calls DataFrame.

    results_calls entry: [chrom_1based, start_bin, end_bin, zscore, ratio_minus1]
      chrom_1based: integer 1-22
      start_bin, end_bin: 0-based bin indices (inclusive)
      zscore:       Stouffer z-score of the segment
      ratio_minus1: median(ratio) - 1 of the segment
    """
    binsize       = wc["binsize"]
    results_calls = wc["results_calls"]
    if not results_calls:
        return pd.DataFrame()

    rows = []
    for entry in results_calls:
        try:
            chrom_1based = int(entry[0])
            start_bin    = int(entry[1])
            end_bin      = int(entry[2])
            zscore       = float(entry[3])
            ratio_m1     = float(entry[4])
        except (IndexError, ValueError, TypeError) as e:
            logger.warning("Skipping malformed call entry %s: %s", entry, e)
            continue

        chrom  = f"chr{chrom_1based}"
        start  = start_bin * binsize
        end    = (end_bin + 1) * binsize
        ratio  = ratio_m1 + 1.0
        log2r  = float(np.log2(ratio)) if ratio > 0 else float("nan")
        call   = "GAIN" if log2r > 0 else "LOSS"

        sub = df_bins[
            (df_bins["chrom"] == chrom) &
            (df_bins["start"] >= start) &
            (df_bins["end"]   <= end)
        ]
        mad_z_mean = (float(np.nanmean(sub["mad_z"]))
                      if len(sub) and "mad_z" in sub.columns
                      else float("nan"))

        rows.append({
            "chrom":           chrom,
            "start":           start,
            "end":             end,
            "type":            call,
            "mean_log2_ratio": log2r,
            "mean_z":          zscore,
            "mean_mad_z":      mad_z_mean,
            "n_bins":          len(sub),
        })

    if not rows:
        return pd.DataFrame()

    df = pd.DataFrame(rows)
    df["_ci"] = df["chrom"].map(CHROM_ORDER).fillna(99).astype(int)
    df = df.sort_values(["_ci", "start"]).drop(columns="_ci").reset_index(drop=True)
    return df


# ── TSV writers (same interface as gxcnv2) ────────────────────────────────────

def _write(df: pd.DataFrame, path: str) -> None:
    with open(path, "w") as f:
        f.write("#" + "\t".join(str(c) for c in df.columns) + "\n")
        df.to_csv(f, sep="\t", index=False, header=False, float_format="%.6g", na_rep="NA")
    logger.info("Written: %s (%d rows)", path, len(df))


def write_bins(df: pd.DataFrame, prefix: str) -> None:
    cols  = ["chrom", "start", "end", "log2_ratio", "z_score", "mad_z"]
    avail = [c for c in cols if c in df.columns]
    _write(df[avail].copy(), f"{prefix}_bins.tsv")


def write_segments(df_calls: pd.DataFrame | None, prefix: str) -> None:
    """
    For gxcnv1, segments = calls.
    Wisecondor does not output separate CBS segments; the Stouffer-z calls
    already represent the segmented result.
    """
    cols = ["chrom", "start", "end", "n_bins", "mean_log2_ratio", "mean_z", "mean_mad_z"]
    if df_calls is None or df_calls.empty:
        with open(f"{prefix}_segments.tsv", "w") as f:
            f.write("#" + "\t".join(cols) + "\n")
        return
    avail = [c for c in cols if c in df_calls.columns]
    _write(df_calls[avail].copy(), f"{prefix}_segments.tsv")


def write_calls(df_calls: pd.DataFrame | None, prefix: str) -> None:
    cols  = ["chrom", "start", "end", "type",
             "mean_log2_ratio", "mean_z", "mean_mad_z", "n_bins"]
    if df_calls is None or df_calls.empty:
        with open(f"{prefix}_calls.tsv", "w") as f:
            f.write("#" + "\t".join(cols) + "\n")
        return
    avail = [c for c in cols if c in df_calls.columns]
    _write(df_calls[avail].copy(), f"{prefix}_calls.tsv")


def write_qcmetrics(metrics: dict, prefix: str) -> None:
    with open(f"{prefix}_qcmetrics.tsv", "w") as f:
        f.write("#metric\tvalue\n")
        for k, v in metrics.items():
            f.write(f"{k}\t{v}\n")
    logger.info("Written: %s_qcmetrics.tsv", prefix)


# ── Main ───────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(
        description="gxcnv1 predict — Wisecondor-based CNV calling with MAD z-score & MAPD"
    )
    p.add_argument("sample_npz", help="wisecondor.py convert output NPZ")
    p.add_argument("out_npz",    help="wisecondor.py test   output NPZ")
    p.add_argument("-o", "--prefix", required=True, help="Output file prefix")
    p.add_argument("--zscore",   type=float, default=5.5,
                   help="Z-score threshold label for qcmetrics (default 5.5)")
    return p.parse_args()


def main():
    args = parse_args()
    logger.info("gxcnv1 predict | sample=%s  out=%s", args.sample_npz, args.out_npz)

    wc = load_wisecondor_result(args.out_npz)
    if wc is None:
        logger.error("Failed to load Wisecondor results — aborting")
        sys.exit(1)

    df_bins = build_bins_df(wc)
    if df_bins.empty:
        logger.error("No bin data extracted from out.npz — aborting")
        sys.exit(1)

    logger.info("Loaded %d bins across %d chromosomes",
                len(df_bins), df_bins["chrom"].nunique())

    # Enhancements
    df_bins  = add_mad_z(df_bins, col="log2_ratio")
    mapd     = compute_mapd(df_bins, col="log2_ratio")
    df_calls = build_calls_df(wc, df_bins)

    logger.info("MAPD = %.4f | calls = %d",
                mapd if not np.isnan(mapd) else -1, len(df_calls))

    # QC metrics
    qc = {
        "n_bins":            len(df_bins),
        "mapd":              f"{mapd:.5g}" if not np.isnan(mapd) else "NA",
        "median_log2_ratio": f"{float(np.nanmedian(df_bins['log2_ratio'])):.5g}",
        "median_z_score":    f"{float(np.nanmedian(df_bins['z_score'])):.5g}",
        "n_calls":           len(df_calls),
        "zscore_threshold":  args.zscore,
        "wc_threshold_z":    f"{wc['threshold_z']:.4f}",
        "reference":         os.path.basename(args.out_npz),
    }

    # Write outputs
    write_bins(df_bins, args.prefix)
    write_segments(df_calls, args.prefix)
    write_calls(df_calls, args.prefix)
    write_qcmetrics(qc, args.prefix)

    logger.info(
        "gxcnv1 predict complete: %d bins / %d calls",
        len(df_bins), len(df_calls),
    )


if __name__ == "__main__":
    main()
