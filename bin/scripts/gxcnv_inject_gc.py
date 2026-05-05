#!/usr/bin/env python3
"""
gxcnv_inject_gc.py
==================
Inject GC fractions from HMMcopy 50kb wig files into gxcnv NPZ sample files,
then reapply GC correction.

gxcnv convert does not compute GC fractions from BAM input (it leaves bin_gc as
NaN), so GC correction is skipped for all BAM-derived NPZ files.  This script
patches the NPZ files in-place by:

  1. Loading the per-sample HMMcopy 50kb wig Normalization file.
  2. Aggregating adjacent 50kb bins → 100kb bin GC fractions (mean of valid bins).
  3. Writing the GC fractions into the bins[:, 3] column of the NPZ.
  4. Re-running polynomial GC correction (same algorithm as gxcnv.convert._gc_correct).
  5. Re-normalising corrected counts to reads-per-million equivalent.
  6. Saving the patched NPZ.

Usage
-----
    python gxcnv_inject_gc.py \\
        --sample-list /path/to/sample_bam_list.txt \\
        --npz-dir     /path/to/npz/ \\
        --wig-suffix  "Output_hmmcopy/{sample_id}.of_fetus.50kb.wig.Normalization.txt"

sample_bam_list.txt format (tab-separated, no header):
    SAMPLE_ID   /path/to/sample_dir
"""

import argparse
import logging
import os
import sys

import numpy as np
import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [gc_inject] %(levelname)s %(message)s",
)
logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# GC correction (mirrors gxcnv.convert._gc_correct)
# ---------------------------------------------------------------------------

def gc_correct(counts: np.ndarray, gc_fractions: np.ndarray,
               mask: np.ndarray, poly_degree: int = 3) -> np.ndarray:
    corrected = np.full(len(counts), np.nan)
    valid = mask & (counts > 0) & np.isfinite(gc_fractions)

    if valid.sum() < poly_degree + 1:
        logger.warning("Too few valid bins for GC correction – using raw counts.")
        corrected[valid] = counts[valid].astype(float)
        return corrected

    gc_bins = np.linspace(0, 1, 101)
    gc_idx  = np.digitize(gc_fractions[valid], gc_bins) - 1
    gc_idx  = np.clip(gc_idx, 0, 99)
    gc_centers = (gc_bins[:-1] + gc_bins[1:]) / 2

    gc_medians = np.full(100, np.nan)
    for i in range(100):
        in_bin = valid.nonzero()[0][gc_idx == i]
        if len(in_bin) >= 3:
            gc_medians[i] = np.median(counts[in_bin])

    # Interpolate missing GC bins
    finite_mask = np.isfinite(gc_medians)
    if finite_mask.sum() < poly_degree + 1:
        corrected[valid] = counts[valid].astype(float)
        return corrected

    coeffs    = np.polyfit(gc_centers[finite_mask], gc_medians[finite_mask], poly_degree)
    predicted = np.polyval(coeffs, gc_fractions[valid])
    predicted = np.where(predicted <= 0, np.nan, predicted)

    global_median = np.nanmedian(gc_medians[finite_mask])
    adj = counts[valid].astype(float) / predicted * global_median
    corrected[valid] = np.where(np.isfinite(adj), adj, np.nan)
    return corrected


# ---------------------------------------------------------------------------
# HMMcopy wig → 100 kb GC map
# ---------------------------------------------------------------------------

def load_gc_from_wig(wig_path: str) -> dict:
    """Return dict {(chrom, bin_100kb_start): gc_mean}."""
    try:
        df = pd.read_csv(
            wig_path, sep="\t",
            names=["chr","start","end","reads","gc","map","valid",
                   "ideal","cor.gc","cor.map","copy"],
            header=0, dtype={"chr": str},
        )
    except Exception as e:
        logger.error("Cannot read wig %s: %s", wig_path, e)
        return {}

    df["gc"] = pd.to_numeric(df["gc"], errors="coerce")
    df = df[df["gc"] > 0].copy()

    # 50 kb → 100 kb bucket
    df["start_100k"] = (df["start"].astype(int) // 100_000) * 100_000
    gc_map = (
        df.groupby(["chr", "start_100k"])["gc"]
        .mean()
        .reset_index()
    )
    result = {
        (row.chr, int(row.start_100k)): float(row.gc)
        for row in gc_map.itertuples(index=False)
    }
    return result


# ---------------------------------------------------------------------------
# Main patch routine
# ---------------------------------------------------------------------------

def patch_npz(npz_path: str, gc_map: dict, dry_run: bool = False) -> bool:
    d = dict(np.load(npz_path, allow_pickle=True))
    bins    = d["bins"]          # (N, 4): chrom_idx, start, end, gc
    chroms  = list(d["chroms"])
    counts  = d["counts"].astype(float)
    mask    = d["mask"].astype(bool)

    # Inject GC fractions
    gc_fractions = np.full(len(bins), np.nan, dtype=np.float32)
    hit = 0
    for i, (chrom_idx, start, end, _) in enumerate(bins):
        chrom = chroms[int(chrom_idx)]
        key   = (chrom, int(start))
        if key in gc_map:
            gc_fractions[i] = gc_map[key]
            hit += 1

    frac = hit / max(len(bins), 1)
    logger.debug("  GC hit rate: %d/%d (%.1f%%)", hit, len(bins), 100 * frac)

    if frac < 0.3:
        logger.warning("  Low GC hit rate %.1f%% — skipping patch.", 100 * frac)
        return False

    # Update bins column 3
    bins[:, 3] = gc_fractions
    d["bins"] = bins

    # Re-run GC correction
    corrected = gc_correct(counts, gc_fractions, mask)

    # Normalise
    valid_sum = np.nansum(corrected[mask])
    if valid_sum > 0:
        corrected = corrected / valid_sum * mask.sum()

    d["corrected"] = corrected

    if not dry_run:
        base = npz_path[:-4] if npz_path.endswith(".npz") else npz_path
        np.savez_compressed(base, **d)

    return True


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sample-list", required=True,
                    help="TSV: SAMPLE_ID<tab>/path/to/sample_dir")
    ap.add_argument("--npz-dir", required=True,
                    help="Directory containing *.gxcnv.npz files")
    ap.add_argument("--wig-suffix", default=
                    "Output_hmmcopy/{sample_id}.of_fetus.50kb.wig.Normalization.txt",
                    help="Relative path template for wig file under sample_dir")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    # Load sample list
    sample_dirs = {}
    with open(args.sample_list) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                sample_dirs[parts[0]] = parts[1]

    logger.info("Loaded %d sample dirs", len(sample_dirs))

    # Find NPZ files
    npz_files = sorted(
        os.path.join(args.npz_dir, fn)
        for fn in os.listdir(args.npz_dir)
        if fn.endswith(".npz")
    )
    logger.info("Found %d NPZ files in %s", len(npz_files), args.npz_dir)

    ok = 0
    for npz_path in npz_files:
        sample_id = os.path.basename(npz_path).replace(".gxcnv.npz", "")
        sample_dir = sample_dirs.get(sample_id)

        if not sample_dir:
            logger.warning("No sample_dir for %s — skipping", sample_id)
            continue

        wig_rel = args.wig_suffix.replace("{sample_id}", sample_id)
        wig_path = os.path.join(sample_dir, wig_rel)

        if not os.path.exists(wig_path):
            logger.warning("Wig not found: %s", wig_path)
            continue

        gc_map = load_gc_from_wig(wig_path)
        if not gc_map:
            logger.warning("Empty GC map for %s", sample_id)
            continue

        success = patch_npz(npz_path, gc_map, dry_run=args.dry_run)
        if success:
            ok += 1
            logger.info("[%d/%d] Patched %s", ok, len(npz_files), sample_id)
        else:
            logger.warning("Patch failed for %s", sample_id)

    logger.info("Done: %d/%d NPZ files patched.", ok, len(npz_files))


if __name__ == "__main__":
    main()
