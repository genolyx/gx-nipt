#!/usr/bin/env python3
"""
build_gxcnv_reference.py
========================
End-to-end gx-cnv reference panel builder for hg19 NIPT samples.

Pipeline steps (all in one run):
  1. Read a sample data TSV (from ken-nipt / gx-nipt analysis outputs).
  2. Filter samples by QC criteria and fetal sex.
  3. For each sample: convert BAM → NPZ (gxcnv convert, parallel).
  4. Inject GC fractions from HMMcopy 50kb wig → 100kb NPZ bins.
  5. Apply polynomial GC correction.
  6. Build reference panel (gxcnv_build_ref_hg19.py / Python API).

Usage
-----
    python build_gxcnv_reference.py \\
        --sample-tsv  /path/to/reference_sample_list_from_json_v3.tsv \\
        --labcode     cordlife \\
        --sex         female \\
        --out-dir     /home/ken/gx-nipt/refs/labs/cordlife/GXCNV/female \\
        --npz-cache   /home/ken/gx-nipt/refs/gxcnv/npz_female \\
        --n-samples   100 \\
        --threads     8

    # Build both female and male references:
    for sex in female male; do
      python build_gxcnv_reference.py --sex $sex --labcode cordlife ...
    done
"""

import argparse
import glob
import logging
import os
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [build_gxcnv] %(levelname)s %(message)s",
)
logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Column name constants (matching reference_sample_list_from_json_v3.tsv)
# ---------------------------------------------------------------------------
COL_SAMPLE_ID  = "sample_id"
COL_GENDER     = "fetal_gender(gd_2)"
COL_RESULT     = "Result"
COL_DUP_RATE   = "duplication_rate(%)"
COL_MAP_RATE   = "mapping_rate(%)"
COL_YFF2       = "YFF_2"
COL_SAMPLE_DIR = "sample_dir"
COL_MONTH      = "month"

# HMMcopy wig path template relative to sample_dir
WIG_TEMPLATE = "Output_hmmcopy/{sample_id}.of_fetus.50kb.wig.Normalization.txt"


# ---------------------------------------------------------------------------
# QC filter
# ---------------------------------------------------------------------------

def filter_samples(df: pd.DataFrame, sex: str,
                   max_dup: float, min_map: float,
                   min_ff: float, n_samples: int,
                   balanced_months: bool) -> pd.DataFrame:
    """Select QC-passing samples of the requested sex."""
    for col in [COL_DUP_RATE, COL_MAP_RATE, COL_YFF2]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df["month_int"] = pd.to_numeric(df[COL_MONTH], errors="coerce")

    sex_code = "XY" if sex.lower() in ("male", "xy", "m") else "XX"

    mask = (
        (df[COL_GENDER] == sex_code) &
        (df[COL_RESULT] == "Low Risk") &
        (df[COL_DUP_RATE] < max_dup) &
        (df[COL_MAP_RATE] > min_map)
    )
    if sex_code == "XY":
        mask &= df[COL_YFF2].notna() & (df[COL_YFF2] >= min_ff)

    candidates = df[mask].copy()
    logger.info("QC-passing %s samples: %d", sex, len(candidates))

    if balanced_months:
        # Take up to 15 per month, newest months first
        selected = (
            candidates
            .sort_values("month_int", ascending=False)
            .groupby(COL_MONTH, group_keys=False)
            .apply(lambda g: g.head(15))
            .reset_index(drop=True)
            .sort_values("month_int", ascending=False)
            .head(n_samples)
        )
    else:
        selected = candidates.sort_values("month_int", ascending=False).head(n_samples)

    logger.info("Selected %d samples (balanced=%s)", len(selected), balanced_months)
    return selected


# ---------------------------------------------------------------------------
# BAM → NPZ conversion (single sample)
# ---------------------------------------------------------------------------

def convert_bam_to_npz(sample_id: str, bam_path: str,
                       npz_path: str, bin_size: int = 100_000,
                       min_mapq: int = 20) -> str | None:
    """Run gxcnv convert for one sample. Returns npz_path on success."""
    if os.path.exists(npz_path):
        logger.debug("NPZ already exists, skipping: %s", sample_id)
        return npz_path
    try:
        from gxcnv.convert import bam_to_npz
        bam_to_npz(bam_path, npz_path, bin_size=bin_size, min_mapq=min_mapq)
        return npz_path
    except Exception as e:
        logger.error("Convert failed for %s: %s", sample_id, e)
        return None


# ---------------------------------------------------------------------------
# GC injection (mirrors gxcnv_inject_gc.py logic)
# ---------------------------------------------------------------------------

def load_gc_from_wig(wig_path: str) -> dict:
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
    df["start_100k"] = (df["start"].astype(int) // 100_000) * 100_000
    gc_map = df.groupby(["chr", "start_100k"])["gc"].mean().reset_index()
    return {(row.chr, int(row.start_100k)): float(row.gc)
            for row in gc_map.itertuples(index=False)}


def gc_correct(counts: np.ndarray, gc_fractions: np.ndarray,
               mask: np.ndarray, poly_degree: int = 3) -> np.ndarray:
    corrected = np.full(len(counts), np.nan)
    valid = mask & (counts > 0) & np.isfinite(gc_fractions)
    if valid.sum() < poly_degree + 1:
        usable = mask & (counts > 0)
        corrected[usable] = counts[usable].astype(float)
        return corrected
    gc_bins    = np.linspace(0, 1, 101)
    gc_idx     = np.clip(np.digitize(gc_fractions[valid], gc_bins) - 1, 0, 99)
    gc_centers = (gc_bins[:-1] + gc_bins[1:]) / 2
    gc_medians = np.full(100, np.nan)
    for i in range(100):
        in_bin = valid.nonzero()[0][gc_idx == i]
        if len(in_bin) >= 3:
            gc_medians[i] = np.median(counts[in_bin])
    finite = np.isfinite(gc_medians)
    if finite.sum() < poly_degree + 1:
        corrected[valid] = counts[valid].astype(float)
        return corrected
    coeffs    = np.polyfit(gc_centers[finite], gc_medians[finite], poly_degree)
    predicted = np.polyval(coeffs, gc_fractions[valid])
    predicted = np.where(predicted <= 0, np.nan, predicted)
    global_median = np.nanmedian(gc_medians[finite])
    adj = counts[valid].astype(float) / predicted * global_median
    corrected[valid] = np.where(np.isfinite(adj), adj, np.nan)
    return corrected


def inject_gc_and_correct(npz_path: str, wig_path: str) -> bool:
    gc_map = load_gc_from_wig(wig_path)
    if not gc_map:
        return False
    d = dict(np.load(npz_path, allow_pickle=True))
    bins   = d["bins"]
    chroms = list(d["chroms"])
    counts = d["counts"].astype(float)
    mask   = d["mask"].astype(bool)

    gc_fractions = np.full(len(bins), np.nan, dtype=np.float32)
    hit = 0
    for i, (chrom_idx, start, end, _) in enumerate(bins):
        key = (chroms[int(chrom_idx)], int(start))
        if key in gc_map:
            gc_fractions[i] = gc_map[key]
            hit += 1

    if hit / max(len(bins), 1) < 0.3:
        logger.warning("Low GC hit rate for %s — skipping GC injection", npz_path)
        return False

    bins[:, 3] = gc_fractions
    d["bins"]  = bins

    corrected = gc_correct(counts, gc_fractions, mask)
    valid_sum = np.nansum(corrected[mask])
    if valid_sum > 0:
        corrected = corrected / valid_sum * mask.sum()
    d["corrected"] = corrected

    base = npz_path[:-4] if npz_path.endswith(".npz") else npz_path
    np.savez_compressed(base, **d)
    return True


# ---------------------------------------------------------------------------
# Reference panel build
# ---------------------------------------------------------------------------

HG19_TARGET_REGIONS = [
    ("DiGeorge_22q11",       "chr22", 18_648_855, 21_802_889),
    ("Williams_7q11",        "chr7",  72_744_454, 74_142_513),
    ("Angelman_15q11",       "chr15", 22_765_628, 28_530_576),
    ("PraderWilli_15q11",    "chr15", 22_765_628, 28_530_576),
    ("Wolf_4p16",            "chr4",       1_000,  4_000_000),
    ("CriDuChat_5p15",       "chr5",       1_000, 11_800_000),
    ("NF1_17q11",            "chr17", 29_421_945, 30_144_979),
    ("Smith_17p11",          "chr17", 16_755_133, 20_519_016),
    ("Langer_Giedion_8q24",  "chr8", 116_600_000, 119_400_000),
    ("Miller_Dieker_17p13",  "chr17",      1_000,  4_000_000),
    ("CHARGE_8q12",          "chr8",  61_720_000,  62_020_000),
    ("Kabuki_12q13",         "chr12", 49_400_000,  50_400_000),
    ("Sotos_5q35",           "chr5", 175_100_000, 177_200_000),
    ("Rubinstein_16p13",     "chr16",  3_800_000,   4_000_000),
    ("Potocki_Lupski_17p11", "chr17", 15_000_000,  20_500_000),
]


def build_reference_panel(npz_paths: list, output_path: str,
                          pca_variance: float = 0.95) -> None:
    from gxcnv.newref import build_reference
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    build_reference(
        npz_paths=npz_paths,
        output_path=output_path,
        target_regions=HG19_TARGET_REGIONS,
        global_pca_variance=pca_variance,
    )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample-tsv",   required=True,
                    help="Path to reference_sample_list TSV")
    ap.add_argument("--labcode",      required=True,
                    help="Lab code (e.g. cordlife)")
    ap.add_argument("--sex",          required=True,
                    choices=["female","male","XX","XY","F","M"],
                    help="Fetal sex of reference samples")
    ap.add_argument("--out-dir",      required=True,
                    help="Output directory for reference.npz")
    ap.add_argument("--npz-cache",    default=None,
                    help="Directory to cache converted NPZ files (default: <out-dir>/npz_cache)")
    ap.add_argument("--n-samples",    type=int, default=100,
                    help="Target number of reference samples (default: 100)")
    ap.add_argument("--balanced-months", action="store_true", default=True,
                    help="Spread selection across months (default: True)")
    ap.add_argument("--max-dup",      type=float, default=15.0,
                    help="Max duplication rate %% (default: 15)")
    ap.add_argument("--min-map",      type=float, default=85.0,
                    help="Min mapping rate %% (default: 85)")
    ap.add_argument("--min-ff",       type=float, default=4.0,
                    help="Min YFF_2 for male samples (default: 4.0)")
    ap.add_argument("--threads",      type=int, default=8,
                    help="Parallel threads for BAM→NPZ conversion (default: 8)")
    ap.add_argument("--pca-variance", type=float, default=0.95,
                    help="PCA variance threshold for reference panel (default: 0.95)")
    ap.add_argument("--bin-size",     type=int, default=100_000,
                    help="Bin size for gxcnv convert (default: 100000)")
    ap.add_argument("--wig-suffix",   default=WIG_TEMPLATE,
                    help="HMMcopy wig path template relative to sample_dir")
    args = ap.parse_args()

    out_dir  = Path(args.out_dir)
    npz_dir  = Path(args.npz_cache) if args.npz_cache else out_dir / "npz_cache"
    npz_dir.mkdir(parents=True, exist_ok=True)
    out_dir.mkdir(parents=True, exist_ok=True)
    ref_path = str(out_dir / "reference.npz")

    # ── Step 1: Load & filter samples ─────────────────────────────────────────
    logger.info("Loading sample TSV: %s", args.sample_tsv)
    df = pd.read_csv(args.sample_tsv, sep="\t", dtype=str)
    df.columns = [c.strip() for c in df.columns]

    selected = filter_samples(
        df, sex=args.sex,
        max_dup=args.max_dup, min_map=args.min_map,
        min_ff=args.min_ff, n_samples=args.n_samples,
        balanced_months=args.balanced_months,
    )

    if selected.empty:
        logger.error("No samples passed QC filters — aborting.")
        sys.exit(1)

    # ── Step 2: BAM → NPZ (parallel) ─────────────────────────────────────────
    logger.info("=== Step 2: BAM → NPZ conversion (%d samples, %d threads) ===",
                len(selected), args.threads)

    convert_tasks = []
    for _, row in selected.iterrows():
        sid     = row[COL_SAMPLE_ID]
        sdir    = str(row.get(COL_SAMPLE_DIR, ""))
        bam     = os.path.join(sdir, f"{sid}.proper_paired.bam")
        npz     = str(npz_dir / f"{sid}.gxcnv.npz")
        convert_tasks.append((sid, bam, npz))

    failed_convert = []
    done = 0
    with ThreadPoolExecutor(max_workers=args.threads) as ex:
        futures = {
            ex.submit(convert_bam_to_npz, sid, bam, npz, args.bin_size): sid
            for sid, bam, npz in convert_tasks
        }
        for fut in as_completed(futures):
            sid = futures[fut]
            result = fut.result()
            done += 1
            if result:
                logger.info("[%d/%d] Converted %s", done, len(convert_tasks), sid)
            else:
                failed_convert.append(sid)
                logger.warning("[%d/%d] FAILED %s", done, len(convert_tasks), sid)

    if failed_convert:
        logger.warning("%d conversions failed: %s", len(failed_convert), failed_convert)

    # ── Step 3: GC injection ──────────────────────────────────────────────────
    logger.info("=== Step 3: GC injection from HMMcopy wig ===")
    gc_ok = 0
    for _, row in selected.iterrows():
        sid  = row[COL_SAMPLE_ID]
        sdir = str(row.get(COL_SAMPLE_DIR, ""))
        npz  = str(npz_dir / f"{sid}.gxcnv.npz")
        if not os.path.exists(npz):
            continue
        wig = os.path.join(sdir, args.wig_suffix.replace("{sample_id}", sid))
        if not os.path.exists(wig):
            logger.warning("Wig not found for %s — skipping GC injection", sid)
            # fall back: raw counts already normalised by patch_npz
            _patch_raw_counts(npz)
            continue
        if inject_gc_and_correct(npz, wig):
            gc_ok += 1
        else:
            _patch_raw_counts(npz)

    logger.info("GC-corrected: %d/%d samples", gc_ok, len(selected))

    # ── Step 4: Build reference panel ─────────────────────────────────────────
    logger.info("=== Step 4: Building reference panel → %s ===", ref_path)
    npz_paths = sorted(glob.glob(str(npz_dir / "*.gxcnv.npz")))
    if not npz_paths:
        logger.error("No NPZ files found in %s", npz_dir)
        sys.exit(1)

    logger.info("Using %d NPZ files for reference panel", len(npz_paths))
    build_reference_panel(npz_paths, ref_path, pca_variance=args.pca_variance)
    logger.info("Reference panel saved: %s", ref_path)
    logger.info("Done.")


def _patch_raw_counts(npz_path: str) -> None:
    """Fallback: use normalised raw counts when GC injection fails."""
    try:
        d = dict(np.load(npz_path, allow_pickle=True))
        if not np.all(np.isnan(d["corrected"])):
            return  # already patched
        counts = d["counts"].astype(float)
        mask   = d["mask"].astype(bool)
        c = np.full(len(counts), np.nan)
        u = mask & (counts > 0)
        c[u] = counts[u]
        s = np.nansum(c[mask])
        if s > 0:
            c = c / s * mask.sum()
        d["corrected"] = c
        base = npz_path[:-4] if npz_path.endswith(".npz") else npz_path
        np.savez_compressed(base, **d)
    except Exception as e:
        logger.warning("Raw-count fallback failed for %s: %s", npz_path, e)


if __name__ == "__main__":
    main()
