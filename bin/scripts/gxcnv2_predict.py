#!/usr/bin/env python3
"""
gxcnv2_predict.py  —  WisecondorX-core CNV prediction with gxcnv2 enhancements.

Core algorithm: WisecondorX predict (identical pipeline, same reference.npz).
Enhancements added on top:
  • MAD-based robust z-score (mad_z column)  — more outlier-resistant than WCX std z-score
  • MAPD noise metric (Median Absolute Pairwise Difference)
  • Per-bin log2(ratio) preserved in output
  • Different TSV output schema (separate bins / segments / calls / qcmetrics)

Two execution modes:
  NPZ mode (original):
    Runs WisecondorX predict as a subprocess, produces BED files, then annotates.
    gxcnv2_predict.py  sample.npz  reference.npz  -o PREFIX  [options]

  Beds mode (preferred — reuses RUN_WCX output, avoids duplicate computation):
    Reads pre-computed WCX BED files produced by RUN_WCX. WC result = gxcnv1, WCX result = gxcnv2.
    gxcnv2_predict.py  --bins-bed PREFIX_bins.bed  --segments-bed PREFIX_segments.bed
                       --aberrations-bed PREFIX_aberrations.bed  -o OUTPUT_PREFIX
"""

import argparse
import os
import subprocess
import sys
import tempfile
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)

import logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [gxcnv2] %(levelname)s %(message)s",
)
logger = logging.getLogger("gxcnv2")

CHROM_ORDER = {f"chr{i}": i for i in range(1, 23)}
CHROM_ORDER["chrX"] = 23
CHROM_ORDER["chrY"] = 24


# ── Run WCX predict ────────────────────────────────────────────────────────────

def run_wcx_predict(sample_npz: str, ref_npz: str, outid: str,
                    alpha: float, zscore: float, seed: int,
                    maskrepeats: int, minrefbins: int,
                    gender: str | None = None) -> bool:
    """
    Run WisecondorX predict as a subprocess, producing BED files.
    Returns True on success, False on failure.
    """
    wcx = _find_wisecondorx()
    if wcx is None:
        logger.error("WisecondorX binary not found on PATH or in conda env")
        return False

    cmd = [
        wcx, "predict",
        sample_npz,
        ref_npz,
        outid,
        "--bed",
        "--alpha",      str(alpha),
        "--zscore",     str(zscore),
        "--seed",       str(seed),
        "--maskrepeats", str(maskrepeats),
        "--minrefbins", str(minrefbins),
    ]
    # NOTE: Do NOT pass --gender to WisecondorX predict.
    # The gender is used only to select the correct reference NPZ (caller's responsibility).
    # Passing --gender triggers WCX sex-chromosome mode which outputs ONLY sex chromosomes.

    logger.info("Running: %s", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        logger.error("WisecondorX predict failed:\n%s", result.stderr[-2000:])
        return False
    logger.info("WisecondorX predict done:\n%s", result.stdout[-500:] or result.stderr[-500:])
    return True


def _find_wisecondorx() -> str | None:
    """Find WisecondorX binary (conda env or PATH)."""
    import shutil
    # Prefer conda nipt env
    candidates = [
        "/opt/conda/envs/nipt/bin/WisecondorX",
        "/opt/conda/envs/nipt/bin/wisecondorx",
    ]
    for c in candidates:
        if os.path.isfile(c):
            return c
    return shutil.which("WisecondorX") or shutil.which("wisecondorx")


# ── Parse WCX BED outputs ─────────────────────────────────────────────────────

def _parse_bed(path: str, numeric_cols: list[str]) -> pd.DataFrame | None:
    """Parse a tab-separated BED file, ignoring comment/header lines."""
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        return None
    try:
        df = pd.read_csv(path, sep="\t", comment="#",
                         header=0, dtype=str, engine="python",
                         on_bad_lines="skip")
    except Exception as e:
        logger.warning("Could not parse %s: %s", path, e)
        return None
    if df.empty:
        return None
    # Normalise column names: strip whitespace, lower
    df.columns = [c.strip().lower() for c in df.columns]
    for col in numeric_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def _normalise_chrom(df: pd.DataFrame) -> pd.DataFrame:
    """Ensure chrom column has 'chr' prefix."""
    if "chrom" in df.columns:
        df = df.copy()
        df["chrom"] = df["chrom"].astype(str).str.strip()
        df["chrom"] = df["chrom"].apply(
            lambda c: c if c.startswith("chr") else f"chr{c}"
        )
    return df


def _ratio_to_log2(df: pd.DataFrame) -> pd.DataFrame:
    """
    Convert WCX 'ratio' column to log2_ratio.

    WCX ratio is a fractional deviation from diploid (0-centred, not 1-centred):
      ratio ≈ 0   → diploid (2 copies)
      ratio ≈ +0.5 → trisomy-like
      ratio ≈ -0.5 → monosomy-like

    log2_ratio = log2(1 + ratio), which gives:
      0     → 0.0  (diploid)
      +0.5  → +0.585 (trisomy)
      -0.5  → -1.0   (monosomy)
    """
    if "ratio" in df.columns:
        ratio = pd.to_numeric(df["ratio"], errors="coerce")
        df = df.copy()
        df["log2_ratio"] = np.where(
            (1 + ratio) > 0, np.log2(1 + ratio), float("nan")
        )
    return df


def load_wcx_bins(outid: str) -> pd.DataFrame | None:
    """Load WCX _bins.bed output.

    Actual columns: chr  start  end  id  ratio  zscore
    """
    path = f"{outid}_bins.bed"
    df = _parse_bed(path, ["start", "end", "ratio", "zscore"])
    if df is None:
        return None
    rename = {"#chr": "chrom", "chr": "chrom", "zscore": "z_score"}
    df = df.rename(columns={k: v for k, v in rename.items() if k in df.columns})
    df = _normalise_chrom(df)
    df = _ratio_to_log2(df)
    return df


def load_wcx_segments(outid: str) -> pd.DataFrame | None:
    """Load WCX _segments.bed output.

    Actual columns: chr  start  end  ratio  zscore
    """
    path = f"{outid}_segments.bed"
    df = _parse_bed(path, ["start", "end", "ratio", "zscore"])
    if df is None:
        return None
    rename = {"#chr": "chrom", "chr": "chrom", "zscore": "z_score"}
    df = df.rename(columns={k: v for k, v in rename.items() if k in df.columns})
    df = _normalise_chrom(df)
    df = _ratio_to_log2(df)
    return df


def load_wcx_aberrations(outid: str) -> pd.DataFrame | None:
    """Load WCX _aberrations.bed output.

    Actual columns: chr  start  end  ratio  zscore  type
    """
    path = f"{outid}_aberrations.bed"
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        return pd.DataFrame()
    df = _parse_bed(path, ["start", "end", "ratio", "zscore"])
    if df is None or df.empty:
        return pd.DataFrame()
    rename = {"#chr": "chrom", "chr": "chrom", "zscore": "z_score"}
    df = df.rename(columns={k: v for k, v in rename.items() if k in df.columns})
    df = _normalise_chrom(df)
    df = _ratio_to_log2(df)
    return df


# ── Direct-path loaders (beds mode: reuse RUN_WCX output) ────────────────────

def _load_wcx_bed_file(path: str) -> pd.DataFrame | None:
    """Load any WCX BED file from a direct file path."""
    df = _parse_bed(path, ["start", "end", "ratio", "zscore"])
    if df is None:
        return None
    rename = {"#chr": "chrom", "chr": "chrom", "zscore": "z_score"}
    df = df.rename(columns={k: v for k, v in rename.items() if k in df.columns})
    df = _normalise_chrom(df)
    df = _ratio_to_log2(df)
    return df


def load_wcx_bins_file(path: str) -> pd.DataFrame | None:
    """Load WCX _bins.bed from a direct file path (beds mode)."""
    return _load_wcx_bed_file(path)


def load_wcx_segments_file(path: str) -> pd.DataFrame | None:
    """Load WCX _segments.bed from a direct file path (beds mode)."""
    return _load_wcx_bed_file(path)


def load_wcx_aberrations_file(path: str) -> pd.DataFrame:
    """Load WCX _aberrations.bed from a direct file path (beds mode)."""
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        return pd.DataFrame()
    df = _load_wcx_bed_file(path)
    return df if df is not None else pd.DataFrame()


# ── Enhancements ───────────────────────────────────────────────────────────────

def add_mad_z(df: pd.DataFrame, col: str = "log2_ratio") -> pd.DataFrame:
    """
    Add MAD-based robust z-score (mad_z) to the bin table.

    More resistant to outlier bins than the WCX within-sample std z-score.
    """
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
    """
    MAPD — Median Absolute Pairwise Difference.

    A reference-free noise metric. Values >0.35 indicate a noisy sample.
    """
    vals = df[col].dropna().values
    if len(vals) < 2:
        return float("nan")
    return float(np.median(np.abs(np.diff(vals))))


def enrich_segments(df_bins: pd.DataFrame, df_segs: pd.DataFrame) -> pd.DataFrame:
    """
    Add mean_mad_z to each CBS segment by joining with the bin table.
    """
    if df_segs is None or df_segs.empty:
        return df_segs
    rows = []
    for _, seg in df_segs.iterrows():
        sub = df_bins[
            (df_bins["chrom"] == seg["chrom"]) &
            (df_bins["start"] >= seg["start"]) &
            (df_bins["end"]   <= seg["end"])
        ]
        mean_madz = float(np.nanmean(sub["mad_z"])) if len(sub) and "mad_z" in sub.columns else float("nan")
        rows.append({**seg.to_dict(), "mean_mad_z": mean_madz, "n_bins": len(sub)})
    return pd.DataFrame(rows)


def build_calls(df_aber: pd.DataFrame, df_segs: pd.DataFrame | None) -> pd.DataFrame:
    """
    Convert WCX aberrations.bed + enriched segments into gxcnv2 calls TSV.
    The type (GAIN / LOSS) is inferred from the segment mean log2_ratio sign.
    """
    if df_aber is None or df_aber.empty:
        return pd.DataFrame()
    rows = []
    for _, row in df_aber.iterrows():
        chrom = row["chrom"]
        start = int(float(row.get("start", 0)))
        end   = int(float(row.get("end",   0)))
        # Match to segment for log2_ratio and z_score
        seg_lr, seg_z, seg_madz, n_bins = float("nan"), float("nan"), float("nan"), 0
        if df_segs is not None and not df_segs.empty:
            match = df_segs[
                (df_segs["chrom"] == chrom) &
                (df_segs["start"] <= start) &
                (df_segs["end"]   >= end)
            ]
            if not match.empty:
                m = match.iloc[0]
                seg_lr   = float(m.get("log2_ratio", float("nan")))
                seg_z    = float(m.get("z_score",    float("nan")))
                seg_madz = float(m.get("mean_mad_z", float("nan")))
                n_bins   = int(m.get("n_bins", 0))
        # Prefer WCX type if available; else infer from log2_ratio sign
        wcx_type = str(row.get("type", "")).strip().upper()
        if wcx_type in ("GAIN", "LOSS"):
            call = wcx_type
        elif wcx_type in ("DUPLICATION", "DUP"):
            call = "GAIN"
        elif wcx_type in ("DELETION", "DEL"):
            call = "LOSS"
        elif np.isfinite(seg_lr):
            call = "GAIN" if seg_lr > 0 else "LOSS"
        else:
            call = "LOSS"
        rows.append({
            "chrom":           chrom,
            "start":           start,
            "end":             end,
            "type":            call,
            "mean_log2_ratio": seg_lr,
            "mean_z":          seg_z,
            "mean_mad_z":      seg_madz,
            "n_bins":          n_bins,
        })
    return pd.DataFrame(rows) if rows else pd.DataFrame()


# ── TSV writers ────────────────────────────────────────────────────────────────

def _sort(df: pd.DataFrame) -> pd.DataFrame:
    if "chrom" not in df.columns:
        return df
    return df.sort_values(
        ["chrom", "start"],
        key=lambda s: s.map(lambda v: CHROM_ORDER.get(v, 99)) if s.name == "chrom" else s,
    )


def _write(df: pd.DataFrame, path: str) -> None:
    with open(path, "w") as f:
        f.write("#" + "\t".join(str(c) for c in df.columns) + "\n")
        df.to_csv(f, sep="\t", index=False, header=False, float_format="%.6g",
                  na_rep="NA")
    logger.info("Written: %s (%d rows)", path, len(df))


def write_bins(df: pd.DataFrame, prefix: str) -> None:
    cols = ["chrom", "start", "end", "log2_ratio", "z_score", "mad_z"]
    avail = [c for c in cols if c in df.columns]
    _write(_sort(df[avail].copy()), f"{prefix}_bins.tsv")


def write_segments(df: pd.DataFrame | None, prefix: str) -> None:
    cols = ["chrom", "start", "end", "n_bins", "log2_ratio", "z_score", "mean_mad_z"]
    if df is None or df.empty:
        with open(f"{prefix}_segments.tsv", "w") as f:
            f.write("#" + "\t".join(cols) + "\n")
        return
    avail = [c for c in cols if c in df.columns]
    _write(_sort(df[avail].copy()), f"{prefix}_segments.tsv")


def write_calls(df: pd.DataFrame | None, prefix: str) -> None:
    cols = ["chrom", "start", "end", "type",
            "mean_log2_ratio", "mean_z", "mean_mad_z", "n_bins"]
    if df is None or df.empty:
        with open(f"{prefix}_calls.tsv", "w") as f:
            f.write("#" + "\t".join(cols) + "\n")
        return
    avail = [c for c in cols if c in df.columns]
    _write(_sort(df[avail].copy()), f"{prefix}_calls.tsv")


def write_qcmetrics(metrics: dict, prefix: str) -> None:
    with open(f"{prefix}_qcmetrics.tsv", "w") as f:
        f.write("#metric\tvalue\n")
        for k, v in metrics.items():
            f.write(f"{k}\t{v}\n")
    logger.info("Written: %s_qcmetrics.tsv", prefix)


# ── Main ───────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(
        description="gxcnv2 predict — WisecondorX-based CNV calling with MAD z-score & MAPD"
    )
    # NPZ mode (original)
    p.add_argument("sample_npz", nargs="?", default=None,
                   help="WisecondorX convert output NPZ (omit when using --bins-bed)")
    p.add_argument("ref_npz",    nargs="?", default=None,
                   help="WisecondorX reference NPZ (omit when using --bins-bed)")
    p.add_argument("-o", "--prefix",  required=True, help="Output file prefix")
    p.add_argument("--gender",        default=None, choices=["M", "F"],
                   help="Force gender (M/F); omit to let WCX auto-detect (NPZ mode only)")
    p.add_argument("--zscore",        type=float, default=6.0,
                   help="WCX z-score aberration threshold (default 6.0; NPZ mode only)")
    p.add_argument("--alpha",         type=float, default=0.01)
    p.add_argument("--seed",          type=int,   default=100)
    p.add_argument("--maskrepeats",   type=int,   default=5)
    p.add_argument("--minrefbins",    type=int,   default=150)
    # Beds mode: reuse pre-computed WCX BED files from RUN_WCX
    p.add_argument("--bins-bed",         default=None,
                   help="Pre-computed WCX _bins.bed (skips WCX predict subprocess)")
    p.add_argument("--segments-bed",     default=None,
                   help="Pre-computed WCX _segments.bed")
    p.add_argument("--aberrations-bed",  default=None,
                   help="Pre-computed WCX _aberrations.bed")
    return p.parse_args()


def _write_empty_outputs(prefix: str, reason: str) -> None:
    """Write empty stub TSV files when sample has no usable data."""
    logger.warning("Writing empty outputs for %s: %s", prefix, reason)
    for suffix in ("_bins.tsv", "_segments.tsv", "_calls.tsv"):
        open(f"{prefix}{suffix}", "w").close()
    write_qcmetrics({"status": reason, "n_bins": 0, "n_calls": 0}, prefix)


def main():
    args = parse_args()

    beds_mode = bool(args.bins_bed)

    if beds_mode:
        # ── Beds mode: reuse RUN_WCX output (WCX result = gxcnv2) ────────────
        logger.info("gxcnv2 predict | beds-mode | bins=%s", args.bins_bed)
        df_bins = load_wcx_bins_file(args.bins_bed)
        df_segs = load_wcx_segments_file(args.segments_bed) if args.segments_bed else None
        df_aber = load_wcx_aberrations_file(args.aberrations_bed) if args.aberrations_bed else pd.DataFrame()
        ref_name = "from_wcx"

        if df_bins is None or df_bins.empty:
            _write_empty_outputs(args.prefix, "empty_bins_bed")
            return

    else:
        # ── NPZ mode: run WCX predict as subprocess (original behaviour) ─────
        if not args.sample_npz or not args.ref_npz:
            logger.error("Provide either --bins-bed OR both sample_npz and ref_npz positional args")
            sys.exit(1)
        logger.info("gxcnv2 predict | npz-mode | sample=%s  ref=%s", args.sample_npz, args.ref_npz)
        ref_name = os.path.basename(args.ref_npz)

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_outid = os.path.join(tmpdir, "wcx_out")

            ok = run_wcx_predict(
                args.sample_npz, args.ref_npz, tmp_outid,
                alpha=args.alpha, zscore=args.zscore, seed=args.seed,
                maskrepeats=args.maskrepeats, minrefbins=args.minrefbins,
                gender=args.gender,
            )
            if not ok:
                logger.error("WisecondorX predict failed — aborting")
                sys.exit(1)

            df_bins = load_wcx_bins(tmp_outid)
            df_segs = load_wcx_segments(tmp_outid)
            df_aber = load_wcx_aberrations(tmp_outid)

        if df_bins is None or df_bins.empty:
            logger.error("No bins output from WisecondorX predict — aborting")
            sys.exit(1)

    logger.info("Loaded %d bins, %d segments",
                len(df_bins), len(df_segs) if df_segs is not None else 0)

    # ── Enhancements (shared by both modes) ───────────────────────────────────
    df_bins  = add_mad_z(df_bins, col="log2_ratio")
    mapd     = compute_mapd(df_bins, col="log2_ratio")

    df_segs  = enrich_segments(df_bins, df_segs)
    df_calls = build_calls(df_aber, df_segs)

    logger.info("MAPD = %.4f | calls = %d", mapd if not np.isnan(mapd) else -1, len(df_calls))

    # ── QC metrics ────────────────────────────────────────────────────────────
    zscore_thresh = args.zscore if not beds_mode else "from_wcx"
    qc = {
        "n_bins":            len(df_bins),
        "mapd":              f"{mapd:.5g}" if not np.isnan(mapd) else "NA",
        "median_log2_ratio": f"{float(np.nanmedian(df_bins['log2_ratio'])):.5g}",
        "median_z_score":    f"{float(np.nanmedian(df_bins['z_score'])):.5g}" if "z_score" in df_bins.columns else "NA",
        "n_segments":        len(df_segs) if df_segs is not None else 0,
        "n_calls":           len(df_calls),
        "zscore_threshold":  zscore_thresh,
        "reference":         ref_name,
    }

    # ── Write outputs ─────────────────────────────────────────────────────────
    write_bins(df_bins, args.prefix)
    write_segments(df_segs, args.prefix)
    write_calls(df_calls, args.prefix)
    write_qcmetrics(qc, args.prefix)

    logger.info(
        "gxcnv2 predict complete: %d bins / %d segments / %d calls",
        len(df_bins),
        len(df_segs) if df_segs is not None else 0,
        len(df_calls),
    )


if __name__ == "__main__":
    main()
