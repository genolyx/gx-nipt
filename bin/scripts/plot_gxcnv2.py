#!/usr/bin/env python3
"""
plot_gxcnv2.py — Genome visualisation for gxcnv2 results.

Visual style is intentionally distinct from both WisecondorX and gxcnv:
  • Primary metric : log2(ratio)  — not Z-score
  • Genome-wide   : filled-area track (not scatter), with ±trisomy guide lines
  • Color scheme  : teal / amber / violet on white background
  • QC plot       : KDE density curve  — not histogram
  • Per-chromosome: compact panels with confidence ribbon + CBS segment line

Outputs:
  {prefix}_genome.png   — genome-wide track plus 4-column chromosome board
  {prefix}_chr{N}.png   — per-chromosome panels (individual files)
  {prefix}_qc.png       — KDE of log2(ratio) distribution
"""

import argparse
import os
import sys
import warnings

import cnv_board

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)

# ── Colour palette ─────────────────────────────────────────────────────────────
C_NORMAL    = "#2EC8D6"   # cyan
C_GAIN      = "#C44BD6"   # purple
C_LOSS      = "#F07D5A"   # salmon
C_FILTERED  = "#C5D0D6"
C_SEG_NORM  = "#607D8B"
C_SEG_GAIN  = "#C44BD6"
C_SEG_LOSS  = "#F07D5A"
C_RIBBON    = "#E7F7F8"
C_GUIDE     = "#9AA5B0"

# Trisomy log2(3/2) ≈ +0.585 and monosomy log2(1/2) ≈ -1.000
LR_TRISOMY   =  0.585
LR_MONOSOMY  = -1.000

DPI_GENOME = 150
DPI_CHR    = 130
DPI_QC     = 130


def _save_png(out: str, dpi: int) -> None:
    try:
        plt.savefig(out, dpi=dpi, bbox_inches="tight", pil_kwargs={"optimize": True})
    except TypeError:
        plt.savefig(out, dpi=dpi, bbox_inches="tight")

CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
CHROM_ORDER = {c: i for i, c in enumerate(CHROMS)}


# ── Data loaders ───────────────────────────────────────────────────────────────

def _load_tsv(path: str) -> pd.DataFrame | None:
    rows, cols = [], None
    try:
        with open(path) as f:
            for line in f:
                line = line.rstrip("\n")
                if line.startswith("##"):
                    continue
                if line.startswith("#"):
                    cols = line.lstrip("#").split("\t")
                    continue
                if cols is None:
                    continue
                rows.append(dict(zip(cols, line.split("\t"))))
    except FileNotFoundError:
        return None
    if not rows:
        return pd.DataFrame()
    return pd.DataFrame(rows)


def load_bins(path: str) -> pd.DataFrame | None:
    df = _load_tsv(path)
    if df is None or df.empty:
        return None
    for col in ("start", "end", "log2_ratio", "z_score", "mad_z"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    df["chrom"] = df["chrom"].str.strip()
    df = df[df["chrom"].isin(CHROM_ORDER)].copy()
    df["chrom_idx"] = df["chrom"].map(CHROM_ORDER)
    df = df.sort_values(["chrom_idx", "start"]).reset_index(drop=True)
    return df


def load_segments(path: str) -> pd.DataFrame | None:
    df = _load_tsv(path)
    if df is None or df.empty:
        return None
    for col in ("start", "end", "mean_log2_ratio", "mean_z"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    df["chrom"] = df["chrom"].str.strip()
    df = df[df["chrom"].isin(CHROM_ORDER)].copy()
    df["chrom_idx"] = df["chrom"].map(CHROM_ORDER)
    return df


def load_calls(path: str) -> pd.DataFrame | None:
    df = _load_tsv(path)
    if df is None or df.empty:
        return pd.DataFrame()
    for col in ("start", "end", "mean_log2_ratio", "mean_z"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    df["chrom"] = df["chrom"].str.strip()
    return df


# ── Genome coordinate helpers ──────────────────────────────────────────────────

def genome_positions(df: pd.DataFrame):
    """
    Compute cumulative genomic x-positions for genome-wide plots.
    Returns (pos_array, chrom_offsets, cumsum_list).
    """
    chrom_lens = df.groupby("chrom_idx")["end"].max().sort_index()
    cumsum = [0]
    for ci in sorted(chrom_lens.index):
        cumsum.append(cumsum[-1] + int(chrom_lens[ci]))
    offsets = {ci: cumsum[i] for i, ci in enumerate(sorted(chrom_lens.index))}
    pos = df["start"].values + df["chrom_idx"].map(offsets).values
    return pos.astype(float), offsets, cumsum


def _bin_color(df: pd.DataFrame, calls: pd.DataFrame | None) -> list[str]:
    """Colour bins by call type (GAIN → violet, LOSS → amber, else teal)."""
    colors = [C_NORMAL if str(f) not in ("nan", "") else C_FILTERED
              for f in df.get("z_score", pd.Series([""] * len(df)))]

    if calls is None or calls.empty:
        return colors

    for _, row in calls.iterrows():
        chrom = row.get("chrom", "")
        s, e  = float(row.get("start", 0)), float(row.get("end", 0))
        typ   = str(row.get("type", ""))
        mask  = (
            (df["chrom"] == chrom) &
            (df["start"] >= s) &
            (df["end"]   <= e)
        )
        c = C_GAIN if typ == "GAIN" else (C_LOSS if typ == "LOSS" else C_NORMAL)
        for i in df.index[mask]:
            colors[i] = c

    return colors


def _seg_color(lr: float) -> str:
    if not np.isfinite(lr):
        return C_SEG_NORM
    if lr > 0.3:
        return C_SEG_GAIN
    if lr < -0.3:
        return C_SEG_LOSS
    return C_SEG_NORM


def _ribbon_mad(df: pd.DataFrame) -> float:
    """Half-width of the ±1.5 MAD confidence ribbon."""
    lr = df["log2_ratio"].dropna().values
    if len(lr) < 2:
        return 0.4
    mad = float(np.median(np.abs(lr - np.median(lr))))
    return max(mad * 1.5, 0.15)


def _clean_lr(lr, lo=-2.5, hi=2.5):
    return np.clip(np.where(np.isfinite(lr), lr, 0), lo, hi)


# ── Genome-wide plot ───────────────────────────────────────────────────────────


def plot_genome(df: pd.DataFrame, calls: pd.DataFrame | None,
                prefix: str, segments: pd.DataFrame | None = None,
                cyto_dict: dict | None = None,
                z_cutoff: float | None = None) -> None:
    """Genome-wide track plus 4-column chromosome board → {prefix}_genome.png."""
    cnv_board.render_board(
        df, calls, prefix, cyto_dict=cyto_dict, log_prefix="plot_gxcnv2",
        dpi=DPI_GENOME, z_cutoff=z_cutoff,
    )


def plot_chromosome(df_chr: pd.DataFrame, calls: pd.DataFrame | None,
                    chrom: str, prefix: str,
                    segments: pd.DataFrame | None = None,
                    cyto_dict: dict | None = None,
                    guides=None) -> None:
    """Single-chromosome panel in the same style as the board."""
    cnv_board.render_chromosome(
        df_chr, calls, chrom, prefix, cyto_dict=cyto_dict,
        log_prefix="plot_gxcnv2", dpi=DPI_CHR, guides=guides,
    )


# ── QC plot — KDE ─────────────────────────────────────────────────────────────

def plot_qc(df: pd.DataFrame, prefix: str) -> None:
    """
    KDE density curve of log2(ratio) distribution.

    KDE emphasises the shape of the distribution (bimodality, heavy tails)
    better than a histogram and looks clearly different from gxcnv's QC plot.
    """
    lr = df["log2_ratio"].dropna().values
    lr = lr[np.isfinite(lr)]
    lr = np.clip(lr, -2.5, 2.5)

    fig, ax = plt.subplots(figsize=(7, 4))
    fig.patch.set_facecolor("white")
    ax.set_facecolor("#FAFBFC")

    if len(lr) >= 4:
        kde = gaussian_kde(lr, bw_method="scott")
        x_grid = np.linspace(-2.5, 2.5, 400)
        y_kde  = kde(x_grid)
        ax.fill_between(x_grid, 0, y_kde, color=C_NORMAL, alpha=0.35, linewidth=0)
        ax.plot(x_grid, y_kde, color=C_NORMAL, lw=2)

        # Shade tails (gain / loss zones)
        ax.fill_between(x_grid, 0, y_kde,
                        where=x_grid >= LR_TRISOMY,
                        color=C_GAIN, alpha=0.35, linewidth=0,
                        label=f"Gain zone (>+{LR_TRISOMY:.2f})")
        ax.fill_between(x_grid, 0, y_kde,
                        where=x_grid <= LR_MONOSOMY,
                        color=C_LOSS, alpha=0.35, linewidth=0,
                        label=f"Loss zone (<{LR_MONOSOMY:.2f})")

    ax.axvline(0,           color=C_GUIDE, lw=1.0)
    ax.axvline(LR_TRISOMY,  color=C_GAIN,  lw=1.0, ls="--", label="Trisomy guide")
    ax.axvline(LR_MONOSOMY, color=C_LOSS,  lw=1.0, ls="--", label="Monosomy guide")

    ax.set_xlabel("log₂(ratio)", fontsize=9)
    ax.set_ylabel("Density", fontsize=9)
    ax.set_title(f"{os.path.basename(prefix)} — log₂(ratio) distribution  [gxcnv2 QC]",
                 fontsize=10)
    ax.legend(fontsize=8, loc="upper right")

    plt.tight_layout(pad=0.5)
    out = f"{prefix}_qc.png"
    _save_png(out, DPI_QC)
    plt.close(fig)
    print(f"[plot_gxcnv2] {out}", flush=True)


# ── Touch helper ──────────────────────────────────────────────────────────────

def _touch(path: str) -> None:
    # Minimal 1x1 PNG so Nextflow process `-s` checks accept skip stubs.
    import base64
    data = base64.b64decode(
        "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR42mP8z8BQDwAEhQGAhKmMIQAAAABJRU5ErkJggg=="
    )
    with open(path, "wb") as fh:
        fh.write(data)


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(
        description="Plot gxcnv2 results (log2-ratio style)"
    )
    ap.add_argument("--bins",     required=True,  help="*_bins.tsv from gxcnv2_predict.py")
    ap.add_argument("--calls",    required=True,  help="*_calls.tsv from gxcnv2_predict.py")
    ap.add_argument("--segments", default=None,   help="*_segments.tsv (auto-detected if omitted)")
    ap.add_argument("-o", "--prefix", required=True)
    ap.add_argument("--chromosomes", nargs="*", default=None,
                    help="Only plot these chromosomes (default: all)")
    args = ap.parse_args()

    df_bins = load_bins(args.bins)
    if df_bins is None or df_bins.empty:
        print("[plot_gxcnv2] Empty bins file — creating stub outputs", flush=True)
        for suf in ("_genome.png", "_qc.png"):
            _touch(f"{args.prefix}{suf}")
        return

    calls = load_calls(args.calls)

    seg_path = args.segments
    if seg_path is None:
        candidate = args.bins.replace("_bins.tsv", "_segments.tsv")
        if os.path.isfile(candidate):
            seg_path = candidate
    segments = load_segments(seg_path) if seg_path and os.path.isfile(seg_path) else None
    if segments is not None:
        print(f"[plot_gxcnv2] Loaded {len(segments)} CBS segments", flush=True)

    cyto_path = cnv_board.resolve_cytoband(getattr(args, "cytoband", None))
    cyto_dict = cnv_board.load_cytobands(cyto_path) if cyto_path else None
    if cyto_dict:
        print(f"[plot_gxcnv2] Loaded cytobands from {cyto_path}", flush=True)

    z_cutoff = cnv_board.read_z_cutoff(args.bins)
    guides = cnv_board.make_guides(df_bins, z_cutoff)
    if z_cutoff is None:
        print("[plot_gxcnv2] No z cutoff found — threshold lines omitted", flush=True)

    # Genome-wide board (top track + 4-column chromosomes)
    plot_genome(df_bins, calls, args.prefix, segments=segments,
                cyto_dict=cyto_dict, z_cutoff=z_cutoff)

    # Per-chromosome
    chrom_list = args.chromosomes or sorted(
        df_bins["chrom"].unique(), key=lambda c: CHROM_ORDER.get(c, 99)
    )
    for chrom in chrom_list:
        df_chr = df_bins[df_bins["chrom"] == chrom]
        if len(df_chr) == 0:
            continue
        plot_chromosome(df_chr, calls, chrom, args.prefix,
                        segments=segments, cyto_dict=cyto_dict, guides=guides)

    # QC
    plot_qc(df_bins, args.prefix)


if __name__ == "__main__":
    main()
