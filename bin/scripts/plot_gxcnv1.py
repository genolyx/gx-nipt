#!/usr/bin/env python3
"""
plot_gxcnv1.py — Genome visualisation for gxcnv1 results.

Visual style is intentionally distinct from both WisecondorX and gxcnv:
  • Primary metric : log2(ratio)  — not Z-score
  • Genome-wide   : filled-area track (not scatter), with ±trisomy guide lines
  • Color scheme  : teal / amber / violet on white background
  • QC plot       : KDE density curve  — not histogram
  • Per-chromosome: compact panels with confidence ribbon + CBS segment line

Outputs:
  {prefix}_genome.png   — genome-wide log2(ratio) track, all chromosomes
  {prefix}_chr{N}.png   — per-chromosome panels (individual files)
  {prefix}_qc.png       — KDE of log2(ratio) distribution
"""

import argparse
import os
import sys
import warnings

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
C_NORMAL    = "#4ECDC4"   # teal
C_GAIN      = "#C94FD8"   # violet
C_LOSS      = "#FF6B35"   # amber-orange
C_FILTERED  = "#C8D6DF"   # muted blue-grey
C_SEG_NORM  = "#607D8B"   # steel grey (segment line)
C_SEG_GAIN  = "#9B27AF"   # deep violet (gain segment)
C_SEG_LOSS  = "#E64A19"   # deep orange (loss segment)
C_RIBBON    = "#CFE8EC"   # pale teal (±MAD ribbon)
C_GUIDE     = "#90A4AE"   # light steel for reference lines

# Trisomy log2(3/2) ≈ +0.585 and monosomy log2(1/2) ≈ -1.000
LR_TRISOMY   =  0.585
LR_MONOSOMY  = -1.000

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
                prefix: str, segments: pd.DataFrame | None = None) -> None:
    """
    Genome-wide log2(ratio) filled-area track.

    Different from gxcnv (scatter Z-score) and WCX built-in (line plots):
    uses a filled-area ribbon to emphasise deviation from baseline.
    """
    pos, offsets, cumsum = genome_positions(df)
    lr = _clean_lr(df["log2_ratio"].values)
    colors = _bin_color(df, calls)

    unique_chroms = sorted(df["chrom_idx"].unique())
    n_chr = len(unique_chroms)

    ribbon_hw = _ribbon_mad(df)

    fig, ax = plt.subplots(figsize=(22, 4))
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    # Alternating chromosome bands (subtle)
    for i, ci in enumerate(unique_chroms):
        x0, x1 = cumsum[i], cumsum[i + 1]
        bg = "#F5F7F8" if i % 2 == 0 else "white"
        ax.axvspan(x0, x1, color=bg, zorder=0, linewidth=0)

    # ±MAD confidence ribbon
    ax.fill_between(
        [0, cumsum[-1]], [-ribbon_hw, -ribbon_hw], [ribbon_hw, ribbon_hw],
        color=C_RIBBON, alpha=0.55, zorder=1, linewidth=0,
    )

    # Trisomy / monosomy guide lines
    ax.axhline(0,           color=C_GUIDE, lw=0.8, zorder=2)
    ax.axhline(LR_TRISOMY,  color=C_GAIN,  lw=0.7, ls="--", alpha=0.45, zorder=2)
    ax.axhline(LR_MONOSOMY, color=C_LOSS,  lw=0.7, ls="--", alpha=0.45, zorder=2)

    # Filled area from 0 to log2_ratio (gain above, loss below)
    lr_arr = np.array(lr, dtype=float)
    ax.fill_between(pos, 0, lr_arr,
                    where=lr_arr >= 0, color=C_GAIN, alpha=0.25,
                    linewidth=0, zorder=3)
    ax.fill_between(pos, 0, lr_arr,
                    where=lr_arr < 0, color=C_LOSS, alpha=0.25,
                    linewidth=0, zorder=3)

    # Bin dots (coloured by call)
    ax.scatter(pos, lr_arr, c=colors, s=1.2, alpha=0.7, linewidths=0, zorder=4)

    # CBS segment overlay
    if segments is not None and not segments.empty:
        for _, seg in segments.iterrows():
            ci = CHROM_ORDER.get(seg["chrom"])
            if ci is None or ci not in unique_chroms:
                continue
            idx = unique_chroms.index(ci)
            off = cumsum[idx]
            slr = _clean_lr(np.array([seg.get("mean_log2_ratio", 0)]))[0]
            x0  = off + seg["start"]
            x1  = off + seg["end"]
            sc  = _seg_color(float(seg.get("mean_log2_ratio", 0)))
            lw  = 3.5 if abs(float(seg.get("mean_log2_ratio", 0))) > 0.3 else 1.2
            ax.hlines(slr, x0, x1, colors=sc, linewidths=lw, zorder=5)

    # Call highlight spans
    if calls is not None and not calls.empty:
        for _, row in calls.iterrows():
            ci = CHROM_ORDER.get(row["chrom"])
            if ci is None or ci not in unique_chroms:
                continue
            idx = unique_chroms.index(ci)
            x0  = cumsum[idx] + float(row["start"])
            x1  = cumsum[idx] + float(row["end"])
            clr = C_GAIN if str(row.get("type", "")) == "GAIN" else C_LOSS
            ax.axvspan(x0, x1, color=clr, alpha=0.22, zorder=6)
            ax.text(
                (x0 + x1) / 2, 1.9,
                str(row.get("type", "")), fontsize=6,
                ha="center", va="top", color=clr, fontweight="bold", zorder=7,
            )

    # X-axis: chromosome labels
    tick_pos = [(cumsum[i] + cumsum[i + 1]) / 2 for i in range(n_chr)]
    tick_lbl = [CHROMS[ci].replace("chr", "") for ci in unique_chroms]
    ax.set_xticks(tick_pos)
    ax.set_xticklabels(tick_lbl, fontsize=7)
    ax.set_xlim(0, cumsum[-1])
    ax.set_ylim(-2.5, 2.5)
    ax.set_ylabel("log₂(ratio)", fontsize=9)
    ax.set_title(
        f"{os.path.basename(prefix)} — Genome-wide CNV  [log₂ ratio]",
        fontsize=10,
    )

    legend_handles = [
        mpatches.Patch(color=C_NORMAL,  label="Normal"),
        mpatches.Patch(color=C_GAIN,    label="Gain (call)"),
        mpatches.Patch(color=C_LOSS,    label="Loss (call)"),
        mpatches.Patch(color=C_RIBBON,  alpha=0.7, label="±1.5 MAD ribbon"),
        plt.Line2D([0], [0], color=C_GAIN, lw=2, ls="--",
                   label=f"Trisomy guide (+{LR_TRISOMY:.2f})"),
        plt.Line2D([0], [0], color=C_LOSS, lw=2, ls="--",
                   label=f"Monosomy guide ({LR_MONOSOMY:.2f})"),
    ]
    ax.legend(handles=legend_handles, fontsize=7, loc="upper right",
              framealpha=0.85, edgecolor="#CCCCCC", ncol=2)

    plt.tight_layout(pad=0.5)
    out = f"{prefix}_genome.png"
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[plot_gxcnv1] {out}", flush=True)


# ── Per-chromosome plot ────────────────────────────────────────────────────────

def plot_chromosome(df_chr: pd.DataFrame, calls: pd.DataFrame | None,
                    chrom: str, prefix: str,
                    segments: pd.DataFrame | None = None) -> None:
    """
    Single-chromosome log2(ratio) panel.

    Layout: filled area + CBS segment line + confidence ribbon.
    Annotated with any GAIN/LOSS calls for this chromosome.
    """
    fig, ax = plt.subplots(figsize=(13, 3.2))
    fig.patch.set_facecolor("white")
    ax.set_facecolor("#FAFBFC")

    x   = (df_chr["start"].values + df_chr["end"].values) / 2
    lr  = _clean_lr(df_chr["log2_ratio"].values)

    chr_calls = (calls[calls["chrom"] == chrom]
                 if calls is not None and not calls.empty else pd.DataFrame())

    colors = _bin_color(df_chr.reset_index(drop=True),
                        chr_calls if not chr_calls.empty else None)

    ribbon_hw = _ribbon_mad(df_chr)
    ax.fill_between(x, -ribbon_hw, ribbon_hw,
                    color=C_RIBBON, alpha=0.55, linewidth=0, zorder=1)

    ax.axhline(0,           color=C_GUIDE, lw=0.8, zorder=2)
    ax.axhline(LR_TRISOMY,  color=C_GAIN,  lw=0.7, ls="--", alpha=0.4, zorder=2)
    ax.axhline(LR_MONOSOMY, color=C_LOSS,  lw=0.7, ls="--", alpha=0.4, zorder=2)

    lr_arr = np.array(lr, dtype=float)
    ax.fill_between(x, 0, lr_arr, where=lr_arr >= 0,
                    color=C_GAIN, alpha=0.22, linewidth=0, zorder=3)
    ax.fill_between(x, 0, lr_arr, where=lr_arr < 0,
                    color=C_LOSS, alpha=0.22, linewidth=0, zorder=3)

    ax.scatter(x, lr_arr, c=colors, s=6, alpha=0.75, linewidths=0, zorder=4)

    # CBS segments for this chromosome
    if segments is not None and not segments.empty:
        chr_segs = segments[segments["chrom"] == chrom]
        for _, seg in chr_segs.iterrows():
            slr = _clean_lr(np.array([seg.get("mean_log2_ratio", 0)]))[0]
            sc  = _seg_color(float(seg.get("mean_log2_ratio", 0)))
            lw  = 3.5 if abs(float(seg.get("mean_log2_ratio", 0))) > 0.3 else 1.5
            ax.hlines(slr, seg["start"], seg["end"],
                      colors=sc, linewidths=lw, zorder=5)

    # Highlight calls
    for _, row in chr_calls.iterrows():
        x0, x1 = float(row["start"]), float(row["end"])
        typ     = str(row.get("type", ""))
        clr     = C_GAIN if typ == "GAIN" else C_LOSS
        ax.axvspan(x0, x1, color=clr, alpha=0.18, zorder=6)
        ax.text((x0 + x1) / 2, 2.1, typ,
                fontsize=7, ha="center", color=clr,
                fontweight="bold", zorder=7)

    chrom_label = chrom.replace("chr", "")
    ax.set_title(f"{os.path.basename(prefix)} — Chr {chrom_label}  [log₂ ratio]",
                 fontsize=9)
    ax.set_xlabel("Genomic position", fontsize=8)
    ax.set_ylabel("log₂(ratio)", fontsize=8)
    ax.set_ylim(-2.5, 2.5)
    ax.xaxis.set_major_formatter(
        ticker.FuncFormatter(lambda v, _: f"{v/1e6:.0f} Mb")
    )

    plt.tight_layout(pad=0.3)
    safe = chrom.replace("/", "_")
    out = f"{prefix}_{safe}.png"
    plt.savefig(out, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"[plot_gxcnv1] {out}", flush=True)


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
    ax.set_title(f"{os.path.basename(prefix)} — log₂(ratio) distribution  [gxcnv1 QC]",
                 fontsize=10)
    ax.legend(fontsize=8, loc="upper right")

    plt.tight_layout(pad=0.5)
    out = f"{prefix}_qc.png"
    plt.savefig(out, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"[plot_gxcnv1] {out}", flush=True)


# ── Touch helper ──────────────────────────────────────────────────────────────

def _touch(path: str) -> None:
    open(path, "a").close()


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(
        description="Plot gxcnv1 results (log2-ratio style)"
    )
    ap.add_argument("--bins",     required=True,  help="*_bins.tsv from gxcnv1_predict.py")
    ap.add_argument("--calls",    required=True,  help="*_calls.tsv from gxcnv1_predict.py")
    ap.add_argument("--segments", default=None,   help="*_segments.tsv (auto-detected if omitted)")
    ap.add_argument("-o", "--prefix", required=True)
    ap.add_argument("--chromosomes", nargs="*", default=None,
                    help="Only plot these chromosomes (default: all)")
    args = ap.parse_args()

    df_bins = load_bins(args.bins)
    if df_bins is None or df_bins.empty:
        print("[plot_gxcnv1] Empty bins file — creating stub outputs", flush=True)
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
        print(f"[plot_gxcnv1] Loaded {len(segments)} CBS segments", flush=True)

    # Genome-wide
    plot_genome(df_bins, calls, args.prefix, segments=segments)

    # Per-chromosome
    chrom_list = args.chromosomes or sorted(
        df_bins["chrom"].unique(), key=lambda c: CHROM_ORDER.get(c, 99)
    )
    for chrom in chrom_list:
        df_chr = df_bins[df_bins["chrom"] == chrom]
        if len(df_chr) == 0:
            continue
        plot_chromosome(df_chr, calls, chrom, args.prefix, segments=segments)

    # QC
    plot_qc(df_bins, args.prefix)


if __name__ == "__main__":
    main()
