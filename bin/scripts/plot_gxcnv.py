#!/usr/bin/env python3
"""
Generate gx-cnv visualisation figures (genome-wide + per-chromosome).

Outputs:
  {prefix}_genome.png       — genome-wide Z-score track with region callouts
  {prefix}_chrXX.png        — per-chromosome Z-score plot for every chr in data
  {prefix}_regions.png      — monitored microdeletion / microduplication regions (bar chart)
  {prefix}_qc.png           — QC histogram of all bin Z-scores

Usage:
  python plot_gxcnv.py --bins SAMPLE_bins.tsv --calls SAMPLE_calls.tsv -o SAMPLE
"""

import argparse
import os
import sys
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore", category=FutureWarning)

# ── Colours ──────────────────────────────────────────────────────────────────
C_NORMAL    = "#3A86FF"
C_DUP       = "#FF006E"
C_DEL       = "#FFBE0B"
C_FILTERED  = "#CCCCCC"
C_HIGH_RISK = "#E63946"
C_LOW_RISK  = "#2DC653"
C_REGION_BG = "#FFE8A1"


CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
CHROM_ORDER = {c: i for i, c in enumerate(CHROMS)}


def _flag_color(flag):
    if flag in ("DUP_CANDIDATE", "DUP"):
        return C_DUP
    if flag in ("DEL_CANDIDATE", "DEL"):
        return C_DEL
    if flag in ("NORMAL",):
        return C_NORMAL
    return C_FILTERED


def assign_aberration_colors(df_bins, calls):
    """
    Color bins based ONLY on HIGH_RISK calls from calls.tsv — equivalent to
    WCX's aberrations.bed coloring logic.

    HIGH_RISK dual_call → bins within that region colored as DUP (gain) or DEL (loss).
    Everything else stays NORMAL (blue) or FILTERED (gray).

    The dual_call type (gain vs loss) is inferred from the track_a_mean_z sign.
    """
    valid_flags = {"NORMAL", "DUP_CANDIDATE", "DEL_CANDIDATE", "DUP", "DEL"}
    colors = np.array(
        [C_FILTERED if f not in valid_flags else C_NORMAL
         for f in df_bins["flag"].fillna("FILTERED")],
        dtype=object
    )

    if calls is None or calls.empty:
        return colors.tolist()

    high_risk = calls[calls.get("dual_call", calls.get("call", "")) == "HIGH_RISK"] \
        if "dual_call" in calls.columns else pd.DataFrame()

    for _, row in high_risk.iterrows():
        chrom   = row.get("chrom", "")
        s_start = int(row.get("start", 0))
        s_end   = int(row.get("end",   0))
        mz      = float(row.get("track_a_mean_z", 0)) if pd.notna(row.get("track_a_mean_z")) else 0.0

        mask = (
            (df_bins["chrom"] == chrom) &
            (df_bins["start"] >= s_start) &
            (df_bins["end"]   <= s_end) &
            df_bins["flag"].isin(valid_flags)
        )
        colors[mask] = C_DUP if mz >= 0 else C_DEL

    return colors.tolist()


def load_segments(path):
    """Load CBS segments from _segments.tsv."""
    rows = []
    with open(path) as f:
        cols = None
        for line in f:
            line = line.rstrip()
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                cols = line.lstrip("#").split("\t")
                continue
            if cols is None:
                continue
            rows.append(dict(zip(cols, line.split("\t"))))
    if not rows:
        return None
    df = pd.DataFrame(rows)
    for c in ("start", "end", "mean_z"):
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    df["chrom"] = df["chrom"].str.strip()
    df = df[df["chrom"].isin(CHROM_ORDER)].copy()
    df["chrom_idx"] = df["chrom"].map(CHROM_ORDER)
    return df


def load_bins(path):
    rows = []
    with open(path) as f:
        cols = None
        for line in f:
            line = line.rstrip()
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                cols = line.lstrip("#").split("\t")
                continue
            if cols is None:
                continue
            rows.append(dict(zip(cols, line.split("\t"))))
    if not rows:
        return None
    df = pd.DataFrame(rows)
    df["z_score"] = pd.to_numeric(df["z_score"], errors="coerce")
    df["start"]   = pd.to_numeric(df["start"],   errors="coerce")
    df["end"]     = pd.to_numeric(df["end"],     errors="coerce")
    df["chrom"]   = df["chrom"].str.strip()
    df = df[df["chrom"].isin(CHROM_ORDER)].copy()
    df["chrom_idx"] = df["chrom"].map(CHROM_ORDER)
    df = df.sort_values(["chrom_idx", "start"]).reset_index(drop=True)
    return df


def load_calls(path):
    rows = []
    with open(path) as f:
        cols = None
        for line in f:
            line = line.rstrip()
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                cols = line.lstrip("#").split("\t")
                continue
            if cols is None:
                continue
            rows.append(dict(zip(cols, line.split("\t"))))
    if not rows:
        return pd.DataFrame()
    df = pd.DataFrame(rows)
    df["start"]       = pd.to_numeric(df["start"],       errors="coerce")
    df["end"]         = pd.to_numeric(df["end"],         errors="coerce")
    df["track_a_mean_z"] = pd.to_numeric(df["track_a_mean_z"], errors="coerce")
    df["risk_pct"]    = pd.to_numeric(df["risk_pct"],    errors="coerce")
    return df


def genome_cumulative_positions(df):
    """Return per-row cumulative genomic position for a genome-wide x-axis."""
    chrom_lens = df.groupby("chrom_idx")["end"].max().sort_index()
    cumsum = [0]
    for c in sorted(chrom_lens.index):
        cumsum.append(cumsum[-1] + int(chrom_lens[c]))
    chrom_offset = {c: cumsum[i] for i, c in enumerate(sorted(chrom_lens.index))}
    pos = df["start"].values + df["chrom_idx"].map(chrom_offset).values
    return pos.astype(float), chrom_offset, cumsum


def _clean_z(z, lo=-10, hi=10):
    return np.clip(np.where(np.isfinite(z), z, 0), lo, hi)


def _seg_color(z):
    if not np.isfinite(z):
        return "#AAAAAA"
    if z >= 3:
        return C_DUP
    if z <= -3:
        return C_DEL
    return "#888888"


def plot_genome(df, calls, prefix, segments=None):
    """Genome-wide Z-score scatter with CBS segment overlay."""
    pos, offsets, cumsum = genome_cumulative_positions(df)
    z = _clean_z(df["z_score"].values, lo=-5, hi=5)
    # Color only HIGH_RISK called regions — mirrors WCX aberrations.bed coloring
    colors = assign_aberration_colors(df, calls)

    unique_chroms = sorted(df["chrom_idx"].unique())
    n_chr         = len(unique_chroms)

    fig, ax = plt.subplots(figsize=(20, 4))
    fig.patch.set_facecolor("#FAFAFA")
    ax.set_facecolor("#FAFAFA")

    # Alternating chromosome backgrounds
    for i, ci in enumerate(unique_chroms):
        x0 = cumsum[i]
        x1 = cumsum[i + 1]
        bg = "#F0F4FF" if i % 2 == 0 else "#FAFAFA"
        ax.axvspan(x0, x1, color=bg, zorder=0)

    ax.scatter(pos, z, c=colors, s=1.5, alpha=0.65, linewidths=0, zorder=2)

    # CBS segment overlay — thick line for significant segments
    if segments is not None:
        for _, seg in segments.iterrows():
            ci = CHROM_ORDER.get(seg["chrom"])
            if ci is None or ci not in unique_chroms:
                continue
            idx   = unique_chroms.index(ci)
            off   = cumsum[idx]
            seg_z = _clean_z(np.array([seg.get("mean_z", 0)]), lo=-5, hi=5)[0]
            x0    = off + seg["start"]
            x1    = off + seg["end"]
            sc    = _seg_color(seg.get("mean_z", 0))
            lw    = 3.5 if abs(seg.get("mean_z", 0)) >= 2.5 else 1.2
            ax.hlines(seg_z, x0, x1, colors=sc, linewidths=lw, zorder=4)

    # Reference bands
    ax.axhline(0,   color="#999999", lw=0.8, zorder=1)
    ax.axhline(-3,  color=C_DEL, lw=0.6, ls="--", alpha=0.5, zorder=1)
    ax.axhline( 3,  color=C_DUP, lw=0.6, ls="--", alpha=0.5, zorder=1)

    # HIGH_RISK region highlights
    high_risk = calls[calls["dual_call"] == "HIGH_RISK"] if "dual_call" in calls.columns else pd.DataFrame()
    for _, row in high_risk.iterrows():
        ci = CHROM_ORDER.get(row["chrom"])
        if ci is None or ci not in unique_chroms:
            continue
        idx = unique_chroms.index(ci)
        x0  = cumsum[idx] + row["start"]
        x1  = cumsum[idx] + row["end"]
        ax.axvspan(x0, x1, color=C_HIGH_RISK, alpha=0.3, zorder=3)
        ax.text((x0 + x1) / 2, 4.3,
                row.get("region_name", ""), fontsize=5.5, ha="center",
                color=C_HIGH_RISK, fontweight="bold", zorder=5)

    # Chromosome labels on x-axis
    tick_pos = [(cumsum[i] + cumsum[i + 1]) / 2 for i in range(n_chr)]
    tick_lbl = [CHROMS[ci].replace("chr", "") for ci in unique_chroms]
    ax.set_xticks(tick_pos)
    ax.set_xticklabels(tick_lbl, fontsize=7)
    ax.set_xlim(0, cumsum[-1])
    ax.set_ylim(-5, 5)
    ax.set_ylabel("Z-score", fontsize=9)
    ax.set_title(f"{os.path.basename(prefix)} — Genome-wide CNV Z-score", fontsize=10)

    # Legend
    legend_handles = [
        mpatches.Patch(color=C_NORMAL,    label="Normal"),
        mpatches.Patch(color=C_DUP,       label="Gain (HIGH RISK call)"),
        mpatches.Patch(color=C_DEL,       label="Loss (HIGH RISK call)"),
        mpatches.Patch(color=C_HIGH_RISK, alpha=0.4, label="HIGH RISK region"),
        plt.Line2D([0], [0], color=C_DUP,    lw=2.5, label="CBS segment (gain)"),
        plt.Line2D([0], [0], color=C_DEL,    lw=2.5, label="CBS segment (loss)"),
        plt.Line2D([0], [0], color="#888888", lw=1.2, label="CBS segment (neutral)"),
    ]
    ax.legend(handles=legend_handles, fontsize=7, loc="upper right",
              framealpha=0.8, edgecolor="#CCCCCC")

    plt.tight_layout(pad=0.5)
    out = f"{prefix}_genome.png"
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[plot_gxcnv] {out}", flush=True)


def plot_chromosome(df_chr, calls, chrom, prefix, segments=None):
    """Z-score scatter + CBS segment overlay for a single chromosome."""
    fig, ax = plt.subplots(figsize=(12, 3))
    fig.patch.set_facecolor("#FAFAFA")
    ax.set_facecolor("#FAFAFA")

    x = df_chr["start"].values + (df_chr["end"].values - df_chr["start"].values) / 2
    z = _clean_z(df_chr["z_score"].values, lo=-5, hi=5)
    # Color only HIGH_RISK called regions (calls.tsv)
    chr_calls = calls[calls["chrom"] == chrom] if not calls.empty else pd.DataFrame()
    c = assign_aberration_colors(df_chr, chr_calls)

    ax.scatter(x, z, c=c, s=7, alpha=0.7, linewidths=0, zorder=2)
    ax.axhline(0,  color="#999999", lw=0.8, zorder=1)
    ax.axhline(-3, color=C_DEL, lw=0.7, ls="--", alpha=0.6, zorder=1)
    ax.axhline( 3, color=C_DUP, lw=0.7, ls="--", alpha=0.6, zorder=1)

    # CBS segment overlay
    if segments is not None:
        chr_segs = segments[segments["chrom"] == chrom]
        for _, seg in chr_segs.iterrows():
            seg_z = _clean_z(np.array([seg.get("mean_z", 0)]), lo=-5, hi=5)[0]
            sc    = _seg_color(seg.get("mean_z", 0))
            lw    = 3.5 if abs(seg.get("mean_z", 0)) >= 2.5 else 1.5
            ax.hlines(seg_z, seg["start"], seg["end"],
                      colors=sc, linewidths=lw, zorder=5)

    # All monitored region highlights (calls.tsv contains ALL regions, not just HIGH_RISK)
    chr_calls = calls[calls["chrom"] == chrom] if not calls.empty else pd.DataFrame()
    ymax = 7.5
    for _, row in chr_calls.iterrows():
        x0, x1 = row["start"], row["end"]
        call    = row.get("dual_call", "LOW_RISK")
        color   = C_HIGH_RISK if call == "HIGH_RISK" else "#CCCCCC"
        alpha   = 0.3 if call == "HIGH_RISK" else 0.12
        ax.axvspan(x0, x1, color=color, alpha=alpha, zorder=3)
        rname = row.get("region_name", "")
        if rname:
            ax.text((x0 + x1) / 2, ymax * 0.88,
                    rname.replace("_", "\n"), fontsize=6, ha="center",
                    color=color if call == "HIGH_RISK" else "#888888",
                    fontweight="bold", zorder=6)

    chrom_name = chrom.replace("chr", "")
    ax.set_title(f"{os.path.basename(prefix)} — Chr {chrom_name}", fontsize=9)
    ax.set_xlabel("Genomic position (bp)", fontsize=8)
    ax.set_ylabel("Z-score", fontsize=8)
    ax.set_ylim(-5, 5)
    ax.xaxis.set_major_formatter(
        matplotlib.ticker.FuncFormatter(lambda v, _: f"{v/1e6:.0f}Mb")
    )

    plt.tight_layout(pad=0.3)
    safe_chr = chrom.replace("/", "_")
    out = f"{prefix}_{safe_chr}.png"
    plt.savefig(out, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"[plot_gxcnv] {out}", flush=True)


def plot_regions(calls, prefix):
    """Bar chart of all monitored regions coloured by risk."""
    if calls.empty or "region_name" not in calls.columns:
        _touch(f"{prefix}_regions.png")
        return

    fig, ax = plt.subplots(figsize=(max(8, len(calls) * 0.55 + 2), 4))
    fig.patch.set_facecolor("#FAFAFA")
    ax.set_facecolor("#FAFAFA")

    risk     = calls["risk_pct"].fillna(0).values
    is_high  = calls["dual_call"].values == "HIGH_RISK"
    colors   = [C_HIGH_RISK if h else C_LOW_RISK for h in is_high]
    regions  = calls["region_name"].values
    x        = np.arange(len(regions))

    bars = ax.bar(x, risk, color=colors, width=0.7, edgecolor="white", linewidth=0.5)
    ax.axhline(50, color="#888888", lw=0.7, ls="--")
    ax.set_xticks(x)
    ax.set_xticklabels(regions, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Risk (%)", fontsize=9)
    ax.set_ylim(0, 105)
    ax.set_title(f"{os.path.basename(prefix)} — Monitored Region Risk", fontsize=10)

    handles = [
        mpatches.Patch(color=C_HIGH_RISK, label="HIGH RISK"),
        mpatches.Patch(color=C_LOW_RISK,  label="Low Risk"),
    ]
    ax.legend(handles=handles, fontsize=8, loc="upper right")

    plt.tight_layout(pad=0.5)
    out = f"{prefix}_regions.png"
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[plot_gxcnv] {out}", flush=True)


def plot_qc(df, prefix):
    """Histogram of per-bin Z-scores for QC."""
    z = df["z_score"].dropna().values
    z = z[np.isfinite(z)]
    z = np.clip(z, -10, 10)

    fig, ax = plt.subplots(figsize=(7, 4))
    fig.patch.set_facecolor("#FAFAFA")
    ax.set_facecolor("#FAFAFA")

    ax.hist(z, bins=100, color=C_NORMAL, edgecolor="white", linewidth=0.3, alpha=0.85)
    ax.axvline(-3, color=C_DEL, lw=1.2, ls="--", label="thresh −3")
    ax.axvline( 3, color=C_DUP, lw=1.2, ls="--", label="thresh +3")
    ax.set_xlabel("Z-score", fontsize=9)
    ax.set_ylabel("Bin count", fontsize=9)
    ax.set_title(f"{os.path.basename(prefix)} — Z-score distribution (QC)", fontsize=10)
    ax.legend(fontsize=8)

    plt.tight_layout(pad=0.5)
    out = f"{prefix}_qc.png"
    plt.savefig(out, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"[plot_gxcnv] {out}", flush=True)


def _touch(path):
    open(path, "a").close()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bins",     required=True)
    ap.add_argument("--calls",    required=True)
    ap.add_argument("--segments", default=None,
                    help="CBS segments TSV (_segments.tsv); auto-detected if omitted")
    ap.add_argument("-o", "--prefix", required=True)
    ap.add_argument("--chromosomes", nargs="*", default=None,
                    help="Only plot these chromosomes (default: all in data)")
    args = ap.parse_args()

    df_bins = load_bins(args.bins)
    if df_bins is None or df_bins.empty:
        print("[plot_gxcnv] Empty bins file — creating stub outputs", flush=True)
        for suf in ("_genome.png", "_regions.png", "_qc.png"):
            _touch(f"{args.prefix}{suf}")
        return

    calls = load_calls(args.calls)

    # Auto-detect segments file next to bins file
    seg_path = args.segments
    if seg_path is None:
        candidate = args.bins.replace("_bins.tsv", "_segments.tsv")
        if os.path.isfile(candidate):
            seg_path = candidate
    segments = load_segments(seg_path) if seg_path and os.path.isfile(seg_path) else None
    if segments is not None:
        print(f"[plot_gxcnv] Loaded {len(segments)} segments from {seg_path}", flush=True)

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

    # Regions bar chart — uses all monitored regions from calls tsv
    plot_regions(calls, args.prefix)

    # QC histogram
    plot_qc(df_bins, args.prefix)


if __name__ == "__main__":
    main()
