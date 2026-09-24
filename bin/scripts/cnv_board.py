"""Shared CNV board for gxcnv1 and gxcnv2.

One PNG: genome-wide log2(ratio) on top, then per-chromosome panels
in a 4-column grid. Colours follow the portal reference figure
(cyan normal, purple gain, salmon loss). Dashed lines sit at the
log2(ratio) height of the caller z-score cutoff, not a fixed ±0.58.
"""

from __future__ import annotations

import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Reference figure palette
C_NORMAL = "#2EC8D6"   # cyan
C_GAIN = "#C44BD6"     # purple
C_LOSS = "#F07D5A"     # salmon
C_FILTERED = "#C5D0D6"
C_ZERO = "#9AA5B0"

Y_LO, Y_HI = -2.0, 2.0

_DEFAULT_CYTO_PATHS = [
    "/opt/gx-nipt/refs/bed/common/cytoBand.txt",
    os.path.join(os.path.dirname(__file__), "..", "..", "refs", "bed", "common", "cytoBand.txt"),
]

CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]


def resolve_cytoband(path: str | None) -> str | None:
    if path and os.path.isfile(path):
        return path
    for candidate in _DEFAULT_CYTO_PATHS:
        candidate = os.path.abspath(candidate)
        if os.path.isfile(candidate):
            return candidate
    return None


def load_cytobands(cyto_file: str) -> dict[str, list]:
    """UCSC cytoBand.txt → {chrom without 'chr': [[start, end, name, stain], ...]}."""
    cyto_dict: dict[str, list] = {}
    cur_chrom = None
    cur_bands: list = []
    with open(cyto_file) as fh:
        for line in fh:
            parts = line.split()
            if len(parts) < 4:
                continue
            chrom_key = parts[0].replace("chr", "")
            if chrom_key != cur_chrom:
                if cur_chrom is not None:
                    cyto_dict[cur_chrom] = cur_bands
                cur_bands = []
                cur_chrom = chrom_key
            cur_bands.append(parts[1:])
    if cur_chrom is not None:
        cyto_dict[cur_chrom] = cur_bands
    return cyto_dict


def _save_png(out: str, dpi: int) -> None:
    try:
        plt.savefig(out, dpi=dpi, bbox_inches="tight", pil_kwargs={"optimize": True})
    except TypeError:
        plt.savefig(out, dpi=dpi, bbox_inches="tight")


def _clean_lr(lr, lo=-2.5, hi=2.5):
    return np.clip(np.where(np.isfinite(lr), lr, 0.0), lo, hi)


def bin_colors(df: pd.DataFrame, calls: pd.DataFrame | None) -> list[str]:
    z = df["z_score"] if "z_score" in df.columns else pd.Series([np.nan] * len(df))
    colors = [
        C_NORMAL if (pd.notna(v) and str(v) not in ("nan", "")) else C_FILTERED
        for v in z
    ]
    if calls is None or calls.empty:
        return colors
    for _, row in calls.iterrows():
        chrom = row.get("chrom", "")
        s, e = float(row.get("start", 0)), float(row.get("end", 0))
        typ = str(row.get("type", "")).upper()
        c = C_GAIN if typ == "GAIN" else (C_LOSS if typ == "LOSS" else C_NORMAL)
        mask = (df["chrom"] == chrom) & (df["start"] >= s) & (df["end"] <= e)
        for i in df.index[mask]:
            colors[i] = c
    return colors


def _fmt_z(z: float) -> str:
    if abs(z - round(z)) < 0.05:
        return f"{z:.0f}"
    return f"{z:.2f}"


def read_z_cutoff(bins_path: str | None) -> float | None:
    """Caller z cutoff stamped on bins (##z_cutoff=) or a sibling qcmetrics file."""
    if not bins_path:
        return None
    try:
        with open(bins_path) as fh:
            for line in fh:
                if not line.startswith("#"):
                    break
                if line.startswith("##") and "z_cutoff=" in line:
                    raw = line.split("z_cutoff=", 1)[1].strip().split()[0]
                    return float(raw)
    except (OSError, ValueError):
        pass

    qc_path = bins_path.replace("_bins.tsv", "_qcmetrics.tsv")
    if qc_path == bins_path or not os.path.isfile(qc_path):
        return None
    metrics: dict[str, str] = {}
    try:
        with open(qc_path) as fh:
            for line in fh:
                if line.startswith("#") or "\t" not in line:
                    continue
                key, val = line.rstrip("\n").split("\t", 1)
                metrics[key] = val
    except OSError:
        return None
    # Wisecondor stores the threshold it actually tested with.
    saw_from_wcx = False
    for key in ("wc_threshold_z", "zscore_threshold", "z_cutoff"):
        raw = metrics.get(key, "")
        if raw == "from_wcx":
            saw_from_wcx = True
            continue
        try:
            return float(raw)
        except ValueError:
            continue
    # Older gxcnv2 beds runs stored "from_wcx". RUN_WCX uses 15 for mom, 6 otherwise.
    if saw_from_wcx:
        name = os.path.basename(bins_path)
        return 15.0 if "_mom_" in name else 6.0
    return None


def log2_at_z(df: pd.DataFrame, z_cutoff: float) -> tuple[float, float] | None:
    """Map a symmetric z cutoff onto this sample's log2(ratio) axis."""
    if "log2_ratio" not in df.columns or "z_score" not in df.columns:
        return None
    lr = pd.to_numeric(df["log2_ratio"], errors="coerce").to_numpy(dtype=float)
    z = pd.to_numeric(df["z_score"], errors="coerce").to_numpy(dtype=float)
    mask = np.isfinite(lr) & np.isfinite(z) & (np.abs(z) >= 0.25)
    if int(mask.sum()) < 30:
        return None
    slope = abs(float(np.median(lr[mask] / z[mask])))
    if not np.isfinite(slope) or slope == 0:
        return None
    hi = slope * float(z_cutoff)
    lo = -slope * float(z_cutoff)
    if not np.isfinite(hi) or not np.isfinite(lo):
        return None
    return hi, lo


def make_guides(df: pd.DataFrame, z_cutoff: float | None):
    """(log2 gain line, log2 loss line, z cutoff), or None when unknown."""
    if z_cutoff is None or not np.isfinite(z_cutoff) or z_cutoff <= 0:
        return None
    mapped = log2_at_z(df, float(z_cutoff))
    if mapped is None:
        return None
    return (mapped[0], mapped[1], float(z_cutoff))


def _legend_handles(guides: tuple[float, float, float] | None):
    handles = [
        plt.Line2D(
            [0], [0], marker="o", color="none", markerfacecolor=C_NORMAL,
            markeredgecolor="none", markersize=7, label="Normal (log₂ ratio)",
        ),
        plt.Line2D(
            [0], [0], marker="s", color="none", markerfacecolor=C_GAIN,
            markeredgecolor="none", markersize=7, label="Gain (call)",
        ),
        plt.Line2D(
            [0], [0], marker="s", color="none", markerfacecolor=C_LOSS,
            markeredgecolor="none", markersize=7, label="Loss (call)",
        ),
    ]
    if guides is not None:
        _hi, _lo, z_cut = guides
        ztxt = _fmt_z(z_cut)
        handles.append(plt.Line2D(
            [0], [0], color=C_GAIN, lw=1.3, ls="--",
            label=f"gain cutoff (z=+{ztxt})",
        ))
        handles.append(plt.Line2D(
            [0], [0], color=C_LOSS, lw=1.3, ls="--",
            label=f"loss cutoff (z=−{ztxt})",
        ))
    return handles


def _style_ax(ax) -> None:
    ax.set_facecolor("white")
    for sp in ax.spines.values():
        sp.set_color("#D0D5DA")
        sp.set_linewidth(0.6)
    ax.tick_params(colors="#4A5560", labelsize=7, length=2.5, width=0.5)


def _guides(ax, guides: tuple[float, float, float] | None, *, strong: bool = False) -> None:
    if strong:
        zero_c, zero_lw = "#5C6770", 1.05
        gain_c, loss_c, guide_lw = "#8E24AA", "#E04A28", 1.45
        dash = (0, (2.4, 1.3))
        alpha = 1.0
    else:
        zero_c, zero_lw = C_ZERO, 0.7
        gain_c, loss_c, guide_lw = C_GAIN, C_LOSS, 0.8
        dash = "--"
        alpha = 0.9
    ax.axhline(0, color=zero_c, lw=zero_lw, zorder=2)
    y_abs = 2.0
    if guides is not None:
        hi, lo, _z = guides
        ax.axhline(hi, color=gain_c, lw=guide_lw, ls=dash, alpha=alpha, zorder=3)
        ax.axhline(lo, color=loss_c, lw=guide_lw, ls=dash, alpha=alpha, zorder=3)
        y_abs = max(2.0, abs(hi), abs(lo))
        if y_abs > 2.0:
            y_abs = min(y_abs * 1.15, 3.5)
    ax.set_ylim(-y_abs, y_abs)
    if y_abs <= 2.05:
        ax.set_yticks([-2, -1, 0, 1, 2])
    else:
        top = int(np.floor(y_abs))
        ax.set_yticks(list(range(-top, top + 1)))


def _shade_calls(ax, calls: pd.DataFrame, x_of) -> None:
    if calls is None or calls.empty:
        return
    for _, row in calls.iterrows():
        typ = str(row.get("type", "")).upper()
        clr = C_GAIN if typ == "GAIN" else C_LOSS
        x0, x1 = x_of(row)
        if x1 <= x0:
            continue
        ax.axvspan(x0, x1, color=clr, alpha=0.35, zorder=1, linewidth=0)


def _draw_cytoband(ax, chrom: str, cyto_dict: dict[str, list] | None, x_max: float) -> None:
    ax.set_xlim(0, x_max)
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.set_xticks([])
    ax.set_facecolor("#F2F2F2")
    for sp in ax.spines.values():
        sp.set_visible(False)
    if not cyto_dict:
        ax.add_patch(mpatches.Rectangle((0, 0.2), x_max, 0.6, facecolor="#C8C8C8", lw=0))
        return
    bands = cyto_dict.get(chrom.replace("chr", ""), [])
    if not bands:
        ax.add_patch(mpatches.Rectangle((0, 0.2), x_max, 0.6, facecolor="#C8C8C8", lw=0))
        return
    for band in bands:
        start, end = float(band[0]), float(band[1])
        stain = band[3] if len(band) > 3 else "gneg"
        if stain.startswith("gpos"):
            try:
                level = float(stain[4:]) / 100.0
            except ValueError:
                level = 0.5
            g = 0.82 - 0.62 * level
        elif stain == "acen":
            g = 0.25
        elif stain in ("gvar", "stalk"):
            g = 0.45
        else:
            g = 0.88
        ax.add_patch(mpatches.Rectangle(
            (start, 0.18), max(end - start, 1.0), 0.64,
            facecolor=(g, g, g), edgecolor="none", lw=0,
        ))


def _attach_cyto(ax, chrom: str, cyto_dict, x_max: float):
    divider = make_axes_locatable(ax)
    ax_cyto = divider.append_axes("bottom", size="14%", pad=0.22)
    _draw_cytoband(ax_cyto, chrom, cyto_dict, x_max)
    return ax_cyto


def genome_positions(df: pd.DataFrame):
    chrom_lens = df.groupby("chrom_idx")["end"].max().sort_index()
    cumsum = [0]
    for ci in sorted(chrom_lens.index):
        cumsum.append(cumsum[-1] + int(chrom_lens[ci]))
    offsets = {ci: cumsum[i] for i, ci in enumerate(sorted(chrom_lens.index))}
    pos = df["start"].values + df["chrom_idx"].map(offsets).values
    return pos.astype(float), offsets, cumsum


def draw_genome(ax, df: pd.DataFrame, calls: pd.DataFrame | None,
                guides: tuple[float, float, float] | None = None) -> None:
    pos, _offsets, cumsum = genome_positions(df)
    lr = _clean_lr(df["log2_ratio"].values)
    colors = bin_colors(df, calls)
    unique = sorted(df["chrom_idx"].unique())

    _style_ax(ax)
    _guides(ax, guides)

    def x_of(row):
        ci = None
        chrom = str(row.get("chrom", ""))
        if chrom in CHROMS:
            ci = CHROMS.index(chrom)
        if ci is None or ci not in unique:
            return 0.0, 0.0
        idx = unique.index(ci)
        return cumsum[idx] + float(row["start"]), cumsum[idx] + float(row["end"])

    _shade_calls(ax, calls if calls is not None else pd.DataFrame(), x_of)
    ax.scatter(pos, lr, c=colors, s=2.0, alpha=0.85, linewidths=0, zorder=4, rasterized=True)

    n = len(unique)
    tick_pos = [(cumsum[i] + cumsum[i + 1]) / 2 for i in range(n)]
    tick_lbl = [CHROMS[ci].replace("chr", "") for ci in unique]
    ax.set_xticks(tick_pos)
    ax.set_xticklabels(tick_lbl, fontsize=8)
    ax.xaxis.tick_top()
    ax.tick_params(axis="x", length=0, pad=2)
    ax.set_xlim(0, cumsum[-1])
    ax.set_ylabel("log₂(ratio)", fontsize=8)
    # Title is a figure suptitle so it sits above these top tick labels.
    ax.legend(
        handles=_legend_handles(guides), fontsize=7.5, loc="upper center",
        bbox_to_anchor=(0.5, -0.08), ncol=5, frameon=True,
        framealpha=0.95, edgecolor="#E2E6EA",
    )


def draw_chromosome(ax, df_chr: pd.DataFrame, calls: pd.DataFrame | None,
                    chrom: str, cyto_dict,
                    guides: tuple[float, float, float] | None = None) -> None:
    df_chr = df_chr.reset_index(drop=True)
    x = (df_chr["start"].values + df_chr["end"].values) / 2
    lr = _clean_lr(df_chr["log2_ratio"].values)
    x_max = float(df_chr["end"].max()) if len(df_chr) else 1.0
    chr_calls = (
        calls[calls["chrom"] == chrom]
        if calls is not None and not calls.empty and "chrom" in calls.columns
        else pd.DataFrame()
    )
    colors = bin_colors(df_chr, chr_calls if not chr_calls.empty else None)

    _style_ax(ax)
    _guides(ax, guides, strong=True)
    _shade_calls(
        ax, chr_calls,
        lambda row: (float(row["start"]), float(row["end"])),
    )
    ax.scatter(x, lr, c=colors, s=3.5, alpha=0.85, linewidths=0, zorder=4, rasterized=True)
    label = chrom.replace("chr", "")
    ax.set_title(f"Chromosome {label}", fontsize=9, pad=3)
    ax.set_xlim(0, x_max)
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda v, _: f"{v / 1e6:.0f}"))
    _attach_cyto(ax, chrom, cyto_dict, x_max)


def render_board(df: pd.DataFrame, calls: pd.DataFrame | None, prefix: str,
                 cyto_dict: dict | None = None, log_prefix: str = "cnv_board",
                 dpi: int = 110, z_cutoff: float | None = None) -> str:
    """Write {prefix}_genome.png as the genome-wide + 4-column board."""
    present = [c for c in CHROMS if c in set(df["chrom"])]
    n = len(present)
    ncols = 4
    nrows = int(np.ceil(n / ncols)) if n else 1

    fig_w = 16.5
    fig_h = 3.6 + nrows * 2.15
    fig = plt.figure(figsize=(fig_w, fig_h), facecolor="white")
    gs = fig.add_gridspec(
        2 + nrows, ncols,
        height_ratios=[2.55, 0.22] + [1.7] * nrows,
        hspace=0.72, wspace=0.32,
        left=0.05, right=0.985, top=0.93, bottom=0.05,
    )

    guides = make_guides(df, z_cutoff)
    if guides is not None:
        print(
            f"[{log_prefix}] z cutoff={guides[2]:.4g} "
            f"→ log2 lines {guides[0]:+.3f} / {guides[1]:+.3f}",
            flush=True,
        )

    sample = os.path.basename(prefix)
    fig.suptitle(
        f"{sample} — Genome-wide CNV [log₂ ratio]",
        fontsize=13, fontweight="bold", y=0.975,
    )
    ax_g = fig.add_subplot(gs[0, :])
    draw_genome(ax_g, df, calls, guides=guides)

    ax_sec = fig.add_subplot(gs[1, :])
    ax_sec.axis("off")
    ax_sec.set_title(
        "Chromosome-wise CNV [log₂ ratio] — 4 chromosomes per row",
        fontsize=12, fontweight="bold", pad=2,
    )

    for i, chrom in enumerate(present):
        r, c = divmod(i, ncols)
        ax = fig.add_subplot(gs[2 + r, c])
        if c != 0:
            ax.set_ylabel("")
        draw_chromosome(ax, df[df["chrom"] == chrom], calls, chrom, cyto_dict, guides=guides)
        if c != 0:
            ax.set_ylabel("")

    chr_axes = [a for a in fig.axes if a.get_title().startswith("Chromosome")]
    for i, ax in enumerate(chr_axes):
        if i % ncols == 0:
            ax.set_ylabel("log₂(ratio)", fontsize=7)
        else:
            ax.set_ylabel("")

    fig.legend(
        handles=_legend_handles(guides), fontsize=8, loc="lower center",
        bbox_to_anchor=(0.5, 0.0), ncol=5, frameon=True,
        framealpha=0.95, edgecolor="#E2E6EA",
    )

    out = f"{prefix}_genome.png"
    _save_png(out, dpi)
    plt.close(fig)
    print(f"[{log_prefix}] {out}", flush=True)
    return out


def render_chromosome(df_chr: pd.DataFrame, calls: pd.DataFrame | None,
                      chrom: str, prefix: str, cyto_dict=None,
                      log_prefix: str = "cnv_board", dpi: int = 110,
                      guides: tuple[float, float, float] | None = None) -> str:
    """Single-chromosome PNG in the same style as a board panel."""
    fig, ax = plt.subplots(figsize=(8.2, 2.6), facecolor="white")
    draw_chromosome(ax, df_chr, calls, chrom, cyto_dict, guides=guides)
    ax.set_ylabel("log₂(ratio)", fontsize=8)
    fig.tight_layout(pad=0.4)
    safe = chrom.replace("/", "_")
    out = f"{prefix}_{safe}.png"
    _save_png(out, dpi)
    plt.close(fig)
    print(f"[{log_prefix}] {out}", flush=True)
    return out
