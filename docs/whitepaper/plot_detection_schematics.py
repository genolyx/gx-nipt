#!/usr/bin/env python3
"""
White-paper schematics for EZD aneuploidy calling.

  1. sca_female_xo_xxx  : chrX UAR vs Z-score  (XO / XXX)
  2. sca_male_xxy_xyy   : UAR[X] vs UAR[Y]     (XXY / XYY)
  3. trisomy_<chr>      : UAR vs Z-score       (e.g. T21)

Decision zones are produced by evaluating the production decision code
(SCADetector.detect_female_sca / detect_male_sca, and the run_decision rule
in ezd_runner.py) so the figure cannot drift from the pipeline logic.
Reference points come from the lab EZD reference tables.

Usage:
  python3 docs/whitepaper/plot_detection_schematics.py \
      --ezd-dir /home/ken/ken-nipt/data/refs/ucl/EZD --group orig \
      --chrom chr21 --out docs/whitepaper/figures
"""

import argparse
import json
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch, Polygon, Rectangle  # noqa: E402

matplotlib.rcParams["hatch.linewidth"] = 0.6

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(REPO, "bin", "scripts", "modules"))
from sca_detector import SCADetector  # noqa: E402

LINE = "#263238"
ND_FILL = "#F4F6F8"
HATCH_EDGE = "#B8C0C8"
NORMAL_PT = "#55606B"

ZONE = {
    "XO Detected": "#F2A184",
    "XO Suspected": "#FBDDD0",
    "XXX Detected": "#F0C94E",
    "XXX Suspected": "#FAEBB8",
    "XXY Detected": "#F5AE62",
    "XXY Suspected": "#FCE3C6",
    "XYY Detected": "#8FD1A0",
    "XYY Suspected": "#D6EFDC",
    "Detected": "#D49BE2",
    "Suspected": "#F1DDF6",
    "Not Detected": ND_FILL,
}
POINT = {
    "XX": NORMAL_PT, "XY": NORMAL_PT, "N": NORMAL_PT,
    "XO": "#D4491F", "XXX": "#B88A00",
    "XXY": "#E07B10", "XYY": "#23874C", "XXYY": "#7B3FA0",
    "P": "#8E2BAA",
}

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 9,
    "axes.edgecolor": "#5F6B76",
    "axes.linewidth": 0.8,
    "xtick.color": "#37474F",
    "ytick.color": "#37474F",
    "svg.fonttype": "none",
})


# ── helpers ─────────────────────────────────────────────────────────────────
def _new_fig(title, subtitle):
    fig = plt.figure(figsize=(8.6, 7.9), facecolor="white")
    ax = fig.add_axes([0.10, 0.315, 0.64, 0.50])
    fig.text(0.10, 0.965, title, fontsize=14, fontweight="bold", color=LINE)
    fig.text(0.10, 0.938, subtitle, fontsize=9.5, color="#5F6B76")
    return fig, ax


def _spread(vals, lo, hi, gap):
    """Positions (same order as vals) at least `gap` apart, centred on the
    originals and kept inside [lo, hi]."""
    order = np.argsort(vals)
    orig = np.array(vals, dtype=float)[order]
    pos = orig.copy()
    for i in range(1, len(pos)):
        pos[i] = max(pos[i], pos[i - 1] + gap)
    if len(pos):
        pos += orig.mean() - pos.mean()
        if pos[0] < lo:
            pos += lo - pos[0]
        if pos[-1] > hi:
            pos -= pos[-1] - hi
    out = np.empty_like(pos)
    out[order] = pos
    return out


def _vlabels(ax, items):
    """items: [(x, name, value_str)] → labels above the axes with leaders."""
    x0, x1 = ax.get_xlim()
    fr = [(x - x0) / (x1 - x0) for x, _, _ in items]
    tx = _spread(fr, 0.06, 0.94, 0.2)
    for (x, name, val), f, t in zip(items, fr, tx):
        ax.annotate(
            f"{name}\n{val}", xy=(f, 1.0), xycoords="axes fraction",
            xytext=(t, 1.06), textcoords="axes fraction",
            ha="center", va="bottom", fontsize=8.2, color=LINE,
            arrowprops=dict(arrowstyle="-", color="#90A4AE", lw=0.7,
                            shrinkA=0, shrinkB=0),
        )


def _hlabels(ax, items):
    """items: [(y, name, value_str)] → labels right of the axes with leaders."""
    y0, y1 = ax.get_ylim()
    fr = [(y - y0) / (y1 - y0) for y, _, _ in items]
    ty = _spread(fr, 0.03, 0.97, 0.085)
    for (y, name, val), f, t in zip(items, fr, ty):
        ax.annotate(
            f"{name} = {val}", xy=(1.0, f), xycoords="axes fraction",
            xytext=(1.04, t), textcoords="axes fraction",
            ha="left", va="center", fontsize=8.2, color=LINE,
            arrowprops=dict(arrowstyle="-", color="#90A4AE", lw=0.7,
                            shrinkA=0, shrinkB=0),
        )


def _zone_text(ax, x, y, s, size=9, weight="bold", color=LINE):
    ax.text(x, y, s, ha="center", va="center", fontsize=size,
            fontweight=weight, color=color, zorder=6,
            bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="none",
                      alpha=0.72))


def _legends(fig, points, zones):
    kw = dict(loc="lower left", fontsize=8, frameon=False, handlelength=1.3,
              columnspacing=1.3, handletextpad=0.5)
    fig.legend(handles=points, bbox_to_anchor=(0.095, 0.222),
               ncol=len(points), **kw)
    fig.legend(handles=zones, bbox_to_anchor=(0.095, 0.197),
               ncol=len(zones), **kw)


def _pt(color, label, ms=5, alpha=1.0, marker="o"):
    return Line2D([], [], ls="", marker=marker, ms=ms, color=color,
                  alpha=alpha, label=label)


def _zp(fill, label, hatch=None):
    return Patch(fc=fill, ec=HATCH_EDGE if hatch else "#C9CED4", lw=0.5,
                 hatch=hatch, label=label)


def _rules(fig, lines, y=0.175):
    fig.text(0.10, y, "Decision rule", fontsize=10, fontweight="bold",
             color=LINE)
    fig.text(0.10, y - 0.02, "\n".join(lines), fontsize=8.6, color="#37474F",
             va="top", linespacing=1.55)


def _save(fig, out_dir, stem):
    os.makedirs(out_dir, exist_ok=True)
    for ext in ("png", "svg", "pdf"):
        path = os.path.join(out_dir, f"{stem}.{ext}")
        fig.savefig(path, dpi=300 if ext == "png" else None, facecolor="white")
        print(path)
    plt.close(fig)


def _fmt(v, nd=2):
    return f"{v:.{nd}f}".replace("-", "\u2212")


# ── 1. female XO / XXX ──────────────────────────────────────────────────────
def plot_female(det, group, ref_csv, out_dir, xlim=(4.5, 6.0), ylim=(-30, 30)):
    p = det.get_female_params(group)
    a, b = p["xo_ur_x_max"], p["ur_x_low"]
    c, d = p["ur_x_high"], p["xxx_ur_x_min"]
    zlo, zhi = p["xo_z_threshold"], p["xxx_z_threshold"]

    fig, ax = _new_fig(
        "Sex-chromosome aneuploidy (female): XO / XXX",
        f"chrX  ·  UAR[X] vs Z-score[X]  ·  {group} reference",
    )
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)

    # Evaluate production logic on every cell (z_normal_high/low included
    # as extra cuts so any effect they have would show up).
    xs = [xlim[0], a, b, c, d, xlim[1]]
    ys = sorted({ylim[0], zlo, p["z_normal_low"], p["z_normal_high"], zhi,
                 ylim[1]})
    for i in range(len(xs) - 1):
        for j in range(len(ys) - 1):
            xm, ym = (xs[i] + xs[i + 1]) / 2, (ys[j] + ys[j + 1]) / 2
            call = det.detect_female_sca(xm, ym, group)
            discordant = call == "Not Detected" and not (b <= xm <= c)
            ax.add_patch(Rectangle(
                (xs[i], ys[j]), xs[i + 1] - xs[i], ys[j + 1] - ys[j],
                fc=ZONE[call], ec=HATCH_EDGE if discordant else "none",
                hatch="////" if discordant else None, lw=0, zorder=0,
            ))

    df = pd.read_csv(ref_csv, sep="\t")
    df["type"] = df["type"].replace({"X0": "XO"})
    counts = {}
    for t, s, al, z in (("XX", 6, 0.30, 2), ("XO", 16, 0.9, 3),
                        ("XXX", 16, 0.9, 3)):
        sub = df[df["type"] == t]
        counts[t] = len(sub)
        ax.scatter(sub["UAR"], sub["Z"], s=s, c=POINT[t], alpha=al,
                   lw=0.3 if t != "XX" else 0, edgecolors="white", zorder=z)

    for x in (a, b, c, d):
        ax.axvline(x, color=LINE, lw=1.2, zorder=4)
    for y in (zlo, zhi):
        ax.axhline(y, color=LINE, lw=1.2, ls=(0, (5, 2.5)), zorder=4)

    _vlabels(ax, [(a, "xo_ur_x_max", _fmt(a)), (b, "ur_x_low", _fmt(b)),
                  (c, "ur_x_high", _fmt(c)), (d, "xxx_ur_x_min", _fmt(d))])
    _hlabels(ax, [(zhi, "xxx_z_threshold", _fmt(zhi, 1)),
                  (zlo, "xo_z_threshold", _fmt(zlo, 1))])

    xl, xr = (xlim[0] + a) / 2, (d + xlim[1]) / 2
    _zone_text(ax, xl, (ylim[0] + zlo) / 2, "XO\nDetected")
    _zone_text(ax, xl, (zlo + zhi) / 2, "XO\nSuspected", weight="normal")
    _zone_text(ax, xr, (zhi + ylim[1]) / 2, "XXX\nDetected")
    _zone_text(ax, xr, (zlo + zhi) / 2, "XXX\nSuspected", weight="normal")
    _zone_text(ax, (b + c) / 2, zhi - 1.4, "Not Detected", size=7.8,
               weight="normal")
    _zone_text(ax, xl, (zhi + ylim[1]) / 2, "Discordant\n\u2192 Not Detected",
               size=8, weight="normal", color="#5F6B76")
    _zone_text(ax, xr, (ylim[0] + zlo) / 2, "Discordant\n\u2192 Not Detected",
               size=8, weight="normal", color="#5F6B76")

    ax.set_xlabel("UAR[X] (%)", fontsize=10)
    ax.set_ylabel("Z-score[X]", fontsize=10)

    _legends(fig, [
        _pt(NORMAL_PT, f"XX reference (n={counts['XX']:,})", 4, 0.6),
        _pt(POINT["XO"], f"Confirmed XO (n={counts['XO']})"),
        _pt(POINT["XXX"], f"Confirmed XXX (n={counts['XXX']})"),
    ], [
        _zp(ZONE["XO Detected"], "XO Detected"),
        _zp(ZONE["XO Suspected"], "XO Suspected"),
        _zp(ZONE["XXX Detected"], "XXX Detected"),
        _zp(ZONE["XXX Suspected"], "XXX Suspected"),
        _zp(ND_FILL, "Not Detected"),
        _zp(ND_FILL, "Discordant", hatch="////"),
    ])

    _rules(fig, [
        f"\u2022 Z-score[X]:  Z < {_fmt(zlo, 1)} \u2192 XO signal;   "
        f"Z > {_fmt(zhi, 1)} \u2192 XXX signal;   otherwise neutral.",
        f"\u2022 UAR[X]:  < {_fmt(a)} XO-detected  |  {_fmt(a)}\u2013{_fmt(b)} "
        f"XO-suspected  |  {_fmt(b)}\u2013{_fmt(c)} normal",
        f"               {_fmt(c)}\u2013{_fmt(d)} XXX-suspected  |  "
        f"\u2265 {_fmt(d)} XXX-detected",
        "\u2022 Detected only when both axes are in the Detected zone; "
        "one axis alone \u2192 Suspected.",
        "\u2022 Opposite signals on the two axes (e.g. low Z with high UAR) "
        "\u2192 Not Detected.",
    ])
    _save(fig, out_dir, f"sca_female_xo_xxx_{group}")


# ── 2. male XXY / XYY ───────────────────────────────────────────────────────
def plot_male(det, group, ref_csv, out_dir, xlim=(4.0, 6.0), ylim=(0.0, 0.16)):
    p = det.get_male_params(group)
    m, k, thr, mg = p["slope"], p["intercept"], p["ur_x_threshold"], p["margin"]

    fig, ax = _new_fig(
        "Sex-chromosome aneuploidy (male): XXY / XYY",
        f"UAR[X] vs UAR[Y]  ·  {group} reference",
    )
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)

    xx = np.linspace(xlim[0], xlim[1], 400)
    b0 = m * xx + k
    b1 = b0 + mg

    def band(x_lo, x_hi, lower, upper):
        sel = (xx >= x_lo) & (xx <= x_hi)
        xs = xx[sel]
        lo = np.clip(lower[sel], *ylim)
        hi = np.clip(upper[sel], *ylim)
        return np.r_[np.c_[xs, lo], np.c_[xs[::-1], hi[::-1]]]

    top = np.full_like(xx, ylim[1])
    bot = np.full_like(xx, ylim[0])
    ax.add_patch(Polygon(band(*xlim, bot, b0), fc=ZONE["Not Detected"],
                         lw=0, zorder=0))
    for x_lo, x_hi, kind in ((xlim[0], thr, "XYY"), (thr, xlim[1], "XXY")):
        xm = (x_lo + x_hi) / 2
        assert det.detect_male_sca(xm, m * xm + k + mg / 2, group) == \
            f"{kind} Suspected"
        assert det.detect_male_sca(xm, m * xm + k + mg * 3, group) == \
            f"{kind} Detected"
        ax.add_patch(Polygon(band(x_lo, x_hi, b0, b1),
                             fc=ZONE[f"{kind} Suspected"], lw=0, zorder=0))
        ax.add_patch(Polygon(band(x_lo, x_hi, b1, top),
                             fc=ZONE[f"{kind} Detected"], lw=0, zorder=0))

    df = pd.read_csv(ref_csv, sep="\t")
    xcol, ycol = "UAR.X", "UAR.Y"
    counts = {}
    for t, s, al, z in (("XY", 6, 0.30, 2), ("XXY", 16, 0.9, 3),
                        ("XYY", 16, 0.9, 3), ("XXYY", 22, 0.95, 3)):
        sub = df[df["type"] == t]
        counts[t] = len(sub)
        if len(sub):
            ax.scatter(sub[xcol], sub[ycol], s=s, c=POINT[t], alpha=al,
                       lw=0.3 if t != "XY" else 0, edgecolors="white",
                       zorder=z, marker="D" if t == "XXYY" else "o")

    ax.plot(xx, b0, color=LINE, lw=1.3, zorder=4)
    ax.plot(xx, b1, color=LINE, lw=1.0, ls=(0, (4, 2.5)), zorder=4)
    yb = m * thr + k
    ax.plot([thr, thr], [yb, ylim[1]], color=LINE, lw=1.2, zorder=4)

    _vlabels(ax, [(thr, "ur_x_threshold", _fmt(thr))])

    def along(x, y, text, dy_pts):
        p0 = ax.transData.transform((x, m * x + k))
        p1 = ax.transData.transform((x + 0.1, m * (x + 0.1) + k))
        ang = np.degrees(np.arctan2(p1[1] - p0[1], p1[0] - p0[0]))
        ax.annotate(text, xy=(x, y), xytext=(0, dy_pts),
                    textcoords="offset points", rotation=ang,
                    rotation_mode="anchor", ha="left", va="center",
                    fontsize=8.2, color=LINE, zorder=6)

    sign = "\u2212" if m < 0 else "+"
    along(4.06, m * 4.06 + k, "boundary", -8)
    along(4.06, m * 4.06 + k + mg, f"boundary + {mg:g}", 8)

    _zone_text(ax, 5.05, 0.135, "XYY Detected")
    _zone_text(ax, 5.68, 0.10, "XXY Detected")
    _zone_text(ax, 4.3, 0.03, "Not Detected", weight="normal")
    xs_ = 5.66
    ax.annotate("Suspected band\n(boundary < UAR[Y] \u2264 boundary + "
                f"{mg:g})", xy=(xs_, m * xs_ + k + mg / 2),
                xytext=(4.72, 0.006), fontsize=7.8, color="#37474F",
                ha="left", va="bottom", zorder=6,
                arrowprops=dict(arrowstyle="->", color="#5F6B76", lw=0.7),
                bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="none",
                          alpha=0.85))

    ax.set_xlabel("UAR[X] (%)", fontsize=10)
    ax.set_ylabel("UAR[Y] (%)", fontsize=10)

    points = [
        _pt(NORMAL_PT, f"XY reference (n={counts['XY']:,})", 4, 0.6),
        _pt(POINT["XXY"], f"Confirmed XXY (n={counts['XXY']})"),
        _pt(POINT["XYY"], f"Confirmed XYY (n={counts['XYY']})"),
    ]
    if counts.get("XXYY"):
        points.append(_pt(POINT["XXYY"], f"XXYY (n={counts['XXYY']})", 4.5,
                          marker="D"))
    _legends(fig, points, [
        _zp(ZONE["XYY Detected"], "XYY Detected"),
        _zp(ZONE["XYY Suspected"], "XYY Suspected"),
        _zp(ZONE["XXY Detected"], "XXY Detected"),
        _zp(ZONE["XXY Suspected"], "XXY Suspected"),
        _zp(ND_FILL, "Not Detected"),
    ])

    _rules(fig, [
        f"\u2022 boundary(x) = {sign}{abs(m):.4f}\u00b7UAR[X] + {k:.4f}   "
        "(fitted to the XY reference cloud)",
        f"\u2022 UAR[Y] > boundary + {mg:g} \u2192 Detected;   "
        f"boundary < UAR[Y] \u2264 boundary + {mg:g} \u2192 Suspected;   "
        "otherwise Not Detected.",
        f"\u2022 Type of an abnormal call:  UAR[X] \u2264 {_fmt(thr)} "
        f"\u2192 XYY;   UAR[X] > {_fmt(thr)} \u2192 XXY.",
    ])
    _save(fig, out_dir, f"sca_male_xxy_xyy_{group}")


# ── 3. trisomy (e.g. chr21) ─────────────────────────────────────────────────
def trisomy_call(uar, z, t):
    """Same rule as ezd_runner.run_decision."""
    if z >= t["Z_max"] and uar >= t["UAR_max"]:
        return "Detected"
    if t["Z_min"] <= z or t["UAR_min"] <= uar:
        return "Suspected"
    return "Not Detected"


def plot_trisomy(chrom, group, ref_csv, thr_tsv, out_dir):
    t = pd.read_csv(thr_tsv, sep="\t").set_index("chr").loc[chrom].to_dict()
    df = pd.read_csv(ref_csv, sep="\t")
    n, pos = df[df["type"] == "N"], df[df["type"] == "P"]
    tri = {"chr21": "T21 (Down)", "chr18": "T18 (Edwards)",
           "chr13": "T13 (Patau)"}.get(chrom, f"Trisomy {chrom[3:]}")

    xpad = (df["UAR"].max() - df["UAR"].min()) * 0.06
    xlim = (df["UAR"].min() - xpad, df["UAR"].max() + xpad)
    ylim = (min(-5.0, df["Z"].min() - 4.5), max(35.0, df["Z"].max() + 2))

    fig, ax = _new_fig(
        f"Autosomal trisomy: {tri}",
        f"{chrom}  ·  UAR vs Z-score  ·  {group} reference",
    )
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)

    xs = sorted({xlim[0], t["UAR_min"], t["UAR_max"], xlim[1]})
    ys = sorted({ylim[0], t["Z_min"], t["Z_max"], ylim[1]})
    for i in range(len(xs) - 1):
        for j in range(len(ys) - 1):
            call = trisomy_call((xs[i] + xs[i + 1]) / 2,
                                (ys[j] + ys[j + 1]) / 2, t)
            ax.add_patch(Rectangle((xs[i], ys[j]), xs[i + 1] - xs[i],
                                   ys[j + 1] - ys[j], fc=ZONE[call], lw=0,
                                   zorder=0))

    ax.scatter(n["UAR"], n["Z"], s=6, c=NORMAL_PT, alpha=0.30, lw=0, zorder=2)
    ax.scatter(pos["UAR"], pos["Z"], s=18, c=POINT["P"], alpha=0.9, lw=0.3,
               edgecolors="white", zorder=3)

    for x in (t["UAR_min"], t["UAR_max"]):
        ax.axvline(x, color=LINE, lw=1.2, zorder=4)
    for y in (t["Z_min"], t["Z_max"]):
        ax.axhline(y, color=LINE, lw=1.2, ls=(0, (5, 2.5)), zorder=4)

    _vlabels(ax, [(t["UAR_min"], "UAR_min", _fmt(t["UAR_min"])),
                  (t["UAR_max"], "UAR_max", _fmt(t["UAR_max"]))])
    _hlabels(ax, [(t["Z_max"], "Z_max", _fmt(t["Z_max"], 1)),
                  (t["Z_min"], "Z_min", _fmt(t["Z_min"], 1))])

    _zone_text(ax, t["UAR_max"] + (xlim[1] - t["UAR_max"]) * 0.3,
               ylim[1] - (ylim[1] - t["Z_max"]) * 0.15,
               f"{tri.split()[0]}\nDetected")
    _zone_text(ax, (xlim[0] + t["UAR_min"]) / 2, (t["Z_max"] + ylim[1]) / 2,
               "Suspected", weight="normal")
    _zone_text(ax, (t["UAR_max"] + xlim[1]) / 2, (ylim[0] + t["Z_min"]) / 2,
               "Suspected", weight="normal")
    _zone_text(ax, (xlim[0] + t["UAR_min"]) / 2, ylim[0] + 1.6,
               "Not Detected", size=8.5, weight="normal")

    ax.set_xlabel(f"UAR[{chrom[3:]}] (%)", fontsize=10)
    ax.set_ylabel(f"Z-score[{chrom[3:]}]", fontsize=10)

    calls = {k: df[df["type"] == k].apply(
        lambda r: trisomy_call(r["UAR"], r["Z"], t), axis=1).value_counts()
        for k in ("N", "P")}
    _legends(fig, [
        _pt(NORMAL_PT, f"Euploid reference (n={len(n):,})", 4, 0.6),
        _pt(POINT["P"], f"Confirmed {tri.split()[0]} (n={len(pos)})"),
    ], [
        _zp(ZONE["Detected"], f"{tri.split()[0]} Detected"),
        _zp(ZONE["Suspected"], "Suspected"),
        _zp(ND_FILL, "Not Detected"),
    ])

    _rules(fig, [
        f"\u2022 Detected:  Z \u2265 Z_max ({_fmt(t['Z_max'], 1)})  AND  "
        f"UAR \u2265 UAR_max ({_fmt(t['UAR_max'])})",
        f"\u2022 Suspected:  Z \u2265 Z_min ({_fmt(t['Z_min'], 1)})  OR  "
        f"UAR \u2265 UAR_min ({_fmt(t['UAR_min'])}),  when not Detected",
        f"\u2022 Not Detected:  Z < {_fmt(t['Z_min'], 1)}  AND  "
        f"UAR < {_fmt(t['UAR_min'])}",
        f"\u2022 Reference check:  euploid max UAR {n['UAR'].max():.3f}, "
        f"max Z {n['Z'].max():.2f}  |  {tri.split()[0]} min UAR "
        f"{pos['UAR'].min():.3f}, min Z {pos['Z'].min():.2f}",
        f"    euploid \u2192 Not Detected {calls['N'].get('Not Detected', 0):,}"
        f"/{len(n):,};   {tri.split()[0]} \u2192 Detected "
        f"{calls['P'].get('Detected', 0)}/{len(pos)}",
    ])
    _save(fig, out_dir, f"trisomy_{chrom}_{group}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--ezd-dir", required=True,
                    help="<ref_dir>/labs/<labcode>/EZD")
    ap.add_argument("--group", default="orig", choices=["orig", "fetus"])
    ap.add_argument("--chrom", default="chr21")
    ap.add_argument("--out", default=os.path.join(REPO, "docs", "whitepaper",
                                                  "figures"))
    args = ap.parse_args()

    gdir = os.path.join(args.ezd_dir, args.group)
    det = SCADetector()
    with open(os.path.join(gdir, "sca_config.json")) as fh:
        det.configs[args.group] = json.load(fh)

    plot_female(det, args.group, os.path.join(gdir, "female.txt"), args.out)
    plot_male(det, args.group, os.path.join(gdir, "male.txt"), args.out)
    plot_trisomy(args.chrom, args.group, os.path.join(gdir, f"{args.chrom}.txt"),
                 os.path.join(gdir, f"{args.group}_thresholds_new.tsv"),
                 args.out)


if __name__ == "__main__":
    main()
