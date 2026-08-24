#!/usr/bin/env python3
"""
Place gxcnv1 / gxcnv2 artefacts beside native WC / WCX (do not overwrite).

Portal needs both image sets:
  - Native WC / WCX plots (Output_WC / Output_WCX flat ken-nipt names)
  - gxcnv1 / gxcnv2 plots (under gxcnv1/ and gxcnv2/, plus optional sidecars)

Sidecar flat names (optional convenience for Portal under Output_WC/WCX):
  gxcnv1/{sample}_{group}_genome.png
      → Output_WC/{sample}.gxcnv1.{group}_z.png
  gxcnv1/{sample}_{group}_calls.tsv
      → Output_WC/{sample}.gxcnv1.{group}.report.txt
  gxcnv2/{sample}_{group}_genome.png
      → Output_WCX/{sample}.gxcnv2.{group}.png
  gxcnv2/{sample}_{group}_calls.tsv
      → Output_WCX/{sample}.gxcnv2.{group}_aberrations.bed

Native flats (``{sample}.wc.{group}_z.png``, ``{sample}.wcx.{group}.png``, …)
are never replaced.
"""

from __future__ import annotations

import argparse
import math
import shutil
import sys
from pathlib import Path


GROUPS = ("orig", "fetus", "mom")


def _strip_chr(chrom: str) -> str:
    chrom = chrom.strip()
    if chrom.lower().startswith("chr"):
        return chrom[3:]
    return chrom


def _read_calls(path: Path) -> list[dict]:
    """Parse gxcnv1/2 calls.tsv (header may be '#chrom\\t...')."""
    default_header = [
        "chrom",
        "start",
        "end",
        "type",
        "mean_log2_ratio",
        "mean_z",
        "mean_mad_z",
        "n_bins",
    ]
    rows: list[dict] = []
    header = None
    with path.open(encoding="utf-8") as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if not line.strip():
                continue
            parts = line.split("\t")
            first = parts[0].lstrip("#").lower()
            if header is None and first == "chrom":
                header = [c.lstrip("#") for c in parts]
                continue
            if header is None:
                header = default_header
            if len(parts) < 6:
                continue
            rec = {header[i]: parts[i] for i in range(min(len(header), len(parts)))}
            if "chrom" not in rec:
                for k in list(rec):
                    if k.lstrip("#") == "chrom":
                        rec["chrom"] = rec[k]
                        break
            rows.append(rec)
    return rows


def calls_to_wc_report(calls_path: Path, out_path: Path) -> None:
    """Write minimal WC-style report.txt from gxcnv1 calls.tsv."""
    rows = _read_calls(calls_path)
    lines = ["# Test results: #", "z-score\teffect\tmbsize\tlocation"]
    for r in rows:
        try:
            log2 = float(r["mean_log2_ratio"])
            z = float(r["mean_z"])
            start = int(float(r["start"]))
            end = int(float(r["end"]))
        except (KeyError, ValueError):
            continue
        effect = (math.pow(2.0, log2) - 1.0) * 100.0
        mbsize = (end - start) / 1e6
        chrom = _strip_chr(r.get("chrom", ""))
        location = f"{chrom}:{start}-{end}"
        lines.append(f"{z:.2f}\t{effect:.2f}\t{mbsize:.2f}\t{location}")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def calls_to_wcx_bed(calls_path: Path, out_path: Path) -> None:
    """Write WCX-style aberrations.bed from gxcnv2 calls.tsv."""
    rows = _read_calls(calls_path)
    lines = ["chr\tstart\tend\tratio\tzscore\ttype"]
    for r in rows:
        try:
            log2 = float(r["mean_log2_ratio"])
            z = float(r["mean_z"])
            start = int(float(r["start"]))
            end = int(float(r["end"]))
        except (KeyError, ValueError):
            continue
        ratio = math.pow(2.0, log2) - 1.0
        chrom = _strip_chr(r.get("chrom", ""))
        ctype = str(r.get("type", "")).strip().lower()
        if ctype not in ("gain", "loss"):
            ctype = ctype.lower()
        lines.append(f"{chrom}\t{start}\t{end}\t{ratio}\t{z}\t{ctype}")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _copy_if_exists(src: Path, dst: Path) -> bool:
    if not src.is_file():
        return False
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)
    return True


def inject_portal_aliases(sample: str, analysis_dir: Path, outdir: Path) -> dict:
    """
    Write gxcnv1/2 sidecars under Output_WC / Output_WCX without touching
    native WC/WCX flat filenames.
    """
    sample_root = analysis_dir / sample
    gxcnv1 = sample_root / "gxcnv1"
    gxcnv2 = sample_root / "gxcnv2"
    if not gxcnv1.is_dir() and (outdir / "gxcnv1").is_dir():
        gxcnv1 = outdir / "gxcnv1"
    if not gxcnv2.is_dir() and (outdir / "gxcnv2").is_dir():
        gxcnv2 = outdir / "gxcnv2"

    summary = {"wc": [], "wcx": [], "skipped": []}

    for group in GROUPS:
        # ── gxcnv1 sidecars (native WC flats untouched) ──────────────
        genome1 = gxcnv1 / f"{sample}_{group}_genome.png"
        calls1 = gxcnv1 / f"{sample}_{group}_calls.tsv"
        gx1_png = outdir / "Output_WC" / f"{sample}.gxcnv1.{group}_z.png"
        gx1_report = outdir / "Output_WC" / f"{sample}.gxcnv1.{group}.report.txt"

        if genome1.is_file() or calls1.is_file():
            if _copy_if_exists(genome1, gx1_png):
                summary["wc"].append(str(gx1_png))
            if calls1.is_file():
                calls_to_wc_report(calls1, gx1_report)
                summary["wc"].append(str(gx1_report))
        else:
            summary["skipped"].append(f"gxcnv1/{group}")

        # ── gxcnv2 sidecars (native WCX flats untouched) ─────────────
        genome2 = gxcnv2 / f"{sample}_{group}_genome.png"
        calls2 = gxcnv2 / f"{sample}_{group}_calls.tsv"
        gx2_png = outdir / "Output_WCX" / f"{sample}.gxcnv2.{group}.png"
        gx2_bed = outdir / "Output_WCX" / f"{sample}.gxcnv2.{group}_aberrations.bed"

        if genome2.is_file() or calls2.is_file():
            if _copy_if_exists(genome2, gx2_png):
                summary["wcx"].append(str(gx2_png))
            if calls2.is_file():
                calls_to_wcx_bed(calls2, gx2_bed)
                summary["wcx"].append(str(gx2_bed))
        else:
            summary["skipped"].append(f"gxcnv2/{group}")

    return summary


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description="Write gxcnv1/2 sidecars beside native WC/WCX (no overwrite)"
    )
    ap.add_argument("--sample", required=True)
    ap.add_argument(
        "--analysis-dir",
        required=True,
        type=Path,
        help="Analysis root containing <sample>/gxcnv1 and <sample>/gxcnv2",
    )
    ap.add_argument(
        "--outdir",
        required=True,
        type=Path,
        help="Portal output directory (contains Output_WC / Output_WCX)",
    )
    args = ap.parse_args(argv)

    summary = inject_portal_aliases(args.sample, args.analysis_dir, args.outdir)
    print(
        f"[portal_cnv_alias] sample={args.sample} "
        f"gxcnv1_sidecar={len(summary['wc'])} gxcnv2_sidecar={len(summary['wcx'])} "
        f"skipped={summary['skipped']}"
    )
    for p in summary["wc"] + summary["wcx"]:
        print(f"  wrote {p}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
