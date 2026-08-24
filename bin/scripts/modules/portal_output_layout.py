#!/usr/bin/env python3
"""
Align gx-nipt output tree with ken-nipt / Portal flat paths.

Portal reads plot paths from result JSON, e.g.:
  Output_EZD/orig_EZD_grid.png
  Output_PRIZM/{sample}_orig_chromosome_line.png
  Output_PRIZM/{sample}_orig_10mb_line.png
  Output_WC/{sample}.wc.orig_z.png
  Output_WCX/{sample}.wcx.orig.png

gx-nipt analysis keeps group subdirs (Output_EZD/orig/...). This script
hoists Portal-facing files to the flat ken-nipt layout, re-syncs gxcnv1/2
from analysis (in case COPY raced ahead of PLOT publishDir), then writes
gxcnv1/2 sidecars beside native WC/WCX flats (never overwriting natives).
"""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

# Same-directory import when invoked as a script from Nextflow / run_nipt.sh
sys.path.insert(0, str(Path(__file__).resolve().parent))
from portal_cnv_alias import inject_portal_aliases  # noqa: E402

GROUPS = ("orig", "fetus", "mom")


def _copy_file(src: Path, dst: Path) -> bool:
    if not src.is_file():
        return False
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists() and src.samefile(dst):
        return False
    shutil.copy2(src, dst)
    return True


def _resolve_gxcnv_dir(analysis_dir: Path, sample: str, name: str) -> Path | None:
    """Prefer nested analysis/<sample>/<sample>/gxcnvN (publishDir layout)."""
    candidates = [
        analysis_dir / sample / name,
        analysis_dir / name,
        # already-copied output tree as last resort
    ]
    for c in candidates:
        if c.is_dir():
            return c
    return None


def sync_gxcnv_trees(analysis_dir: Path, sample: str, outdir: Path) -> list[str]:
    """Copy any missing gxcnv1/2 files from analysis into outdir."""
    written: list[str] = []
    for name in ("gxcnv1", "gxcnv2"):
        src_dir = _resolve_gxcnv_dir(analysis_dir, sample, name)
        if src_dir is None:
            continue
        dst_dir = outdir / name
        dst_dir.mkdir(parents=True, exist_ok=True)
        for src in src_dir.iterdir():
            if not src.is_file():
                continue
            dst = dst_dir / src.name
            if dst.is_file() and dst.stat().st_size > 0:
                # Prefer newer / larger analysis artefact when output is stale stub
                if dst.stat().st_mtime >= src.stat().st_mtime and dst.stat().st_size >= src.stat().st_size:
                    continue
            if _copy_file(src, dst):
                written.append(str(dst))
    return written


def sync_wcx_plot_dirs(analysis_dir: Path, sample: str, outdir: Path) -> list[str]:
    """Ensure WCX --plot dirs land under outdir/Output_WCX/{group}/ before flatten.

    publishDir + COPY_TO_OUTPUT usually already copy them; this covers races /
    partial copies so chr_plots packing still works.
    """
    written: list[str] = []
    for group in GROUPS:
        plots_name = f"{sample}.wcx.{group}.plots"
        candidates = [
            analysis_dir / sample / "Output_WCX" / group / plots_name,
            analysis_dir / "Output_WCX" / group / plots_name,
            outdir / "Output_WCX" / group / plots_name,
        ]
        src_dir = next((c for c in candidates if c.is_dir()), None)
        if src_dir is None:
            continue
        dst_dir = outdir / "Output_WCX" / group / plots_name
        if src_dir.resolve() != dst_dir.resolve():
            dst_dir.mkdir(parents=True, exist_ok=True)
            for src in src_dir.iterdir():
                if src.is_file() and _copy_file(src, dst_dir / src.name):
                    written.append(str(dst_dir / src.name))
    return written


def flatten_portal_plots(outdir: Path, sample: str) -> list[str]:
    """Hoist Portal plot/report files to ken-nipt flat paths."""
    written: list[str] = []

    # ── EZD ──────────────────────────────────────────────────────────────
    for group in GROUPS:
        src = outdir / "Output_EZD" / group / f"{group}_EZD_grid.png"
        dst = outdir / "Output_EZD" / f"{group}_EZD_grid.png"
        if _copy_file(src, dst):
            written.append(str(dst))

    # ── PRIZM (all portal PNGs for the group) ────────────────────────────
    for group in GROUPS:
        src_dir = outdir / "Output_PRIZM" / group
        if not src_dir.is_dir():
            continue
        for src in sorted(src_dir.glob(f"{sample}_{group}_*.png")):
            dst = outdir / "Output_PRIZM" / src.name
            if _copy_file(src, dst):
                written.append(str(dst))

    # ── WC flat hoist (native WC only; gxcnv1 uses sidecars / gxcnv1/) ─
    for group in GROUPS:
        for fname in (
            f"{sample}.wc.{group}_z.png",
            f"{sample}.wc.{group}.report.txt",
        ):
            src = outdir / "Output_WC" / group / fname
            dst = outdir / "Output_WC" / fname
            if _copy_file(src, dst):
                written.append(str(dst))

    # ── WCX flat hoist + native genome_wide rename ───────────────────────
    for group in GROUPS:
        for fname in (
            f"{sample}.wcx.{group}.png",
            f"{sample}.wcx.{group}_aberrations.bed",
        ):
            src = outdir / "Output_WCX" / group / fname
            dst = outdir / "Output_WCX" / fname
            if _copy_file(src, dst):
                written.append(str(dst))

        flat_png = outdir / "Output_WCX" / f"{sample}.wcx.{group}.png"
        if not flat_png.is_file():
            candidates = [
                outdir
                / "Output_WCX"
                / group
                / f"{sample}.wcx.{group}.plots"
                / "genome_wide.png",
                outdir / "Output_WCX" / group / "genome_wide.png",
                outdir / "Output_WCX" / "chr_plots" / group / "genome_wide.png",
            ]
            for src in candidates:
                if _copy_file(src, flat_png):
                    written.append(str(flat_png))
                    break

        # Microdeletion Zoom-in: Portal JSON → Output_WCX/chr_plots/{group}/chrN.png
        # Source is WisecondorX --plot output: {sample}.wcx.{group}.plots/chr*.png
        plots_src = (
            outdir / "Output_WCX" / group / f"{sample}.wcx.{group}.plots"
        )
        plots_dst = outdir / "Output_WCX" / "chr_plots" / group
        if plots_src.is_dir():
            for src in sorted(plots_src.glob("chr*.png")):
                dst = plots_dst / src.name
                if _copy_file(src, dst):
                    written.append(str(dst))

    return written


def apply_portal_layout(sample: str, analysis_dir: Path, outdir: Path) -> dict:
    for d in (
        outdir / "Output_EZD",
        outdir / "Output_PRIZM",
        outdir / "Output_WC",
        outdir / "Output_WCX",
        outdir / "gxcnv1",
        outdir / "gxcnv2",
    ):
        d.mkdir(parents=True, exist_ok=True)

    synced = sync_gxcnv_trees(analysis_dir, sample, outdir)
    synced.extend(sync_wcx_plot_dirs(analysis_dir, sample, outdir))
    flattened = flatten_portal_plots(outdir, sample)
    alias = inject_portal_aliases(sample, analysis_dir, outdir)

    return {
        "synced": synced,
        "flattened": flattened,
        "alias_wc": alias.get("wc", []),
        "alias_wcx": alias.get("wcx", []),
        "alias_skipped": alias.get("skipped", []),
    }


def verify_portal_paths(outdir: Path, sample: str) -> list[str]:
    """Return missing Portal plot paths (relative)."""
    missing: list[str] = []
    for group in GROUPS:
        expected = [
            f"Output_EZD/{group}_EZD_grid.png",
            f"Output_PRIZM/{sample}_{group}_chromosome_line.png",
            f"Output_PRIZM/{sample}_{group}_10mb_line.png",
            f"Output_PRIZM/{sample}_{group}_chromosome_heatmap.png",
            f"Output_PRIZM/{sample}_{group}_10mb_heatmap.png",
            f"Output_WC/{sample}.wc.{group}_z.png",
            f"Output_WCX/{sample}.wcx.{group}.png",
            f"gxcnv1/{sample}_{group}_genome.png",
            f"gxcnv2/{sample}_{group}_genome.png",
        ]
        for rel in expected:
            if not (outdir / rel).is_file():
                missing.append(rel)
    return missing


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description="Flatten gx-nipt outputs to ken-nipt / Portal paths"
    )
    ap.add_argument("--sample", required=True)
    ap.add_argument(
        "--analysis-dir",
        required=True,
        type=Path,
        help="Per-sample analysis dir (…/analysis/<work>/<order_id>)",
    )
    ap.add_argument(
        "--outdir",
        required=True,
        type=Path,
        help="Portal output directory (…/output/<work>/<order_id>)",
    )
    args = ap.parse_args(argv)

    if not args.outdir.is_dir():
        print(f"[portal_output_layout] ERROR: outdir missing: {args.outdir}", file=sys.stderr)
        return 2

    summary = apply_portal_layout(args.sample, args.analysis_dir, args.outdir)
    missing = verify_portal_paths(args.outdir, args.sample)
    print(
        f"[portal_output_layout] sample={args.sample} "
        f"synced={len(summary['synced'])} flattened={len(summary['flattened'])} "
        f"alias_wc={len(summary['alias_wc'])} alias_wcx={len(summary['alias_wcx'])} "
        f"missing={len(missing)}"
    )
    for p in summary["flattened"][:20]:
        print(f"  flat {p}")
    if len(summary["flattened"]) > 20:
        print(f"  ... +{len(summary['flattened']) - 20} more")
    for p in missing:
        print(f"  MISSING {p}", file=sys.stderr)
    # Missing plots are warnings — pipeline may have skipped a group.
    return 0


if __name__ == "__main__":
    sys.exit(main())
