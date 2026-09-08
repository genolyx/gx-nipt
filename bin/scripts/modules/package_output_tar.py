#!/usr/bin/env python3
"""
Build a Portal-facing <order_id>.output.tar.

Goals vs naive ``tar Output_* gxcnv*``:
  - Omit Output_hmmcopy (HMMcopy intermediate; analysis/ keeps it for PRIZM).
  - Prefer ken-nipt flat layout: top-level files under Output_EZD / PRIZM / WC
    (nested orig|fetus|mom copies are duplicates after portal_output_layout).
  - Output_WCX: flat files + chr_plots/ (Portal microdeletion zoom).
  - gxcnv1 / gxcnv2: PNG plots + *_calls.tsv (JSON result paths; drop bins/segments).
"""

from __future__ import annotations

import argparse
import os
import tarfile
from pathlib import Path

# Never ship to Portal archive
_EXCLUDE_TOP_DIRS = frozenset({"Output_hmmcopy"})

# Only top-level files (no nested group trees)
_FLAT_ONLY_DIRS = frozenset(
    {
        "Output_EZD",
        "Output_PRIZM",
        "Output_WC",
        "Output_QC",
        "Output_FF",
        "Output_MD",
        "Output_Result",
    }
)

_GXCNV_DIRS = frozenset({"gxcnv1", "gxcnv2"})
_GROUPS = ("orig", "fetus", "mom")


def _same_bytes(a: Path, b: Path) -> bool:
    if not a.is_file() or not b.is_file():
        return False
    if a.stat().st_size != b.stat().st_size:
        return False
    return a.read_bytes() == b.read_bytes()


def _is_gxcnv_clone(path: Path, gxcnv: Path) -> bool:
    return gxcnv.is_file() and _same_bytes(path, gxcnv)


def _native_wc_png(outdir: Path, order_id: str, group: str) -> Path | None:
    gx1 = outdir / "gxcnv1" / f"{order_id}_{group}_genome.png"
    nested = outdir / "Output_WC" / group / f"{order_id}.wc.{group}_z.png"
    flat = outdir / "Output_WC" / f"{order_id}.wc.{group}_z.png"
    for p in (nested, flat):
        if p.is_file() and p.stat().st_size > 0 and not _is_gxcnv_clone(p, gx1):
            return p
    return None


def _native_wcx_png(outdir: Path, order_id: str, group: str) -> Path | None:
    gx2 = outdir / "gxcnv2" / f"{order_id}_{group}_genome.png"
    candidates = [
        outdir / "Output_WCX" / group / f"{order_id}.wcx.{group}.plots" / "genome_wide.png",
        outdir / "Output_WCX" / group / "genome_wide.png",
        outdir / "Output_WCX" / "chr_plots" / group / "genome_wide.png",
        outdir / "Output_WCX" / group / f"{order_id}.wcx.{group}.png",
        outdir / "Output_WCX" / f"{order_id}.wcx.{group}.png",
    ]
    for p in candidates:
        if p.is_file() and p.stat().st_size > 0 and not _is_gxcnv_clone(p, gx2):
            return p
    return None


def _iter_archive_members(outdir: Path, order_id: str) -> list[tuple[Path, str]]:
    """Return (realpath, archive name) pairs."""
    members: list[tuple[Path, str]] = []
    skip_flats: set[str] = set()

    json_path = outdir / f"{order_id}.json"
    if json_path.is_file():
        members.append((json_path, json_path.name))

    # Prefer native Wisecondor / WisecondorX images under Portal flat names.
    for group in _GROUPS:
        wc_png = _native_wc_png(outdir, order_id, group)
        wc_arc = f"Output_WC/{order_id}.wc.{group}_z.png"
        if wc_png is not None:
            members.append((wc_png, wc_arc))
        skip_flats.add(wc_arc)

        wcx_png = _native_wcx_png(outdir, order_id, group)
        wcx_arc = f"Output_WCX/{order_id}.wcx.{group}.png"
        if wcx_png is not None:
            members.append((wcx_png, wcx_arc))
        skip_flats.add(wcx_arc)

    for child in sorted(outdir.iterdir()):
        name = child.name
        if name in _EXCLUDE_TOP_DIRS:
            continue

        if child.is_file():
            continue

        if not child.is_dir():
            continue

        if name in _FLAT_ONLY_DIRS:
            for f in sorted(child.iterdir()):
                if not f.is_file():
                    continue
                arc = f"{name}/{f.name}"
                if arc in skip_flats:
                    continue
                members.append((f, arc))
            continue

        if name == "Output_WCX":
            for f in sorted(child.iterdir()):
                if not f.is_file():
                    continue
                arc = f"{name}/{f.name}"
                if arc in skip_flats:
                    continue
                members.append((f, arc))
            chr_plots = child / "chr_plots"
            if chr_plots.is_dir():
                for f in sorted(chr_plots.rglob("*")):
                    if f.is_file():
                        members.append((f, f.relative_to(outdir).as_posix()))
            continue

        if name in _GXCNV_DIRS:
            for f in sorted(child.iterdir()):
                if not f.is_file():
                    continue
                # Portal JSON points at *_calls.tsv; keep PNGs + calls only.
                if f.suffix.lower() == ".png" or f.name.endswith("_calls.tsv"):
                    members.append((f, f"{name}/{f.name}"))
            continue

        if name.startswith("Output_"):
            for f in sorted(child.iterdir()):
                if f.is_file():
                    members.append((f, f"{name}/{f.name}"))

    seen_arc: set[str] = set()
    seen_path: set[Path] = set()
    out: list[tuple[Path, str]] = []
    for path, arc in members:
        rp = path.resolve()
        if arc in seen_arc or rp in seen_path:
            continue
        seen_arc.add(arc)
        seen_path.add(rp)
        out.append((path, arc))
    return out


def build_tar(outdir: Path, order_id: str, tar_path: Path) -> dict:
    members = _iter_archive_members(outdir, order_id)
    if not members:
        raise RuntimeError(f"no files to archive under {outdir}")

    tar_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = tar_path.with_suffix(tar_path.suffix + ".tmp")
    if tmp.exists():
        tmp.unlink()

    bytes_sum = 0
    with tarfile.open(tmp, "w") as tf:
        for path, arcname in members:
            tf.add(path, arcname=arcname)
            bytes_sum += path.stat().st_size

    os.replace(tmp, tar_path)
    tar_size = tar_path.stat().st_size
    return {
        "members": len(members),
        "payload_bytes": bytes_sum,
        "tar_bytes": tar_size,
        "tar_path": str(tar_path),
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--outdir", required=True, help="output/<work>/<order> directory")
    ap.add_argument("--order-id", required=True)
    ap.add_argument(
        "--tar-path",
        default="",
        help="default: <outdir>/<order-id>.output.tar",
    )
    args = ap.parse_args()

    outdir = Path(args.outdir).resolve()
    order_id = args.order_id.strip()
    tar_path = (
        Path(args.tar_path).resolve()
        if args.tar_path.strip()
        else outdir / f"{order_id}.output.tar"
    )

    if not outdir.is_dir():
        raise SystemExit(f"outdir not found: {outdir}")

    summary = build_tar(outdir, order_id, tar_path)
    print(
        f"[package_output_tar] {summary['tar_path']} "
        f"members={summary['members']} "
        f"payload={summary['payload_bytes'] / (1024 * 1024):.1f}MB "
        f"tar={summary['tar_bytes'] / (1024 * 1024):.1f}MB"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
