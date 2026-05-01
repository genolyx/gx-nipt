#!/usr/bin/env python3
"""
gxcnv_build_ref_hg19.py
=======================
Build a gxcnv reference panel with hg19 clinical target region coordinates.

The gxcnv CLI `gxcnv newref` uses hg38 default regions and does not expose a
--target-regions option.  This wrapper calls the Python API directly with hg19
liftover coordinates so the reference matches the hg19-aligned BAMs used in
gx-nipt.

Usage
-----
    python gxcnv_build_ref_hg19.py \
        /path/to/normal_samples/*.gxcnv.npz \
        -o /path/to/refs/labs/cordlife/GXCNV/reference.npz

    # Or point at a directory of NPZ files:
    python gxcnv_build_ref_hg19.py \
        /path/to/normal_npz_dir/ \
        -o /path/to/output/reference.npz

Step 0 (convert BAMs to NPZ) — run once per normal sample:
    gxcnv convert sample.bam sample.gxcnv.npz --bin-size 100000

Recommended: ≥ 50 normal female samples for a robust reference panel.
Mixed-sex panels are supported; sex is predicted automatically via GMM.
"""

import argparse
import glob
import logging
import os
import sys

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [gxcnv_build_ref_hg19] %(levelname)s %(message)s",
)
logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# hg19 (GRCh37) clinical target regions
# Approximate UCSC liftover of the default hg38 regions in gxcnv newref.py.
# Each entry: (name, chrom, start, end)
# ---------------------------------------------------------------------------
HG19_TARGET_REGIONS = [
    ("DiGeorge_22q11",      "chr22", 18_648_855, 21_802_889),
    ("Williams_7q11",       "chr7",  72_744_454, 74_142_513),
    ("Angelman_15q11",      "chr15", 22_765_628, 28_530_576),
    ("PraderWilli_15q11",   "chr15", 22_765_628, 28_530_576),
    ("Wolf_4p16",           "chr4",  1_000,      4_000_000),
    ("CriDuChat_5p15",      "chr5",  1_000,      11_800_000),
    ("NF1_17q11",           "chr17", 29_421_945, 30_144_979),
    ("Smith_17p11",         "chr17", 16_755_133, 20_519_016),
    ("Langer_Giedion_8q24", "chr8",  116_600_000, 119_400_000),
    ("Miller_Dieker_17p13", "chr17", 1_000,       4_000_000),
    ("CHARGE_8q12",         "chr8",  61_720_000,  62_020_000),
    ("Kabuki_12q13",        "chr12", 49_400_000,  50_400_000),
    ("Sotos_5q35",          "chr5",  175_100_000, 177_200_000),
    ("Rubinstein_16p13",    "chr16", 3_800_000,   4_000_000),
    ("Potocki_Lupski_17p11","chr17", 15_000_000,  20_500_000),
]


def collect_npz_paths(inputs: list[str]) -> list[str]:
    paths = []
    for item in inputs:
        if os.path.isdir(item):
            found = sorted(glob.glob(os.path.join(item, "*.npz")))
            paths.extend(found)
        elif os.path.isfile(item):
            paths.append(item)
        else:
            expanded = sorted(glob.glob(item))
            if expanded:
                paths.extend(expanded)
            else:
                logger.warning("No files matched: %s", item)
    return paths


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build a gxcnv reference panel with hg19 clinical target regions."
    )
    parser.add_argument(
        "inputs",
        nargs="+",
        help="NPZ sample files or directories containing NPZ files (glob patterns ok)",
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        metavar="OUTPUT_NPZ",
        help="Output reference NPZ path (e.g. refs/labs/cordlife/GXCNV/reference.npz)",
    )
    parser.add_argument(
        "--pca-variance",
        type=float,
        default=0.95,
        metavar="F",
        help="Cumulative variance threshold for global PCA (default: 0.95)",
    )
    args = parser.parse_args()

    npz_paths = collect_npz_paths(args.inputs)
    if not npz_paths:
        logger.error("No NPZ files found from: %s", args.inputs)
        sys.exit(1)

    logger.info("Found %d NPZ files to use for reference panel.", len(npz_paths))
    for p in npz_paths:
        logger.info("  %s", p)

    try:
        from gxcnv.newref import build_reference
    except ImportError as e:
        logger.error("gxcnv is not installed: %s", e)
        logger.error("Install with: pip install git+https://github.com/genolyx/gx-cnv")
        sys.exit(1)

    out_path = args.output
    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)

    logger.info("Building reference panel with hg19 target regions -> %s", out_path)
    build_reference(
        npz_paths=npz_paths,
        output_path=out_path,
        target_regions=HG19_TARGET_REGIONS,
        global_pca_variance=args.pca_variance,
    )
    logger.info("Reference panel built successfully: %s", out_path)


if __name__ == "__main__":
    main()
