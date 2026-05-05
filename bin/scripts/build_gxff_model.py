#!/usr/bin/env python3
"""
build_gxff_model.py
===================
End-to-end gx-FF model builder for hg19 NIPT samples.

Pipeline steps (all in one run):
  1. Read sample data TSV (reference_sample_list_from_json_v3.tsv).
  2. Filter male-fetus samples by QC criteria (YFF_2 as ground-truth FF).
  3. Stratified sampling across FF bins and months.
  4. Generate training config TSV (FILEPATH=wig, BAM_PATH=BAM).
  5. Train gx-FF ensemble (coverage from wig + fragment from BAM).
  6. Save model + cross-validation metrics.

Usage
-----
    python build_gxff_model.py \\
        --sample-tsv  /path/to/reference_sample_list_from_json_v3.tsv \\
        --out-dir     /home/ken/gx-nipt/refs/gxff/model \\
        --n-samples   170 \\
        --cv-folds    5 \\
        --threads     8

Prerequisites
-------------
  • gxff must be installed (pip install git+https://github.com/genolyx/gx-FF).
  • The patched pipeline.py (supporting BAM_PATH column) must be in place:
      /opt/conda/envs/nipt/lib/python3.11/site-packages/gxff/core/pipeline.py
    This patch is located at: <repo>/bin/gxff_patch/pipeline.py
    To apply inside Docker, mount it as a bind-override before training.
"""

import argparse
import logging
import os
import sys
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [build_gxff] %(levelname)s %(message)s",
)
logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Column name constants
# ---------------------------------------------------------------------------
COL_SAMPLE_ID  = "sample_id"
COL_GENDER     = "fetal_gender(gd_2)"
COL_RESULT     = "Result"
COL_DUP_RATE   = "duplication_rate(%)"
COL_MAP_RATE   = "mapping_rate(%)"
COL_YFF2       = "YFF_2"
COL_SAMPLE_DIR = "sample_dir"
COL_MONTH      = "month"

WIG_TEMPLATE = "Output_hmmcopy/{sample_id}.of_fetus.50kb.wig.Normalization.txt"


# ---------------------------------------------------------------------------
# Sample selection
# ---------------------------------------------------------------------------

def select_samples(df: pd.DataFrame, n_samples: int,
                   max_dup: float, min_map: float, min_ff: float) -> pd.DataFrame:
    """Select male-fetus samples with stratified FF bins."""
    for col in [COL_DUP_RATE, COL_MAP_RATE, COL_YFF2]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df["month_int"] = pd.to_numeric(df[COL_MONTH], errors="coerce")

    mask = (
        (df[COL_GENDER] == "XY") &
        (df[COL_RESULT] == "Low Risk") &
        (df[COL_YFF2] >= min_ff) &
        (df[COL_DUP_RATE] < max_dup) &
        (df[COL_MAP_RATE] > min_map)
    )
    candidates = df[mask].copy()
    logger.info("QC-passing male samples: %d", len(candidates))

    # Stratified sampling across FF bins (low-FF over-represented)
    bins   = [(min_ff, 6), (6, 8), (8, 10), (10, 15), (15, 999)]
    # Proportional targets summing to n_samples
    counts = [len(candidates[(candidates[COL_YFF2] >= lo) & (candidates[COL_YFF2] < hi)])
              for lo, hi in bins]
    total  = max(sum(counts), 1)
    targets = [max(1, round(n_samples * c / total)) for c in counts]
    # Ensure exact total
    diff = n_samples - sum(targets)
    targets[-1] += diff

    selected = []
    for (lo, hi), n in zip(bins, targets):
        pool = candidates[(candidates[COL_YFF2] >= lo) & (candidates[COL_YFF2] < hi)]
        pool = pool.sort_values("month_int", ascending=False)
        take = min(n, len(pool))
        selected.append(pool.head(take))
        logger.info("  FF %4.1f–%4.1f%%: %d/%d selected", lo, hi, take, len(pool))

    result = pd.concat(selected).drop_duplicates(COL_SAMPLE_ID)
    logger.info("Total selected: %d samples", len(result))
    return result


# ---------------------------------------------------------------------------
# Training config generation
# ---------------------------------------------------------------------------

def make_training_config(selected: pd.DataFrame, wig_template: str) -> pd.DataFrame:
    """Build training config TSV rows (FILEPATH=wig, BAM_PATH=bam)."""
    rows = []
    missing_wig = []
    missing_bam = []

    for _, row in selected.iterrows():
        sid  = row[COL_SAMPLE_ID]
        sdir = str(row.get(COL_SAMPLE_DIR, ""))
        yff2 = float(row[COL_YFF2])

        wig = os.path.join(sdir, wig_template.replace("{sample_id}", sid))
        bam = os.path.join(sdir, f"{sid}.proper_paired.bam")

        if not os.path.exists(wig):
            missing_wig.append(sid)
            continue
        if not os.path.exists(bam):
            missing_bam.append(sid)
            # Still include — fragment features will be skipped gracefully
            bam = ""

        rows.append({
            "SAMPLE_ID":    sid,
            "FILEPATH":     wig,   # coverage features from HMMcopy wig
            "SEX_FETUS":    "M",
            "FF_REFERENCE": round(yff2, 4),
            "BAM_PATH":     bam,   # fragment features from BAM (patched pipeline)
        })

    if missing_wig:
        logger.warning("%d samples missing wig → excluded: %s", len(missing_wig), missing_wig[:5])
    if missing_bam:
        logger.warning("%d samples missing BAM → fragment features will be skipped: %s",
                       len(missing_bam), missing_bam[:5])

    logger.info("Training config rows: %d", len(rows))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Training
# ---------------------------------------------------------------------------

def train_model(config_path: str, out_dir: str,
                genome: str = "hg19",
                features: list = None,
                augment: bool = True,
                low_ff_weight: float = 3.0,
                threads: int = 8,
                cv_folds: int = 5) -> None:
    """Call gxff train (requires patched pipeline.py for BAM_PATH support)."""
    if features is None:
        features = ["coverage", "fragment"]

    from gxff.core.pipeline import GxFFPipeline
    from gxff.core.config import GxFFConfig

    config = GxFFConfig(
        genome=genome,
        use_coverage="coverage" in features,
        use_fragment="fragment" in features,
        use_nucleosome="nucleosome" in features,
        augment=augment,
        low_ff_weight=low_ff_weight,
        n_threads=threads,
        cv_folds=cv_folds,
    )
    pipeline = GxFFPipeline(config)
    pipeline.train(config_path=config_path, output_dir=out_dir)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample-tsv",    required=True,
                    help="Path to reference_sample_list TSV")
    ap.add_argument("--out-dir",       required=True,
                    help="Output directory for model files")
    ap.add_argument("--n-samples",     type=int, default=170,
                    help="Target training sample count (default: 170)")
    ap.add_argument("--max-dup",       type=float, default=15.0)
    ap.add_argument("--min-map",       type=float, default=85.0)
    ap.add_argument("--min-ff",        type=float, default=4.0,
                    help="Min YFF_2 %% for male samples (default: 4.0)")
    ap.add_argument("--genome",        default="hg19")
    ap.add_argument("--features",      nargs="+",
                    default=["coverage", "fragment"],
                    choices=["coverage", "fragment", "nucleosome"])
    ap.add_argument("--no-augment",    action="store_true")
    ap.add_argument("--low-ff-weight", type=float, default=3.0)
    ap.add_argument("--cv-folds",      type=int, default=5)
    ap.add_argument("--threads",       type=int, default=8)
    ap.add_argument("--wig-suffix",    default=WIG_TEMPLATE)
    ap.add_argument("--config-only",   action="store_true",
                    help="Only generate training config TSV, do not train")
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # ── Step 1: Load & filter ─────────────────────────────────────────────────
    logger.info("Loading sample TSV: %s", args.sample_tsv)
    df = pd.read_csv(args.sample_tsv, sep="\t", dtype=str)
    df.columns = [c.strip() for c in df.columns]

    selected = select_samples(df, n_samples=args.n_samples,
                              max_dup=args.max_dup, min_map=args.min_map,
                              min_ff=args.min_ff)
    if selected.empty:
        logger.error("No samples passed filters — aborting.")
        sys.exit(1)

    # ── Step 2: Generate training config ──────────────────────────────────────
    logger.info("=== Step 2: Generating training config ===")
    config_df = make_training_config(selected, wig_template=args.wig_suffix)

    if config_df.empty:
        logger.error("Training config is empty — aborting.")
        sys.exit(1)

    config_path = str(out_dir / "training_config.tsv")
    config_df.to_csv(config_path, sep="\t", index=False)
    logger.info("Training config saved: %s (%d samples)", config_path, len(config_df))

    if args.config_only:
        logger.info("--config-only flag set — skipping training.")
        return

    # ── Step 3: Train model ───────────────────────────────────────────────────
    logger.info("=== Step 3: Training gx-FF model ===")
    logger.info("Features: %s | Samples: %d | CV folds: %d | Threads: %d",
                args.features, len(config_df), args.cv_folds, args.threads)
    train_model(
        config_path=config_path,
        out_dir=str(out_dir),
        genome=args.genome,
        features=args.features,
        augment=not args.no_augment,
        low_ff_weight=args.low_ff_weight,
        threads=args.threads,
        cv_folds=args.cv_folds,
    )

    # ── Done ──────────────────────────────────────────────────────────────────
    metrics_path = out_dir / "cv_metrics.tsv"
    if metrics_path.exists():
        metrics = pd.read_csv(metrics_path, sep="\t")
        logger.info("CV metrics: %s", metrics.to_dict(orient="records"))
    logger.info("Model saved to: %s", out_dir)
    logger.info("Done.")


if __name__ == "__main__":
    main()
