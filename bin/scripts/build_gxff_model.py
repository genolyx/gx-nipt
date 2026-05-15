#!/usr/bin/env python3
"""
build_gxff_model.py
===================
End-to-end gx-FF model builder for hg19 NIPT samples.

Pipeline steps (all in one run):
  1. Load reference_sample_list TSV.
  2. Filter & select samples:
       - Male  (XY): YFF_2 as ground-truth FF
       - Female (XX): M-SeqFF as ground-truth FF
       Stratified FF bins, with forced minimum for Low FF bucket.
  3. Split into train / test sets; save test set JSON.
  4. Generate training config TSV (FILEPATH=of_orig wig, BAM_PATH=BAM).
  5. Train gx-FF ensemble via patched pipeline (parallel feature extraction).
  6. Post-training calibration: fit M-gxFF linear correction on train set.
  7. Validate on held-out test set; save metrics TSV.

Usage
-----
    python build_gxff_model.py \\
        --sample-tsv  /path/to/reference_sample_list_from_json_v3.tsv \\
        --out-dir     /home/ken/gx-nipt/refs/gxff/model_v11 \\
        --n-samples   900 \\
        --threads     32

    # config-only (no training):
    python build_gxff_model.py --sample-tsv ... --out-dir ... --config-only

    # validation only on existing model:
    python build_gxff_model.py --sample-tsv ... --out-dir ... --validate-only

Prerequisites
-------------
  • Docker image gx-nipt:latest (or direct Python env with gxff installed).
  • patched pipeline.py must be in place inside the Python env.
    Copy: <repo>/bin/gxff_patch/pipeline.py
      → /opt/conda/envs/nipt/lib/python3.11/site-packages/gxff/core/pipeline.py
"""

from __future__ import annotations

import argparse
import json
import logging
import multiprocessing as _mp
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy import stats

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [build_gxff] %(levelname)s %(message)s",
    stream=sys.stderr,
    force=True,
)
logger = logging.getLogger(__name__)

# ── Column constants ──────────────────────────────────────────────────────────
COL_SAMPLE_ID  = "sample_id"
COL_GENDER     = "fetal_gender(gd_2)"
COL_RESULT     = "Result"
COL_DUP_RATE   = "duplication_rate(%)"
COL_MAP_RATE   = "mapping_rate(%)"
COL_YFF2       = "YFF_2"
COL_MSEQFF     = "M-SeqFF"
COL_SEQFF      = "SeqFF"
COL_SAMPLE_DIR = "sample_dir"
COL_MONTH      = "month"

# WIG suffix – must match pipeline prediction (of_orig, not of_fetus)
WIG_TEMPLATE = "Output_hmmcopy/{sample_id}.of_orig.50kb.wig.Normalization.txt"


# ─────────────────────────────────────────────────────────────────────────────
# 1. Sample selection
# ─────────────────────────────────────────────────────────────────────────────

def _stratified_select(
    candidates: pd.DataFrame,
    ff_col: str,
    n_total: int,
    bins: List[Tuple[float, float]],
    low_ff_min_count: int,
) -> pd.DataFrame:
    """
    Proportional stratified sampling across FF bins.
    The first bin (lowest FF) is guaranteed at least `low_ff_min_count` samples.
    """
    counts = [
        len(candidates[(candidates[ff_col] >= lo) & (candidates[ff_col] < hi)])
        for lo, hi in bins
    ]
    total = max(sum(counts), 1)
    targets = [max(1, round(n_total * c / total)) for c in counts]

    # Guarantee minimum low-FF count
    if targets[0] < low_ff_min_count and counts[0] >= low_ff_min_count:
        diff = low_ff_min_count - targets[0]
        targets[0] = low_ff_min_count
        # Take proportionally from higher bins
        for i in range(len(targets) - 1, 0, -1):
            take = min(diff, max(0, targets[i] - 1))
            targets[i] -= take
            diff -= take
            if diff <= 0:
                break

    # Exact total adjustment
    diff = n_total - sum(targets)
    targets[-1] = max(1, targets[-1] + diff)

    selected = []
    for (lo, hi), n in zip(bins, targets):
        pool = candidates[(candidates[ff_col] >= lo) & (candidates[ff_col] < hi)].copy()
        pool = pool.sort_values(COL_MONTH, ascending=False)
        take = min(n, len(pool))
        selected.append(pool.head(take))
        logger.info("  FF %4.1f–%4.1f%%: %d/%d selected", lo, hi, take, len(pool))

    return pd.concat(selected).drop_duplicates(COL_SAMPLE_ID)


def select_samples(
    df: pd.DataFrame,
    n_samples: int,
    male_fraction: float,
    max_dup: float,
    min_map: float,
    min_ff: float,
    low_ff_min_count: int,
) -> pd.DataFrame:
    """
    Select male + female samples with stratified FF bins.

    Male  ground truth : YFF_2
    Female ground truth: M-SeqFF
    """
    for col in [COL_DUP_RATE, COL_MAP_RATE, COL_YFF2, COL_MSEQFF, COL_MONTH]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    qc_mask = (
        (df[COL_RESULT] == "Low Risk") &
        (df[COL_DUP_RATE] < max_dup) &
        (df[COL_MAP_RATE] > min_map)
    )

    # ── Male ────────────────────────────────────────────────────────────────
    n_male = round(n_samples * male_fraction)
    male_mask = qc_mask & (df[COL_GENDER] == "XY") & (df[COL_YFF2] >= min_ff)
    male_cand = df[male_mask].copy()
    logger.info("QC-passing male candidates: %d", len(male_cand))

    bins_male = [(min_ff, 6), (6, 8), (8, 10), (10, 15), (15, 999)]
    logger.info("=== Male sample selection ===")
    male_sel = _stratified_select(
        male_cand, COL_YFF2, n_male, bins_male, low_ff_min_count
    )
    male_sel = male_sel.copy()
    male_sel["_ff_ref"] = male_sel[COL_YFF2]
    male_sel["_sex"]    = "M"
    logger.info("Male selected: %d", len(male_sel))

    # ── Female ───────────────────────────────────────────────────────────────
    n_female = n_samples - n_male
    female_mask = (
        qc_mask &
        (df[COL_GENDER] == "XX") &
        (df[COL_MSEQFF] >= min_ff) &
        df[COL_MSEQFF].notna()
    )
    female_cand = df[female_mask].copy()
    logger.info("QC-passing female candidates: %d", len(female_cand))

    bins_female = [(min_ff, 6), (6, 8), (8, 10), (10, 15), (15, 999)]
    logger.info("=== Female sample selection ===")
    female_sel = _stratified_select(
        female_cand, COL_MSEQFF, n_female, bins_female, low_ff_min_count
    )
    female_sel = female_sel.copy()
    female_sel["_ff_ref"] = female_sel[COL_MSEQFF]
    female_sel["_sex"]    = "F"
    logger.info("Female selected: %d", len(female_sel))

    result = pd.concat([male_sel, female_sel]).drop_duplicates(COL_SAMPLE_ID)
    logger.info("Total selected: %d (male=%d, female=%d)",
                len(result), (result["_sex"] == "M").sum(), (result["_sex"] == "F").sum())
    return result


# ─────────────────────────────────────────────────────────────────────────────
# 2. Train / test split
# ─────────────────────────────────────────────────────────────────────────────

def train_test_split(
    selected: pd.DataFrame,
    test_fraction: float,
    random_seed: int = 42,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Stratified train/test split preserving gender and FF-bin distribution.
    Saves test set as JSON for validate_gxff.py.
    """
    rng = np.random.default_rng(random_seed)
    train_rows, test_rows = [], []

    for sex, grp in selected.groupby("_sex"):
        ff_col = COL_YFF2 if sex == "M" else COL_MSEQFF
        grp = grp.copy()
        # FF bins for stratification
        bins_e = [0, 6, 8, 10, 15, 999]
        grp["_bin"] = pd.cut(grp[ff_col], bins=bins_e, right=True)

        for _, bin_grp in grp.groupby("_bin", observed=True):
            if len(bin_grp) == 0:
                continue
            n_test = max(1, round(len(bin_grp) * test_fraction))
            idx = rng.choice(len(bin_grp), size=n_test, replace=False)
            test_rows.append(bin_grp.iloc[idx])
            train_rows.append(bin_grp.drop(bin_grp.index[idx]))

    train_df = pd.concat(train_rows).drop_duplicates(COL_SAMPLE_ID)
    test_df  = pd.concat(test_rows).drop_duplicates(COL_SAMPLE_ID)
    logger.info("Train: %d | Test: %d", len(train_df), len(test_df))
    return train_df, test_df


def save_test_set_json(test_df: pd.DataFrame, wig_template: str, out_path: str) -> None:
    records = []
    for _, row in test_df.iterrows():
        sid  = row[COL_SAMPLE_ID]
        sdir = str(row[COL_SAMPLE_DIR])
        wig  = os.path.join(sdir, wig_template.replace("{sample_id}", sid))
        ff   = float(row["_ff_ref"])
        sex  = str(row["_sex"])
        records.append({"id": sid, "wig": wig, "ff": ff, "gender": sex})
    with open(out_path, "w") as f:
        json.dump(records, f, indent=2)
    logger.info("Test set JSON saved: %s (%d samples)", out_path, len(records))


# ─────────────────────────────────────────────────────────────────────────────
# 3. Training config
# ─────────────────────────────────────────────────────────────────────────────

def make_training_config(
    selected: pd.DataFrame,
    wig_template: str,
) -> pd.DataFrame:
    """Build training config TSV rows (FILEPATH=wig, BAM_PATH=bam)."""
    rows = []
    missing_wig, missing_bam = [], []

    for _, row in selected.iterrows():
        sid  = row[COL_SAMPLE_ID]
        sdir = str(row[COL_SAMPLE_DIR])
        sex  = str(row["_sex"])
        ff   = float(row["_ff_ref"])

        wig = os.path.join(sdir, wig_template.replace("{sample_id}", sid))
        bam = os.path.join(sdir, f"{sid}.proper_paired.bam")

        if not os.path.exists(wig):
            missing_wig.append(sid)
            continue

        bam_val = bam if os.path.exists(bam) else ""
        if not bam_val:
            missing_bam.append(sid)

        rows.append({
            "SAMPLE_ID":    sid,
            "FILEPATH":     wig,                   # coverage features from HMMcopy wig
            "SEX_FETUS":    sex,
            "FF_REFERENCE": round(ff / 100.0, 6),  # convert % → fraction [0,1]
            "BAM_PATH":     bam_val,               # fragment features from BAM
        })

    if missing_wig:
        logger.warning("%d samples missing WIG → excluded: %s",
                       len(missing_wig), missing_wig[:5])
    if missing_bam:
        logger.warning("%d samples missing BAM → fragment features skipped: %s",
                       len(missing_bam), missing_bam[:5])

    logger.info("Training config rows: %d (male=%d, female=%d, low_ff_lt6=%d)",
                len(rows),
                sum(1 for r in rows if r["SEX_FETUS"] == "M"),
                sum(1 for r in rows if r["SEX_FETUS"] == "F"),
                sum(1 for r in rows if r["FF_REFERENCE"] < 0.06))
    return pd.DataFrame(rows)


# ─────────────────────────────────────────────────────────────────────────────
# 4. Training
# ─────────────────────────────────────────────────────────────────────────────

def train_model(
    config_path: str,
    out_dir: str,
    genome: str = "hg19",
    features: Optional[List[str]] = None,
    augment: bool = True,
    low_ff_weight: float = 5.0,
    threads: int = 8,
    cv_folds: int = 5,
) -> None:
    """Call gxff train (requires patched pipeline.py for BAM_PATH support)."""
    if features is None:
        features = ["coverage", "fragment"]

    from gxff.core.pipeline import GxFFPipeline
    from gxff.core.config import GxFFConfig

    config = GxFFConfig(
        genome=genome,
        features=features,
        augment=augment,
        low_ff_weight=low_ff_weight,
        threads=threads,
        cv_folds=cv_folds,
    )
    pipeline = GxFFPipeline(config)
    pipeline.train(config_path=config_path, output_dir=out_dir)


# ─────────────────────────────────────────────────────────────────────────────
# Module-level worker: predict a single sample (picklable for ProcessPoolExecutor)
# ─────────────────────────────────────────────────────────────────────────────

def _predict_worker(args: Tuple) -> Tuple:
    """
    Returns (sample_id, gxff_pct, truth_pct, error_msg)
    Loads model fresh per process to avoid shared-state issues.
    """
    model_path, genome, features, sid, wig, bam, truth = args
    try:
        from gxff.core.pipeline import GxFFPipeline
        from gxff.core.config import GxFFConfig
        import numpy as np
        config   = GxFFConfig(genome=genome, features=features)
        pipeline = GxFFPipeline(config=config, model_path=model_path)
        feat = pipeline._extract_features(
            bincount_path=wig,
            bam_path=bam if (bam and os.path.exists(bam)) else None,
        )
        ff_arr, lgbm_arr, dnn_arr = pipeline.model.predict(feat.reshape(1, -1))
        gxff_pct  = float(np.squeeze(ff_arr))   * 100.0
        lgbm_pct  = float(np.squeeze(lgbm_arr)) * 100.0
        dnn_pct   = float(np.squeeze(dnn_arr))  * 100.0
        return (sid, gxff_pct, lgbm_pct, dnn_pct, truth, None)
    except Exception as exc:
        return (sid, None, None, None, truth, str(exc))


def _parallel_predict(
    model_path: str,
    samples: List[Tuple],   # list of (sid, wig, bam, truth)
    genome: str,
    features: List[str],
    n_workers: int,
    label: str = "Predicting",
) -> List[Tuple]:
    """Predict a list of samples in parallel. Returns list of worker results."""
    worker_args = [
        (model_path, genome, features, sid, wig, bam, truth)
        for sid, wig, bam, truth in samples
    ]
    results = []
    completed = failed = 0
    total = len(worker_args)
    mp_ctx = _mp.get_context("forkserver")
    with ProcessPoolExecutor(max_workers=n_workers, mp_context=mp_ctx) as ex:
        futures = {ex.submit(_predict_worker, a): a[3] for a in worker_args}
        for fut in as_completed(futures):
            completed += 1
            res = fut.result()
            if res[1] is None:
                failed += 1
                logger.warning("%s failed %s: %s", label, res[0], res[5])
            results.append(res)
            if completed % 50 == 0 or completed == total:
                logger.info("%s progress: %d/%d done, %d failed",
                            label, completed, total, failed)
    return results


# ─────────────────────────────────────────────────────────────────────────────
# 5. M-gxFF calibration (linear correction on training set)
# ─────────────────────────────────────────────────────────────────────────────

def calibrate_mgxff(
    model_path: str,
    train_df: pd.DataFrame,
    wig_template: str,
    genome: str = "hg19",
    features: Optional[List[str]] = None,
    n_workers: int = 32,
) -> Tuple[float, float]:
    """
    Fit linear calibration: M-gxFF = slope * gxFF + intercept
    using male training samples (YFF_2 as ground truth). Parallel prediction.
    """
    if features is None:
        features = ["coverage", "fragment"]

    male_train = train_df[train_df["_sex"] == "M"].copy()
    logger.info("Calibration: predicting %d male training samples (workers=%d)...",
                len(male_train), n_workers)

    samples = []
    for _, row in male_train.iterrows():
        sid  = row[COL_SAMPLE_ID]
        sdir = str(row[COL_SAMPLE_DIR])
        wig  = os.path.join(sdir, wig_template.replace("{sample_id}", sid))
        bam  = os.path.join(sdir, f"{sid}.proper_paired.bam")
        ff   = float(row["_ff_ref"])
        if os.path.exists(wig):
            samples.append((sid, wig, bam, ff))

    results = _parallel_predict(model_path, samples, genome, features, n_workers,
                                label="Calibration")

    gxff_preds, yff2_truths = [], []
    for sid, gxff_pct, _, _, truth, err in results:
        if gxff_pct is not None:
            gxff_preds.append(gxff_pct)
            yff2_truths.append(truth)

    if len(gxff_preds) < 10:
        logger.warning("Too few calibration samples (%d) — skipping M-gxFF fit",
                       len(gxff_preds))
        return 1.0, 0.0

    slope, intercept, r, _, _ = stats.linregress(gxff_preds, yff2_truths)
    logger.info(
        "M-gxFF calibration (N=%d): M-gxFF = %.4f * gxFF + %.4f  (r=%.4f)",
        len(gxff_preds), slope, intercept, r,
    )
    return float(slope), float(intercept)


# ─────────────────────────────────────────────────────────────────────────────
# 6. Validation
# ─────────────────────────────────────────────────────────────────────────────

def validate_model(
    model_path: str,
    test_json: str,
    out_dir: str,
    slope: float,
    intercept: float,
    wig_template: str,
    genome: str = "hg19",
    features: Optional[List[str]] = None,
    n_workers: int = 32,
) -> pd.DataFrame:
    """
    Predict on held-out test set (parallel), compute metrics for raw gxFF and M-gxFF.
    Saves validation_per_sample.tsv and validation_summary.tsv.
    """
    if features is None:
        features = ["coverage", "fragment"]

    with open(test_json) as f:
        test_samples = json.load(f)

    logger.info("Validating on %d test samples (workers=%d)...",
                len(test_samples), n_workers)

    samples = []
    gender_map = {}
    for s in test_samples:
        sid   = s["id"]
        wig   = s["wig"]
        truth = s["ff"]
        sex   = s.get("gender", "?")
        bam   = str(Path(wig).parent.parent / f"{sid}.proper_paired.bam")
        samples.append((sid, wig, bam, truth))
        gender_map[sid] = sex

    results = _parallel_predict(model_path, samples, genome, features, n_workers,
                                label="Validation")

    rows = []
    for sid, gxff_pct, lgbm_pct, dnn_pct, truth, err in results:
        sex = gender_map.get(sid, "?")
        if gxff_pct is not None:
            mgxff_pct = slope * gxff_pct + intercept
            rows.append({
                "sample_id": sid, "gender": sex, "ff_true": truth,
                "gxff":  round(gxff_pct, 3),
                "mgxff": round(mgxff_pct, 3),
                "lgbm":  round(lgbm_pct, 3),
                "dnn":   round(dnn_pct, 3),
                "error_gxff":  round(gxff_pct - truth, 3),
                "error_mgxff": round(mgxff_pct - truth, 3),
            })
        else:
            rows.append({
                "sample_id": sid, "gender": sex, "ff_true": truth,
                "gxff": float("nan"), "mgxff": float("nan"),
                "lgbm": float("nan"), "dnn": float("nan"),
                "error_gxff": float("nan"), "error_mgxff": float("nan"),
            })

    df = pd.DataFrame(rows)
    per_sample_path = os.path.join(out_dir, "validation_per_sample.tsv")
    df.to_csv(per_sample_path, sep="\t", index=False)

    # Summary metrics
    def met(true, pred, label):
        t, p = np.array(true), np.array(pred)
        ok = ~np.isnan(p)
        t, p = t[ok], p[ok]
        r   = np.corrcoef(t, p)[0, 1] if len(t) > 1 else float("nan")
        mae  = np.mean(np.abs(t - p))
        rmse = np.sqrt(np.mean((t - p) ** 2))
        return {"algorithm": label, "n": len(t),
                "pearson_r": round(r, 4), "mae_pct": round(mae, 3), "rmse_pct": round(rmse, 3)}

    summary_rows = []
    for grp, mask in [("All", df.gender.notna()),
                      ("Male", df.gender == "M"),
                      ("Female", df.gender == "F")]:
        sub = df[mask]
        if len(sub) == 0:
            continue
        for alg, col in [("gxFF", "gxff"), ("M-gxFF", "mgxff")]:
            m = met(sub.ff_true, sub[col], alg)
            m["group"] = grp
            summary_rows.append(m)

    summary = pd.DataFrame(summary_rows)[["algorithm","group","n","pearson_r","mae_pct","rmse_pct"]]
    summary_path = os.path.join(out_dir, "validation_summary.tsv")
    summary.to_csv(summary_path, sep="\t", index=False)

    print("\n" + "=" * 70)
    print("VALIDATION SUMMARY")
    print("=" * 70)
    print(summary.to_string(index=False))
    print("=" * 70)

    logger.info("Validation saved: %s", out_dir)
    return summary


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--sample-tsv",       required=True,
                    help="reference_sample_list TSV")
    ap.add_argument("--out-dir",          required=True,
                    help="Output directory for model + logs")
    ap.add_argument("--n-samples",        type=int,   default=900,
                    help="Target total training samples (default: 900)")
    ap.add_argument("--male-fraction",    type=float, default=0.5,
                    help="Fraction of male samples (default: 0.5)")
    ap.add_argument("--max-dup",          type=float, default=15.0)
    ap.add_argument("--min-map",          type=float, default=85.0)
    ap.add_argument("--min-ff",           type=float, default=3.0,
                    help="Min FF %% for inclusion (default: 3.0, includes Low FF)")
    ap.add_argument("--low-ff-min-count", type=int,   default=30,
                    help="Min samples in lowest FF bin (<6%%) per gender (default: 30)")
    ap.add_argument("--test-fraction",    type=float, default=0.2,
                    help="Fraction held out for test set (default: 0.2)")
    ap.add_argument("--genome",           default="hg19")
    ap.add_argument("--features",         nargs="+",
                    default=["coverage", "fragment"],
                    choices=["coverage", "fragment", "nucleosome"])
    ap.add_argument("--no-augment",       action="store_true")
    ap.add_argument("--low-ff-weight",    type=float, default=5.0,
                    help="Sample weight for low-FF (<5%%) samples (default: 5.0)")
    ap.add_argument("--cv-folds",         type=int,   default=5)
    ap.add_argument("--threads",          type=int,   default=32)
    ap.add_argument("--wig-suffix",       default=WIG_TEMPLATE)
    ap.add_argument("--config-only",      action="store_true",
                    help="Generate training config only, skip training")
    ap.add_argument("--skip-calibration", action="store_true",
                    help="Skip M-gxFF calibration step (fast re-validation)")
    ap.add_argument("--validate-only",    action="store_true",
                    help="Only run validation on existing model+test JSON")
    ap.add_argument("--random-seed",      type=int,   default=42)
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    model_path    = str(out_dir / "gxff_model.pkl")
    config_path   = str(out_dir / "training_config.tsv")
    test_json     = str(out_dir / "test_set.json")
    calib_path    = str(out_dir / "mgxff_calibration.json")

    # ── validate-only mode ───────────────────────────────────────────────────
    if args.validate_only:
        if not os.path.exists(model_path):
            logger.error("Model not found: %s", model_path)
            sys.exit(1)
        if not os.path.exists(test_json):
            logger.error("Test set JSON not found: %s", test_json)
            sys.exit(1)
        slope, intercept = 1.0, 0.0
        if os.path.exists(calib_path):
            with open(calib_path) as f:
                calib = json.load(f)
            slope, intercept = calib["slope"], calib["intercept"]
            logger.info("Loaded calibration: slope=%.4f, intercept=%.4f", slope, intercept)
        validate_model(model_path, test_json, str(out_dir),
                       slope, intercept, args.wig_suffix, args.genome, args.features,
                       n_workers=args.threads)
        return

    # ── Step 1: Load & select ────────────────────────────────────────────────
    logger.info("Loading sample TSV: %s", args.sample_tsv)
    df = pd.read_csv(args.sample_tsv, sep="\t", dtype=str)
    df.columns = [c.strip() for c in df.columns]

    selected = select_samples(
        df,
        n_samples=args.n_samples,
        male_fraction=args.male_fraction,
        max_dup=args.max_dup,
        min_map=args.min_map,
        min_ff=args.min_ff,
        low_ff_min_count=args.low_ff_min_count,
    )
    if selected.empty:
        logger.error("No samples passed filters — aborting.")
        sys.exit(1)

    # ── Step 2: Train / test split ───────────────────────────────────────────
    logger.info("=== Step 2: Train/test split (test_fraction=%.0f%%) ===",
                args.test_fraction * 100)
    train_df, test_df = train_test_split(selected, args.test_fraction, args.random_seed)
    save_test_set_json(test_df, args.wig_suffix, test_json)

    # ── Step 3: Training config ───────────────────────────────────────────────
    logger.info("=== Step 3: Generating training config ===")
    config_df = make_training_config(train_df, args.wig_suffix)
    if config_df.empty:
        logger.error("Training config is empty — aborting.")
        sys.exit(1)
    config_df.to_csv(config_path, sep="\t", index=False)
    logger.info("Training config saved: %s (%d rows)", config_path, len(config_df))

    if args.config_only:
        logger.info("--config-only: stopping before training.")
        return

    # ── Step 4: Train ─────────────────────────────────────────────────────────
    logger.info("=== Step 4: Training gx-FF model ===")
    logger.info("Features: %s | Train N: %d | low_ff_weight: %.1f | threads: %d",
                args.features, len(config_df), args.low_ff_weight, args.threads)
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

    if not os.path.exists(model_path):
        logger.error("Training failed — model not found at %s", model_path)
        sys.exit(1)

    # ── Step 5: M-gxFF calibration ────────────────────────────────────────────
    slope, intercept = 1.0, 0.0
    if not args.skip_calibration:
        logger.info("=== Step 5: M-gxFF calibration ===")
        slope, intercept = calibrate_mgxff(
            model_path, train_df, args.wig_suffix, args.genome, args.features,
            n_workers=args.threads,
        )
        with open(calib_path, "w") as f:
            json.dump({"slope": slope, "intercept": intercept}, f, indent=2)
        logger.info("Calibration saved: %s", calib_path)

    # ── Step 6: Validation ────────────────────────────────────────────────────
    logger.info("=== Step 6: Test set validation ===")
    validate_model(
        model_path, test_json, str(out_dir),
        slope, intercept, args.wig_suffix, args.genome, args.features,
        n_workers=args.threads,
    )

    logger.info("=== All done. Model: %s ===", model_path)


if __name__ == "__main__":
    main()
