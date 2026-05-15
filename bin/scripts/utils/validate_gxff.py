#!/usr/bin/env python3
"""
validate_gxff.py
================
test_set_v6.json 기반으로 model_v6 / model_v10 성능을 검증합니다.

사용법 (Docker 내부):
  python validate_gxff.py \
      --test-json  /home/ken/gx-nipt/refs/gxff/test_set_v6.json \
      --models     v6=/home/ken/gx-nipt/refs/gxff/model_v6/gxff_model.pkl \
                   v10=/home/ken/gx-nipt/refs/gxff/model_v10/gxff_model.pkl \
      --out        /home/ken/gx-nipt/refs/gxff/validation_results.tsv
"""

import argparse
import json
import logging
import math
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [validate] %(levelname)s %(message)s",
    stream=sys.stderr,
    force=True,
)
logger = logging.getLogger(__name__)


def pearson_r(a, b):
    a, b = np.array(a, dtype=float), np.array(b, dtype=float)
    if len(a) < 2:
        return float("nan")
    return float(np.corrcoef(a, b)[0, 1])


def metrics(true_pct, pred_pct):
    """Compute MAE, RMSE, Pearson r (all in % scale)."""
    t = np.array(true_pct, dtype=float)
    p = np.array(pred_pct, dtype=float)
    mae  = float(np.mean(np.abs(t - p)))
    rmse = float(np.sqrt(np.mean((t - p) ** 2)))
    r    = pearson_r(t, p)
    return {"n": len(t), "pearson_r": r, "mae": mae, "rmse": rmse}


def derive_bam(wig_path: str, sample_id: str) -> str:
    """WIG path → proper_paired BAM path."""
    sample_dir = str(Path(wig_path).parent.parent)
    return os.path.join(sample_dir, f"{sample_id}.proper_paired.bam")


def predict_samples(model_path, test_samples, features, genome="hg19"):
    """
    Batch-predict all test samples using the given model.

    Parameters
    ----------
    features : list[str]  e.g. ["coverage"] or ["coverage","fragment"]
    """
    import pickle
    from gxff.core.config import GxFFConfig
    from gxff.core.pipeline import GxFFPipeline

    config = GxFFConfig(genome=genome, features=features)
    pipeline = GxFFPipeline(config=config, model_path=model_path)

    results = []
    for i, s in enumerate(test_samples, 1):
        sid   = s["id"]
        wig   = s["wig"]
        bam   = derive_bam(wig, sid) if "fragment" in features else None
        truth = s["ff"]        # percentage scale
        gender = s.get("gender", "?")

        try:
            # _extract_features supports both bincount + bam simultaneously
            feat = pipeline._extract_features(
                bincount_path=wig,
                bam_path=bam if (bam and os.path.exists(bam)) else None,
            )
            ff_arr, lgbm_arr, dnn_arr = pipeline.model.predict(feat.reshape(1, -1))
            ff_frac  = float(np.squeeze(ff_arr))
            ff_pct   = ff_frac * 100.0
            lgbm_pct = float(np.squeeze(lgbm_arr)) * 100.0
            dnn_pct  = float(np.squeeze(dnn_arr))  * 100.0
            err = ff_pct - truth
            results.append({
                "sample_id": sid,
                "gender": gender,
                "ff_true": truth,
                "ff_pred": round(ff_pct, 3),
                "ff_lgbm": round(lgbm_pct, 3),
                "ff_dnn":  round(dnn_pct, 3),
                "error":   round(err, 3),
                "abs_error": round(abs(err), 3),
            })
        except Exception as exc:
            logger.warning("Failed %s: %s", sid, exc)
            results.append({
                "sample_id": sid,
                "gender": gender,
                "ff_true": truth,
                "ff_pred": float("nan"),
                "ff_lgbm": float("nan"),
                "ff_dnn":  float("nan"),
                "error":   float("nan"),
                "abs_error": float("nan"),
            })

        if i % 20 == 0 or i == len(test_samples):
            logger.info("Predicted %d/%d", i, len(test_samples))

    return pd.DataFrame(results)


def guess_features(name, path):
    """Detect which features the model was trained on from its internal scaler."""
    try:
        import pickle
        with open(path, "rb") as f:
            model = pickle.load(f)
        # GxFFEnsemble stores the scaler inside the LightGBM sub-model
        n_feat = model._lgbm._scaler.n_features_in_
        logger.info("%s: lgbm scaler n_features_in_=%d", name, n_feat)
        # coverage-only: 49346 features; coverage+fragment: 49801 features
        return ["coverage", "fragment"] if n_feat > 49400 else ["coverage"]
    except Exception as exc:
        logger.warning("%s: could not detect features (%s) — defaulting to coverage-only", name, exc)
        return ["coverage"]


def print_metrics_table(label, df):
    groups = {"All": df, "Male (M)": df[df.gender=="M"], "Female (F)": df[df.gender=="F"]}
    rows = []
    for gname, gdf in groups.items():
        valid = gdf.dropna(subset=["ff_pred"])
        if len(valid) == 0:
            continue
        m = metrics(valid.ff_true, valid.ff_pred)
        rows.append({
            "Model": label, "Group": gname,
            "N": m["n"],
            "Pearson_r": round(m["pearson_r"], 4),
            "MAE_%": round(m["mae"], 3),
            "RMSE_%": round(m["rmse"], 3),
        })
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--test-json", required=True)
    ap.add_argument("--models", nargs="+", required=True,
                    help="name=path pairs, e.g. v6=/path/to/model_v6/gxff_model.pkl")
    ap.add_argument("--out",   default="/home/ken/gx-nipt/refs/gxff/validation_results.tsv")
    ap.add_argument("--genome", default="hg19")
    args = ap.parse_args()

    # Parse model specs
    model_specs = {}
    for spec in args.models:
        name, path = spec.split("=", 1)
        model_specs[name] = path
    logger.info("model_specs parsed: %s", list(model_specs.keys()))

    # Load test set
    with open(args.test_json) as f:
        test_samples = json.load(f)
    logger.info("Test samples: %d", len(test_samples))

    all_metrics = []
    all_predictions = []

    for model_name, model_path in model_specs.items():
        if not os.path.exists(model_path):
            logger.error("Model not found: %s", model_path)
            continue

        features = guess_features(model_name, model_path)
        logger.info("=== Predicting with %s (features: %s) ===", model_name, features)

        df = predict_samples(model_path, test_samples, features, args.genome)
        df.insert(0, "model", model_name)
        all_predictions.append(df)

        mt = print_metrics_table(model_name, df)
        all_metrics.append(mt)

    # Save per-sample predictions
    pred_path = args.out.replace(".tsv", "_per_sample.tsv")
    pd.concat(all_predictions).to_csv(pred_path, sep="\t", index=False)
    logger.info("Per-sample predictions saved: %s", pred_path)

    # Print & save summary metrics
    summary = pd.concat(all_metrics)
    print("\n" + "="*70)
    print("VALIDATION SUMMARY")
    print("="*70)
    print(summary.to_string(index=False))
    print("="*70)

    summary.to_csv(args.out, sep="\t", index=False)
    logger.info("Summary saved: %s", args.out)


if __name__ == "__main__":
    main()
