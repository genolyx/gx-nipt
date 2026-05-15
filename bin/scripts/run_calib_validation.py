#!/usr/bin/env python3
"""
run_calib_validation.py
=======================
model_v11처럼 이미 gxff_model.pkl이 저장된 경우,
calibration (M-gxFF) + test set validation만 수행.

Usage:
    python run_calib_validation.py \
        --model-dir  /home/ken/gx-nipt/refs/gxff/model_v11 \
        --sample-tsv /home/ken/ken-nipt/bin/scripts/utils/reference/reference_sample_list_from_json_v3.tsv \
        --threads    120
"""
import argparse
import json
import logging
import os
import sys

import numpy as np
import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [calib] %(levelname)s %(message)s",
    stream=sys.stderr,
    force=True,
)
log = logging.getLogger()

sys.path.insert(0, "/home/ken/gx-nipt/bin/scripts")
from build_gxff_model import (
    calibrate_mgxff,
    validate_model,
    COL_SAMPLE_ID,
    COL_GENDER,
    COL_YFF2,
    COL_MSEQFF,
    WIG_TEMPLATE,
)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model-dir",  required=True)
    ap.add_argument("--sample-tsv", required=True)
    ap.add_argument("--threads",    type=int, default=120)
    ap.add_argument("--genome",     default="hg19")
    ap.add_argument("--features",   nargs="+", default=["coverage", "fragment"])
    args = ap.parse_args()

    model_path  = os.path.join(args.model_dir, "gxff_model.pkl")
    test_json   = os.path.join(args.model_dir, "test_set.json")
    config_path = os.path.join(args.model_dir, "training_config.tsv")
    calib_path  = os.path.join(args.model_dir, "mgxff_calibration.json")

    if not os.path.exists(model_path):
        log.error("Model not found: %s", model_path)
        sys.exit(1)
    if not os.path.exists(test_json):
        log.error("test_set.json not found: %s", test_json)
        sys.exit(1)

    # Reconstruct train_df from training_config.tsv + sample_tsv
    log.info("Loading training_config.tsv: %s", config_path)
    cfg = pd.read_csv(config_path, sep="\t")
    train_ids = set(cfg["SAMPLE_ID"].tolist())

    log.info("Loading sample TSV: %s", args.sample_tsv)
    ref = pd.read_csv(args.sample_tsv, sep="\t", dtype=str)
    ref.columns = [c.strip() for c in ref.columns]
    ref[COL_YFF2]   = pd.to_numeric(ref[COL_YFF2],   errors="coerce")
    ref[COL_MSEQFF] = pd.to_numeric(ref[COL_MSEQFF], errors="coerce")

    merged = ref[ref[COL_SAMPLE_ID].isin(train_ids)].copy()
    merged["_sex"] = merged[COL_GENDER].map({"XY": "M", "XX": "F"})
    merged["_ff_ref"] = np.where(
        merged["_sex"] == "M", merged[COL_YFF2], merged[COL_MSEQFF]
    )
    log.info(
        "Train samples reconstructed: %d (male=%d, female=%d)",
        len(merged),
        (merged["_sex"] == "M").sum(),
        (merged["_sex"] == "F").sum(),
    )

    # Step 5: M-gxFF calibration (parallel)
    log.info("=== Step 5: M-gxFF calibration (workers=%d) ===", args.threads)
    slope, intercept = calibrate_mgxff(
        model_path, merged, WIG_TEMPLATE,
        args.genome, args.features, n_workers=args.threads,
    )
    with open(calib_path, "w") as f:
        json.dump({"slope": slope, "intercept": intercept}, f, indent=2)
    log.info("Calibration saved: slope=%.4f, intercept=%.4f", slope, intercept)

    # Step 6: Test set validation (parallel)
    log.info("=== Step 6: Test set validation (workers=%d) ===", args.threads)
    validate_model(
        model_path, test_json, args.model_dir,
        slope, intercept, WIG_TEMPLATE,
        args.genome, args.features, n_workers=args.threads,
    )
    log.info("=== Done ===")


if __name__ == "__main__":
    main()
