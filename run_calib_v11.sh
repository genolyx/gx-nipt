#!/usr/bin/env bash
# Calibration + validation only for model_v11 (model already trained)
GXNIPT_DIR="/home/ken/gx-nipt"
KENNIPT_DIR="/home/ken/ken-nipt"
PATCH_SRC="${GXNIPT_DIR}/bin/gxff_patch/pipeline.py"
PATCH_DEST="/opt/conda/envs/nipt/lib/python3.11/site-packages/gxff/core/pipeline.py"
MODEL_DIR="${GXNIPT_DIR}/refs/gxff/model_v11"
LOG="${MODEL_DIR}/calib_validation.log"

nohup docker run --rm \
    --cpus=120 --memory=200g \
    -v "${GXNIPT_DIR}:${GXNIPT_DIR}" \
    -v "${KENNIPT_DIR}:${KENNIPT_DIR}" \
    -v "${PATCH_SRC}:${PATCH_DEST}:ro" \
    -e HOME=/root \
    -e OMP_NUM_THREADS=120 \
    -e OPENBLAS_NUM_THREADS=1 \
    -e MKL_NUM_THREADS=120 \
    gx-nipt:latest \
    /opt/conda/envs/nipt/bin/python - \
        "${MODEL_DIR}" \
        "${KENNIPT_DIR}/bin/scripts/utils/reference/reference_sample_list_from_json_v3.tsv" \
        "120" \
<<'PYEOF'
import sys, json, os, logging
import numpy as np, pandas as pd
from scipy import stats
from pathlib import Path

logging.basicConfig(level=logging.INFO,
    format="%(asctime)s [calib] %(levelname)s %(message)s",
    stream=sys.stderr, force=True)
log = logging.getLogger()

model_dir  = sys.argv[1]
sample_tsv = sys.argv[2]
n_workers  = int(sys.argv[3])

sys.path.insert(0, '/home/ken/gx-nipt/bin/scripts')
from build_gxff_model import (
    calibrate_mgxff, validate_model,
    COL_SAMPLE_ID, COL_GENDER, COL_SAMPLE_DIR, COL_YFF2, COL_MSEQFF,
    WIG_TEMPLATE
)

model_path = os.path.join(model_dir, 'gxff_model.pkl')
test_json  = os.path.join(model_dir, 'test_set.json')
calib_path = os.path.join(model_dir, 'mgxff_calibration.json')
config_path = os.path.join(model_dir, 'training_config.tsv')

# Reconstruct train_df from training_config.tsv + sample_tsv
log.info("Loading training_config.tsv ...")
cfg = pd.read_csv(config_path, sep='\t')
log.info("Loading sample TSV ...")
ref = pd.read_csv(sample_tsv, sep='\t', dtype=str)
ref.columns = [c.strip() for c in ref.columns]
ref[COL_YFF2]   = pd.to_numeric(ref[COL_YFF2],   errors='coerce')
ref[COL_MSEQFF] = pd.to_numeric(ref[COL_MSEQFF], errors='coerce')

# Merge to get sample_dir and ff_ref
train_ids = set(cfg['SAMPLE_ID'].tolist())
merged = ref[ref[COL_SAMPLE_ID].isin(train_ids)].copy()
merged['_sex'] = merged[COL_GENDER].map({'XY': 'M', 'XX': 'F'})
merged['_ff_ref'] = np.where(
    merged['_sex'] == 'M',
    merged[COL_YFF2],
    merged[COL_MSEQFF]
)
log.info("Train samples reconstructed: %d (male=%d, female=%d)",
    len(merged), (merged['_sex']=='M').sum(), (merged['_sex']=='F').sum())

# Step 5: Calibration
log.info("=== Step 5: M-gxFF calibration (parallel, workers=%d) ===", n_workers)
slope, intercept = calibrate_mgxff(
    model_path, merged, WIG_TEMPLATE, 'hg19',
    ['coverage', 'fragment'], n_workers=n_workers
)
with open(calib_path, 'w') as f:
    json.dump({'slope': slope, 'intercept': intercept}, f, indent=2)
log.info("Calibration saved: slope=%.4f, intercept=%.4f", slope, intercept)

# Step 6: Validation
log.info("=== Step 6: Test set validation (parallel, workers=%d) ===", n_workers)
validate_model(
    model_path, test_json, model_dir,
    slope, intercept, WIG_TEMPLATE, 'hg19',
    ['coverage', 'fragment'], n_workers=n_workers
)
log.info("=== Done ===")
PYEOF
> "$LOG" 2>&1 &
echo "PID: $! | Log: $LOG"
