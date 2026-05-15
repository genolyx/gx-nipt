#!/usr/bin/env bash
# model_v10 only re-validation (with correct coverage+fragment features)
set -e

GXNIPT_DIR="/home/ken/gx-nipt"
KENNIPT_DIR="/home/ken/ken-nipt"

nohup docker run --rm \
  -v "${GXNIPT_DIR}:${GXNIPT_DIR}" \
  -v "${KENNIPT_DIR}:${KENNIPT_DIR}" \
  -e HOME=/root \
  gx-nipt:latest \
  /opt/conda/envs/nipt/bin/python "${GXNIPT_DIR}/bin/scripts/validate_gxff.py" \
    --test-json "${GXNIPT_DIR}/refs/gxff/test_set_v6.json" \
    --models \
      "v10=${GXNIPT_DIR}/refs/gxff/model_v10/gxff_model.pkl" \
    --out "${GXNIPT_DIR}/refs/gxff/validation_v10_results.tsv" \
    --genome hg19 \
  > "${GXNIPT_DIR}/refs/gxff/validation_v10.log" 2>&1 &

echo "v10 validation PID: $!"
echo "Log: ${GXNIPT_DIR}/refs/gxff/validation_v10.log"
