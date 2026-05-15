#!/usr/bin/env bash
# test set validation for model_v6 and model_v10
set -e

GXNIPT_DIR="/home/ken/gx-nipt"
KENNIPT_DIR="/home/ken/ken-nipt"

docker run --rm \
  -v "${GXNIPT_DIR}:${GXNIPT_DIR}" \
  -v "${KENNIPT_DIR}:${KENNIPT_DIR}" \
  -e HOME=/root \
  gx-nipt:latest \
  /opt/conda/envs/nipt/bin/python "${GXNIPT_DIR}/bin/scripts/validate_gxff.py" \
    --test-json "${GXNIPT_DIR}/refs/gxff/test_set_v6.json" \
    --models \
      "v6=${GXNIPT_DIR}/refs/gxff/model_v6/gxff_model.pkl" \
      "v10=${GXNIPT_DIR}/refs/gxff/model_v10/gxff_model.pkl" \
    --out "${GXNIPT_DIR}/refs/gxff/validation_results.tsv" \
    --genome hg19

echo ""
echo "Done! Results:"
echo "  Summary  : ${GXNIPT_DIR}/refs/gxff/validation_results.tsv"
echo "  Per-sample: ${GXNIPT_DIR}/refs/gxff/validation_results_per_sample.tsv"
