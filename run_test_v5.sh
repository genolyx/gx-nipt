#!/usr/bin/env bash
# Quick test runner for 3 samples with model_v5
set -euo pipefail

GXFF_MODEL="/home/ken/gx-nipt/refs/gxff/model_v5/gxff_model.pkl"
REF_DIR="/home/ken/gx-nipt/refs"
NF_BIN="/home/ken/gx-exome/nextflow"
PIDS=()

for SAMPLE in GNCI26030002 GNCI26030001 GNCI26030005; do
  BAM="/home/ken/ken-nipt/analysis/2603/${SAMPLE}/${SAMPLE}.proper_paired.bam"
  echo "[$(date +%T)] Starting $SAMPLE ..."
  bash /home/ken/gx-nipt/bin/run_nipt.sh \
    --sample-name  "$SAMPLE" \
    --order-id     "$SAMPLE" \
    --work-dir     "2603" \
    --root-dir     "/home/ken/gx-nipt" \
    --labcode      "cordlife" \
    --ref-dir      "$REF_DIR" \
    --from-bam     "$BAM" \
    --nextflow     "$NF_BIN" \
    --force \
    --fresh \
    --gxff-model   "$GXFF_MODEL" \
    --run-gxcnv \
    > /tmp/gxnipt_${SAMPLE}.log 2>&1 &
  PIDS+=($!)
  echo "  PID=${PIDS[-1]} -> /tmp/gxnipt_${SAMPLE}.log"
done

echo "[$(date +%T)] Waiting for all 3 samples..."
for PID in "${PIDS[@]}"; do
  wait "$PID" || true
done

echo ""
echo "=== RESULTS ==="
for SAMPLE in GNCI26030002 GNCI26030001 GNCI26030005; do
  STATUS=$(tail -1 /tmp/gxnipt_${SAMPLE}.log 2>/dev/null)
  echo "$SAMPLE: $STATUS"
done
