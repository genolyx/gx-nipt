#!/usr/bin/env bash
# train_gxff.sh
# ═══════════════════════════════════════════════════════════════════════
# gx-FF end-to-end training wrapper
#   1. Sample selection (male + female, stratified, low-FF guaranteed)
#   2. Train/test split (20%)
#   3. gxff train (coverage + fragment, parallel workers)
#   4. M-gxFF linear calibration
#   5. Test set validation (gxFF + M-gxFF metrics)
#
# Usage:
#   ./train_gxff.sh --out-dir refs/gxff/model_v11 [options]
#
# Options (passed through to build_gxff_model.py):
#   --n-samples        N          Total training samples (default: 900)
#   --threads          N          Parallel workers (default: 32)
#   --low-ff-weight    W          Weight for low-FF samples (default: 5.0)
#   --low-ff-min-count N          Min samples in lowest FF bin (default: 30)
#   --features         ...        coverage, fragment, nucleosome
#   --config-only                 Generate config only, no training
#   --validate-only               Validate existing model (no training)
#   --skip-calibration            Skip M-gxFF calibration
# ═══════════════════════════════════════════════════════════════════════
set -euo pipefail

GXNIPT_DIR="/home/ken/gx-nipt"
KENNIPT_DIR="/home/ken/ken-nipt"
SAMPLE_TSV="${KENNIPT_DIR}/bin/scripts/utils/reference/reference_sample_list_from_json_v3.tsv"
PATCH_SRC="${GXNIPT_DIR}/bin/gxff_patch/pipeline.py"
PATCH_DEST="/opt/conda/envs/nipt/lib/python3.11/site-packages/gxff/core/pipeline.py"
BUILD_SCRIPT="${GXNIPT_DIR}/bin/scripts/build_gxff_model.py"

# ── Parse --out-dir (required) ──────────────────────────────────────────────
OUT_DIR=""
EXTRA_ARGS=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --out-dir) OUT_DIR="$2"; shift 2 ;;
        *)         EXTRA_ARGS+=("$1"); shift ;;
    esac
done

if [[ -z "$OUT_DIR" ]]; then
    echo "Usage: $0 --out-dir <path> [options]"
    echo "  e.g. $0 --out-dir /home/ken/gx-nipt/refs/gxff/model_v11 --threads 32"
    exit 1
fi

# Make absolute
OUT_DIR="$(realpath -m "$OUT_DIR")"
mkdir -p "$OUT_DIR"

LOG_FILE="${OUT_DIR}/build_gxff.log"

echo "============================================================"
echo "  gx-FF training: $OUT_DIR"
echo "  Log: $LOG_FILE"
echo "  Extra args: ${EXTRA_ARGS[*]:-<none>}"
echo "============================================================"

nohup docker run --rm \
    --cpus=120 \
    --memory=200g \
    -v "${GXNIPT_DIR}:${GXNIPT_DIR}" \
    -v "${KENNIPT_DIR}:${KENNIPT_DIR}" \
    -v "${PATCH_SRC}:${PATCH_DEST}:ro" \
    -e HOME=/root \
    -e OMP_NUM_THREADS=120 \
    -e OPENBLAS_NUM_THREADS=1 \
    -e MKL_NUM_THREADS=120 \
    -e NUMEXPR_NUM_THREADS=120 \
    gx-nipt:latest \
    /opt/conda/envs/nipt/bin/python "${BUILD_SCRIPT}" \
        --sample-tsv  "${SAMPLE_TSV}" \
        --out-dir     "${OUT_DIR}" \
        --n-samples   900 \
        --threads     120 \
        --low-ff-weight    5.0 \
        --low-ff-min-count 30 \
        --features coverage fragment \
        "${EXTRA_ARGS[@]}" \
    > "${LOG_FILE}" 2>&1 &

PID=$!
echo "Docker container PID (nohup): $PID"
echo ""
echo "Monitor with:"
echo "  tail -f ${LOG_FILE}"
echo ""
echo "After training:"
echo "  grep -E 'VALIDATION|mae_pct|pearson' ${LOG_FILE}"
