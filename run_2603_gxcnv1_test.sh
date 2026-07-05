#!/usr/bin/env bash
# 2603 GNCI 샘플 3개 — gxcnv1/gxcnv2 통합 테스트
set -euo pipefail

NEXTFLOW=/home/ken/gx-exome/nextflow
REPO=/home/ken/gx-nipt
REF_DIR=${REPO}/refs
LABCODE=cordlife
WORK_DIR=test_2603_gxcnv1

SAMPLES=(
    "GNCI26030001:30"
    "GNCI26030002:30"
    "GNCI26030005:30"
)

mkdir -p ${REPO}/log/${WORK_DIR}

for entry in "${SAMPLES[@]}"; do
    SID=$(echo $entry | cut -d: -f1)
    AGE=$(echo $entry | cut -d: -f2)
    BAM=/home/ken/ken-nipt/analysis/2603/${SID}/${SID}.proper_paired.bam

    if [ ! -f "$BAM" ]; then
        echo "[SKIP] $SID - BAM not found: $BAM"
        continue
    fi

    OUTDIR=${REPO}/output/${WORK_DIR}/${SID}
    ANALYSISDIR=${REPO}/analysis/${WORK_DIR}/${SID}
    mkdir -p ${OUTDIR} ${ANALYSISDIR} ${REPO}/log/${WORK_DIR}/${SID}

    echo ""
    echo "=========================================="
    echo " Running: $SID"
    echo "=========================================="

    ${NEXTFLOW} run ${REPO}/main.nf \
        --sample_name  ${SID} \
        --labcode      ${LABCODE} \
        --root_dir     ${REPO} \
        --work_dir     ${WORK_DIR} \
        --outdir       ${OUTDIR} \
        --analysisdir  ${ANALYSISDIR} \
        --from_bam     ${BAM} \
        --age          ${AGE} \
        --force        true \
        --gxff_model   ${REF_DIR}/gxff/model_v5/gxff_model.pkl \
        --run_gxcnv    false \
        --run_gxcnv2   true \
        --run_gxcnv1   true \
        --run_wcx      true \
        --ref_dir      ${REF_DIR} \
        --tracedir     ${REPO}/log/${WORK_DIR}/${SID}/pipeline_info \
        2>&1 | tee ${REPO}/log/${WORK_DIR}/${SID}/pipeline.log

    echo "[DONE] $SID"
done

echo ""
echo "=========================================="
echo " 모든 샘플 완료"
echo "=========================================="
