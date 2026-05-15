#!/usr/bin/env bash
# UCL High Risk 샘플 gxcnv2 통합 테스트 (새 reference ucl_batch_0511)
set -euo pipefail

NEXTFLOW=/home/ken/gx-exome/nextflow
REPO=/home/ken/gx-nipt
REF_DIR=${REPO}/refs
LABCODE=ucl_batch_0511
WORK_DIR=test_ucl_gxcnv2

# 선별된 High Risk 샘플 (sample_id:month:age:disease)
SAMPLES=(
    "GNMF26020033:2602:35:T21_XX"
    "GNMF26040103:2604:38:T21_XX"
    "GNMF25120084:2512:32:T18_XY"
    "GNMF26020042:2602:40:T18+XYY_XY"
    "GNMF25090053:2509:29:T13_XY"
    "GNMF26010029:2601:36:T22_XY"
    "GNMF25070004:2507:33:XO_XX"
)

mkdir -p ${REPO}/log/${WORK_DIR}

for entry in "${SAMPLES[@]}"; do
    SID=$(echo $entry | cut -d: -f1)
    MON=$(echo $entry | cut -d: -f2)
    AGE=$(echo $entry | cut -d: -f3)
    DISEASE=$(echo $entry | cut -d: -f4)
    BAM=/home/ken/ken-nipt/analysis/${MON}/${SID}/${SID}.proper_paired.bam

    if [ ! -f "$BAM" ]; then
        echo "[SKIP] $SID - BAM not found: $BAM"
        continue
    fi

    OUTDIR=${REPO}/output/${WORK_DIR}/${SID}
    ANALYSISDIR=${REPO}/analysis/${WORK_DIR}/${SID}
    mkdir -p ${OUTDIR} ${ANALYSISDIR} ${REPO}/log/${WORK_DIR}/${SID}

    echo ""
    echo "=========================================="
    echo " Running: $SID  ($DISEASE)"
    echo "=========================================="

    nextflow run ${REPO}/main.nf \
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
        --run_gxcnv    true \
        --run_gxcnv2   true \
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
