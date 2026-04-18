#!/bin/bash
# WC/WCX Reference 재생성 스크립트 (버그 수정 후)

set -e  # 에러 발생시 중단

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../../../" && pwd)"

SAMPLE_LIST="${REPO_ROOT}/data/refs/ucl/reference_make/reference_sample_list_UCL_filtered.tsv"
LABCODE="ucl"
OUTPUT_DIR="${REPO_ROOT}/data/refs/ucl_new"
SOURCE_DIR="${REPO_ROOT}/data/refs/ucl"

echo "=========================================="
echo "WC/WCX Reference 재생성 (버그 수정 버전)"
echo "=========================================="
echo ""
echo "Sample List: $SAMPLE_LIST"
echo "Output Dir: $OUTPUT_DIR"
echo "Source Dir: $SOURCE_DIR"
echo ""

# 1. 기존 파일 백업
echo "Step 1: 기존 파일 백업..."
mkdir -p ${OUTPUT_DIR}/WCX_buggy_backup
mkdir -p ${OUTPUT_DIR}/WC_buggy_backup

if ls ${OUTPUT_DIR}/WCX/*.npz 1> /dev/null 2>&1; then
    mv ${OUTPUT_DIR}/WCX/*.npz ${OUTPUT_DIR}/WCX_buggy_backup/
    echo "  ✓ WCX 백업 완료: $(ls ${OUTPUT_DIR}/WCX_buggy_backup/*.npz | wc -l) 파일"
else
    echo "  - WCX 파일 없음 (skip)"
fi

if ls ${OUTPUT_DIR}/WC/*.npz 1> /dev/null 2>&1; then
    mv ${OUTPUT_DIR}/WC/*.npz ${OUTPUT_DIR}/WC_buggy_backup/
    echo "  ✓ WC 백업 완료: $(ls ${OUTPUT_DIR}/WC_buggy_backup/*.npz | wc -l) 파일"
else
    echo "  - WC 파일 없음 (skip)"
fi

echo ""

# 2. WCX 재생성
echo "Step 2: WCX Reference 재생성..."
echo "----------------------------------------"
python3 "${SCRIPT_DIR}/create_reference.py" \
    --sample-list ${SAMPLE_LIST} \
    --labcode ${LABCODE} \
    --output-dir ${OUTPUT_DIR} \
    --reference-source ${SOURCE_DIR} \
    --ref-type wcx

echo ""
echo "✓ WCX 생성 완료"
echo "생성된 파일:"
ls -lh ${OUTPUT_DIR}/WCX/*.npz | awk '{print "  - " $9 " (" $5 ")"}'
echo ""

# 3. WC 재생성
echo "Step 3: WC Reference 재생성..."
echo "----------------------------------------"
python3 "${SCRIPT_DIR}/create_reference.py" \
    --sample-list ${SAMPLE_LIST} \
    --labcode ${LABCODE} \
    --output-dir ${OUTPUT_DIR} \
    --reference-source ${SOURCE_DIR} \
    --ref-type wc

echo ""
echo "✓ WC 생성 완료"
echo "생성된 파일:"
ls -lh ${OUTPUT_DIR}/WC/*.npz | awk '{print "  - " $9 " (" $5 ")"}'
echo ""

# 4. 최종 요약
echo "=========================================="
echo "✅ WC/WCX Reference 재생성 완료!"
echo "=========================================="
echo ""
echo "📊 생성된 파일 요약:"
echo ""
echo "WCX (7개 파일):"
ls -lh ${OUTPUT_DIR}/WCX/*.npz | wc -l | awk '{print "  Total: " $1 " files"}'
echo ""
echo "WC (3개 파일):"
ls -lh ${OUTPUT_DIR}/WC/*.npz | wc -l | awk '{print "  Total: " $1 " files"}'
echo ""
echo "💾 백업 위치:"
echo "  - ${OUTPUT_DIR}/WCX_buggy_backup/"
echo "  - ${OUTPUT_DIR}/WC_buggy_backup/"
echo ""
echo "🎯 다음 단계: 테스트 샘플로 품질 검증"
echo ""
