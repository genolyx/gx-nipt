#!/bin/bash
###############################################################################
# WisecondorX Male Reference 재생성 스크립트
# 
# 용도: orig_M_200k_proper_paired.npz, fetus_M_200k_of.npz만 재생성
# 
# 사용법:
#   bash recreate_wcx_male_refs.sh
###############################################################################

set -e

# 색상 정의
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# 경로 설정
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKSPACE_ROOT="$(cd "${SCRIPT_DIR}/../../../../" && pwd)"
DATA_DIR="${WORKSPACE_ROOT}/data"
ANALYSIS_DIR="${WORKSPACE_ROOT}/analysis"

# Docker 설정
DOCKER_IMAGE="nipt_docker_v1.3"
SAMPLE_LIST="${DATA_DIR}/refs/ucl/reference_make/reference_sample_list_UCL_filtered.tsv"
OUTPUT_DIR="${DATA_DIR}/refs/ucl_new/WCX"

echo -e "${BLUE}==========================================================${NC}"
echo -e "${BLUE}WisecondorX Male Reference 재생성${NC}"
echo -e "${BLUE}==========================================================${NC}"
echo ""
echo "Docker Image: ${DOCKER_IMAGE}"
echo "Sample List: ${SAMPLE_LIST}"
echo "Output Dir: ${OUTPUT_DIR}"
echo ""

# 샘플 리스트 확인
if [ ! -f "${SAMPLE_LIST}" ]; then
    echo -e "${RED}Error: Sample list not found: ${SAMPLE_LIST}${NC}"
    exit 1
fi

# 출력 디렉토리 확인/생성
mkdir -p "${OUTPUT_DIR}"

# 기존 Male reference 백업
echo -e "${YELLOW}[1/4] 기존 Male reference 백업${NC}"
BACKUP_DIR="${OUTPUT_DIR}/backup_$(date +%Y%m%d_%H%M%S)"
mkdir -p "${BACKUP_DIR}"

if [ -f "${OUTPUT_DIR}/orig_M_200k_proper_paired.npz" ]; then
    cp -v "${OUTPUT_DIR}/orig_M_200k_proper_paired.npz" "${BACKUP_DIR}/"
    echo -e "${GREEN}✓ orig_M_200k_proper_paired.npz 백업 완료${NC}"
else
    echo -e "${YELLOW}⚠ orig_M_200k_proper_paired.npz 파일 없음${NC}"
fi

if [ -f "${OUTPUT_DIR}/fetus_M_200k_of.npz" ]; then
    cp -v "${OUTPUT_DIR}/fetus_M_200k_of.npz" "${BACKUP_DIR}/"
    echo -e "${GREEN}✓ fetus_M_200k_of.npz 백업 완료${NC}"
else
    echo -e "${YELLOW}⚠ fetus_M_200k_of.npz 파일 없음${NC}"
fi
echo ""

# Docker 명령 준비
SAMPLE_LIST_DOCKER="/refs/ucl/reference_make/reference_sample_list_UCL_filtered.tsv"
OUTPUT_DIR_DOCKER="/refs/ucl_new/WCX"

echo -e "${YELLOW}[2/4] Docker 컨테이너에서 Male reference 재생성${NC}"
echo "Running: python3 /scripts/recreate_wcx_male_refs.py ..."
echo ""

# Docker 실행
docker run --rm \
    --entrypoint bash \
    -v "${ANALYSIS_DIR}:/analysis:ro" \
    -v "${DATA_DIR}/refs:/refs" \
    -v "${SCRIPT_DIR}:/scripts:ro" \
    "${DOCKER_IMAGE}" \
    -c "python3 /scripts/recreate_wcx_male_refs.py ${SAMPLE_LIST_DOCKER} ${OUTPUT_DIR_DOCKER} orig fetus"

DOCKER_EXIT_CODE=$?

if [ ${DOCKER_EXIT_CODE} -ne 0 ]; then
    echo -e "${RED}✗ Docker 실행 실패 (exit code: ${DOCKER_EXIT_CODE})${NC}"
    exit 1
fi

echo ""
echo -e "${YELLOW}[3/4] 생성된 파일 확인${NC}"
echo ""

# 생성된 파일 확인
if [ -f "${OUTPUT_DIR}/orig_M_200k_proper_paired.npz" ]; then
    echo -e "${GREEN}✓ orig_M_200k_proper_paired.npz${NC}"
    ls -lh "${OUTPUT_DIR}/orig_M_200k_proper_paired.npz"
else
    echo -e "${RED}✗ orig_M_200k_proper_paired.npz 생성 실패${NC}"
fi

if [ -f "${OUTPUT_DIR}/fetus_M_200k_of.npz" ]; then
    echo -e "${GREEN}✓ fetus_M_200k_of.npz${NC}"
    ls -lh "${OUTPUT_DIR}/fetus_M_200k_of.npz"
else
    echo -e "${RED}✗ fetus_M_200k_of.npz 생성 실패${NC}"
fi

echo ""
echo -e "${YELLOW}[4/4] 전/후 비교${NC}"
echo ""

echo "=== 백업 (이전) ==="
ls -lh "${BACKUP_DIR}/" 2>/dev/null || echo "백업 파일 없음"

echo ""
echo "=== 현재 (새로 생성) ==="
ls -lh "${OUTPUT_DIR}"/orig_M_200k_proper_paired.npz "${OUTPUT_DIR}"/fetus_M_200k_of.npz 2>/dev/null || echo "파일 없음"

echo ""
echo -e "${GREEN}==========================================================${NC}"
echo -e "${GREEN}완료!${NC}"
echo -e "${GREEN}==========================================================${NC}"
echo ""
echo "백업 위치: ${BACKUP_DIR}"
echo ""
echo "복원이 필요한 경우:"
echo "  cp ${BACKUP_DIR}/*.npz ${OUTPUT_DIR}/"
echo ""
