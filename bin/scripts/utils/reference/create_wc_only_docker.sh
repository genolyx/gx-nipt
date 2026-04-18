#!/bin/bash
# 
# Wisecondor Reference 생성 스크립트 (WC만)
# Docker 컨테이너 기반 실행
#

set -e  # 에러 발생시 중단

REF_UTILS="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${REF_UTILS}/../../../../" && pwd)"

# 설정
SAMPLE_LIST="/refs/ucl/reference_make/reference_sample_list_UCL_filtered.tsv"
WC_OUTPUT_DIR="/refs/ucl_new/WC"
LOG_DIR="${REPO_ROOT}/data/refs/ucl_new/logs"
DOCKER_IMAGE="nipt_docker_v1.3"

# 로그 디렉토리 생성
mkdir -p "$LOG_DIR"

# 시작 시간
START_TIME=$(date +%s)
echo "=========================================="
echo "Wisecondor Reference 생성 시작: $(date)"
echo "=========================================="

# 함수: Docker 컨테이너에서 스크립트 실행
run_in_docker() {
    local script=$1
    local group=$2
    local output_dir=$3
    local log_file=$4
    
    echo ""
    echo "[$group] 시작: $(date)"
    
    docker run --rm \
        --entrypoint bash \
        -v "${REPO_ROOT}/analysis:/analysis:ro" \
        -v "${REPO_ROOT}/data/refs:/refs" \
        -v "${REF_UTILS}:/scripts:ro" \
        ${DOCKER_IMAGE} -c "python3 /scripts/${script} ${SAMPLE_LIST} ${group} ${output_dir}" \
        2>&1 | tee "$log_file"
    
    local exit_code=${PIPESTATUS[0]}
    
    if [ $exit_code -eq 0 ]; then
        echo "[$group] 완료: $(date) ✓"
    else
        echo "[$group] 실패: $(date) ✗"
        return $exit_code
    fi
}

# ===========================================
# Wisecondor Reference 생성
# ===========================================
echo ""
echo "==========================================="
echo "Wisecondor Reference 생성"
echo "==========================================="

run_in_docker "create_wc_docker.py" "orig" "${WC_OUTPUT_DIR}" "${LOG_DIR}/wc_orig.log"
run_in_docker "create_wc_docker.py" "fetus" "${WC_OUTPUT_DIR}" "${LOG_DIR}/wc_fetus.log"
run_in_docker "create_wc_docker.py" "mom" "${WC_OUTPUT_DIR}" "${LOG_DIR}/wc_mom.log"

echo ""
echo "Wisecondor 생성 파일:"
ls -lh "${REPO_ROOT}/data/refs/ucl_new/WC/"*.npz 2>/dev/null || echo "파일 없음"

# ===========================================
# 완료
# ===========================================
END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))
HOURS=$((ELAPSED / 3600))
MINUTES=$(((ELAPSED % 3600) / 60))
SECONDS=$((ELAPSED % 60))

echo ""
echo "=========================================="
echo "Wisecondor Reference 생성 완료!"
echo "소요 시간: ${HOURS}시간 ${MINUTES}분 ${SECONDS}초"
echo "완료 시각: $(date)"
echo "=========================================="
echo ""
echo "생성된 파일:"
ls -lh "${REPO_ROOT}/data/refs/ucl_new/WC/"*.npz 2>/dev/null
