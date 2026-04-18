#!/bin/bash
# WisecondorX Reference 생성 스크립트 (Docker 컨테이너 내부용)

set -e

SAMPLE_LIST=$1
GROUP=$2
OUTPUT_DIR=$3

echo "=== WisecondorX Reference 생성 ==="
echo "Group: $GROUP"
echo "Output: $OUTPUT_DIR"

# 출력 디렉토리 생성
mkdir -p "$OUTPUT_DIR"

# 샘플 리스트에서 Male/Female NPZ 파일 수집
MALE_FILES=""
FEMALE_FILES=""

# 헤더 제외하고 읽기
tail -n +2 "$SAMPLE_LIST" | while IFS=$'\t' read -r month sample_id sample_dir rest; do
    gender=$(echo "$rest" | cut -f9)
    
    # NPZ 파일 경로 구성
    if [ "$GROUP" == "orig" ]; then
        npz_file="${sample_dir}/Output_WCX/${sample_id}.wcx.of_orig.npz"
    elif [ "$GROUP" == "fetus" ]; then
        npz_file="${sample_dir}/Output_WCX/${sample_id}.wcx.of_fetus.npz"
    elif [ "$GROUP" == "mom" ]; then
        npz_file="${sample_dir}/Output_WCX/${sample_id}.wcx.of_mom.npz"
    fi
    
    if [ -f "$npz_file" ]; then
        if [ "$gender" == "XY" ]; then
            MALE_FILES="$MALE_FILES $npz_file"
        elif [ "$gender" == "XX" ]; then
            FEMALE_FILES="$FEMALE_FILES $npz_file"
        fi
    fi
done

# NPZ 파일 개수 계산
MALE_COUNT=$(echo $MALE_FILES | wc -w)
FEMALE_COUNT=$(echo $FEMALE_FILES | wc -w)

echo "Male NPZ files: $MALE_COUNT"
echo "Female NPZ files: $FEMALE_COUNT"

# WisecondorX newref 실행
WCX_BIN="/opt/conda/envs/nipt/bin/WisecondorX"

if [ "$GROUP" == "mom" ]; then
    # mom은 combined reference (female 중심)
    OUTPUT_NPZ="${OUTPUT_DIR}/${GROUP}_200k_of.npz"
    echo "Creating combined reference for mom..."
    $WCX_BIN newref $FEMALE_FILES $MALE_FILES "$OUTPUT_NPZ" --binsize 200000 --yfrac 0
    echo "Created: $OUTPUT_NPZ"
else
    # orig, fetus는 Male/Female 각각 생성
    if [ $MALE_COUNT -gt 0 ]; then
        OUTPUT_NPZ="${OUTPUT_DIR}/${GROUP}_M_200k_of.npz"
        echo "Creating male reference..."
        $WCX_BIN newref $MALE_FILES "$OUTPUT_NPZ" --binsize 200000 --nipt --yfrac 0
        echo "Created: $OUTPUT_NPZ"
    fi
    
    if [ $FEMALE_COUNT -gt 0 ]; then
        OUTPUT_NPZ="${OUTPUT_DIR}/${GROUP}_F_200k_of.npz"
        echo "Creating female reference..."
        $WCX_BIN newref $FEMALE_FILES "$OUTPUT_NPZ" --binsize 200000 --nipt --yfrac 0
        echo "Created: $OUTPUT_NPZ"
    fi
fi

echo "=== WCX Reference 생성 완료 ==="
