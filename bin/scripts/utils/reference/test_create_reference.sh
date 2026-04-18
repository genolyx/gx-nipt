#!/bin/bash

# NIPT Reference Creation Pipeline Test Script

echo "========================================="
echo "NIPT Reference Creation Pipeline Test"
echo "========================================="
echo ""

# 설정
SAMPLE_LIST="/home/ken/ken-nipt/reference_sample_list.tsv"
LABCODE="ucl"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# 1. 스크립트 존재 확인
if [ ! -f "$SCRIPT_DIR/create_reference.py" ]; then
    echo "Error: create_reference.py not found in $SCRIPT_DIR"
    exit 1
fi

# 2. 샘플 리스트 존재 확인
if [ ! -f "$SAMPLE_LIST" ]; then
    echo "Error: Sample list not found: $SAMPLE_LIST"
    exit 1
fi

# 3. Preview mode 테스트
echo "Step 1: Preview sample filtering..."
echo "-----------------------------------"
python3 $SCRIPT_DIR/create_reference.py \
    --sample-list $SAMPLE_LIST \
    --preview-only

if [ $? -ne 0 ]; then
    echo "Error: Preview failed"
    exit 1
fi

echo ""
echo "Preview completed successfully!"
echo ""

# 4. 사용자 확인
read -p "Do you want to create references? (y/N): " answer

if [ "$answer" != "y" ] && [ "$answer" != "Y" ]; then
    echo "Reference creation cancelled."
    exit 0
fi

# 5. Reference 생성 (EZD만 우선)
echo ""
echo "Step 2: Creating EZD reference..."
echo "-----------------------------------"
python3 $SCRIPT_DIR/create_reference.py \
    --sample-list $SAMPLE_LIST \
    --labcode $LABCODE \
    --ref-type ezd \
    --groups orig

if [ $? -ne 0 ]; then
    echo "Error: EZD reference creation failed"
    exit 1
fi

echo ""
echo "========================================="
echo "Test completed successfully!"
echo "========================================="
echo ""
echo "Next steps:"
echo "1. Check filtered sample list: ${SAMPLE_LIST%.tsv}_filtered.tsv"
echo "2. Check reference files: /Work/NIPT/data/refs/$LABCODE/EZD/orig/"
echo "3. Run full pipeline: python3 $SCRIPT_DIR/create_reference.py --sample-list $SAMPLE_LIST --labcode $LABCODE --ref-type all"

