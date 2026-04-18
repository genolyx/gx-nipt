#!/bin/bash
# Test script for integrated sample list generation

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../../../" && pwd)"

echo "=========================================="
echo "Test: Integrated Sample List Generation"
echo "=========================================="
echo ""

# Test 1: Generate sample list only
echo "Test 1: Generate sample list from analysis dirs"
python3 "${SCRIPT_DIR}/create_reference.py" \
    --generate-sample-list \
    --analysis-dirs "${REPO_ROOT}/analysis/2507" "${REPO_ROOT}/analysis/2508" \
    --prefix GNMF \
    --sample-list /tmp/test_sample_list.tsv

if [ $? -eq 0 ]; then
    echo "✓ Sample list generated successfully"
    echo "  Lines: $(wc -l < /tmp/test_sample_list.tsv)"
    echo "  Preview:"
    head -3 /tmp/test_sample_list.tsv
else
    echo "✗ Sample list generation failed"
    exit 1
fi

echo ""
echo "Test 2: Generate sample list AND create reference (preview)"
python3 "${SCRIPT_DIR}/create_reference.py" \
    --analysis-dirs "${REPO_ROOT}/analysis/2507" \
    --prefix GNMF \
    --labcode ucl \
    --preview-only

if [ $? -eq 0 ]; then
    echo "✓ Integrated workflow successful"
else
    echo "✗ Integrated workflow failed"
    exit 1
fi

echo ""
echo "=========================================="
echo "All tests passed!"
echo "=========================================="
