#!/usr/bin/env bash
# =========================================================
# test_compare.sh
# gx-nipt vs ken-nipt comparison test (3 cordlife samples)
#
# Runs each sample through the full gx-nipt pipeline
# (alignment → FF/gender → WCX → gx-FF → gx-cnv),
# then prints a side-by-side comparison with ken-nipt results.
# =========================================================
set -euo pipefail

LABCODE="cordlife"
WORK_DIR="2603"
AGE=29

GXNIPT_ROOT="/home/ken/gx-nipt"
KENMIPT_ROOT="/home/ken/ken-nipt"
GXFF_MODEL="${GXNIPT_ROOT}/refs/gxff/model_v2/gxff_model.pkl"

# Sample list: ORDER_ID  R1_PATH  R2_PATH
declare -A R1 R2
R1[GNCI26030001]="${KENMIPT_ROOT}/fastq/2603/GNCI26030001/2602270008_S8_R1_001.fastq.gz"
R2[GNCI26030001]="${KENMIPT_ROOT}/fastq/2603/GNCI26030001/2602270008_S8_R2_001.fastq.gz"
R1[GNCI26030002]="${KENMIPT_ROOT}/fastq/2603/GNCI26030002/2603020001_S5_R1_001.fastq.gz"
R2[GNCI26030002]="${KENMIPT_ROOT}/fastq/2603/GNCI26030002/2603020001_S5_R2_001.fastq.gz"
R1[GNCI26030005]="${KENMIPT_ROOT}/fastq/2603/GNCI26030005/2603030005_S21_R1_001.fastq.gz"
R2[GNCI26030005]="${KENMIPT_ROOT}/fastq/2603/GNCI26030005/2603030005_S21_R2_001.fastq.gz"

SAMPLES=(GNCI26030001 GNCI26030002 GNCI26030005)

run_sample() {
    local SAMPLE="$1"
    echo ""
    echo "======================================================"
    echo "  Running: $SAMPLE"
    echo "======================================================"

    mkdir -p "${GXNIPT_ROOT}/log/${WORK_DIR}/${SAMPLE}"
    mkdir -p "${GXNIPT_ROOT}/analysis/${WORK_DIR}/${SAMPLE}"
    mkdir -p "${GXNIPT_ROOT}/output/${WORK_DIR}/${SAMPLE}"

    local LOG="${GXNIPT_ROOT}/log/${WORK_DIR}/${SAMPLE}/pipeline_stdout.log"

    docker run --rm \
      -u 1002:1000 --group-add 988 \
      -v "${GXNIPT_ROOT}:/home/ken/gx-nipt" \
      -v "${GXNIPT_ROOT}:/app" \
      -v "${KENMIPT_ROOT}:/home/ken/ken-nipt" \
      -v /var/run/docker.sock:/var/run/docker.sock \
      -e HOME=/tmp \
      -e NXF_HOME=/tmp/.nextflow \
      gx-daemon:latest \
      nextflow \
        -log "${GXNIPT_ROOT}/log/${WORK_DIR}/${SAMPLE}/nextflow.log" \
        run /app/main.nf \
        -profile standard \
        -ansi-log false \
        -work-dir "${GXNIPT_ROOT}/analysis/${WORK_DIR}/${SAMPLE}/work" \
        --sample_name "${SAMPLE}" \
        --labcode "${LABCODE}" \
        --fastq_r1 "${R1[$SAMPLE]}" \
        --fastq_r2 "${R2[$SAMPLE]}" \
        --root_dir "${GXNIPT_ROOT}" \
        --work_dir "${WORK_DIR}" \
        --analysisdir "${GXNIPT_ROOT}/analysis/${WORK_DIR}/${SAMPLE}" \
        --outdir "${GXNIPT_ROOT}/output/${WORK_DIR}/${SAMPLE}" \
        --ref_dir "${GXNIPT_ROOT}/refs" \
        --scratch_dir /tmp/nipt_scratch \
        --tracedir "${GXNIPT_ROOT}/log/${WORK_DIR}/${SAMPLE}" \
        --age "${AGE}" \
        --gxff_model "${GXFF_MODEL}" \
        --run_gxcnv true \
        >"${LOG}" 2>&1
    local EXIT=$?

    if [[ $EXIT -eq 0 ]]; then
        echo "  ✅ DONE (log: $LOG)"
    else
        echo "  ❌ FAILED (exit=$EXIT, log: $LOG)"
        tail -20 "$LOG"
    fi
    return $EXIT
}

print_comparison() {
    echo ""
    echo "======================================================"
    echo "  COMPARISON: gx-nipt vs ken-nipt"
    echo "======================================================"
    printf "%-18s %-10s %-10s %-10s %-10s %-10s %-10s %-10s\n" \
        "SAMPLE" "TOOL" "GENDER" "SeqFF" "gxFF" "YFF_2" "M-SeqFF" "QC"
    echo "--------------------------------------------------------------------------------------------------------------------------------------"

    for SAMPLE in "${SAMPLES[@]}"; do
        # ── ken-nipt results ──────────────────────────────────────────────
        local KEN_GENDER KEN_SEQFF KEN_YFF2
        KEN_GENDER=$(awk 'NR==3{print $3}' \
            "${KENMIPT_ROOT}/analysis/${WORK_DIR}/${SAMPLE}/Output_FF/${SAMPLE}.gender.txt" 2>/dev/null || echo "N/A")
        KEN_SEQFF=$(awk -F'[,]' 'NR==2{gsub(/"/,"",$2); printf "%.2f%%", $2*100}' \
            "${KENMIPT_ROOT}/analysis/${WORK_DIR}/${SAMPLE}/Output_FF/${SAMPLE}.seqff.txt" 2>/dev/null || echo "N/A")

        printf "%-18s %-10s %-10s %-10s %-10s %-10s %-10s %-10s\n" \
            "$SAMPLE" "ken-nipt" "$KEN_GENDER" "$KEN_SEQFF" "-" "-" "-" "-"

        # ── gx-nipt results ───────────────────────────────────────────────
        local GXN_DIR="${GXNIPT_ROOT}/output/${WORK_DIR}/${SAMPLE}"
        local GXN_JSON="${GXN_DIR}/${SAMPLE}.json"

        if [[ -f "$GXN_JSON" ]]; then
            python3 - <<PYEOF
import json, sys
try:
    d = json.load(open("$GXN_JSON"))["final_results"]
    gender   = d.get("fetal_gender",         "N/A")
    seqff    = d.get("fetal_fraction_seqff",  "N/A")
    gxff     = d.get("fetal_fraction_gxff",   "N/A")
    yff2     = d.get("fetal_fraction_yff",    "N/A")
    mseqff   = d.get("fetal_fraction_mseqff", "N/A")
    qc       = d.get("QC_result",             "N/A")
    print(f"{'$SAMPLE':<18} {'gx-nipt':<10} {gender:<10} {seqff:<10} {gxff:<10} {yff2:<10} {mseqff:<10} {qc:<10}")
except Exception as e:
    print(f"{'$SAMPLE':<18} {'gx-nipt':<10} ERROR: {e}")
PYEOF
        else
            # Fall back to individual output files
            local GENDER_FILE="${GXNIPT_ROOT}/analysis/${WORK_DIR}/${SAMPLE}/${SAMPLE}/Output_ff_gender/gender.txt"
            local GXN_GENDER GXN_SEQFF GXN_GXFF GXN_YFF2
            GXN_GENDER=$(grep -i "final_gender" "$GENDER_FILE" 2>/dev/null | awk '{print $2}' || echo "N/A")
            printf "%-18s %-10s %-10s %-10s %-10s %-10s %-10s %-10s\n" \
                "$SAMPLE" "gx-nipt" "$GXN_GENDER" "?" "?" "?" "?" "?"
        fi

        echo ""
    done

    # Also show gx-cnv calls per sample
    echo ""
    echo "gx-cnv HIGH_RISK calls:"
    for SAMPLE in "${SAMPLES[@]}"; do
        local CALLS_DIR="${GXNIPT_ROOT}/analysis/${WORK_DIR}/${SAMPLE}/${SAMPLE}/gxcnv"
        local CALLS_TSV=$(ls "${CALLS_DIR}/${SAMPLE}_calls.tsv" 2>/dev/null | head -1)
        if [[ -n "$CALLS_TSV" ]]; then
            local N_CALLS
            N_CALLS=$(grep -v "^#" "$CALLS_TSV" 2>/dev/null | wc -l)
            echo "  $SAMPLE: $N_CALLS HIGH_RISK region(s)"
        else
            echo "  $SAMPLE: (gx-cnv output not found)"
        fi
    done
}

# ── Main ──────────────────────────────────────────────────────────────────────
echo "gx-nipt comparison test"
echo "  Samples : ${SAMPLES[*]}"
echo "  gx-FF   : $GXFF_MODEL"
echo "  gx-cnv  : ENABLED (auto gender-aware reference)"
echo ""

FAILED=0
for SAMPLE in "${SAMPLES[@]}"; do
    run_sample "$SAMPLE" || FAILED=$((FAILED+1))
done

print_comparison

echo ""
if [[ $FAILED -eq 0 ]]; then
    echo "All ${#SAMPLES[@]} samples completed successfully."
else
    echo "$FAILED / ${#SAMPLES[@]} samples FAILED."
    exit 1
fi
