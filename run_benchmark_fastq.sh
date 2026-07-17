#!/usr/bin/env bash
# FASTQ → full pipeline benchmark for GNCI26030001/002/003
set -euo pipefail

NEXTFLOW=/home/ken/gx-exome/nextflow
REPO=/home/ken/gx-nipt
REF_DIR=${REPO}/refs
export PATH="${REPO}/bin/nxf-docker:${PATH}"
LABCODE=cordlife
WORK_DIR=benchmark_fastq_2603
RESULTS=${REPO}/tmp/benchmark_fastq_results.tsv

mkdir -p "${REPO}/tmp" "${REPO}/log/${WORK_DIR}"

run_sample() {
  local SID=$1 AGE=$2 R1=$3 R2=$4
  local OUTDIR=${REPO}/output/${WORK_DIR}/${SID}
  local ANALYSISDIR=${REPO}/analysis/${WORK_DIR}/${SID}
  local LOGDIR=${REPO}/log/${WORK_DIR}/${SID}
  mkdir -p "${OUTDIR}" "${ANALYSISDIR}" "${LOGDIR}"

  echo ""
  echo "=========================================="
  echo " FASTQ benchmark: ${SID} (age=${AGE})"
  echo " R1: ${R1}"
  echo " R2: ${R2}"
  echo "=========================================="

  local START_TS START_SEC END_TS END_SEC WALL_SEC WALL_MIN EXIT_CODE STATUS NXF_DUR
  START_TS=$(date -Iseconds)
  START_SEC=$(date +%s)

  set +e
  ${NEXTFLOW} run ${REPO}/main.nf \
    --sample_name  "${SID}" \
    --labcode      "${LABCODE}" \
    --root_dir     "${REPO}" \
    --work_dir     "${WORK_DIR}" \
    --outdir       "${OUTDIR}" \
    --analysisdir  "${ANALYSISDIR}" \
    --fastq_r1     "${R1}" \
    --fastq_r2     "${R2}" \
    --age          "${AGE}" \
    --force        true \
    --gxff_model   "${REF_DIR}/gxff/model_v5/gxff_model.pkl" \
    --run_gxcnv    false \
    --run_gxcnv2   true \
    --run_gxcnv1   true \
    --run_wcx      true \
    --ref_dir      "${REF_DIR}" \
    --tracedir     "${LOGDIR}/pipeline_info" \
    2>&1 | tee "${LOGDIR}/pipeline_fastq.log"
  EXIT_CODE=${PIPESTATUS[0]}
  set -e

  END_TS=$(date -Iseconds)
  END_SEC=$(date +%s)
  WALL_SEC=$((END_SEC - START_SEC))
  WALL_MIN=$(awk "BEGIN {printf \"%.2f\", ${WALL_SEC}/60}")

  if [[ ${EXIT_CODE} -eq 0 ]]; then STATUS="SUCCESS"; else STATUS="FAILED"; fi
  NXF_DUR=$(grep 'Duration  :' "${LOGDIR}/pipeline_fastq.log" | tail -1 | sed 's/.*Duration  : //' || echo "N/A")

  echo -e "${SID}\t${START_TS}\t${END_TS}\t${WALL_SEC}\t${WALL_MIN}\t${STATUS}\t${NXF_DUR}" >> "${RESULTS}"
  echo "[DONE] ${SID}: wall=${WALL_MIN}min (${WALL_SEC}s) status=${STATUS} nxf=${NXF_DUR}"
}

echo -e "sample\tstart\tend\twall_sec\twall_min\tstatus\tnxf_duration" > "${RESULTS}"
TOTAL_START=$(date +%s)

run_sample GNCI26030001 26 \
  /home/ken/ken-nipt/fastq/2603/GNCI26030001/2602270008_S8_R1_001.fastq.gz \
  /home/ken/ken-nipt/fastq/2603/GNCI26030001/2602270008_S8_R2_001.fastq.gz

run_sample GNCI26030002 29 \
  /home/ken/ken-nipt/fastq/2603/GNCI26030002/2603020001_S5_R1_001.fastq.gz \
  /home/ken/ken-nipt/fastq/2603/GNCI26030002/2603020001_S5_R2_001.fastq.gz

run_sample GNCI26030003 26 \
  /home/ken/ken-nipt/fastq/2603/GNCI26030003/2602250001_S13_R1_001.fastq.gz \
  /home/ken/ken-nipt/fastq/2603/GNCI26030003/2602250001_S13_R2_001.fastq.gz

TOTAL_END=$(date +%s)
TOTAL_WALL=$((TOTAL_END - TOTAL_START))
TOTAL_MIN=$(awk "BEGIN {printf \"%.2f\", ${TOTAL_WALL}/60}")

echo ""
echo "=========================================="
echo " Benchmark complete — total wall: ${TOTAL_MIN} min"
echo " Results: ${RESULTS}"
echo "=========================================="
cat "${RESULTS}"
