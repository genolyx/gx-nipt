#!/usr/bin/env bash
# 8-sample benchmark — max 4 concurrent Nextflow runs (slot pool)
# When one sample finishes, launch the next so active count stays <= MAX_SLOTS.
# max_cpus=16 per sample (gx default throughput profile)
set -euo pipefail

NEXTFLOW=/home/ken/gx-exome/nextflow
REPO=/home/ken/gx-nipt
REF_DIR=${REPO}/refs
export PATH="${REPO}/bin/nxf-docker:${PATH}"
LABCODE=cordlife
WORK_DIR=benchmark_slot4_8_2603
MAX_SLOTS=4
MAX_CPUS=16
RESULTS=${REPO}/tmp/benchmark_slot4_8_results.tsv

mkdir -p "${REPO}/tmp" "${REPO}/log/${WORK_DIR}"

resolve_fastq_pair() {
  local sid=$1
  local dir="/home/ken/ken-nipt/fastq/2603/${sid}"
  local r1 r2
  r1=$(find "${dir}" -maxdepth 1 -name '*R1*.fastq.gz' | head -1)
  r2=$(find "${dir}" -maxdepth 1 -name '*R2*.fastq.gz' | head -1)
  if [[ -z "${r1}" || -z "${r2}" ]]; then
    echo "ERROR: missing FASTQ pair for ${sid} in ${dir}" >&2
    exit 1
  fi
  echo "${r1}|${r2}"
}

run_sample_bg() {
  local SID=$1
  local AGE=$2
  local R1=$3
  local R2=$4
  local OUTDIR=${REPO}/output/${WORK_DIR}/${SID}
  local ANALYSISDIR=${REPO}/analysis/${WORK_DIR}/${SID}
  local LOGDIR=${REPO}/log/${WORK_DIR}/${SID}
  local NXF_WORK=${REPO}/work-${WORK_DIR}-${SID}
  mkdir -p "${OUTDIR}" "${ANALYSISDIR}" "${LOGDIR}" "${NXF_WORK}"

  (
    local START_TS START_SEC END_TS END_SEC WALL_SEC WALL_MIN EXIT_CODE STATUS NXF_DUR
    START_TS=$(date -Iseconds)
    START_SEC=$(date +%s)

    set +e
    ${NEXTFLOW} run ${REPO}/main.nf \
      -work-dir "${NXF_WORK}" \
      --sample_name  "${SID}" \
      --labcode      "${LABCODE}" \
      --root_dir     "${REPO}" \
      --work_dir     "${WORK_DIR}" \
      --outdir       "${OUTDIR}" \
      --analysisdir  "${ANALYSISDIR}" \
      --fastq_r1     "${R1}" \
      --fastq_r2     "${R2}" \
      --age          "${AGE}" \
      --max_cpus     "${MAX_CPUS}" \
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

    echo -e "${SID}\t${START_TS}\t${END_TS}\t${WALL_SEC}\t${WALL_MIN}\t${STATUS}\t${NXF_DUR}" >> "${RESULTS}.${SID}.tmp"
    mv "${RESULTS}.${SID}.tmp" "${RESULTS}.${SID}.done"
    echo "[DONE] ${SID}: wall=${WALL_MIN}min (${WALL_SEC}s) status=${STATUS} nxf=${NXF_DUR}" | tee -a "${REPO}/log/${WORK_DIR}/summary.log"
  ) &
}

launch_sample() {
  local SID=$1
  local pair
  pair=$(resolve_fastq_pair "${SID}")
  run_sample_bg "${SID}" 30 "${pair%%|*}" "${pair##*|}"
  echo "[LAUNCH] ${SID} pid=$! slot_pool active=$((RUNNING + 1))/${MAX_SLOTS}"
}

echo -e "sample\tstart\tend\twall_sec\twall_min\tstatus\tnxf_duration" > "${RESULTS}"
rm -f "${REPO}/tmp/benchmark_slot4_8_"*.done 2>/dev/null || true

SAMPLES=(
  GNCI26030001 GNCI26030002 GNCI26030003 GNCI26030004
  GNCI26030005 GNCI26030006 GNCI26030007 GNCI26030008
)
TOTAL=${#SAMPLES[@]}
NEXT_IDX=0
RUNNING=0
FAIL=0

BATCH_START=$(date -Iseconds)
BATCH_START_SEC=$(date +%s)

echo "=========================================="
echo " 8-sample slot-pool benchmark"
echo " max_slots=${MAX_SLOTS}  max_cpus=${MAX_CPUS}  work_dir=${WORK_DIR}"
echo " batch_start=${BATCH_START}"
echo "=========================================="

while (( NEXT_IDX < TOTAL || RUNNING > 0 )); do
  while (( RUNNING < MAX_SLOTS && NEXT_IDX < TOTAL )); do
    launch_sample "${SAMPLES[NEXT_IDX]}"
    NEXT_IDX=$((NEXT_IDX + 1))
    RUNNING=$((RUNNING + 1))
  done
  if (( RUNNING > 0 )); then
    if ! wait -n; then FAIL=$((FAIL + 1)); fi
    RUNNING=$((RUNNING - 1))
    echo "[SLOT] freed — running=${RUNNING} queued=$((TOTAL - NEXT_IDX))"
  fi
done

BATCH_END_SEC=$(date +%s)
BATCH_WALL=$((BATCH_END_SEC - BATCH_START_SEC))
BATCH_MIN=$(awk "BEGIN {printf \"%.2f\", ${BATCH_WALL}/60}")

for SID in "${SAMPLES[@]}"; do
  if [[ -f "${RESULTS}.${SID}.done" ]]; then
    cat "${RESULTS}.${SID}.done" >> "${RESULTS}"
  else
    echo -e "${SID}\tN/A\tN/A\tN/A\tN/A\tMISSING\tN/A" >> "${RESULTS}"
  fi
done

echo ""
echo "=========================================="
echo " Batch wall (first launch → last finish): ${BATCH_MIN} min (${BATCH_WALL}s)"
echo " Failed wrappers: ${FAIL}"
echo " Results: ${RESULTS}"
echo "=========================================="
cat "${RESULTS}"
