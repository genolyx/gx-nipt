#!/usr/bin/env bash
# 8-sample parallel FASTQ benchmark with --max_cpus 16
# Uses GNCI26030001..008 from ken-nipt/fastq/2603
set -euo pipefail

NEXTFLOW=/home/ken/gx-exome/nextflow
REPO=/home/ken/gx-nipt
REF_DIR=${REPO}/refs
export PATH="${REPO}/bin/nxf-docker:${PATH}"
LABCODE=cordlife
WORK_DIR=benchmark_parallel8_2603
MAX_CPUS=16
RESULTS=${REPO}/tmp/benchmark_parallel8_results.tsv
PIDS=()

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

    # flock-free append via temp file
    echo -e "${SID}\t${START_TS}\t${END_TS}\t${WALL_SEC}\t${WALL_MIN}\t${STATUS}\t${NXF_DUR}" >> "${RESULTS}.${SID}.tmp"
    mv "${RESULTS}.${SID}.tmp" "${RESULTS}.${SID}.done"
    echo "[DONE] ${SID}: wall=${WALL_MIN}min (${WALL_SEC}s) status=${STATUS} nxf=${NXF_DUR}" | tee -a "${LOGDIR}/summary.log"
  ) &
  PIDS+=($!)
  echo "[LAUNCH] ${SID} pid=$! work=${NXF_WORK}"
}

echo -e "sample\tstart\tend\twall_sec\twall_min\tstatus\tnxf_duration" > "${RESULTS}"
rm -f "${REPO}/tmp/benchmark_parallel8_"*.done 2>/dev/null || true

BATCH_START=$(date -Iseconds)
BATCH_START_SEC=$(date +%s)

echo "=========================================="
echo " 8-sample parallel benchmark"
echo " max_cpus=${MAX_CPUS}  work_dir=${WORK_DIR}"
echo " batch_start=${BATCH_START}"
echo "=========================================="

for i in $(seq 1 8); do
  SID=$(printf 'GNCI260300%02d' "${i}")
  pair=$(resolve_fastq_pair "${SID}")
  R1="${pair%%|*}"
  R2="${pair##*|}"
  run_sample_bg "${SID}" 30 "${R1}" "${R2}"
done

echo ""
echo "Waiting for ${#PIDS[@]} Nextflow runs..."
FAIL=0
for pid in "${PIDS[@]}"; do
  if ! wait "${pid}"; then FAIL=$((FAIL + 1)); fi
done

BATCH_END=$(date -Iseconds)
BATCH_END_SEC=$(date +%s)
BATCH_WALL=$((BATCH_END_SEC - BATCH_START_SEC))
BATCH_MIN=$(awk "BEGIN {printf \"%.2f\", ${BATCH_WALL}/60}")

# Merge per-sample results
for i in $(seq 1 8); do
  SID=$(printf 'GNCI260300%02d' "${i}")
  if [[ -f "${RESULTS}.${SID}.done" ]]; then
    cat "${RESULTS}.${SID}.done" >> "${RESULTS}"
  else
    echo -e "${SID}\tN/A\tN/A\tN/A\tN/A\tMISSING\tN/A" >> "${RESULTS}"
  fi
done

echo ""
echo "=========================================="
echo " Batch wall (launch → last finish): ${BATCH_MIN} min (${BATCH_WALL}s)"
echo " Failed wrappers: ${FAIL}"
echo " Results: ${RESULTS}"
echo "=========================================="
cat "${RESULTS}"
