#!/usr/bin/env bash
# =========================================================
#  gx-nipt: Daemon-facing pipeline launcher
# =========================================================
#  gx-daemon (nipt plugin) calls this script. It owns the
#  host-side directory layout contract, then delegates the
#  actual analysis to the Nextflow pipeline (main.nf).
#
#  Contract with gx-daemon (same style as ken-nipt + sgnipt):
#
#    $ROOT_DIR/
#      ├── fastq/<work_dir>/<order_id>/{R1,R2}.fastq.gz    # inputs
#      ├── analysis/<work_dir>/<order_id>/...              # intermediates
#      ├── log/<work_dir>/<order_id>/pipeline.log          # logs
#      ├── output/<work_dir>/<order_id>/                   # daemon reads this
#      │    ├── <order_id>.json                            # required (result)
#      │    ├── <order_id>.output.tar                      # required (archive)
#      │    ├── <order_id>_progress.txt                    # progress/status
#      │    ├── <order_id>.completed | .failed
#      │    └── Output_*/                                  # detailed artefacts
#      └── config/<labcode>/pipeline_config.json           # lab config
#
#  SSD scratch (optional):
#    --use-ssd --scratch-dir /fast/nvme
#      * Nextflow workDir and intermediate BAMs live on SSD
#      * Only <order_id>.proper_paired.bam + the final output tree
#        end up on persistent storage
#      * Per-sample SSD subtree is removed on success/failure
#
# =========================================================

set -euo pipefail

# -----------------------------------------------------------
# Defaults
# -----------------------------------------------------------
SAMPLE_NAME=""
ORDER_ID=""
WORK_DIR=""
ROOT_DIR=""
FASTQ_R1=""
FASTQ_R2=""
LABCODE=""
AGE=""
FROM_BAM=""

# Optional NIPT algorithm knobs
ALGORITHM_ONLY="false"
FORCE="false"
FRESH="false"
GXFF_MODEL=""
GXCNV_REFERENCE=""
RUN_WCX="true"

# SSD scratch
USE_SSD="false"
SCRATCH_DIR=""
SSD_MAX_USAGE_GB="200"

# Reference data root (bind-mounted into the Docker container at the same path)
REF_DIR=""

# Repo (where main.nf lives) - resolvable relative to this script
REPO_DIR="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )/.." &>/dev/null && pwd )"

# Nextflow binary (override via env or --nextflow)
NEXTFLOW_BIN="${NEXTFLOW_BIN:-nextflow}"

usage() {
    cat <<'EOF'
gx-nipt daemon launcher

Required:
  --order-id <id>            Order id (used as sample name by default)
  --work-dir <name>          Batch work directory (e.g. 250430_01)
  --root-dir <path>          Host root (e.g. /data/gx-nipt or /home/ken/gx-nipt-data)
  --labcode  <code>          Lab code (cordlife, ucl, vn, ...)
  --fastq-r1 <name|path>     R1 file (basename under <root>/fastq/<work>/<order>/ or abs path)
  --fastq-r2 <name|path>     R2 file (basename under <root>/fastq/<work>/<order>/ or abs path)
  --age      <int>           Maternal age

Optional:
  --sample-name <name>       Sample name (defaults to <order-id>)
  --from-bam <path>          Start from existing proper_paired BAM
  --algorithm-only           Skip alignment; reuse existing BAM
  --force                    Re-run even if <order_id>.completed marker exists
  --fresh                    Clear any .nextflow/ resume cache for this sample
  --gxff-model <.pkl>        Enable gx-FF (LightGBM + DNN) ensemble
                             (deferred: no default model yet; leave unset
                             to fall back to seqFF-only FF estimation)
  --gxcnv-reference <.npz>   Enable gx-cnv parallel CNV track
  --gxcnv-model <.npz|.pkl>  Alias for --gxcnv-reference (deferred: no
                             default reference yet; leave unset to fall
                             back to WisecondorX-only CNV calling)
  --no-wcx                   Disable WisecondorX (once gx-cnv is validated)
  --use-ssd                  Enable SSD scratch (Strategy B)
  --scratch-dir <path>       SSD mount point (default: /tmp/nipt_scratch)
  --ssd-max-usage-gb <n>     Abort if SSD usage exceeds this (default: 200)
  --ref-dir <path>           Reference data root on host (default: /data/reference)
                             Layout expected:
                               <ref-dir>/genomes/hg19/hg19.fa(+BWA index)
                               <ref-dir>/hmmcopy/hg19.{50kb,10mb}.{gc,map}.wig
                               <ref-dir>/models/seqff_model.pkl
                               <ref-dir>/labs/<labcode>/{WC,WCX,EZD,PRIZM}/<group>/...
                               <ref-dir>/labs/<labcode>/bed/...
  --nextflow <bin>           Path to nextflow binary

Environment overrides:
  NEXTFLOW_BIN=/usr/local/bin/nextflow
  NXF_OPTS="-Xms1g -Xmx4g"
EOF
    exit "${1:-0}"
}

# -----------------------------------------------------------
# Argument parsing (GNU-style long options)
# -----------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        --order-id)          ORDER_ID="$2"; shift 2 ;;
        --sample-name)       SAMPLE_NAME="$2"; shift 2 ;;
        --work-dir)          WORK_DIR="$2"; shift 2 ;;
        --root-dir)          ROOT_DIR="$2"; shift 2 ;;
        --labcode)           LABCODE="$2"; shift 2 ;;
        --fastq-r1)          FASTQ_R1="$2"; shift 2 ;;
        --fastq-r2)          FASTQ_R2="$2"; shift 2 ;;
        --age)               AGE="$2"; shift 2 ;;
        --from-bam)          FROM_BAM="$2"; shift 2 ;;
        --algorithm-only)    ALGORITHM_ONLY="true"; shift ;;
        --force|-f)          FORCE="true"; shift ;;
        --fresh)             FRESH="true"; shift ;;
        --gxff-model)        GXFF_MODEL="$2"; shift 2 ;;
        --gxcnv-reference|--gxcnv-model)
                             GXCNV_REFERENCE="$2"; shift 2 ;;
        --no-wcx)            RUN_WCX="false"; shift ;;
        --use-ssd)           USE_SSD="true"; shift ;;
        --scratch-dir)       SCRATCH_DIR="$2"; shift 2 ;;
        --ssd-max-usage-gb)  SSD_MAX_USAGE_GB="$2"; shift 2 ;;
        --ref-dir)           REF_DIR="$2"; shift 2 ;;
        --nextflow)          NEXTFLOW_BIN="$2"; shift 2 ;;
        -h|--help)           usage 0 ;;
        *)
            echo "[run_nipt] Unknown option: $1" >&2
            usage 1
            ;;
    esac
done

# Defaults derived from arguments
SAMPLE_NAME="${SAMPLE_NAME:-$ORDER_ID}"
SCRATCH_DIR="${SCRATCH_DIR:-/tmp/nipt_scratch}"
# Allow env overrides (NIPT_REF_DIR from gx-daemon) before falling back
# to the nextflow.config default.
REF_DIR="${REF_DIR:-${NIPT_REF_DIR:-}}"

# -----------------------------------------------------------
# Validate required args
# -----------------------------------------------------------
missing=()
[[ -z "$ORDER_ID" ]] && missing+=("--order-id")
[[ -z "$WORK_DIR" ]] && missing+=("--work-dir")
[[ -z "$ROOT_DIR" ]] && missing+=("--root-dir")
[[ -z "$LABCODE"  ]] && missing+=("--labcode")
if [[ -z "$FROM_BAM" && "$ALGORITHM_ONLY" == "false" ]]; then
    [[ -z "$FASTQ_R1" ]] && missing+=("--fastq-r1")
    [[ -z "$FASTQ_R2" ]] && missing+=("--fastq-r2")
    [[ -z "$AGE"      ]] && missing+=("--age")
fi
if (( ${#missing[@]} > 0 )); then
    echo "[run_nipt] Missing required options: ${missing[*]}" >&2
    usage 1
fi

# -----------------------------------------------------------
# Resolved paths (mirrors ken-nipt / gx-daemon layout)
# -----------------------------------------------------------
HOST_FASTQ_DIR="${ROOT_DIR}/fastq/${WORK_DIR}/${ORDER_ID}"
HOST_ANALYSIS_DIR="${ROOT_DIR}/analysis/${WORK_DIR}/${ORDER_ID}"
HOST_OUTPUT_DIR="${ROOT_DIR}/output/${WORK_DIR}/${ORDER_ID}"
HOST_LOG_DIR="${ROOT_DIR}/log/${WORK_DIR}/${ORDER_ID}"
HOST_CONFIG_DIR="${ROOT_DIR}/config"
HOST_LAB_CONFIG="${HOST_CONFIG_DIR}/${LABCODE}/pipeline_config.json"

PROGRESS_FILE="${HOST_OUTPUT_DIR}/${ORDER_ID}_progress.txt"
COMPLETED_FILE="${HOST_OUTPUT_DIR}/${ORDER_ID}.completed"
FAILED_FILE="${HOST_OUTPUT_DIR}/${ORDER_ID}.failed"
PIPELINE_LOG="${HOST_LOG_DIR}/pipeline.log"

mkdir -p "$HOST_FASTQ_DIR" "$HOST_ANALYSIS_DIR" "$HOST_OUTPUT_DIR" "$HOST_LOG_DIR" \
         "${HOST_CONFIG_DIR}/${LABCODE}"

# Tee all output to pipeline.log while still echoing to daemon's stdout
exec > >(tee -a "$PIPELINE_LOG") 2>&1

_ts() { date '+%Y-%m-%d %H:%M:%S%z'; }

progress() {
    local stage="$1"
    local pct="${2:-}"
    local msg="${3:-}"
    printf '[%s] %-12s %s%s\n' "$(_ts)" "$stage" "${pct:+(${pct}%) }" "$msg" \
        | tee -a "$PROGRESS_FILE" > /dev/null
}

progress "START" "0" "order=${ORDER_ID} work=${WORK_DIR} labcode=${LABCODE}"

# -----------------------------------------------------------
# Completion-marker guard (fast-return if already done)
# -----------------------------------------------------------
if [[ "$FORCE" == "false" && -f "$COMPLETED_FILE" ]]; then
    echo "[run_nipt] Already completed: $COMPLETED_FILE (use --force to re-run)"
    progress "COMPLETED" "100" "already-completed marker present"
    exit 0
fi

# Clear failure marker from any previous attempt
rm -f "$FAILED_FILE"

# -----------------------------------------------------------
# Ensure lab config is visible under <root>/config/<lab>/
# main.nf reads from exactly that path. Copy from repo's
# conf/labs/<lab>/pipeline_config.json when missing so the
# daemon can run without manual setup.
# -----------------------------------------------------------
if [[ ! -f "$HOST_LAB_CONFIG" ]]; then
    REPO_LAB_CONFIG="${REPO_DIR}/conf/labs/${LABCODE}/pipeline_config.json"
    if [[ -f "$REPO_LAB_CONFIG" ]]; then
        cp -f "$REPO_LAB_CONFIG" "$HOST_LAB_CONFIG"
        echo "[run_nipt] Seeded lab config: $HOST_LAB_CONFIG (from $REPO_LAB_CONFIG)"
    else
        echo "[run_nipt] ERROR: lab config not found: $HOST_LAB_CONFIG" >&2
        echo "[run_nipt]        (also missing in repo: $REPO_LAB_CONFIG)" >&2
        echo "FAILED: lab config missing (${LABCODE})" > "$FAILED_FILE"
        exit 2
    fi
fi

# -----------------------------------------------------------
# FASTQ resolution
#   - Daemon usually puts absolute paths in --fastq-r1/r2.
#   - main.nf expects basenames under <root>/fastq/<work>/<sample>/
#     so we symlink when needed.
# -----------------------------------------------------------
link_into_fastq_dir() {
    local src="$1"
    if [[ -z "$src" ]]; then return 0; fi
    if [[ "$src" != /* ]]; then
        # Already a bare filename — assume it lives in HOST_FASTQ_DIR
        if [[ ! -f "${HOST_FASTQ_DIR}/${src}" ]]; then
            echo "[run_nipt] ERROR: FASTQ missing: ${HOST_FASTQ_DIR}/${src}" >&2
            return 1
        fi
        return 0
    fi
    if [[ ! -f "$src" ]]; then
        echo "[run_nipt] ERROR: FASTQ not found at absolute path: $src" >&2
        return 1
    fi
    local base
    base="$(basename "$src")"
    local dest="${HOST_FASTQ_DIR}/${base}"
    if [[ ! -e "$dest" ]]; then
        ln -sf "$src" "$dest"
        echo "[run_nipt] Linked $src -> $dest"
    fi
}

if [[ -z "$FROM_BAM" ]]; then
    link_into_fastq_dir "$FASTQ_R1" || { echo "FAILED: R1 missing" > "$FAILED_FILE"; exit 2; }
    link_into_fastq_dir "$FASTQ_R2" || { echo "FAILED: R2 missing" > "$FAILED_FILE"; exit 2; }
    FASTQ_R1_NAME="$(basename "$FASTQ_R1")"
    FASTQ_R2_NAME="$(basename "$FASTQ_R2")"
else
    FASTQ_R1_NAME=""
    FASTQ_R2_NAME=""
fi

progress "PREPARED" "5" "inputs ready under ${HOST_FASTQ_DIR}"

# -----------------------------------------------------------
# Fresh run: purge Nextflow resume state for this sample
# -----------------------------------------------------------
if [[ "$FRESH" == "true" ]]; then
    rm -rf "${REPO_DIR}/.nextflow" "${REPO_DIR}/work" 2>/dev/null || true
    echo "[run_nipt] Fresh: cleared .nextflow/ and work/"
fi

# -----------------------------------------------------------
# Build Nextflow command
# -----------------------------------------------------------
NF_ARGS=(
    run "${REPO_DIR}/main.nf"
    --sample_name "$SAMPLE_NAME"
    --labcode     "$LABCODE"
    --root_dir    "$ROOT_DIR"
    --work_dir    "$WORK_DIR"
    --outdir      "$HOST_OUTPUT_DIR"
    --analysisdir "$HOST_ANALYSIS_DIR"
)

if [[ -n "$AGE" ]];       then NF_ARGS+=( --age "$AGE" ); fi
if [[ -n "$FROM_BAM" ]];  then NF_ARGS+=( --from_bam "$FROM_BAM" ); fi
if [[ -n "$FASTQ_R1_NAME" ]]; then NF_ARGS+=( --fastq_r1 "$FASTQ_R1_NAME" --fastq_r2 "$FASTQ_R2_NAME" ); fi
if [[ "$ALGORITHM_ONLY" == "true" ]]; then NF_ARGS+=( --algorithm_only true ); fi
if [[ "$FORCE" == "true" ]];          then NF_ARGS+=( --force true ); fi
if [[ -n "$GXFF_MODEL" ]];       then NF_ARGS+=( --gxff_model "$GXFF_MODEL" ); fi
if [[ -n "$GXCNV_REFERENCE" ]];  then NF_ARGS+=( --gxcnv_reference "$GXCNV_REFERENCE" ); fi
if [[ "$RUN_WCX" == "false" ]];  then NF_ARGS+=( --run_wcx false ); fi
if [[ -n "$REF_DIR" ]];          then NF_ARGS+=( --ref_dir "$REF_DIR" ); fi

# -----------------------------------------------------------
# Reference-directory preflight
# -----------------------------------------------------------
# Use the CLI override if provided, else inherit nextflow.config default.
EFFECTIVE_REF_DIR="${REF_DIR:-/data/reference}"
REQUIRED_REF_PATHS=(
    "${EFFECTIVE_REF_DIR}/genomes/hg19/hg19.fa"
    "${EFFECTIVE_REF_DIR}/hmmcopy/hg19.50kb.gc.wig"
    "${EFFECTIVE_REF_DIR}/hmmcopy/hg19.50kb.map.wig"
    "${EFFECTIVE_REF_DIR}/labs/${LABCODE}"
)
missing_refs=()
for p in "${REQUIRED_REF_PATHS[@]}"; do
    [[ -e "$p" ]] || missing_refs+=( "$p" )
done
if (( ${#missing_refs[@]} > 0 )); then
    echo "[run_nipt] ERROR: required reference files/dirs are missing:" >&2
    for p in "${missing_refs[@]}"; do echo "         - $p" >&2; done
    echo "[run_nipt] Hint: see docs/reference_generation.md for the expected layout," >&2
    echo "           or override the root with --ref-dir (or NIPT_REF_DIR env)." >&2
    {
        echo "Pipeline aborted: missing reference data under ${EFFECTIVE_REF_DIR}"
        printf '  - %s\n' "${missing_refs[@]}"
    } > "$FAILED_FILE"
    progress "FAILED" "" "missing reference data under ${EFFECTIVE_REF_DIR}"
    exit 2
fi

if [[ "$USE_SSD" == "true" ]]; then
    NF_ARGS=( -profile ssd_scratch "${NF_ARGS[@]}" )
    NF_ARGS+=( --scratch_dir "$SCRATCH_DIR" --ssd_max_usage_gb "$SSD_MAX_USAGE_GB" )
fi

# Nextflow house-keeping files go under a per-run tracedir inside log/
TRACE_DIR="${HOST_LOG_DIR}/pipeline_info"
mkdir -p "$TRACE_DIR"
NF_ARGS+=( --tracedir "$TRACE_DIR" )

progress "RUN" "10" "nextflow run main.nf"
echo "[run_nipt] cwd=$REPO_DIR"
echo "[run_nipt] cmd=${NEXTFLOW_BIN} ${NF_ARGS[*]}"

pushd "$REPO_DIR" >/dev/null

# -----------------------------------------------------------
# Execute Nextflow
# -----------------------------------------------------------
set +e
"$NEXTFLOW_BIN" "${NF_ARGS[@]}"
NF_EXIT=$?
set -e

popd >/dev/null

if [[ "$NF_EXIT" -ne 0 ]]; then
    echo "[run_nipt] Nextflow failed with exit code $NF_EXIT" >&2
    {
        echo "Pipeline failed for ${ORDER_ID} at $(_ts)"
        echo "Reason: nextflow exit ${NF_EXIT}"
    } >> "$FAILED_FILE"
    progress "FAILED" "" "nextflow exit ${NF_EXIT}"
    exit "$NF_EXIT"
fi

progress "PIPELINE_DONE" "90" "nextflow completed"

# -----------------------------------------------------------
# Package daemon-expected artefacts
#
# gx-daemon looks for:
#   output/<work>/<order>/<order>.json         (result summary)
#   output/<work>/<order>/<order>.output.tar   (all Output_*)
#   output/<work>/<order>/<order>_progress.txt (we keep writing to it)
# -----------------------------------------------------------
RESULT_JSON_SRC="${HOST_ANALYSIS_DIR}/${SAMPLE_NAME}/Output_Result/${SAMPLE_NAME}.result.json"
if [[ ! -f "$RESULT_JSON_SRC" ]]; then
    # REPORT_WORKFLOW also publishes the JSON under output/Output_Result/
    RESULT_JSON_SRC="${HOST_OUTPUT_DIR}/Output_Result/${SAMPLE_NAME}.result.json"
fi

if [[ ! -f "$RESULT_JSON_SRC" ]]; then
    echo "[run_nipt] ERROR: report JSON not produced: $RESULT_JSON_SRC" >&2
    {
        echo "Pipeline completed but report JSON missing."
        echo "Expected: $RESULT_JSON_SRC"
    } >> "$FAILED_FILE"
    progress "FAILED" "" "missing report JSON"
    exit 3
fi

cp -f "$RESULT_JSON_SRC" "${HOST_OUTPUT_DIR}/${ORDER_ID}.json"
echo "[run_nipt] Copied result JSON -> ${HOST_OUTPUT_DIR}/${ORDER_ID}.json"

# Build tar (deterministic relative paths; skip if Output_* dirs missing)
TAR_FILE="${HOST_OUTPUT_DIR}/${ORDER_ID}.output.tar"
(
    cd "$HOST_OUTPUT_DIR"
    # Collect Output_* subdirs + top-level JSON/HTML
    TAR_ITEMS=()
    while IFS= read -r -d '' dir; do
        TAR_ITEMS+=( "$(basename "$dir")" )
    done < <(find . -maxdepth 1 -type d -name 'Output_*' -print0 2>/dev/null || true)
    # Always include the JSON; include review HTML if present.
    [[ -f "${ORDER_ID}.json" ]] && TAR_ITEMS+=( "${ORDER_ID}.json" )
    [[ -f "Output_Result/${SAMPLE_NAME}.review.html" ]] && true  # already under Output_Result/
    if (( ${#TAR_ITEMS[@]} == 0 )); then
        echo "[run_nipt] WARNING: no Output_* dirs to archive" >&2
    else
        tar -cf "$TAR_FILE" "${TAR_ITEMS[@]}"
        echo "[run_nipt] Built archive: $TAR_FILE ($(du -h "$TAR_FILE" | cut -f1))"
    fi
)

# -----------------------------------------------------------
# Success markers
# -----------------------------------------------------------
echo "$(_ts)" > "$COMPLETED_FILE"
progress "COMPLETED" "100" "artefacts ready for daemon pickup"
echo "[run_nipt] SUCCESS: ${ORDER_ID}"
exit 0
