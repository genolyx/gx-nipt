/**
 * gx-FF Module
 *
 * Improved Fetal Fraction estimation using gx-FF (LightGBM + DNN ensemble).
 * Runs in parallel with seqFF and produces an ensemble FF value.
 *
 * Inputs:
 *   - bam: aligned BAM file (with index)
 *   - bincount: pre-computed 50kb bin count file (optional, from HMMcopy)
 *   - model_pkl: trained gx-FF model file (.pkl)
 *   - sample_id: sample identifier
 *
 * Outputs:
 *   - gxff_tsv: TSV with FF_GXFF, FF_LGBM, FF_DNN, QC_FLAGS
 */

process GXFF_PREDICT {
    tag "${sample_id}"
    label 'process_medium'

    container 'genolyx/gx-ff:latest'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path bincount
    path model_pkl

    output:
    tuple val(sample_id), path("${sample_id}.gxff.tsv"), emit: gxff_tsv
    path "${sample_id}.gxff.log",                        emit: log

    script:
    def bincount_arg = bincount.name != 'NO_FILE' ? "--bincount ${bincount}" : "--bam ${bam}"
    """
    set -euo pipefail

    python -m gxff predict \\
        ${bincount_arg} \\
        --model ${model_pkl} \\
        --sample-id ${sample_id} \\
        --out ${sample_id}.gxff.tsv \\
        2>&1 | tee ${sample_id}.gxff.log

    # Validate output
    if [ ! -s "${sample_id}.gxff.tsv" ]; then
        echo "ERROR: gx-FF produced empty output for ${sample_id}" >&2
        exit 1
    fi
    """

    stub:
    """
    echo -e "SAMPLE_ID\tFF_GXFF\tFF_LGBM\tFF_DNN\tQC_FLAGS" > ${sample_id}.gxff.tsv
    echo -e "${sample_id}\t0.1234\t0.1210\t0.1258\tPASS"    >> ${sample_id}.gxff.tsv
    touch ${sample_id}.gxff.log
    """
}


/**
 * Ensemble FF calculation: combine seqFF and gx-FF results.
 *
 * Strategy:
 *   - If FF_GXFF QC_FLAGS == "PASS": weighted average (gxFF 0.6, seqFF 0.4)
 *   - If FF_GXFF has low_ff flag (FF < 5%): use gxFF exclusively (better at low FF)
 *   - If FF_GXFF fails QC: fall back to seqFF only
 *
 * Output columns:
 *   FF_FINAL, FF_SEQFF, FF_GXFF, FF_LGBM, FF_DNN, FF_METHOD, QC_FLAGS
 */
process GXFF_ENSEMBLE {
    tag "${sample_id}"
    label 'process_low'

    container 'genolyx/gx-nipt:latest'

    input:
    tuple val(sample_id), path(seqff_tsv), path(gxff_tsv)

    output:
    tuple val(sample_id), path("${sample_id}.ff_ensemble.tsv"), emit: ff_tsv
    path "${sample_id}.ff_ensemble.log",                        emit: log

    script:
    """
    set -euo pipefail

    python3 - <<'PYEOF'
import sys
import csv
import math
import logging

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [GXFF_ENSEMBLE] %(levelname)s %(message)s",
    handlers=[
        logging.FileHandler("${sample_id}.ff_ensemble.log"),
        logging.StreamHandler(sys.stdout),
    ],
)
logger = logging.getLogger(__name__)

sample_id = "${sample_id}"

# ── Read seqFF ────────────────────────────────────────────────────────────────
# seqff.txt is written by ff_gender_improved.py write_yff_txt() as a
# key-value file (tab-separated, first line is "value" header):
#   value
#   SeqFF<TAB>15.3968
#   status<TAB>OK
seqff_val = None
try:
    with open("${seqff_tsv}") as f:
        for line in f:
            line = line.strip()
            if not line or line == "value":
                continue
            parts = line.split("\\t")
            if len(parts) >= 2 and parts[0] == "SeqFF":
                try:
                    seqff_val = float(parts[1])
                except (ValueError, TypeError):
                    pass
                break
    logger.info("seqFF value: %s", seqff_val)
except Exception as e:
    logger.warning("Could not read seqFF txt: %s", e)

# ── Read gx-FF ────────────────────────────────────────────────────────────────
gxff_val = None
lgbm_val = None
dnn_val  = None
qc_flags = "UNKNOWN"
try:
    with open("${gxff_tsv}") as f:
        reader = csv.DictReader(f, delimiter="\\t")
        for row in reader:
            gxff_val = float(row.get("FF_GXFF", "nan"))
            lgbm_val = float(row.get("FF_LGBM", "nan"))
            dnn_val  = float(row.get("FF_DNN",  "nan"))
            qc_flags = row.get("QC_FLAGS", "UNKNOWN")
            break
    logger.info("gx-FF value: %.4f  QC: %s", gxff_val or float("nan"), qc_flags)
except Exception as e:
    logger.warning("Could not read gx-FF TSV: %s", e)

# ── Ensemble decision ─────────────────────────────────────────────────────────
gxff_ok  = gxff_val is not None and not math.isnan(gxff_val) and "FAIL" not in qc_flags.upper()
seqff_ok = seqff_val is not None and not math.isnan(seqff_val)

if gxff_ok and seqff_ok:
    is_low_ff = gxff_val < 0.05 or "low_ff" in qc_flags.lower()
    if is_low_ff:
        # At low FF, gx-FF is more accurate (trained with oversampling)
        ff_final = gxff_val
        method   = "gxff_only_low_ff"
        logger.info("Low-FF regime: using gx-FF exclusively (FF=%.4f)", ff_final)
    else:
        # Standard regime: weighted ensemble
        ff_final = 0.6 * gxff_val + 0.4 * seqff_val
        method   = "ensemble_0.6gxff_0.4seqff"
        logger.info(
            "Ensemble: gxFF=%.4f seqFF=%.4f -> final=%.4f",
            gxff_val, seqff_val, ff_final,
        )
elif gxff_ok:
    ff_final = gxff_val
    method   = "gxff_only"
    logger.warning("seqFF unavailable; using gx-FF only")
elif seqff_ok:
    ff_final = seqff_val
    method   = "seqff_only"
    logger.warning("gx-FF unavailable/failed; falling back to seqFF")
else:
    ff_final = float("nan")
    method   = "FAILED"
    logger.error("Both FF methods failed for sample %s", sample_id)

# ── Write output ──────────────────────────────────────────────────────────────
header = ["SAMPLE_ID", "FF_FINAL", "FF_SEQFF", "FF_GXFF", "FF_LGBM", "FF_DNN", "FF_METHOD", "QC_FLAGS"]
row = [
    sample_id,
    f"{ff_final:.4f}" if not math.isnan(ff_final) else "NA",
    f"{seqff_val:.4f}" if seqff_ok else "NA",
    f"{gxff_val:.4f}" if gxff_ok else "NA",
    f"{lgbm_val:.4f}" if lgbm_val is not None and not math.isnan(lgbm_val) else "NA",
    f"{dnn_val:.4f}"  if dnn_val  is not None and not math.isnan(dnn_val)  else "NA",
    method,
    qc_flags,
]
with open("${sample_id}.ff_ensemble.tsv", "w") as out:
    out.write("\\t".join(header) + "\\n")
    out.write("\\t".join(row)    + "\\n")

logger.info("FF ensemble complete: FF_FINAL=%.4f  method=%s", ff_final if not math.isnan(ff_final) else -1, method)
PYEOF
    """

    stub:
    """
    echo -e "SAMPLE_ID\tFF_FINAL\tFF_SEQFF\tFF_GXFF\tFF_LGBM\tFF_DNN\tFF_METHOD\tQC_FLAGS" > ${sample_id}.ff_ensemble.tsv
    echo -e "${sample_id}\t0.1240\t0.1220\t0.1234\t0.1210\t0.1258\tensemble_0.6gxff_0.4seqff\tPASS" >> ${sample_id}.ff_ensemble.tsv
    touch ${sample_id}.ff_ensemble.log
    """
}
