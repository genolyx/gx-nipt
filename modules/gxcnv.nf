/**
 * gx-cnv Module
 *
 * Hybrid dual-track CNV detection engine for sWGS NIPT data.
 * Runs in parallel with WisecondorX. Results are compared and
 * both are included in the final report until gx-cnv is validated.
 *
 * Algorithm:
 *   Track A: Whole-genome Z-score normalisation + CBS segmentation (WisecondorX-inspired)
 *   Track B: Regional PCA denoising + Laplace-smoothed Mahalanobis distance (BinDel-inspired)
 *   Decision: AND-gate — HIGH_RISK only when both tracks independently exceed thresholds
 *
 * Commands:
 *   gxcnv convert  BAM → NPZ (GC-corrected bin counts)
 *   gxcnv newref   NPZ files → reference panel
 *   gxcnv predict  sample NPZ + reference → CNV calls
 */

// ─────────────────────────────────────────────────────────────────────────────
// Step 1: Convert BAM → NPZ
// ─────────────────────────────────────────────────────────────────────────────
process GXCNV_CONVERT {
    tag "${sample_id}"
    label 'process_medium'
    label 'nipt_docker'

    input:
    tuple val(sample_id), path(bam), path(bai)
    val   bin_size
    path  blacklist_bed

    output:
    tuple val(sample_id), path("${sample_id}.gxcnv.npz"), emit: npz
    path  "${sample_id}.gxcnv_convert.log",               emit: log

    script:
    def blacklist_arg = blacklist_bed.name != 'NO_FILE' ? "--blacklist ${blacklist_bed}" : ""
    """
    set -euo pipefail

    gxcnv convert \\
        ${bam} \\
        ${sample_id}.gxcnv.npz \\
        --bin-size ${bin_size} \\
        --min-mapq 1 \\
        ${blacklist_arg} \\
        2>&1 | tee ${sample_id}.gxcnv_convert.log

    if [ ! -s "${sample_id}.gxcnv.npz" ]; then
        echo "ERROR: gxcnv convert produced empty NPZ for ${sample_id}" >&2
        exit 1
    fi
    """

    stub:
    """
    python3 -c "import numpy as np; np.savez_compressed('${sample_id}.gxcnv.npz', bins=np.zeros(30000))"
    touch ${sample_id}.gxcnv_convert.log
    """
}


// ─────────────────────────────────────────────────────────────────────────────
// Step 2: Build reference panel (run once, not per-sample)
// ─────────────────────────────────────────────────────────────────────────────
process GXCNV_NEWREF {
    tag "gxcnv_newref"
    label 'process_high'
    label 'nipt_docker'

    input:
    path npz_files   // collection of normal sample NPZ files
    val  output_name
    val  pca_variance

    output:
    path "${output_name}.npz", emit: reference
    path "${output_name}.newref.log", emit: log

    script:
    """
    set -euo pipefail

    gxcnv newref \\
        ${npz_files.join(' ')} \\
        -o ${output_name}.npz \\
        --pca-variance ${pca_variance} \\
        2>&1 | tee ${output_name}.newref.log

    if [ ! -s "${output_name}.npz" ]; then
        echo "ERROR: gxcnv newref produced empty reference" >&2
        exit 1
    fi
    """

    stub:
    """
    python3 -c "import numpy as np; np.savez_compressed('${output_name}.npz', ref=np.zeros((50,30000)))"
    touch ${output_name}.newref.log
    """
}


// ─────────────────────────────────────────────────────────────────────────────
// Step 3: Predict CNVs (per-sample)
// ─────────────────────────────────────────────────────────────────────────────
process GXCNV_PREDICT {
    tag "${sample_id}"
    label 'process_medium'
    label 'nipt_docker'

    publishDir "${analysisdir}/${sample_id}/gxcnv", mode: 'copy', pattern: "*.tsv", overwrite: true
    publishDir "${analysisdir}/${sample_id}/gxcnv", mode: 'copy', pattern: "*.txt", overwrite: true

    input:
    tuple val(sample_id), path(sample_npz)
    path  reference_npz
    val   thresh_z
    val   thresh_p
    val   fetal_fraction   // from FF ensemble step; "NA" if unavailable
    val   analysisdir

    output:
    tuple val(sample_id), path("${sample_id}_bins.tsv"),     emit: bins_tsv
    tuple val(sample_id), path("${sample_id}_segments.tsv"), emit: segments_tsv
    tuple val(sample_id), path("${sample_id}_calls.tsv"),    emit: calls_tsv
    tuple val(sample_id), path("${sample_id}_regions.tsv"),  emit: regions_tsv
    tuple val(sample_id), path("${sample_id}_qcmetrics.tsv"),emit: qcmetrics_tsv
    tuple val(sample_id), path("${sample_id}_sex.txt"),      emit: sex_txt
    path  "${sample_id}.gxcnv_predict.log",                  emit: log

    script:
    def ff_arg = fetal_fraction != 'NA' ? "--fetal-fraction ${fetal_fraction}" : ""
    """
    set -euo pipefail

    gxcnv predict \\
        ${sample_npz} \\
        ${reference_npz} \\
        -o ${sample_id} \\
        --thresh-z ${thresh_z} \\
        --thresh-p ${thresh_p} \\
        ${ff_arg} \\
        2>&1 | tee ${sample_id}.gxcnv_predict.log

    # Validate required outputs
    for f in ${sample_id}_bins.tsv ${sample_id}_calls.tsv ${sample_id}_regions.tsv; do
        if [ ! -f "\$f" ]; then
            echo "ERROR: gxcnv predict missing output: \$f" >&2
            exit 1
        fi
    done
    """

    stub:
    """
    for suffix in bins segments calls regions qcmetrics; do
        echo -e "##gxcnv_version=0.1.0\\n#chrom\tstart\tend\tdual_call" > ${sample_id}_\${suffix}.tsv
    done
    echo -e "##gxcnv_version=0.1.0\\n#metric\tvalue\\npredicted_sex\tF" > ${sample_id}_sex.txt
    touch ${sample_id}.gxcnv_predict.log
    """
}


// ─────────────────────────────────────────────────────────────────────────────
// Step 4: Compare gx-cnv vs WisecondorX results
//         Produces a concordance summary for validation purposes.
// ─────────────────────────────────────────────────────────────────────────────
process GXCNV_COMPARE {
    tag "${sample_id}"
    label 'process_low'
    label 'nipt_docker'

    publishDir "${analysisdir}/${sample_id}/gxcnv", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id),
          path(gxcnv_calls),
          path(gxcnv_regions),
          path(wcx_aberrations)   // WisecondorX aberrations.bed or equivalent
    val   analysisdir

    output:
    tuple val(sample_id), path("${sample_id}.cnv_comparison.tsv"), emit: comparison
    path  "${sample_id}.cnv_comparison.log",                       emit: log

    script:
    """
    set -euo pipefail

    python3 - <<'PYEOF'
import sys
import csv
import logging

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [GXCNV_COMPARE] %(levelname)s %(message)s",
    handlers=[
        logging.FileHandler("${sample_id}.cnv_comparison.log"),
        logging.StreamHandler(sys.stdout),
    ],
)
logger = logging.getLogger(__name__)

sample_id = "${sample_id}"

# ── Parse gx-cnv HIGH_RISK calls ─────────────────────────────────────────────
gxcnv_calls = set()
try:
    with open("${gxcnv_calls}") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\\t")
            if len(parts) >= 4:
                region_name = parts[3] if len(parts) > 3 else f"{parts[0]}:{parts[1]}-{parts[2]}"
                gxcnv_calls.add(region_name)
    logger.info("gx-cnv HIGH_RISK calls: %s", gxcnv_calls)
except Exception as e:
    logger.warning("Could not parse gx-cnv calls: %s", e)

# ── Parse gx-cnv all regions ─────────────────────────────────────────────────
gxcnv_regions = {}
try:
    with open("${gxcnv_regions}") as f:
        reader = None
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                cols = line.lstrip("#").strip().split("\\t")
                reader = cols
                continue
            if reader is None:
                continue
            parts = line.strip().split("\\t")
            row = dict(zip(reader, parts))
            rname = row.get("region_name", "")
            if rname:
                gxcnv_regions[rname] = {
                    "gxcnv_dual_call": row.get("dual_call", ""),
                    "gxcnv_track_a_z": row.get("track_a_mean_z", "NA"),
                    "gxcnv_track_b_p": row.get("track_b_pvalue", "NA"),
                    "gxcnv_risk_pct":  row.get("risk_pct", "NA"),
                }
except Exception as e:
    logger.warning("Could not parse gx-cnv regions: %s", e)

# ── Parse WisecondorX aberrations ─────────────────────────────────────────────
wcx_calls = set()
try:
    with open("${wcx_aberrations}") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\\t")
            if len(parts) >= 4:
                wcx_calls.add(parts[3])
    logger.info("WisecondorX calls: %s", wcx_calls)
except Exception as e:
    logger.warning("Could not parse WisecondorX aberrations: %s", e)

# ── Concordance analysis ──────────────────────────────────────────────────────
all_regions = set(gxcnv_regions.keys()) | gxcnv_calls | wcx_calls

rows = []
for region in sorted(all_regions):
    gxcnv_call = "HIGH_RISK" if region in gxcnv_calls else "LOW_RISK"
    wcx_call   = "ABERRANT"  if region in wcx_calls   else "NORMAL"
    concordant = (gxcnv_call == "HIGH_RISK") == (wcx_call == "ABERRANT")
    status     = "CONCORDANT" if concordant else "DISCORDANT"

    region_info = gxcnv_regions.get(region, {})
    rows.append({
        "sample_id":       sample_id,
        "region":          region,
        "gxcnv_call":      gxcnv_call,
        "wcx_call":        wcx_call,
        "concordance":     status,
        "gxcnv_track_a_z": region_info.get("gxcnv_track_a_z", "NA"),
        "gxcnv_track_b_p": region_info.get("gxcnv_track_b_p", "NA"),
        "gxcnv_risk_pct":  region_info.get("gxcnv_risk_pct",  "NA"),
    })

# Summary stats
n_total       = len(rows)
n_concordant  = sum(1 for r in rows if r["concordance"] == "CONCORDANT")
n_gxcnv_only  = sum(1 for r in rows if r["gxcnv_call"] == "HIGH_RISK" and r["wcx_call"] == "NORMAL")
n_wcx_only    = sum(1 for r in rows if r["gxcnv_call"] == "LOW_RISK"  and r["wcx_call"] == "ABERRANT")

logger.info(
    "Concordance summary: total=%d concordant=%d gxcnv_only=%d wcx_only=%d",
    n_total, n_concordant, n_gxcnv_only, n_wcx_only,
)

# ── Write output ──────────────────────────────────────────────────────────────
fieldnames = ["sample_id", "region", "gxcnv_call", "wcx_call", "concordance",
              "gxcnv_track_a_z", "gxcnv_track_b_p", "gxcnv_risk_pct"]
with open("${sample_id}.cnv_comparison.tsv", "w", newline="") as out:
    out.write(f"## sample={sample_id}\\n")
    out.write(f"## n_regions={n_total}\\n")
    out.write(f"## n_concordant={n_concordant}\\n")
    out.write(f"## n_gxcnv_only={n_gxcnv_only}\\n")
    out.write(f"## n_wcx_only={n_wcx_only}\\n")
    writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\\t")
    writer.writeheader()
    writer.writerows(rows)

logger.info("Comparison written to ${sample_id}.cnv_comparison.tsv")
PYEOF
    """

    stub:
    """
    echo -e "## sample=${sample_id}\\n#sample_id\tregion\tgxcnv_call\twcx_call\tconcordance" > ${sample_id}.cnv_comparison.tsv
    touch ${sample_id}.cnv_comparison.log
    """
}
