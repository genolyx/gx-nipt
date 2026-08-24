/**
 * gxcnv2 Module
 *
 * WisecondorX-core CNV detection engine, independent of gxcnv.
 * Reuses BED files from RUN_WCX — no duplicate convert/predict.
 *
 * Algorithm (gxcnv2_predict.py --bins-bed mode):
 *   Reads WCX BED files: _bins.bed, _segments.bed, _aberrations.bed
 *   Enhancements: MAD-based robust z-score · MAPD · log2(ratio)
 *
 * Visualisation (plot_gxcnv2.py):
 *   log2(ratio) genome-wide track · teal/amber/violet palette · KDE QC
 *
 * Data flow:
 *   RUN_WCX  →  (bins.bed, segments.bed, aberrations.bed)  →  GXCNV2_PREDICT  →  TSV
 *   GXCNV2_PREDICT  →  (bins.tsv, calls.tsv)  →  GXCNV2_PLOT  →  PNG
 *
 * WCX result = gxcnv2  (identical underlying calls, gxcnv2 visuals)
 */


// ─────────────────────────────────────────────────────────────────────────────
// Parse WCX BED files from RUN_WCX into gxcnv2 TSV
// ─────────────────────────────────────────────────────────────────────────────
process GXCNV2_PREDICT {
    tag "${sample_id}"
    label 'process_medium'
    label 'nipt_docker'

    publishDir {
        def _parts = sample_id.toString().tokenize('_')
        "${analysisdir}/${_parts[0..-2].join('_')}/gxcnv2"
    }, mode: 'copy', pattern: "*.tsv", overwrite: true

    input:
    tuple val(sample_id), path(bins_bed), path(segs_bed), path(aber_bed)   // from RUN_WCX.out.wcx_beds
    val   analysisdir

    output:
    tuple val(sample_id), path("${sample_id}_bins.tsv"),     emit: bins_tsv
    tuple val(sample_id), path("${sample_id}_segments.tsv"), emit: segments_tsv
    tuple val(sample_id), path("${sample_id}_calls.tsv"),    emit: calls_tsv
    tuple val(sample_id), path("${sample_id}_qcmetrics.tsv"),emit: qcmetrics_tsv
    path  "${sample_id}.gxcnv2_predict.log",                 emit: log

    script:
    """
    set -euo pipefail

    export TMPDIR="\${NXF_TASK_WORKDIR}"
    export MPLCONFIGDIR="\${NXF_TASK_WORKDIR}"

    # ── Skip if bins.bed is empty (empty BAM → empty beds from RUN_WCX) ─────────
    if [ ! -s "${bins_bed}" ]; then
        echo "[GXCNV2] Empty bins.bed from RUN_WCX — creating stub outputs." >&2
        for suffix in bins segments calls qcmetrics; do
            printf "#chrom\\tstart\\tend\\n" > ${sample_id}_\${suffix}.tsv
        done
        touch ${sample_id}.gxcnv2_predict.log
        exit 0
    fi

    # ── Run gxcnv2 predict in beds mode (reuse RUN_WCX output) ─────────────────
    python3 /opt/gx-nipt/bin/scripts/gxcnv2_predict.py \\
        --bins-bed        ${bins_bed} \\
        --segments-bed    ${segs_bed} \\
        --aberrations-bed ${aber_bed} \\
        -o                ${sample_id} \\
        2>&1 | tee ${sample_id}.gxcnv2_predict.log

    for f in ${sample_id}_bins.tsv ${sample_id}_calls.tsv; do
        if [ ! -f "\$f" ]; then
            echo "ERROR: gxcnv2 predict missing output: \$f" >&2
            exit 1
        fi
    done
    """

    stub:
    """
    for suffix in bins segments calls qcmetrics; do
        echo -e "#chrom\\tstart\\tend" > ${sample_id}_\${suffix}.tsv
    done
    touch ${sample_id}.gxcnv2_predict.log
    """
}


// ─────────────────────────────────────────────────────────────────────────────
// ─────────────────────────────────────────────────────────────────────────────
// Plot — genome-wide + per-chromosome log2(ratio) figures
// ─────────────────────────────────────────────────────────────────────────────
process GXCNV2_PLOT {
    tag "${sample_id}"
    label 'process_low'
    label 'nipt_docker'

    publishDir {
        def _p = sample_id.toString().tokenize('_')
        "${analysisdir}/${_p[0..-2].join('_')}/gxcnv2"
    }, mode: 'copy', overwrite: true, pattern: "*.png"

    input:
    tuple val(sample_id), path(bins_tsv), path(calls_tsv)
    val   analysisdir

    output:
    tuple val(sample_id), path("${sample_id}_genome.png"), emit: genome_png
    tuple val(sample_id), path("${sample_id}_qc.png"),     emit: qc_png
    tuple val(sample_id), path("${sample_id}_chr*.png"),   emit: chr_pngs, optional: true

    script:
    def _p    = sample_id.toString().tokenize('_')
    def base_id = _p[0..-2].join('_')
    """
    set -euo pipefail

    export MPLCONFIGDIR="\${NXF_TASK_WORKDIR}"

    python3 /opt/gx-nipt/bin/scripts/plot_gxcnv2.py \\
        --bins  ${bins_tsv} \\
        --calls ${calls_tsv} \\
        -o      ${sample_id}

    for f in ${sample_id}_genome.png ${sample_id}_qc.png; do
        if [ ! -s "\$f" ]; then
            echo "[GXCNV2_PLOT] ERROR: missing or empty plot: \$f" >&2
            echo "[GXCNV2_PLOT] FAIL_REASON=gxcnv2 plot failed for ${sample_id} (\$f)" >&2
            exit 1
        fi
    done
    """

    stub:
    """
    touch ${sample_id}_genome.png ${sample_id}_qc.png
    """
}
