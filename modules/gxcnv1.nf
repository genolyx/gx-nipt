/**
 * gxcnv1 Module
 *
 * Wisecondor-core CNV detection engine, independent of gxcnv and gxcnv2.
 * Reuses RUN_WC sample.npz + out.npz — no duplicate convert/test.
 *
 * Algorithm (gxcnv1_predict.py):
 *   Reads wisecondor test out.npz: PCA normalisation + Stouffer z-score
 *   Enhancements: MAD-based robust z-score · MAPD · log2(ratio)
 *
 * Visualisation (plot_gxcnv1.py):
 *   log2(ratio) genome-wide track · teal/amber/violet palette · KDE QC
 *
 * Data flow:
 *   RUN_WC  →  (sample.npz, out.npz)  →  GXCNV1_PREDICT  →  TSV
 *   GXCNV1_PREDICT  →  (bins.tsv, calls.tsv)  →  GXCNV1_PLOT  →  PNG
 *
 * WC result = gxcnv1  (identical underlying calls, gxcnv1 visuals)
 */


// ─────────────────────────────────────────────────────────────────────────────
// Parse out.npz from RUN_WC into gxcnv2-compatible TSV
// ─────────────────────────────────────────────────────────────────────────────
process GXCNV1_PREDICT {
    tag "${sample_id}"
    label 'process_medium'
    label 'nipt_docker'

    publishDir {
        def _parts = sample_id.toString().tokenize('_')
        "${analysisdir}/${_parts[0..-2].join('_')}/gxcnv1"
    }, mode: 'copy', pattern: "*.tsv", overwrite: true

    input:
    tuple val(sample_id), path(sample_npz), path(out_npz)   // both from RUN_WC.out
    val   analysisdir

    output:
    tuple val(sample_id), path("${sample_id}_bins.tsv"),     emit: bins_tsv
    tuple val(sample_id), path("${sample_id}_segments.tsv"), emit: segments_tsv
    tuple val(sample_id), path("${sample_id}_calls.tsv"),    emit: calls_tsv
    tuple val(sample_id), path("${sample_id}_qcmetrics.tsv"),emit: qcmetrics_tsv
    path  "${sample_id}.gxcnv1_predict.log",                 emit: log

    script:
    """
    set -euo pipefail

    export LOGNAME="\${LOGNAME:-nipt}"
    export MPLCONFIGDIR="\${NXF_TASK_WORKDIR}"

    # ── Skip if NPZ is empty (empty BAM → empty out.npz from RUN_WC) ───────────
    if [ ! -s "${sample_npz}" ] || [ ! -s "${out_npz}" ]; then
        echo "[GXCNV1] Empty NPZ from RUN_WC — creating stub outputs." >&2
        for suffix in bins segments calls qcmetrics; do
            printf "#chrom\\tstart\\tend\\n" > ${sample_id}_\${suffix}.tsv
        done
        touch ${sample_id}.gxcnv1_predict.log
        exit 0
    fi

    # ── Parse out.npz → gxcnv2-compatible TSV (Python3) ─────────────────────
    # out.npz was already produced by RUN_WC (wisecondor test), not re-run here.
    python3 /opt/gx-nipt/bin/scripts/gxcnv1_predict.py \\
        ${sample_npz} \\
        ${out_npz} \\
        -o ${sample_id} \\
        2>&1 | tee ${sample_id}.gxcnv1_predict.log

    for f in ${sample_id}_bins.tsv ${sample_id}_calls.tsv; do
        if [ ! -f "\$f" ]; then
            echo "ERROR: gxcnv1 predict missing output: \$f" >&2
            exit 1
        fi
    done
    """

    stub:
    """
    for suffix in bins segments calls qcmetrics; do
        echo -e "#chrom\\tstart\\tend" > ${sample_id}_\${suffix}.tsv
    done
    touch ${sample_id}.gxcnv1_predict.log
    """
}


// ─────────────────────────────────────────────────────────────────────────────
// Plot — genome-wide + QC figures (same visual style as gxcnv2)
// ─────────────────────────────────────────────────────────────────────────────
process GXCNV1_PLOT {
    tag "${sample_id}"
    label 'process_low'
    label 'nipt_docker'

    publishDir {
        def _p = sample_id.toString().tokenize('_')
        "${analysisdir}/${_p[0..-2].join('_')}/gxcnv1"
    }, mode: 'copy', overwrite: true, pattern: "*.png"

    input:
    tuple val(sample_id), path(bins_tsv), path(calls_tsv)
    val   analysisdir

    output:
    tuple val(sample_id), path("${sample_id}_genome.png"), emit: genome_png
    tuple val(sample_id), path("${sample_id}_qc.png"),     emit: qc_png
    tuple val(sample_id), path("${sample_id}_chr*.png"),   emit: chr_pngs, optional: true

    script:
    def cyto_file = "${params.ref_dir}/bed/common/cytoBand.txt"
    """
    set -euo pipefail

    export MPLCONFIGDIR="\${NXF_TASK_WORKDIR}"

    python3 /opt/gx-nipt/bin/scripts/plot_gxcnv1.py \\
        --bins     ${bins_tsv} \\
        --calls    ${calls_tsv} \\
        --cytoband ${cyto_file} \\
        -o         ${sample_id}

    for f in ${sample_id}_genome.png ${sample_id}_qc.png; do
        if [ ! -s "\$f" ]; then
            echo "[GXCNV1_PLOT] ERROR: missing or empty plot: \$f" >&2
            echo "[GXCNV1_PLOT] FAIL_REASON=gxcnv1 plot failed for ${sample_id} (\$f)" >&2
            exit 1
        fi
    done
    """

    stub:
    """
    touch ${sample_id}_genome.png ${sample_id}_qc.png
    """
}
