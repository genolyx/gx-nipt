/**
 * gxcnv1 Module
 *
 * Wisecondor-core CNV detection engine, independent of gxcnv and gxcnv2.
 * Uses the WC reference.npz files already used by RUN_WC, and reuses
 * RUN_WC's own out.npz directly — no duplicate test step needed.
 *
 * Algorithm (gxcnv1_predict.py):
 *   Reads wisecondor test out.npz: PCA normalisation + Stouffer z-score (already computed)
 *   Enhancements: MAD-based robust z-score · MAPD metric · log2(ratio) output
 *
 * Visualisation (plot_gxcnv1.py):
 *   log2(ratio) filled-area genome-wide track · teal/amber/violet palette
 *   KDE QC plot · call segment overlay  (same visual style as gxcnv2)
 *
 * Data flow:
 *   RUN_WC  →  (sample.npz, out.npz)  →  GXCNV1_PREDICT  →  TSV
 *   GXCNV1_PREDICT  →  (bins.tsv, calls.tsv)  →  GXCNV1_PLOT  →  PNG
 *
 * WC result = gxcnv1  (identical underlying calls, gxcnv2-style visuals)
 */


// ─────────────────────────────────────────────────────────────────────────────
// Step 1: BAM → sample NPZ  (same Wisecondor convert as RUN_WC)
// ─────────────────────────────────────────────────────────────────────────────
process GXCNV1_CONVERT {
    tag "${sample_id}"
    label 'process_low'
    label 'nipt_docker'

    input:
    tuple val(sample_id), path(bam), path(bai)
    val   binsize   // must match the reference (200000)

    output:
    tuple val(sample_id), path("${sample_id}.gxcnv1.npz"), emit: npz
    path  "${sample_id}.gxcnv1_convert.log",               emit: log

    script:
    """
    set -euo pipefail

    export LOGNAME="\${LOGNAME:-nipt}"
    export MPLCONFIGDIR="\${NXF_TASK_WORKDIR}"

    READ_COUNT=\$(samtools view -c ${bam})
    if [ "\${READ_COUNT}" -eq 0 ]; then
        echo "[GXCNV1] ${sample_id} BAM is empty — skipping convert." >&2
        touch ${sample_id}.gxcnv1.npz
        touch ${sample_id}.gxcnv1_convert.log
        exit 0
    fi

    # wisecondor.py uses Python2; -binsize uses a single dash
    python2 /opt/wisecondor/wisecondor.py convert \\
        ${bam} \\
        ${sample_id}.gxcnv1.npz \\
        -binsize ${binsize} \\
        2>&1 | tee ${sample_id}.gxcnv1_convert.log

    if [ ! -s "${sample_id}.gxcnv1.npz" ]; then
        echo "ERROR: Wisecondor convert produced empty NPZ for ${sample_id}" >&2
        exit 1
    fi
    """

    stub:
    """
    python3 -c "import numpy as np; np.savez_compressed('${sample_id}.gxcnv1', sample={})"
    touch ${sample_id}.gxcnv1_convert.log
    """
}


// ─────────────────────────────────────────────────────────────────────────────
// Step 2: Parse out.npz from RUN_WC into gxcnv2-compatible TSV
//         Receives both sample.npz and out.npz produced by RUN_WC,
//         so wisecondor.py test is NOT run again here.
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
// Step 3: Plot — genome-wide + per-chromosome log2(ratio) figures
//         Same visual style as gxcnv2 (filled-area, teal/amber/violet, KDE QC)
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
    """

    stub:
    """
    touch ${sample_id}_genome.png ${sample_id}_qc.png
    """
}
