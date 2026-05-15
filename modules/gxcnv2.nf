/**
 * gxcnv2 Module
 *
 * WisecondorX-core CNV detection engine, independent of gxcnv.
 * Uses the same WCX reference.npz files already present for RUN_WCX,
 * so no additional reference building is required.
 *
 * Algorithm (gxcnv2_predict.py):
 *   Coverage normalise → PCA correction → Within-sample normalise → CBS
 *   Enhancements: MAD-based robust z-score · MAPD metric · log2(ratio) output
 *
 * Visualisation (plot_gxcnv2.py):
 *   log2(ratio) filled-area genome-wide track · teal/amber/violet palette
 *   KDE QC plot · CBS segment overlay (distinct from gxcnv and WCX visuals)
 *
 * Commands:
 *   WisecondorX convert  BAM → sample NPZ        (GXCNV2_CONVERT)
 *   gxcnv2_predict.py    sample NPZ + ref → TSV  (GXCNV2_PREDICT)
 *   plot_gxcnv2.py       TSV → PNG               (GXCNV2_PLOT)
 *
 * Reference layout (identical to WCX, no new files needed):
 *   {ref_dir}/labs/{labcode}/WCX/orig_{M|F}_200k_proper_paired.npz
 *   {ref_dir}/labs/{labcode}/WCX/fetus_{M|F}_200k_of.npz
 *   {ref_dir}/labs/{labcode}/WCX/mom_200k_of.npz
 */


// ─────────────────────────────────────────────────────────────────────────────
// Step 1: BAM → sample NPZ  (same WisecondorX convert as RUN_WCX)
// ─────────────────────────────────────────────────────────────────────────────
process GXCNV2_CONVERT {
    tag "${sample_id}"
    label 'process_low'
    label 'nipt_docker'

    input:
    tuple val(sample_id), path(bam), path(bai)
    val   binsize   // must match the reference (200000)

    output:
    tuple val(sample_id), path("${sample_id}.gxcnv2.npz"), emit: npz
    path  "${sample_id}.gxcnv2_convert.log",               emit: log

    script:
    """
    set -euo pipefail

    export TMPDIR="\${NXF_TASK_WORKDIR}"
    export MPLCONFIGDIR="\${NXF_TASK_WORKDIR}"

    READ_COUNT=\$(samtools view -c ${bam})
    if [ "\${READ_COUNT}" -eq 0 ]; then
        echo "[GXCNV2] ${sample_id} BAM is empty — skipping convert." >&2
        touch ${sample_id}.gxcnv2.npz
        touch ${sample_id}.gxcnv2_convert.log
        exit 0
    fi

    # WisecondorX convert appends .npz automatically; pass prefix only
    WisecondorX convert \\
        ${bam} \\
        ${sample_id}.gxcnv2 \\
        --binsize ${binsize} \\
        2>&1 | tee ${sample_id}.gxcnv2_convert.log

    if [ ! -s "${sample_id}.gxcnv2.npz" ]; then
        echo "ERROR: WisecondorX convert produced empty NPZ for ${sample_id}" >&2
        exit 1
    fi
    """

    stub:
    """
    python3 -c "import numpy as np; np.savez_compressed('${sample_id}.gxcnv2', **{'1': np.zeros(100, dtype='int32')})"
    touch ${sample_id}.gxcnv2_convert.log
    """
}


// ─────────────────────────────────────────────────────────────────────────────
// Step 2: Predict CNVs using gxcnv2_predict.py + WCX reference
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
    tuple val(sample_id), path(sample_npz)
    path  gender_txt     // <sample>.gender.txt — final_gender line
    val   labcode
    val   zscore         // aberration z-score threshold (e.g. 3.5)
    val   analysisdir

    output:
    tuple val(sample_id), path("${sample_id}_bins.tsv"),     emit: bins_tsv
    tuple val(sample_id), path("${sample_id}_segments.tsv"), emit: segments_tsv
    tuple val(sample_id), path("${sample_id}_calls.tsv"),    emit: calls_tsv
    tuple val(sample_id), path("${sample_id}_qcmetrics.tsv"),emit: qcmetrics_tsv
    path  "${sample_id}.gxcnv2_predict.log",                 emit: log

    script:
    def _parts   = sample_id.toString().tokenize('_')
    def grp      = _parts[-1]   // orig / fetus / mom
    def base_sid = _parts[0..-2].join('_')
    """
    set -euo pipefail

    export TMPDIR="\${NXF_TASK_WORKDIR}"
    export MPLCONFIGDIR="\${NXF_TASK_WORKDIR}"
    export R_TEMPDIR="\${NXF_TASK_WORKDIR}"

    # ── Resolve gender ──────────────────────────────────────────────────────────
    GENDER=\$(awk -F'[\\t ]+' 'tolower(\$1) == "final_gender" {
        g = toupper(\$2)
        if (g == "XY" || g == "MALE"   || g == "M") { print "M"; exit }
        if (g == "XX" || g == "FEMALE" || g == "F") { print "F"; exit }
    }' ${gender_txt})

    if [ -z "\${GENDER}" ]; then
        echo "[GXCNV2] Could not parse final_gender; defaulting to F." >&2
        GENDER="F"
    fi
    echo "[GXCNV2] group=${grp}  gender=\${GENDER}"

    # ── Resolve WCX reference (same files as RUN_WCX) ──────────────────────────
    # Naming mirrors modules/wisecondor.nf RUN_WCX exactly.
    if [ "${grp}" = "mom" ]; then
        REF_NPZ="${params.ref_dir}/labs/${labcode}/WCX/mom_200k_of.npz"
    elif [ "\${GENDER}" = "M" ]; then
        case "${grp}" in
            orig)  REF_NPZ="${params.ref_dir}/labs/${labcode}/WCX/orig_M_200k_proper_paired.npz" ;;
            fetus) REF_NPZ="${params.ref_dir}/labs/${labcode}/WCX/fetus_M_200k_of.npz" ;;
        esac
    else
        case "${grp}" in
            orig)  REF_NPZ="${params.ref_dir}/labs/${labcode}/WCX/orig_F_200k_proper_paired.npz" ;;
            fetus) REF_NPZ="${params.ref_dir}/labs/${labcode}/WCX/fetus_F_200k_of.npz" ;;
        esac
    fi
    echo "[GXCNV2] ref=\${REF_NPZ}"

    # ── Skip if NPZ is empty (empty BAM from convert step) ─────────────────────
    if [ ! -s "${sample_npz}" ]; then
        echo "[GXCNV2] Empty sample NPZ — creating stub outputs." >&2
        for suffix in bins segments calls qcmetrics; do
            printf "#chrom\\tstart\\tend\\n" > ${sample_id}_\${suffix}.tsv
        done
        touch ${sample_id}.gxcnv2_predict.log
        exit 0
    fi

    # ── Skip if reference not found ─────────────────────────────────────────────
    if [ ! -f "\${REF_NPZ}" ]; then
        echo "[GXCNV2] WARNING: reference not found: \${REF_NPZ}" >&2
        echo "[GXCNV2] Skipping — creating stub outputs." >&2
        for suffix in bins segments calls qcmetrics; do
            printf "##gxcnv2_skip=true\\n#chrom\\tstart\\tend\\n" > ${sample_id}_\${suffix}.tsv
        done
        touch ${sample_id}.gxcnv2_predict.log
        exit 0
    fi

    # ── Run gxcnv2 predict ──────────────────────────────────────────────────────
    python3 /opt/gx-nipt/bin/scripts/gxcnv2_predict.py \\
        ${sample_npz} \\
        "\${REF_NPZ}" \\
        -o ${sample_id} \\
        --gender \${GENDER} \\
        --zscore ${zscore} \\
        --alpha   0.01 \\
        --seed    100 \\
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
// Step 3: Plot — genome-wide + per-chromosome log2(ratio) figures
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
    """

    stub:
    """
    touch ${sample_id}_genome.png ${sample_id}_qc.png
    """
}
