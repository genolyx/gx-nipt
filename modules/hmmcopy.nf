/*
 * =========================================================
 *  Module: HMMcopy
 *  readCounter → HMMcopy.R normalization
 *
 *  Each process operates on a single (sample, group) trio
 *  entry. Parallelisation across orig / fetus / mom is driven
 *  by the caller workflow (workflows/hmmcopy.nf).
 *
 *  Reference wig tracks are read from
 *    ${params.ref_dir}/hmmcopy/hg19.{50kb,10mb}.{gc,map}.wig
 * =========================================================
 */

// ---------------------------------------------------------
//  readCounter: BAM → wig (per group, per resolution)
// ---------------------------------------------------------
process READCOUNTER {
    tag "${sample_name}:${group}:${res}"
    label 'process_low'
    label 'nipt_docker'

    publishDir "${analysisdir}/${sample_name}/Output_hmmcopy", mode: 'copy', overwrite: true

    input:
        tuple val(sample_name), val(group), path(bam), path(bai)
        val  res        // '50kb' or '10mb'
        val  analysisdir

    output:
        tuple val(sample_name), val(group),
              path("${sample_name}.of_${group}.${res}.wig"),
              emit: wig

    script:
        def window = (res == '50kb') ? 50000 : 10000000
        """
        set -euo pipefail

        readCounter \\
            --window ${window} \\
            --quality 20 \\
            --chromosome "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY" \\
            ${bam} \\
            > ${sample_name}.of_${group}.${res}.wig

        echo "[READCOUNTER ${res}] ${sample_name}/${group} complete"
        """
}

// ---------------------------------------------------------
//  HMMcopy.R normalisation (per group, per resolution)
// ---------------------------------------------------------
process HMMCOPY_R {
    tag "${sample_name}:${group}:${res}"
    label 'process_low'
    label 'nipt_docker'

    publishDir "${analysisdir}/${sample_name}/Output_hmmcopy", mode: 'copy', overwrite: true

    input:
        tuple val(sample_name), val(group), path(wig)
        val  res        // '50kb' or '10mb'
        val  analysisdir

    output:
        tuple val(sample_name), val(group),
              path("${sample_name}.of_${group}.${res}.wig.Normalization.txt"),
              emit: norm

    script:
        def gc_wig   = "${params.ref_dir}/hmmcopy/hg19.${res}.gc.wig"
        def map_wig  = "${params.ref_dir}/hmmcopy/hg19.${res}.map.wig"
        def r_script = "/opt/gx-nipt/bin/scripts/HMMcopy.R"
        """
        set -euo pipefail

        export R_TEMPDIR="\${NXF_TASK_WORKDIR}"
        export TMPDIR="\${NXF_TASK_WORKDIR}"

        Rscript ${r_script} \\
            ${wig} \\
            ${gc_wig} \\
            ${map_wig} \\
            ${res}

        echo "[HMMCOPY ${res}] ${sample_name}/${group} complete"
        """
}
