/*
 * =========================================================
 *  Module: Picard MarkDuplicates
 * =========================================================
 */

process PICARD_MARKDUP {
    tag "${sample_name}"
    label 'process_medium'
    label 'nipt_docker'

    publishDir "${analysisdir}/${sample_name}", mode: 'copy', overwrite: true,
               pattern: "${sample_name}.dedup.*"

    input:
        val  sample_name
        path bam
        path config_json
        val  analysisdir

    output:
        path "${sample_name}.dedup.bam",     emit: bam
        path "${sample_name}.dedup.bam.bai", emit: bai
        path "${sample_name}.dedup.metrics", emit: metrics

    script:
        def picard_mem = "20G"
        def picard_jar = "/Work/NIPT/bin/picard/picard.jar"
        """
        set -euo pipefail


        java -Xmx${picard_mem} -Djava.io.tmpdir="\$(realpath .)" \\
            -jar ${picard_jar} MarkDuplicates \\
            INPUT=${bam} \\
            OUTPUT=${sample_name}.dedup.bam \\
            METRICS_FILE=${sample_name}.dedup.metrics \\
            REMOVE_DUPLICATES=true \\
            ASSUME_SORTED=true \\
            VALIDATION_STRINGENCY=LENIENT \\
            2> ${sample_name}.picard.log

        samtools index ${sample_name}.dedup.bam

        echo "[PICARD] MarkDuplicates complete: ${sample_name}.dedup.bam"
        """
}
