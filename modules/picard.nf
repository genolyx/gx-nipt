/*
 * =========================================================
 *  Module: Picard MarkDuplicates
 * =========================================================
 */

process PICARD_MARKDUP {
    tag "${sample_name}"
    label 'process_medium'
    label 'nipt_docker'

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
        """
        set -euo pipefail

        picard -Xmx${picard_mem} MarkDuplicates \\
            INPUT=${bam} \\
            OUTPUT=${sample_name}.dedup.bam \\
            METRICS_FILE=${sample_name}.dedup.metrics \\
            REMOVE_DUPLICATES=false \\
            ASSUME_SORTED=true \\
            VALIDATION_STRINGENCY=LENIENT \\
            2> ${sample_name}.picard.log

        samtools index ${sample_name}.dedup.bam

        echo "[PICARD] MarkDuplicates complete: ${sample_name}.dedup.bam"
        """
}
