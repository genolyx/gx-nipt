/*
 * =========================================================
 *  Module: Samtools operations
 * =========================================================
 */

process SAMTOOLS_SORT_INDEX {
    tag "${sample_name}"
    label 'process_medium'
    label 'nipt_docker'

    input:
        val  sample_name
        path bam
        val  suffix        // 'sorted', 'unique', etc.
        val  analysisdir

    output:
        path "${sample_name}.${suffix}.bam",     emit: bam
        path "${sample_name}.${suffix}.bam.bai", emit: bai

    script:
        def threads = task.cpus
        def mem     = task.memory.toGiga().intValue()
        """
        set -euo pipefail

        samtools sort \\
            -@ ${threads} \\
            -m ${mem}G \\
            -o ${sample_name}.${suffix}.bam \\
            ${bam}

        samtools index \\
            -@ ${threads} \\
            ${sample_name}.${suffix}.bam

        echo "[SAMTOOLS] Sort+Index complete: ${sample_name}.${suffix}.bam"
        """
}

process SAMTOOLS_FILTER_UNIQUE {
    tag "${sample_name}"
    label 'process_medium'
    label 'nipt_docker'

    input:
        val  sample_name
        path bam
        val  analysisdir

    output:
        path "${sample_name}.unique.bam",     emit: bam
        path "${sample_name}.unique.bam.bai", emit: bai

    script:
        def threads = task.cpus
        """
        set -euo pipefail

        # Filter: unique reads only (MQ >= 1, no secondary/supplementary)
        samtools view \\
            -@ ${threads} \\
            -b -q 1 \\
            -F 0x100 -F 0x800 \\
            ${bam} \\
        | samtools sort -@ ${threads} -o ${sample_name}.unique.bam -

        samtools index -@ ${threads} ${sample_name}.unique.bam

        echo "[SAMTOOLS] Unique filter complete: ${sample_name}.unique.bam"
        """
}

process SAMTOOLS_PROPER_PAIRED {
    tag "${sample_name}"
    label 'process_medium'
    label 'nipt_docker'

    input:
        val  sample_name
        path bam
        val  analysisdir

    output:
        path "${sample_name}.proper_paired.bam",     emit: bam
        path "${sample_name}.proper_paired.bam.bai", emit: bai

    script:
        def threads = task.cpus
        """
        set -euo pipefail

        # Filter: properly paired reads only
        samtools view \\
            -@ ${threads} \\
            -b -f 0x2 \\
            ${bam} \\
        | samtools sort -@ ${threads} -o ${sample_name}.proper_paired.bam -

        samtools index -@ ${threads} ${sample_name}.proper_paired.bam

        echo "[SAMTOOLS] Proper-paired filter complete: ${sample_name}.proper_paired.bam"
        """
}
