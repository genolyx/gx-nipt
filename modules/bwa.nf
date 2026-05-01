/*
 * =========================================================
 *  Module: BWA-MEM2 Alignment
 * =========================================================
 */

process BWA_ALIGN {
    tag "${sample_name}"
    label 'process_high'
    label 'nipt_docker'

    publishDir "${analysisdir}/${sample_name}", mode: 'copy', overwrite: true,
               pattern: "${sample_name}.sorted.bam*"

    input:
        val  sample_name
        path fastq_pair    // [r1, r2]
        path config_json
        val  labcode
        val  analysisdir

    output:
        // Output a coordinate-sorted BAM directly (same as ken-nipt: bwa | samtools sort pipe).
        // This avoids writing an intermediate raw.bam and ensures identical read ordering
        // for downstream duplicate marking.
        path "${sample_name}.sorted.bam",     emit: bam
        path "${sample_name}.sorted.bam.bai", emit: bai

    script:
        def r1 = fastq_pair[0]
        def r2 = fastq_pair[1]
        def threads     = params.max_cpus ?: 24
        def sort_threads = Math.max(1, threads.toInteger() - 1)
        def ref_genome  = "${params.ref_dir}/genomes/hg19/hg19.fa"
        """
        set -euo pipefail

        echo "[BWA] Aligning ${sample_name} ..."

        bwa-mem2 mem \\
            -t ${threads} \\
            ${ref_genome} \\
            ${r1} ${r2} \\
            2> ${sample_name}.bwa.log \\
        | samtools sort \\
            -@ ${sort_threads} \\
            -m 2G \\
            -o ${sample_name}.sorted.bam \\
            -T ${sample_name}.sort_tmp -

        samtools index -@ ${threads} ${sample_name}.sorted.bam

        echo "[BWA] Alignment + sort complete: ${sample_name}.sorted.bam"
        """
}
