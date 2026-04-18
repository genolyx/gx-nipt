/*
 * =========================================================
 *  Module: BWA-MEM2 Alignment
 * =========================================================
 */

process BWA_ALIGN {
    tag "${sample_name}"
    label 'process_high'
    label 'nipt_docker'

    input:
        val  sample_name
        path fastq_pair    // [r1, r2]
        path config_json
        val  labcode
        val  analysisdir

    output:
        path "${sample_name}.raw.bam", emit: bam

    script:
        def r1 = fastq_pair[0]
        def r2 = fastq_pair[1]
        def bwa_threads = params.max_cpus ?: 24
        def ref_genome  = "/data/refs/genome/hg19/hg19.fa"
        """
        set -euo pipefail

        echo "[BWA] Aligning ${sample_name} ..."

        bwa-mem2 mem \\
            -t ${bwa_threads} \\
            -R "@RG\\tID:${sample_name}\\tSM:${sample_name}\\tPL:ILLUMINA\\tLB:${sample_name}" \\
            ${ref_genome} \\
            ${r1} ${r2} \\
            2> ${sample_name}.bwa.log \\
        | samtools view -bS -o ${sample_name}.raw.bam -

        echo "[BWA] Alignment complete: ${sample_name}.raw.bam"
        """
}
