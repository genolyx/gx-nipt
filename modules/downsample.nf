/*
 * =========================================================
 *  Module: FASTQ Downsampling
 *  If FASTQ size > max_fq_size, downsample to downsample_size reads
 * =========================================================
 */

process DOWNSAMPLE_FASTQ {
    tag "${sample_name}"
    label 'process_low'
    label 'nipt_docker'

    input:
        val  sample_name
        path fastq_pair    // [r1, r2]
        path config_json

    output:
        path "*.fastq.gz", emit: fastq

    script:
        def r1 = fastq_pair[0]
        def r2 = fastq_pair[1]
        """
        set -euo pipefail

        python3 /opt/gx-nipt/bin/scripts/utils/downsample_fastq.py \\
            --r1 ${r1} \\
            --r2 ${r2} \\
            --config ${config_json} \\
            --sample ${sample_name} \\
            --outdir .

        echo "[DOWNSAMPLE] Complete for ${sample_name}"
        """
}
