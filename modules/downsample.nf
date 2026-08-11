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

    // Distinct names from staged inputs — Nextflow excludes input basenames
    // from output matching even if the files are deleted and recreated.
    output:
        tuple path("${sample_name}_ds_R1.fastq.gz"),
              path("${sample_name}_ds_R2.fastq.gz"),
              emit: fastq

    script:
        def r1 = fastq_pair[0]
        def r2 = fastq_pair[1]
        """
        set -euo pipefail

        # Dereference staged inputs (may be host symlinks). Never write through
        # them — shell redirects via symlink truncate the host-side source.
        cp -L "${r1}" __in_R1.fastq.gz
        cp -L "${r2}" __in_R2.fastq.gz

        python3 /opt/gx-nipt/bin/scripts/utils/downsample_fastq.py \\
            --r1 __in_R1.fastq.gz \\
            --r2 __in_R2.fastq.gz \\
            --config ${config_json} \\
            --sample ${sample_name} \\
            --outdir .

        # Publish under names that cannot collide with staged input basenames
        mv -f "${sample_name}_R1.fastq.gz" "${sample_name}_ds_R1.fastq.gz"
        mv -f "${sample_name}_R2.fastq.gz" "${sample_name}_ds_R2.fastq.gz"

        rm -f __in_R1.fastq.gz __in_R2.fastq.gz

        echo "[DOWNSAMPLE] Complete for ${sample_name}"
        """
}
