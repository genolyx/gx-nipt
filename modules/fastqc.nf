/*
 * =========================================================
 *  Module: FastQC
 * =========================================================
 */

process FASTQC {
    tag "${sample_name}"
    label 'process_low'
    label 'nipt_docker'

    publishDir "${analysisdir}/${sample_name}/Output_QC", mode: 'copy', overwrite: true

    input:
        val  sample_name
        path fastq_pair
        val  analysisdir

    output:
        path "*.html", emit: html
        path "*.zip",  emit: zip

    script:
        def r1 = fastq_pair[0]
        def r2 = fastq_pair[1]
        """
        set -euo pipefail

        # javax.imageio uses java.io.tmpdir; redirect to writable workdir.
        export JAVA_TOOL_OPTIONS="-Djava.io.tmpdir=\${NXF_TASK_WORKDIR}"

        fastqc -d . -o . -t 2 ${r1} ${r2}

        echo "[FASTQC] Complete for ${sample_name}"
        """
}
