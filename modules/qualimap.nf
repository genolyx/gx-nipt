/*
 * =========================================================
 *  Module: Qualimap + QC Filter
 * =========================================================
 */

process QUALIMAP {
    tag "${sample_name}"
    label 'process_medium'
    label 'nipt_docker'

    publishDir "${analysisdir}/${sample_name}/Output_QC", mode: 'copy', overwrite: true

    input:
        val  sample_name
        path bam
        path config_json
        val  analysisdir

    output:
        path "genome_results.txt",        emit: genome_results
        path "${sample_name}.qc.txt",     emit: qc_txt
        path "qualimapReport.html",       emit: html, optional: true

    script:
        // Qualimap writes output into a subdirectory named after the BAM file:
        // e.g. GNCI26030001.proper_paired_stats/genome_results.txt
        // Strip the .bam suffix and append _stats to get the subdir name.
        def bam_base    = bam.name.replaceAll(/\.bam$/, '')
        def stats_dir   = "${bam_base}_stats"
        """
        set -euo pipefail

        # Qualimap (Java) needs a writable tmpdir; /tmp is root-owned in --user containers.
        export JAVA_TOOL_OPTIONS="-Djava.io.tmpdir=\${NXF_TASK_WORKDIR}"

        # Run Qualimap (writes to ./${stats_dir}/)
        qualimap bamqc \\
            -bam ${bam} \\
            -outdir . \\
            -outformat HTML \\
            --java-mem-size=8G \\
            -nt ${task.cpus} \\
            2> qualimap.log

        # Hoist results to workdir root so Nextflow can stage them
        cp ${stats_dir}/genome_results.txt  genome_results.txt
        cp ${stats_dir}/qualimapReport.html qualimapReport.html || true

        # Parse genome_results.txt -> sample-level QC summary
        python3 /opt/gx-nipt/bin/scripts/utils/parse_qualimap.py \\
            --input genome_results.txt \\
            --sample ${sample_name} \\
            --output ${sample_name}.qc.txt

        echo "[QUALIMAP] Complete for ${sample_name}"
        """
}

process QC_FILTER {
    tag "${sample_name}"
    label 'process_low'
    label 'nipt_docker'

    publishDir "${analysisdir}/${sample_name}/Output_QC", mode: 'copy', overwrite: true

    input:
        val  sample_name
        path qc_txt
        path config_json
        val  analysisdir

    output:
        path "${sample_name}.qc.filter.txt", emit: qc_result

    script:
        """
        set -euo pipefail

        python3 /opt/gx-nipt/bin/scripts/utils/qc_filter.py \\
            --qc-txt ${qc_txt} \\
            --config ${config_json} \\
            --sample ${sample_name} \\
            --output ${sample_name}.qc.filter.txt

        echo "[QC_FILTER] Complete for ${sample_name}"
        """
}
