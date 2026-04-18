/*
 * =========================================================
 *  Module: Microdeletion Detection
 * =========================================================
 */

process RUN_MD_DETECTION {
    tag "${sample_name}"
    label 'process_low'
    label 'nipt_docker'

    publishDir "${analysisdir}/${sample_name}/Output_MD", mode: 'copy', overwrite: true

    input:
        val  sample_name
        path wc_results   // collected WC + WCX result files
        path config_json
        val  labcode
        val  analysisdir

    output:
        path "${sample_name}.md_result.tsv", emit: md_result

    script:
        def bed_dir = "/data/refs/${labcode}/bed"
        """
        set -euo pipefail

        python3 /opt/gx-nipt/bin/scripts/modules/process_md_result.py \\
            --sample ${sample_name} \\
            --labcode ${labcode} \\
            --config ${config_json} \\
            --analysis-dir ${analysisdir} \\
            --bed-dir ${bed_dir} \\
            --output ${sample_name}.md_result.tsv

        echo "[MD] Microdeletion detection complete for ${sample_name}"
        """
}
