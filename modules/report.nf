/*
 * =========================================================
 *  Module: Report Generation & Output Copy
 *
 *  Portal-compatible output directory structure:
 *    output/{work_dir}/{sample_name}/
 *      ├── Output_QC/
 *      ├── Output_FF/
 *      ├── Output_EZD/{orig,fetus,mom}/
 *      ├── Output_PRIZM/{orig,fetus,mom}/
 *      ├── Output_WC/{orig,fetus,mom}/
 *      ├── Output_WCX/{orig,fetus,mom}/
 *      ├── Output_hmmcopy/
 *      ├── Output_MD/
 *      └── Output_Result/
 *           ├── {sample_name}.result.json
 *           └── {sample_name}.review.html
 * =========================================================
 */

process GENERATE_JSON {
    tag "${sample_name}"
    label 'process_low'
    label 'nipt_docker'

    publishDir "${analysisdir}/${sample_name}/Output_Result", mode: 'copy', overwrite: true

    input:
        val  sample_name
        path all_results   // collected upstream outputs
        path config_json
        val  labcode
        val  analysisdir

    output:
        path "${sample_name}.result.json", emit: json_file

    script:
        """
        set -euo pipefail

        python3 /opt/gx-nipt/bin/scripts/modules/generate_json_output.py \\
            --sample ${sample_name} \\
            --labcode ${labcode} \\
            --config ${config_json} \\
            --analysis-dir ${analysisdir} \\
            --output ${sample_name}.result.json

        echo "[REPORT] JSON generated for ${sample_name}"
        """
}

process GENERATE_HTML {
    tag "${sample_name}"
    label 'process_low'
    label 'nipt_docker'

    publishDir "${analysisdir}/${sample_name}/Output_Result", mode: 'copy', overwrite: true

    input:
        val  sample_name
        path json_file
        path config_json
        val  labcode
        val  analysisdir

    output:
        path "${sample_name}.review.html", emit: html_file

    script:
        """
        set -euo pipefail

        python3 /opt/gx-nipt/bin/scripts/modules/html_review_page.py \\
            --sample ${sample_name} \\
            --json ${json_file} \\
            --config ${config_json} \\
            --labcode ${labcode} \\
            --analysis-dir ${analysisdir} \\
            --output ${sample_name}.review.html

        echo "[REPORT] HTML generated for ${sample_name}"
        """
}

process COPY_TO_OUTPUT {
    tag "${sample_name}"
    label 'process_low'

    input:
        val  sample_name
        path json_file
        path html_file
        val  analysisdir
        val  outdir

    output:
        path ".copy_done", emit: done

    script:
        """
        set -euo pipefail

        # Create portal-compatible output directory structure
        mkdir -p ${outdir}/Output_QC
        mkdir -p ${outdir}/Output_FF
        mkdir -p ${outdir}/Output_EZD/orig
        mkdir -p ${outdir}/Output_EZD/fetus
        mkdir -p ${outdir}/Output_EZD/mom
        mkdir -p ${outdir}/Output_PRIZM/orig
        mkdir -p ${outdir}/Output_PRIZM/fetus
        mkdir -p ${outdir}/Output_PRIZM/mom
        mkdir -p ${outdir}/Output_WC/orig
        mkdir -p ${outdir}/Output_WC/fetus
        mkdir -p ${outdir}/Output_WC/mom
        mkdir -p ${outdir}/Output_WCX/orig
        mkdir -p ${outdir}/Output_WCX/fetus
        mkdir -p ${outdir}/Output_WCX/mom
        mkdir -p ${outdir}/Output_hmmcopy
        mkdir -p ${outdir}/Output_MD
        mkdir -p ${outdir}/Output_Result

        # Copy analysis results to output directory
        rsync -a --ignore-missing-args \\
            ${analysisdir}/${sample_name}/Output_QC/ \\
            ${outdir}/Output_QC/

        rsync -a --ignore-missing-args \\
            ${analysisdir}/${sample_name}/Output_FF/ \\
            ${outdir}/Output_FF/

        rsync -a --ignore-missing-args \\
            ${analysisdir}/${sample_name}/Output_EZD/ \\
            ${outdir}/Output_EZD/

        rsync -a --ignore-missing-args \\
            ${analysisdir}/${sample_name}/Output_PRIZM/ \\
            ${outdir}/Output_PRIZM/

        rsync -a --ignore-missing-args \\
            ${analysisdir}/${sample_name}/Output_WC/ \\
            ${outdir}/Output_WC/

        rsync -a --ignore-missing-args \\
            ${analysisdir}/${sample_name}/Output_WCX/ \\
            ${outdir}/Output_WCX/

        rsync -a --ignore-missing-args \\
            ${analysisdir}/${sample_name}/Output_hmmcopy/ \\
            ${outdir}/Output_hmmcopy/

        rsync -a --ignore-missing-args \\
            ${analysisdir}/${sample_name}/Output_MD/ \\
            ${outdir}/Output_MD/

        # Copy final report
        cp ${json_file} ${outdir}/Output_Result/
        cp ${html_file} ${outdir}/Output_Result/

        # Write completion marker (portal compatibility)
        echo "\$(date -Iseconds)" > ${outdir}/${sample_name}.completed

        touch .copy_done
        echo "[COPY] Output copied to ${outdir}"
        """
}
