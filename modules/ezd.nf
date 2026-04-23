/*
 * =========================================================
 *  Module: EZD (Empirical Z-score Distribution)
 *
 *  Drives ezd_runner.py's Nextflow adapter CLI. The adapter
 *  expects a single (sample, group) pair per invocation and
 *  produces a flat output directory consumed by Nextflow.
 *
 *  Reference files:
 *    ${params.ref_dir}/labs/<labcode>/EZD/<group>/
 *      sca_config.json
 *      <group>_thresholds_new.tsv
 *      <chromosome tables for plotting>
 * =========================================================
 */

process RUN_EZD {
    tag "${sample_name}:${group}"
    label 'process_low'
    label 'nipt_docker'

    publishDir "${analysisdir}/${sample_name}/Output_EZD/${group}", mode: 'copy', overwrite: true

    input:
        tuple val(sample_name), val(group), path(norm_50kb)
        path config_json
        val  labcode
        val  analysisdir

    output:
        tuple val(sample_name), val(group),
              path("${group}_ezd_results.tsv"),
              emit: ezd_result
        path "Trisomy_detect_result_${group}_with_SCA.tsv", emit: trisomy_result
        path "*.png",                                       emit: plots, optional: true

    script:
        """
        set -euo pipefail

        export TMPDIR="\${NXF_TASK_WORKDIR}"
        export MPLCONFIGDIR="\${NXF_TASK_WORKDIR}"

        python3 /opt/gx-nipt/bin/scripts/modules/ezd_runner.py \\
            --sample   ${sample_name} \\
            --group    ${group} \\
            --norm-file ${norm_50kb} \\
            --ref-dir  ${params.ref_dir} \\
            --labcode  ${labcode} \\
            --config   ${config_json} \\
            --outdir   .

        echo "[EZD] ${group} complete for ${sample_name}"
        """
}
