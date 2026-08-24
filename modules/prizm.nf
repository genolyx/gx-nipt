/*
 * =========================================================
 *  Module: PRIZM (Z-score based trisomy detection)
 *
 *  Improvements over ken-nipt:
 *   - orig / fetus / mom 3개 그룹 병렬 실행
 *   - FDR correction 옵션 추가
 * =========================================================
 */

process RUN_PRIZM {
    tag "${sample_name}:${group}"
    label 'process_low'
    label 'nipt_docker'

    publishDir { "${analysisdir}/${sample_name}/Output_PRIZM/${group}" }, mode: 'copy', overwrite: true

    input:
        tuple val(sample_name), val(group), path(count_10mb)
        path gender_txt
        path config_json
        val  labcode
        val  analysisdir

    output:
        tuple val(sample_name), val(group),
              path("${sample_name}.of_${group}.prizm.qc.txt"), emit: prizm_result
        path "*.png",                                          emit: plots, optional: true
        path "*.enhanced_qc.txt",                              emit: enhanced_qc, optional: true
        path "*.trisomy_detection.*",                          emit: trisomy, optional: true
        path "*.prizm_summary.txt",                            emit: summary, optional: true

    script:
        def ref_dir = "${params.ref_dir}/labs/${labcode}/PRIZM/${group}"
        """
        set -euo pipefail

        export TMPDIR="\${NXF_TASK_WORKDIR}"
        export MPLCONFIGDIR="\${NXF_TASK_WORKDIR}"

        python3 /opt/gx-nipt/bin/scripts/modules/prizm_runner.py \\
            --sample ${sample_name} \\
            --group ${group} \\
            --count-file ${count_10mb} \\
            --gender-file ${gender_txt} \\
            --config ${config_json} \\
            --ref-dir ${ref_dir} \\
            --outdir . \\
            --labcode ${labcode}

        # Required Portal review plots (heatmap + chr line plot)
        for f in \\
            ${sample_name}_${group}_chromosome_heatmap.png \\
            ${sample_name}_${group}_chromosome_line.png
        do
            if [ ! -s "\$f" ]; then
                echo "[PRIZM] ERROR: missing or empty plot: \$f" >&2
                echo "[PRIZM] FAIL_REASON=PRIZM plot failed for ${sample_name}/${group} (\$f)" >&2
                exit 1
            fi
        done

        echo "[PRIZM] ${group} complete for ${sample_name}"
        """
}
