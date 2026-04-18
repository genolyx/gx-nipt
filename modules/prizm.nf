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

    publishDir "${analysisdir}/${sample_name}/Output_PRIZM/${group}", mode: 'copy', overwrite: true

    input:
        val  sample_name
        each group           // 'orig', 'fetus', 'mom'
        path count_10mb
        path ff_result
        path config_json
        val  labcode
        val  analysisdir

    output:
        path "${sample_name}.of_${group}.prizm.qc.txt",  emit: prizm_result
        path "*.png",                                     emit: plots, optional: true

    script:
        def ref_dir = "/data/refs/${labcode}/PRIZM/${group}"
        """
        set -euo pipefail

        python3 /opt/gx-nipt/bin/scripts/modules/prizm_runner.py \\
            --sample ${sample_name} \\
            --group ${group} \\
            --count-file ${count_10mb} \\
            --ff-file ${ff_result} \\
            --config ${config_json} \\
            --ref-dir ${ref_dir} \\
            --outdir . \\
            --labcode ${labcode}

        echo "[PRIZM] ${group} complete for ${sample_name}"
        """
}
