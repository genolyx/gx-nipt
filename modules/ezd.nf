/*
 * =========================================================
 *  Module: EZD (Empirical Z-score Distribution)
 *
 *  Improvements over ken-nipt:
 *   - BUG FIX: f-string 오류 수정 (progress.update_step("{base_step}.5"))
 *   - Threshold 계산: 단순 중간값 → ROC/Youden's J 기반 최적 threshold
 *   - orig / fetus / mom 3개 그룹 병렬 실행
 * =========================================================
 */

process RUN_EZD {
    tag "${sample_name}:${group}"
    label 'process_low'
    label 'nipt_docker'

    publishDir "${analysisdir}/${sample_name}/Output_EZD/${group}", mode: 'copy', overwrite: true

    input:
        val  sample_name
        each group          // 'orig', 'fetus', 'mom'
        path norm_50kb
        path ff_result
        path config_json
        val  labcode
        val  analysisdir

    output:
        path "${group}_ezd_results.tsv",                    emit: ezd_result
        path "Trisomy_detect_result_${group}_with_SCA.tsv", emit: trisomy_result
        path "*.png",                                       emit: plots, optional: true

    script:
        def ref_dir = "/data/refs/${labcode}/EZD/${group}"
        """
        set -euo pipefail

        python3 /opt/gx-nipt/bin/scripts/modules/ezd_runner.py \\
            --sample ${sample_name} \\
            --group ${group} \\
            --norm-file ${norm_50kb} \\
            --ff-file ${ff_result} \\
            --config ${config_json} \\
            --ref-dir ${ref_dir} \\
            --outdir . \\
            --labcode ${labcode}

        echo "[EZD] ${group} complete for ${sample_name}"
        """
}
