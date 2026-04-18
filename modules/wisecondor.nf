/*
 * =========================================================
 *  Module: Wisecondor & WisecondorX
 * =========================================================
 */

process RUN_WC {
    tag "${sample_name}:${group}"
    label 'process_low'
    label 'nipt_docker'

    publishDir "${analysisdir}/${sample_name}/Output_WC/${group}", mode: 'copy', overwrite: true

    input:
        val  sample_name
        each group
        path bam
        path ff_result
        path config_json
        val  labcode
        val  analysisdir

    output:
        path "${sample_name}.wc.${group}.report.txt", emit: wc_result
        path "*.npz",                                 emit: npz, optional: true

    script:
        def ref_npz = "/data/refs/${labcode}/WC/${group}/${group}_200k.npz"
        def binsize  = 200000
        """
        set -euo pipefail

        # Convert BAM to npz
        WisecondorX convert \\
            ${bam} \\
            ${sample_name}.wc.${group}.npz \\
            --binsize ${binsize}

        # Run WisecondorX predict
        WisecondorX predict \\
            ${sample_name}.wc.${group}.npz \\
            ${ref_npz} \\
            ${sample_name}.wc.${group} \\
            --zscore 6.0 \\
            --bed

        # Rename output for compatibility
        mv ${sample_name}.wc.${group}_bins.bed ${sample_name}.wc.${group}.report.txt 2>/dev/null || true

        echo "[WC] ${group} complete for ${sample_name}"
        """
}

process RUN_WCX {
    tag "${sample_name}:${group}"
    label 'process_low'
    label 'nipt_docker'

    publishDir "${analysisdir}/${sample_name}/Output_WCX/${group}", mode: 'copy', overwrite: true

    input:
        val  sample_name
        each group
        path bam
        path ff_result
        path config_json
        val  labcode
        val  analysisdir

    output:
        path "${sample_name}.wcx.${group}.report.txt", emit: wcx_result
        path "*.npz",                                  emit: npz, optional: true

    script:
        def binsize = 200000
        """
        set -euo pipefail

        # Read gender from ff_result to select gender-specific reference
        GENDER=\$(python3 -c "
import sys
with open('${ff_result}') as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) >= 2 and parts[0] == 'gd_2':
            gender = parts[2] if len(parts) > 2 else parts[1]
            print('male' if 'XY' in gender else 'female')
            sys.exit(0)
print('female')
")

        if [ "\$GENDER" = "male" ]; then
            REF_NPZ="/data/refs/${labcode}/WCX/${group}/M_200k.npz"
        else
            REF_NPZ="/data/refs/${labcode}/WCX/${group}/F_200k.npz"
        fi

        # Convert BAM to npz
        WisecondorX convert \\
            ${bam} \\
            ${sample_name}.wcx.${group}.npz \\
            --binsize ${binsize}

        # Run WisecondorX predict with gender-aware reference
        WisecondorX predict \\
            ${sample_name}.wcx.${group}.npz \\
            \${REF_NPZ} \\
            ${sample_name}.wcx.${group} \\
            --zscore 6.0 \\
            --bed

        mv ${sample_name}.wcx.${group}_bins.bed ${sample_name}.wcx.${group}.report.txt 2>/dev/null || true

        echo "[WCX] ${group} complete for ${sample_name}"
        """
}
