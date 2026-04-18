/*
 * =========================================================
 *  Module: Fetal Fraction & Gender Detection
 *
 *  Improvements over ken-nipt:
 *   - YFF1: Y-chromosome coverage ratio (gd_1) - BUG FIX: use config BED paths
 *   - YFF2: Adjusted YFF via wig normalization (gd_2)
 *   - SeqFF: Sequence-based FF (GC-corrected coverage)
 *   - Fragment FF: Fragment size distribution-based FF (gd_4)
 *   - GENDER_DECISION: Weighted ensemble of all FF estimates
 * =========================================================
 */

process CALCULATE_YFF {
    tag "${sample_name}"
    label 'process_low'
    label 'nipt_docker'

    publishDir "${analysisdir}/${sample_name}/Output_FF", mode: 'copy', overwrite: true

    input:
        val  sample_name
        path bam
        path config_json
        val  labcode
        val  analysisdir

    output:
        path "${sample_name}.yff1.txt", emit: yff_txt

    script:
        """
        set -euo pipefail

        python3 /opt/gx-nipt/bin/scripts/modules/ff_gender_improved.py \\
            --mode yff1 \\
            --bam ${bam} \\
            --config ${config_json} \\
            --labcode ${labcode} \\
            --sample ${sample_name} \\
            --output ${sample_name}.yff1.txt

        echo "[YFF1] Complete for ${sample_name}"
        """
}

process CALCULATE_YFF2 {
    tag "${sample_name}"
    label 'process_low'
    label 'nipt_docker'

    publishDir "${analysisdir}/${sample_name}/Output_FF", mode: 'copy', overwrite: true

    input:
        val  sample_name
        path config_json
        val  labcode
        val  analysisdir

    output:
        path "${sample_name}.yff2.txt", emit: yff2_txt

    script:
        def wig_norm = "${analysisdir}/${sample_name}/Output_hmmcopy/${sample_name}.proper_paired.50kb.wig.Normalization.txt"
        """
        set -euo pipefail

        python3 /opt/gx-nipt/bin/scripts/modules/ff_gender_improved.py \\
            --mode yff2 \\
            --wig-norm ${wig_norm} \\
            --config ${config_json} \\
            --labcode ${labcode} \\
            --sample ${sample_name} \\
            --output ${sample_name}.yff2.txt

        echo "[YFF2] Complete for ${sample_name}"
        """
}

process CALCULATE_SEQFF {
    tag "${sample_name}"
    label 'process_low'
    label 'nipt_docker'

    publishDir "${analysisdir}/${sample_name}/Output_FF", mode: 'copy', overwrite: true

    input:
        val  sample_name
        path bam
        path config_json
        val  labcode
        val  analysisdir

    output:
        path "${sample_name}.seqff.txt", emit: seqff_txt

    script:
        def seqff_model = "/data/refs/models/seqff_model.pkl"
        """
        set -euo pipefail

        python3 /opt/gx-nipt/bin/scripts/modules/ff_gender_improved.py \\
            --mode seqff \\
            --bam ${bam} \\
            --config ${config_json} \\
            --labcode ${labcode} \\
            --sample ${sample_name} \\
            --seqff-model ${seqff_model} \\
            --output ${sample_name}.seqff.txt

        echo "[SeqFF] Complete for ${sample_name}"
        """
}

process CALCULATE_FRAGMENT_FF {
    tag "${sample_name}"
    label 'process_low'
    label 'nipt_docker'

    publishDir "${analysisdir}/${sample_name}/Output_FF", mode: 'copy', overwrite: true

    input:
        val  sample_name
        path bam
        path config_json
        val  analysisdir

    output:
        path "${sample_name}.fragment_ff.txt", emit: frag_ff_txt

    script:
        def fragment_cutoff = 160
        """
        set -euo pipefail

        # Fragment-size based FF:
        # Fetal cfDNA is shorter than maternal cfDNA.
        # Ratio of short fragments (<160bp) to total is used as FF proxy.
        python3 /opt/gx-nipt/bin/scripts/modules/ff_gender_improved.py \\
            --mode fragment_ff \\
            --bam ${bam} \\
            --config ${config_json} \\
            --fragment-cutoff ${fragment_cutoff} \\
            --sample ${sample_name} \\
            --output ${sample_name}.fragment_ff.txt

        echo "[Fragment FF] Complete for ${sample_name}"
        """
}

process GENDER_DECISION {
    tag "${sample_name}"
    label 'process_low'
    label 'nipt_docker'

    publishDir "${analysisdir}/${sample_name}/Output_FF", mode: 'copy', overwrite: true

    input:
        val  sample_name
        path yff1_txt
        path yff2_txt
        path seqff_txt
        path frag_ff_txt
        path config_json
        val  analysisdir

    output:
        path "${sample_name}.fetal_fraction.txt", emit: ff_result
        path "${sample_name}.gender.txt",         emit: gender_txt

    script:
        """
        set -euo pipefail

        # Weighted ensemble: YFF2 (primary) > SeqFF > Fragment FF > YFF1
        # Gender decision: gd_2 (YFF2 ratio) is primary indicator
        python3 /opt/gx-nipt/bin/scripts/modules/ff_gender_improved.py \\
            --mode gender_decision \\
            --yff1 ${yff1_txt} \\
            --yff2 ${yff2_txt} \\
            --seqff ${seqff_txt} \\
            --frag-ff ${frag_ff_txt} \\
            --config ${config_json} \\
            --sample ${sample_name} \\
            --ff-output ${sample_name}.fetal_fraction.txt \\
            --gender-output ${sample_name}.gender.txt

        echo "[GENDER_DECISION] Complete for ${sample_name}"
        """
}
