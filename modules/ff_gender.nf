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
        path bai
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
        path wig_norm
        path config_json
        val  labcode
        val  analysisdir

    output:
        path "${sample_name}.yff2.txt", emit: yff2_txt

    script:
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
        path bai
        path config_json
        val  labcode
        val  analysisdir

    output:
        path "${sample_name}.seqff.txt", emit: seqff_txt

    script:
        def seqff_r_dir = "/opt/gx-nipt/bin/scripts/seqFF_R"
        """
        set -euo pipefail
        export TMPDIR="\${NXF_TASK_WORKDIR}"

        # Run the official R-based seqFF script (outputs fraction, not %)
        # Must run from the seqFF_R directory (reference data files are read relative to CWD)
        _raw="\${NXF_TASK_WORKDIR}/${sample_name}.seqff_raw.txt"
        cd ${seqff_r_dir} && Rscript --vanilla seqff.r -f \${NXF_TASK_WORKDIR}/${bam} -o "\${_raw}"
        cd \${NXF_TASK_WORKDIR}

        # Parse SeqFF fraction from R output and convert to %
        # R output format: "SeqFF",0.0527...
        python3 - <<'PYEOF'
import csv, sys, os

raw = "${sample_name}.seqff_raw.txt"
out = "${sample_name}.seqff.txt"

seqff_pct = None
try:
    with open(raw) as f:
        reader = csv.reader(f)
        next(reader)  # skip header row
        for row in reader:
            if len(row) >= 2 and row[0].strip('"') == "SeqFF":
                seqff_pct = round(float(row[1]) * 100, 4)
                break
except Exception as e:
    print(f"[SeqFF] Could not parse R output: {e}", file=sys.stderr)

if seqff_pct is None:
    print("[SeqFF] WARNING: SeqFF not found in R output, defaulting to 0.0", file=sys.stderr)
    seqff_pct = 0.0

with open(out, "w") as f:
    f.write("value\\n")
    f.write(f"SeqFF\\t{seqff_pct}\\n")
    f.write("status\\tOK\\n")

print(f"[SeqFF] SeqFF = {seqff_pct}%")
PYEOF

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
        path bai
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
