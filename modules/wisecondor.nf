/*
 * =========================================================
 *  Module: Wisecondor & WisecondorX (trio version)
 *
 *  Each process operates on one (sample, group, bam) tuple.
 *  Parallelisation across orig / fetus / mom is driven by
 *  the caller workflow (workflows/wisecondor.nf).
 *
 *  Gender for WCX is resolved from the GENDER_DECISION output
 *  (<sample>.gender.txt), which carries a line:
 *      final_gender\tMALE|FEMALE
 * =========================================================
 */

process RUN_WC {
    tag "${sample_name}:${group}"
    label 'process_low'
    label 'nipt_docker'

    publishDir "${analysisdir}/${sample_name}/Output_WC/${group}", mode: 'copy', overwrite: true

    input:
        tuple val(sample_name), val(group), path(bam), path(bai)
        path config_json
        val  labcode
        val  analysisdir

    output:
        tuple val(sample_name), val(group),
              path("${sample_name}.wc.${group}.report.txt"),
              emit: wc_result
        path "*.npz",                                 emit: npz, optional: true

    script:
        // Actual NPZ naming convention in refs:
        //   orig  → orig_200k_proper_paired.npz
        //   fetus → fetus_200k_of.npz
        //   mom   → mom_200k_of.npz
        def wc_ref_names = [orig: 'orig_200k_proper_paired.npz',
                            fetus: 'fetus_200k_of.npz',
                            mom:   'mom_200k_of.npz']
        def ref_npz = "${params.ref_dir}/labs/${labcode}/WC/${wc_ref_names[group]}"
        def binsize  = 200000
        """
        set -euo pipefail

        # python2 getpass.getuser() fails when uid has no /etc/passwd entry
        # (Docker --user 1002:1000).  Set LOGNAME so it returns early.
        export LOGNAME="\${LOGNAME:-nipt}"
        # wisecondor.py uses matplotlib for report plots; redirect font cache.
        export MPLCONFIGDIR="\${NXF_TASK_WORKDIR}"

        # Skip WC gracefully when the BAM has no mapped reads (e.g. empty
        # fetus/mom BAM from TLEN-split of a very-low-FF sample).
        READ_COUNT=\$(samtools view -c ${bam})
        if [ "\${READ_COUNT}" -eq 0 ]; then
            echo "[WC] ${group} BAM is empty (0 reads) — writing empty report." >&2
            printf "WC\\t${group}\\tSKIPPED (empty BAM)\\n" > ${sample_name}.wc.${group}.report.txt
            exit 0
        fi

        # ── Step 1: BAM → NPZ (old Python-2 Wisecondor) ──────────────────────
        python2 /opt/wisecondor/wisecondor.py convert \\
            ${bam} \\
            ${sample_name}.wc.${group}.npz \\
            -binsize ${binsize}

        # ── Step 2: Test (predict against reference) ──────────────────────────
        python2 /opt/wisecondor/wisecondor.py test \\
            ${sample_name}.wc.${group}.npz \\
            ${sample_name}.wc.${group}.out.npz \\
            ${ref_npz}

        # ── Step 3: Report → report.txt ───────────────────────────────────────
        python2 /opt/wisecondor/wisecondor.py report \\
            ${sample_name}.wc.${group}.npz \\
            ${sample_name}.wc.${group}.out.npz \\
            > ${sample_name}.wc.${group}.report.txt

        echo "[WC] ${group} complete for ${sample_name}"
        """
}

process RUN_WCX {
    tag "${sample_name}:${group}"
    label 'process_low'
    label 'nipt_docker'

    publishDir "${analysisdir}/${sample_name}/Output_WCX/${group}", mode: 'copy', overwrite: true

    input:
        tuple val(sample_name), val(group), path(bam), path(bai)
        path gender_txt
        path config_json
        val  labcode
        val  analysisdir

    output:
        tuple val(sample_name), val(group),
              path("${sample_name}.wcx.${group}_aberrations.bed"),
              emit: wcx_result
        path "${sample_name}.wcx.${group}.plots", emit: plots_dir, optional: true
        path "*.npz",                             emit: npz,      optional: true

    script:
        def binsize = 200000
        """
        set -euo pipefail

        # WisecondorX calls CBS.R (R) and uses matplotlib; both need writable temp dirs.
        export TMPDIR="\${NXF_TASK_WORKDIR}"
        export R_TEMPDIR="\${NXF_TASK_WORKDIR}"
        export MPLCONFIGDIR="\${NXF_TASK_WORKDIR}"

        # Parse gender from GENDER_DECISION output (final_gender line)
        # Accepts MALE/FEMALE (canonical) and XY/XX (written by ff_gender_improved.py)
        GENDER=\$(awk -F'[\\t ]+' 'tolower(\$1) == "final_gender" {
            g = toupper(\$2);
            if (g == "MALE" || g == "M" || g == "XY") { print "male"; exit }
            if (g == "FEMALE" || g == "F" || g == "XX") { print "female"; exit }
        }' ${gender_txt})
        if [ -z "\${GENDER}" ]; then
            echo "[WCX] Could not parse final_gender from ${gender_txt}; defaulting to female." >&2
            GENDER="female"
        fi

        # Actual NPZ naming convention in refs:
        #   orig  → orig_{M|F}_200k_proper_paired.npz  (gender-split)
        #   fetus → fetus_{M|F}_200k_of.npz            (gender-split)
        #   mom   → mom_200k_of.npz                    (no gender split)
        if [ "${group}" = "mom" ]; then
            REF_NPZ="${params.ref_dir}/labs/${labcode}/WCX/mom_200k_of.npz"
        elif [ "\${GENDER}" = "male" ]; then
            case "${group}" in
                orig)  REF_NPZ="${params.ref_dir}/labs/${labcode}/WCX/orig_M_200k_proper_paired.npz" ;;
                fetus) REF_NPZ="${params.ref_dir}/labs/${labcode}/WCX/fetus_M_200k_of.npz" ;;
            esac
        else
            case "${group}" in
                orig)  REF_NPZ="${params.ref_dir}/labs/${labcode}/WCX/orig_F_200k_proper_paired.npz" ;;
                fetus) REF_NPZ="${params.ref_dir}/labs/${labcode}/WCX/fetus_F_200k_of.npz" ;;
            esac
        fi

        # Skip WCX gracefully when the BAM has no mapped reads.
        READ_COUNT=\$(samtools view -c ${bam})
        if [ "\${READ_COUNT}" -eq 0 ]; then
            echo "[WCX] ${group} BAM is empty (0 reads) — writing empty aberrations file." >&2
            touch ${sample_name}.wcx.${group}_aberrations.bed
            exit 0
        fi

        # Convert BAM to npz
        # NOTE: WisecondorX convert appends .npz automatically,
        # so pass the prefix without extension.
        WisecondorX convert \\
            ${bam} \\
            ${sample_name}.wcx.${group} \\
            --binsize ${binsize}

        # Run WisecondorX predict with gender-aware reference
        WisecondorX predict \\
            ${sample_name}.wcx.${group}.npz \\
            \${REF_NPZ} \\
            ${sample_name}.wcx.${group} \\
            --zscore 6.0 \\
            --bed

        echo "[WCX] ${group} complete for ${sample_name}"
        """
}
