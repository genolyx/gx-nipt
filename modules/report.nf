/*
 * =========================================================
 *  Module: Report Generation & Output Copy
 *
 *  Portal-compatible output directory structure (ken-nipt flat plots):
 *    output/{work_dir}/{sample_name}/
 *      ├── Output_QC/
 *      ├── Output_FF/
 *      ├── Output_EZD/{orig,fetus,mom}_EZD_grid.png   (+ nested group dirs kept)
 *      ├── Output_PRIZM/{sample}_{group}_*.png
 *      ├── Output_WC/{sample}.wc.{group}_z.png
 *      ├── Output_WCX/{sample}.wcx.{group}.png
 *      ├── Output_hmmcopy/
 *      ├── Output_MD/
 *      ├── gxcnv1/ gxcnv2/
 *      └── Output_Result/
 * =========================================================
 */

process GENERATE_JSON {
    tag "${sample_name}"
    label 'process_low'
    label 'nipt_docker'

    publishDir { "${analysisdir}/${sample_name}/Output_Result" }, mode: 'copy', overwrite: true

    input:
        val  sample_name
        path all_results   // collected upstream outputs (staging only)
        path config_json
        val  labcode
        val  analysisdir

    output:
        path "${sample_name}.result.json", emit: json_file

    script:
        // generate_json_output.py uses single-dash flags (ken-nipt legacy CLI)
        // and writes to ``<output_dir>/<sample_name>/<sample_name>.json``.
        // We run it with output_dir=. and rename to the Nextflow-expected name.
        def bed_dir = "${params.ref_dir}/labs/${labcode}/bed"
        // Maternal age drives risk_before tables; daemon/run_nipt pass --age into
        // params.age. Fall back to 35 only when from_bam/algorithm_only omit it.
        def maternal_age = (params.age != null && params.age.toString().trim()) ? params.age : 35
        """
        set -euo pipefail

        export TMPDIR="\${NXF_TASK_WORKDIR}"
        export MPLCONFIGDIR="\${NXF_TASK_WORKDIR}"

        python3 /opt/gx-nipt/bin/scripts/modules/generate_json_output.py \\
            -analysis_dir   ${analysisdir} \\
            -sample_name    ${sample_name} \\
            -output_dir     . \\
            -target_bed_dir ${bed_dir} \\
            -ref_dir        ${params.ref_dir}/labs/${labcode} \\
            -config_file    ${config_json} \\
            -version        "gx-nipt-1.0" \\
            -age            ${maternal_age}

        # The script writes <cwd>/<sample_name>/<sample_name>.json — hoist
        # it up so it matches the declared output.
        if [ -f "./${sample_name}/${sample_name}.json" ]; then
            mv "./${sample_name}/${sample_name}.json" "./${sample_name}.result.json"
        else
            echo "[REPORT] ERROR: expected JSON not produced at ./${sample_name}/${sample_name}.json" >&2
            ls -la ./${sample_name} || true
            exit 1
        fi

        echo "[REPORT] JSON generated for ${sample_name} (age=${maternal_age})"
        """
}

process GENERATE_HTML {
    tag "${sample_name}"
    label 'process_low'
    label 'nipt_docker'

    publishDir { "${analysisdir}/${sample_name}/Output_Result" }, mode: 'copy', overwrite: true

    input:
        val  sample_name
        path json_file
        path config_json
        val  labcode
        val  analysisdir

    output:
        path "${sample_name}.review.html", emit: html_file

    script:
        // html_review_page.py CLI: -json_file / -output_dir / -sample_name.
        // It writes ``<output_dir>/<order_id>_report.html`` where order_id is
        // parsed from the JSON's final_results. We glob the produced file and
        // rename it to the Nextflow-expected name.
        """
        set -euo pipefail

        export TMPDIR="\${NXF_TASK_WORKDIR}"
        export MPLCONFIGDIR="\${NXF_TASK_WORKDIR}"

        python3 /opt/gx-nipt/bin/scripts/modules/html_review_page.py \\
            -json_file   ${json_file} \\
            -output_dir  . \\
            -sample_name ${sample_name}

        # The script names the file <order_id>_report.html. Rename the first
        # match to the Nextflow-declared output.
        html_src=\$(ls -1 ./*_report.html 2>/dev/null | head -n1 || true)
        if [ -z "\${html_src}" ]; then
            echo "[REPORT] ERROR: HTML was not produced." >&2
            ls -la . >&2
            exit 1
        fi
        mv "\${html_src}" "./${sample_name}.review.html"

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
        // Barrier only: count of gxcnv1/2 genome plot basenames (no path staging —
        // gxcnv1 & gxcnv2 share *_{group}_genome.png names).
        val  gxcnv_plots_ready

    output:
        path ".copy_done", emit: done

    script:
        """
        set -euo pipefail

        export TMPDIR="\${NXF_TASK_WORKDIR}"
        export MPLCONFIGDIR="\${NXF_TASK_WORKDIR}"

        echo "[COPY] gxcnv plot barrier count=${gxcnv_plots_ready}"

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
        mkdir -p ${outdir}/gxcnv1
        mkdir -p ${outdir}/gxcnv2

        # Copy analysis results to output directory
        _copy() { [ -d "\$1" ] && cp -a "\$1/." "\$2/" || true; }
        _copy ${analysisdir}/${sample_name}/Output_QC      ${outdir}/Output_QC
        _copy ${analysisdir}/${sample_name}/Output_FF      ${outdir}/Output_FF
        _copy ${analysisdir}/${sample_name}/Output_EZD     ${outdir}/Output_EZD
        _copy ${analysisdir}/${sample_name}/Output_PRIZM   ${outdir}/Output_PRIZM
        _copy ${analysisdir}/${sample_name}/Output_WC      ${outdir}/Output_WC
        _copy ${analysisdir}/${sample_name}/Output_WCX     ${outdir}/Output_WCX
        _copy ${analysisdir}/${sample_name}/Output_hmmcopy ${outdir}/Output_hmmcopy
        _copy ${analysisdir}/${sample_name}/Output_MD      ${outdir}/Output_MD
        _copy ${analysisdir}/${sample_name}/gxcnv1         ${outdir}/gxcnv1
        _copy ${analysisdir}/${sample_name}/gxcnv2         ${outdir}/gxcnv2

        # Flatten to ken-nipt / Portal paths + re-sync gxcnv1/2 + WC/WCX aliases.
        # Nested Output_*/{group}/ trees are kept for MD / internal use.
        python3 ${projectDir}/bin/scripts/modules/portal_output_layout.py \\
            --sample ${sample_name} \\
            --analysis-dir ${analysisdir} \\
            --outdir ${outdir}

        # Copy final report
        cp ${json_file} ${outdir}/Output_Result/
        cp ${html_file} ${outdir}/Output_Result/

        # Do NOT write <sample>.completed here. When sample_name == order_id
        # that file is the daemon success marker; writing it before
        # run_nipt.sh publishes <order_id>.json causes false COMPLETED
        # recovery. run_nipt.sh writes the marker after JSON + tar are ready.

        touch .copy_done
        echo "[COPY] Output copied to ${outdir}"
        """
}
