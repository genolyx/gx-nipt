/*
 * =========================================================
 *  Report Workflow
 *  All results → JSON + HTML report → Output directory
 * =========================================================
 */

include { GENERATE_JSON }    from '../modules/report'
include { GENERATE_HTML }    from '../modules/report'
include { COPY_TO_OUTPUT }   from '../modules/report'

workflow REPORT_WORKFLOW {
    take:
        sample_name    // string
        ch_ezd_result  // channel
        ch_prizm_result// channel
        ch_ff_result   // channel
        ch_md_result   // channel
        ch_qc_result   // channel
        ch_config      // channel: pipeline_config.json
        labcode        // string
        analysisdir    // string
        outdir         // string

    main:
        // Collect all upstream results before generating report
        ch_all_results = ch_ezd_result
            .mix(ch_prizm_result)
            .mix(ch_ff_result)
            .mix(ch_md_result)
            .mix(ch_qc_result)
            .collect()

        // Generate final JSON output
        GENERATE_JSON(
            sample_name,
            ch_all_results,
            ch_config,
            labcode,
            analysisdir
        )

        // Generate HTML review page
        GENERATE_HTML(
            sample_name,
            GENERATE_JSON.out.json_file,
            ch_config,
            labcode,
            analysisdir
        )

        // Copy results to portal-compatible output directory
        COPY_TO_OUTPUT(
            sample_name,
            GENERATE_JSON.out.json_file,
            GENERATE_HTML.out.html_file,
            analysisdir,
            outdir
        )

    emit:
        report_json = GENERATE_JSON.out.json_file
        report_html = GENERATE_HTML.out.html_file
}
