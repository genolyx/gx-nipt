/*
 * =========================================================
 *  Report Workflow
 *  All results → JSON + HTML report → Output directory
 *
 *  Upstream channels that carry per-group tuples
 *  (sample, group, path) — EZD / PRIZM — are first
 *  reduced to their path component before being collected
 *  into a single staging bundle for GENERATE_JSON.
 * =========================================================
 */

include { GENERATE_JSON }    from '../modules/report'
include { GENERATE_HTML }    from '../modules/report'
include { COPY_TO_OUTPUT }   from '../modules/report'

workflow REPORT_WORKFLOW {
    take:
        sample_name      // string
        ch_ezd_result    // channel: tuple(sample, group, path)
        ch_prizm_result  // channel: tuple(sample, group, path)
        ch_ff_result     // channel: path
        ch_md_result     // channel: path
        ch_qc_result     // channel: path
        ch_config        // channel: pipeline_config.json
        labcode          // string
        analysisdir      // string
        outdir           // string
        ch_ff_ensemble   // channel: tuple(sample_id, path) — ensures GXFF_ENSEMBLE finishes before report

    main:
        // Normalise per-group tuple channels to plain path channels
        // so `.mix().collect()` produces a clean list of files.
        ch_ezd_paths    = ch_ezd_result.map   { sid, grp, p -> p }
        ch_prizm_paths  = ch_prizm_result.map { sid, grp, p -> p }

        // Strip sample_id key from ff_ensemble tuple; mix in so report waits
        // for GXFF_ENSEMBLE to publish ff_ensemble.tsv before generating JSON.
        ch_ff_ensemble_path = ch_ff_ensemble.map { _sid, p -> p }

        ch_all_results = ch_ezd_paths
            .mix(ch_prizm_paths)
            .mix(ch_ff_result)
            .mix(ch_md_result)
            .mix(ch_qc_result)
            .mix(ch_ff_ensemble_path)
            .collect()

        GENERATE_JSON(
            sample_name,
            ch_all_results,
            ch_config,
            labcode,
            analysisdir
        )

        GENERATE_HTML(
            sample_name,
            GENERATE_JSON.out.json_file,
            ch_config,
            labcode,
            analysisdir
        )

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
