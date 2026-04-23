/*
 * =========================================================
 *  Microdeletion Workflow
 *  WC/WCX results → MD target comparison → MD result TSV
 *
 *  The upstream WC_WORKFLOW emits per-group tuples
 *  (sample, group, result_path).  The MD adapter script
 *  (process_md_result.py) scans ${analysisdir} for the
 *  published per-group outputs, so here we only need to
 *  stage the result files as a single collected batch to
 *  fire the MD job exactly once per sample.
 * =========================================================
 */

include { RUN_MD_DETECTION }  from '../modules/microdeletion'

workflow MD_WORKFLOW {
    take:
        sample_name   // string
        ch_wc_result  // channel: tuple(sample, group, result_path)  [WC ∪ WCX]
        ch_config     // channel: pipeline_config.json
        labcode       // string
        analysisdir   // string

    main:
        // Reduce per-group tuples to a single list of result files so
        // RUN_MD_DETECTION is invoked exactly once per sample with all
        // upstream CNV results staged.
        ch_wc_files = ch_wc_result
            .map { sid, grp, p -> p }
            .collect()

        RUN_MD_DETECTION(
            sample_name,
            ch_wc_files,
            ch_config,
            labcode,
            analysisdir
        )

    emit:
        md_result = RUN_MD_DETECTION.out.md_result
}
