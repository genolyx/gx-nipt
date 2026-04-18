/*
 * =========================================================
 *  Microdeletion Workflow
 *  WC/WCX results → MD target comparison → MD result TSV
 * =========================================================
 */

include { RUN_MD_DETECTION }  from '../modules/microdeletion'

workflow MD_WORKFLOW {
    take:
        sample_name   // string
        ch_wc_result  // channel: WC + WCX result files
        ch_config     // channel: pipeline_config.json
        labcode       // string
        analysisdir   // string

    main:
        RUN_MD_DETECTION(
            sample_name,
            ch_wc_result,
            ch_config,
            labcode,
            analysisdir
        )

    emit:
        md_result = RUN_MD_DETECTION.out.md_result
}
