/*
 * =========================================================
 *  EZD Workflow (trio)
 *
 *  Input:
 *    norm_50kb : channel of (sample, group, norm.txt) tuples
 *                emitted by HMMCOPY_WORKFLOW for every group.
 *
 *  Runs RUN_EZD once per (sample, group) — Nextflow
 *  parallelises across orig / fetus / mom automatically.
 * =========================================================
 */

include { RUN_EZD } from '../modules/ezd'

workflow EZD_WORKFLOW {
    take:
        norm_50kb     // channel: tuple(sample, group, norm.txt)
        ch_config     // path channel: pipeline_config.json
        labcode       // string
        analysisdir   // string

    main:
        RUN_EZD(
            norm_50kb,
            ch_config,
            labcode,
            analysisdir,
        )

    emit:
        // (sample, group, ezd_results.tsv)
        ezd_result = RUN_EZD.out.ezd_result
}
