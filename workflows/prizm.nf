/*
 * =========================================================
 *  PRIZM Workflow (trio)
 *
 *  Input:
 *    norm_10mb : channel of (sample, group, norm.txt) tuples
 *                emitted by HMMCOPY_WORKFLOW for every group.
 *                PRIZM ingests this 10mb normalisation TXT
 *                directly as its "count_file_10mb".
 *    gender_txt: path channel — <sample>.gender.txt from
 *                GENDER_DECISION, used by PRIZM adapter to
 *                pick male/female reference CSVs.
 * =========================================================
 */

include { RUN_PRIZM } from '../modules/prizm'

workflow PRIZM_WORKFLOW {
    take:
        norm_10mb       // channel: tuple(sample, group, norm.txt)
        gender_txt      // path channel: <sample>.gender.txt
        ch_config       // path channel: pipeline_config.json
        labcode         // string
        analysisdir     // string

    main:
        RUN_PRIZM(
            norm_10mb,
            gender_txt,
            ch_config,
            labcode,
            analysisdir,
        )

    emit:
        // (sample, group, prizm.qc.txt)
        prizm_result = RUN_PRIZM.out.prizm_result
}
