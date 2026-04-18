/*
 * =========================================================
 *  PRIZM Workflow
 *  HMMcopy 10mb count → PRIZM (orig / fetus / mom)
 * =========================================================
 */

include { RUN_PRIZM }  from '../modules/prizm'

workflow PRIZM_WORKFLOW {
    take:
        sample_name    // string
        ch_count_10mb  // channel: 10mb count file
        ch_ff_result   // channel: fetal_fraction.txt
        ch_config      // channel: pipeline_config.json
        labcode        // string
        analysisdir    // string

    main:
        ch_groups = Channel.of('orig', 'fetus', 'mom')

        RUN_PRIZM(
            sample_name,
            ch_groups,
            ch_count_10mb,
            ch_ff_result,
            ch_config,
            labcode,
            analysisdir
        )

    emit:
        prizm_result = RUN_PRIZM.out.prizm_result
}
