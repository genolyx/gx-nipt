/*
 * =========================================================
 *  EZD Workflow
 *  HMMcopy 50kb norm → EZD (orig / fetus / mom)
 * =========================================================
 */

include { RUN_EZD }  from '../modules/ezd'

workflow EZD_WORKFLOW {
    take:
        sample_name   // string
        ch_norm_50kb  // channel: .50kb.wig.Normalization.txt
        ch_ff_result  // channel: fetal_fraction.txt (for fetus/mom split)
        ch_config     // channel: pipeline_config.json
        labcode       // string
        analysisdir   // string

    main:
        // Run EZD for all three groups in parallel
        ch_groups = Channel.of('orig', 'fetus', 'mom')

        RUN_EZD(
            sample_name,
            ch_groups,
            ch_norm_50kb,
            ch_ff_result,
            ch_config,
            labcode,
            analysisdir
        )

    emit:
        ezd_result = RUN_EZD.out.ezd_result
}
