/*
 * =========================================================
 *  Wisecondor Workflow
 *  BAM → WC (Wisecondor) + WCX (WisecondorX) → orig/fetus/mom
 * =========================================================
 */

include { RUN_WC }   from '../modules/wisecondor'
include { RUN_WCX }  from '../modules/wisecondor'

workflow WC_WORKFLOW {
    take:
        sample_name   // string
        ch_bam        // channel: proper_paired.bam
        ch_ff_result  // channel: fetal_fraction.txt
        ch_config     // channel: pipeline_config.json
        labcode       // string
        analysisdir   // string

    main:
        ch_groups = Channel.of('orig', 'fetus', 'mom')

        // WC (Wisecondor) for each group
        RUN_WC(
            sample_name,
            ch_groups,
            ch_bam,
            ch_ff_result,
            ch_config,
            labcode,
            analysisdir
        )

        // WCX (WisecondorX) for each group - gender-aware reference
        RUN_WCX(
            sample_name,
            ch_groups,
            ch_bam,
            ch_ff_result,
            ch_config,
            labcode,
            analysisdir
        )

    emit:
        wc_result = RUN_WC.out.wc_result.mix(RUN_WCX.out.wcx_result)
}
