/*
 * =========================================================
 *  Fetal Fraction & Gender Detection Workflow
 *  YFF1 + YFF2 + SeqFF + Fragment FF → Gender Decision
 * =========================================================
 */

include { CALCULATE_YFF }        from '../modules/ff_gender'
include { CALCULATE_YFF2 }       from '../modules/ff_gender'
include { CALCULATE_SEQFF }      from '../modules/ff_gender'
include { CALCULATE_FRAGMENT_FF} from '../modules/ff_gender'
include { GENDER_DECISION }      from '../modules/ff_gender'

workflow FF_GENDER_WORKFLOW {
    take:
        sample_name   // string
        ch_bam        // channel: proper_paired.bam
        ch_config     // channel: pipeline_config.json
        labcode       // string
        analysisdir   // string

    main:
        // YFF1: Y-chromosome coverage ratio (gd_1)
        CALCULATE_YFF(
            sample_name,
            ch_bam,
            ch_config,
            labcode,
            analysisdir
        )

        // YFF2: Adjusted YFF using wig normalization (gd_2)
        // Depends on HMMcopy 50kb wig (produced in parallel)
        // Uses a deferred channel join - wig file is passed via analysisdir path
        CALCULATE_YFF2(
            sample_name,
            ch_config,
            labcode,
            analysisdir
        )

        // SeqFF: Sequence-based FF from GC-corrected coverage
        CALCULATE_SEQFF(
            sample_name,
            ch_bam,
            ch_config,
            labcode,
            analysisdir
        )

        // Fragment-size based FF (gd_4)
        CALCULATE_FRAGMENT_FF(
            sample_name,
            ch_bam,
            ch_config,
            analysisdir
        )

        // Final gender decision: combine all FF estimates
        GENDER_DECISION(
            sample_name,
            CALCULATE_YFF.out.yff_txt,
            CALCULATE_YFF2.out.yff2_txt,
            CALCULATE_SEQFF.out.seqff_txt,
            CALCULATE_FRAGMENT_FF.out.frag_ff_txt,
            ch_config,
            analysisdir
        )

    emit:
        ff_result  = GENDER_DECISION.out.ff_result
        gender_txt = GENDER_DECISION.out.gender_txt
}
