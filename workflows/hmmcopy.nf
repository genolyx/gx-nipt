/*
 * =========================================================
 *  HMMcopy Workflow
 *  BAM → readCounter (50kb + 10mb) → HMMcopy.R normalization
 * =========================================================
 */

include { READCOUNTER_50KB }  from '../modules/hmmcopy'
include { READCOUNTER_10MB }  from '../modules/hmmcopy'
include { HMMCOPY_R_50KB }    from '../modules/hmmcopy'
include { HMMCOPY_R_10MB }    from '../modules/hmmcopy'
include { COUNT_10MB }        from '../modules/hmmcopy'

workflow HMMCOPY_WORKFLOW {
    take:
        sample_name   // string
        ch_bam        // channel: proper_paired.bam
        labcode       // string
        analysisdir   // string

    main:
        // readCounter at 50kb resolution (for EZD)
        READCOUNTER_50KB(
            sample_name,
            ch_bam,
            analysisdir
        )

        // readCounter at 10mb resolution (for PRIZM)
        READCOUNTER_10MB(
            sample_name,
            ch_bam,
            analysisdir
        )

        // HMMcopy normalization at 50kb
        HMMCOPY_R_50KB(
            sample_name,
            READCOUNTER_50KB.out.wig,
            analysisdir
        )

        // HMMcopy normalization at 10mb
        HMMCOPY_R_10MB(
            sample_name,
            READCOUNTER_10MB.out.wig,
            analysisdir
        )

        // Count matrix for PRIZM (10mb bins)
        COUNT_10MB(
            sample_name,
            READCOUNTER_10MB.out.wig,
            analysisdir
        )

    emit:
        norm_50kb  = HMMCOPY_R_50KB.out.norm_txt
        norm_10mb  = HMMCOPY_R_10MB.out.norm_txt
        count_10mb = COUNT_10MB.out.count_txt
}
