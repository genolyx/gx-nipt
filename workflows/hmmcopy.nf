/*
 * =========================================================
 *  HMMcopy Workflow (trio version)
 *
 *  Inputs:
 *    - bam_trio : tuple (sample, orig_bam, orig_bai,
 *                               fetus_bam, fetus_bai,
 *                               mom_bam, mom_bai)
 *                 produced by SAMTOOLS_SPLIT_FETUS_MOM.
 *
 *  Flow (per group x resolution):
 *    BAM → readCounter → HMMcopy.R normalisation
 *
 *  Emits per-resolution channels keyed on (sample, group):
 *    norm_50kb : [sample, group, norm.txt]  - used by EZD
 *    norm_10mb : [sample, group, norm.txt]  - used by PRIZM
 *
 *  Note: PRIZM ingests the 10mb normalisation TXT directly;
 *  no separate wig→count conversion is needed (kept in sync
 *  with prizm_runner.run_prizm_analysis()).
 * =========================================================
 */

// Nextflow DSL2 forbids calling the same process twice in one
// workflow context — alias the shared processes once per resolution.
include { READCOUNTER as READCOUNTER_50K } from '../modules/hmmcopy'
include { READCOUNTER as READCOUNTER_10M } from '../modules/hmmcopy'
include { HMMCOPY_R   as HMMCOPY_R_50K   } from '../modules/hmmcopy'
include { HMMCOPY_R   as HMMCOPY_R_10M   } from '../modules/hmmcopy'

workflow HMMCOPY_WORKFLOW {
    take:
        bam_trio     // channel: tuple(sample, o_bam, o_bai, f_bam, f_bai, m_bam, m_bai)
        labcode      // string (kept for symmetry; unused here)
        analysisdir  // string

    main:
        // Flatten the trio tuple into one (sample, group, bam, bai)
        // channel entry per group so Nextflow parallelises across
        // orig / fetus / mom naturally.
        ch_bam_by_group = bam_trio.flatMap { t ->
            def sample = t[0]
            [
                tuple(sample, 'orig',  t[1], t[2]),
                tuple(sample, 'fetus', t[3], t[4]),
                tuple(sample, 'mom',   t[5], t[6]),
            ]
        }

        // ── readCounter at 50kb (for EZD) ─────────────────────
        READCOUNTER_50K(ch_bam_by_group, '50kb', analysisdir)

        // ── readCounter at 10mb (for PRIZM) ───────────────────
        READCOUNTER_10M(ch_bam_by_group, '10mb', analysisdir)

        // ── HMMcopy.R normalisation @ 50kb ────────────────────
        HMMCOPY_R_50K(READCOUNTER_50K.out.wig, '50kb', analysisdir)

        // ── HMMcopy.R normalisation @ 10mb ────────────────────
        HMMCOPY_R_10M(READCOUNTER_10M.out.wig, '10mb', analysisdir)

    emit:
        // Each channel element: (sample, group, path)
        norm_50kb = HMMCOPY_R_50K.out.norm
        norm_10mb = HMMCOPY_R_10M.out.norm
}
