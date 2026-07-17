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
 *    norm_50kb : [sample, group, norm.txt]  - used by EZD / FF
 *    norm_10mb : [sample, group, norm.txt]  - used by PRIZM
 *
 *  ken-nipt input BAM per group (must stay in sync):
 *    orig 50kb  → proper_paired.bam   (max reads for EZD / FF)
 *    orig 10mb  → of_orig.bam         (Uniform_2017_allY.bed filter; PRIZM orig)
 *    fetus/mom  → of_fetus / of_mom   (both resolutions)
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
        bam_of_orig  // channel: tuple(sample, of_orig_bam, of_orig_bai) — PRIZM orig 10mb only
        labcode      // string (kept for symmetry; unused here)
        analysisdir  // string

    main:
        // 50kb: all groups from trio (orig slot = proper_paired.bam, ken-nipt EZD / FF path).
        ch_bam_by_group_50kb = bam_trio.flatMap { t ->
            def sample = t[0]
            [
                tuple(sample, 'orig',  t[1], t[2]),
                tuple(sample, 'fetus', t[3], t[4]),
                tuple(sample, 'mom',   t[5], t[6]),
            ]
        }

        // 10mb: fetus/mom from trio; orig from bed-filtered of_orig.bam (ken-nipt PRIZM path).
        ch_bam_10mb_fetus_mom = bam_trio.flatMap { t ->
            def sample = t[0]
            [
                tuple(sample, 'fetus', t[3], t[4]),
                tuple(sample, 'mom',   t[5], t[6]),
            ]
        }
        ch_bam_10mb_orig = bam_of_orig.map { t ->
            tuple(t[0], 'orig', t[1], t[2])
        }
        ch_bam_by_group_10mb = ch_bam_10mb_fetus_mom.mix(ch_bam_10mb_orig)

        // ── readCounter at 50kb (for EZD) ─────────────────────
        READCOUNTER_50K(ch_bam_by_group_50kb, '50kb', analysisdir)

        // ── readCounter at 10mb (for PRIZM) ───────────────────
        READCOUNTER_10M(ch_bam_by_group_10mb, '10mb', analysisdir)

        // ── HMMcopy.R normalisation @ 50kb ────────────────────
        HMMCOPY_R_50K(READCOUNTER_50K.out.wig, '50kb', analysisdir)

        // ── HMMcopy.R normalisation @ 10mb ────────────────────
        HMMCOPY_R_10M(READCOUNTER_10M.out.wig, '10mb', analysisdir)

    emit:
        // Each channel element: (sample, group, path)
        norm_50kb = HMMCOPY_R_50K.out.norm
        norm_10mb = HMMCOPY_R_10M.out.norm
}
