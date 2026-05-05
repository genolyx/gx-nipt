/*
 * =========================================================
 *  Wisecondor / gx-cnv Workflow (trio version)
 *
 *  CNV detection stack (parallel execution per group):
 *    WC   – Wisecondor    (orig / fetus / mom)
 *    WCX  – WisecondorX   (orig / fetus / mom, gender-aware)
 *    gx-cnv – Hybrid dual-track (orig BAM only, genome-wide)
 *
 *  Inputs:
 *    bam_trio    – tuple emitted by SAMTOOLS_SPLIT_FETUS_MOM
 *                  (sample, o_bam, o_bai, f_bam, f_bai, m_bam, m_bai)
 *    gender_txt  – path to <sample>.gender.txt produced by
 *                  GENDER_DECISION (carries final_gender line).
 *
 *  Emits per-group tuples:
 *    wc_result          – (sample, group, result_path) for WC ∪ WCX
 *    gxcnv_calls        – (sample, calls_tsv) — orig only
 *    gxcnv_comparison   – (sample, comparison_tsv) — orig only
 * =========================================================
 */

include { RUN_WC }            from '../modules/wisecondor'
include { RUN_WCX }           from '../modules/wisecondor'
include { GXCNV_CONVERT }     from '../modules/gxcnv'
include { GXCNV_GC_INJECT }   from '../modules/gxcnv'
include { GXCNV_PREDICT }     from '../modules/gxcnv'
include { GXCNV_COMPARE }     from '../modules/gxcnv'

workflow WC_WORKFLOW {
    take:
        bam_trio         // channel: tuple(sample, o_bam, o_bai, f_bam, f_bai, m_bam, m_bai)
        gender_txt       // path channel: <sample>.gender.txt
        ch_config        // path: pipeline_config.json
        labcode          // val: lab identifier
        analysisdir      // val: output directory root
        run_gxcnv        // val: boolean — enable gx-cnv (reference auto-resolved from ref_dir/labcode/GXCNV/{sex})
        gxcnv_bin_size   // val: bin size for gxcnv convert (default: 100000)
        gxcnv_thresh_z   // val: Track A Z-score threshold (default: -3.0)
        gxcnv_thresh_p   // val: Track B p-value threshold (default: 0.05)
        wig_norm_orig    // path: HMMcopy 50kb normalization (orig group) for GC injection

    main:
        // ── Flatten trio → per-group (sample, group, bam, bai) tuples ──
        ch_bam_by_group = bam_trio.flatMap { t ->
            def sample = t[0]
            [
                tuple(sample, 'orig',  t[1], t[2]),
                tuple(sample, 'fetus', t[3], t[4]),
                tuple(sample, 'mom',   t[5], t[6]),
            ]
        }

        // ── WC (Wisecondor) ─────────────────────────────────────────────
        RUN_WC(
            ch_bam_by_group,
            ch_config,
            labcode,
            analysisdir
        )

        // ── WCX (WisecondorX) ───────────────────────────────────────────
        ch_wcx_result = Channel.empty()
        if ( params.run_wcx ) {
            RUN_WCX(
                ch_bam_by_group,
                gender_txt,
                ch_config,
                labcode,
                analysisdir
            )
            ch_wcx_result = RUN_WCX.out.wcx_result
        }

        // ── gx-cnv (orig BAM only, validation mode) ─────────────────────
        // Reference is gender-aware: ref_dir/labs/{labcode}/GXCNV/{female|male}/reference.npz
        ch_gxcnv_calls      = Channel.empty()
        ch_gxcnv_comparison = Channel.empty()

        if ( run_gxcnv ) {
            // gxcnv operates on the maternal-plasma (orig) BAM
            ch_bam_orig = bam_trio.map { t -> tuple(t[0], t[1], t[2]) }

            // Step 1: BAM → NPZ (GC fractions will be NaN at this point)
            GXCNV_CONVERT(
                ch_bam_orig,
                gxcnv_bin_size,
                file('NO_FILE')   // optional blacklist BED — not used by default
            )

            // Step 1b: Inject GC fractions from HMMcopy WIG + re-apply GC correction
            // The reference panel was built with GC-corrected NPZ files; without this
            // step the sample Z-scores are systematically biased.
            GXCNV_GC_INJECT(
                GXCNV_CONVERT.out.npz,
                wig_norm_orig
            )

            GXCNV_PREDICT(
                GXCNV_GC_INJECT.out.npz,
                gender_txt,       // parsed in process script to select female/male reference
                labcode,
                gxcnv_thresh_z,
                gxcnv_thresh_p,
                Channel.value('NA'),  // FF passed as placeholder (FF join TODO)
                analysisdir
            )
            ch_gxcnv_calls = GXCNV_PREDICT.out.calls_tsv

            // Compare against WCX orig track (if WCX enabled)
            if ( params.run_wcx ) {
                ch_wcx_orig = ch_wcx_result
                    .filter { sid, grp, p -> grp == 'orig' }
                    .map    { sid, grp, p -> tuple(sid, p) }

                ch_compare_input = GXCNV_PREDICT.out.calls_tsv
                    .join( GXCNV_PREDICT.out.regions_tsv )
                    .join( ch_wcx_orig, remainder: true )
                    .map { sid, calls, regions, aber ->
                        def aber_file = (aber == null) ? file('NO_FILE') : aber
                        tuple( sid, calls, regions, aber_file )
                    }
                GXCNV_COMPARE( ch_compare_input, analysisdir )
                ch_gxcnv_comparison = GXCNV_COMPARE.out.comparison
            }
        }

    emit:
        // Each element: (sample, group, result_path) — WC ∪ WCX
        wc_result        = RUN_WC.out.wc_result.mix( ch_wcx_result )
        gxcnv_calls      = ch_gxcnv_calls
        gxcnv_comparison = ch_gxcnv_comparison
}
