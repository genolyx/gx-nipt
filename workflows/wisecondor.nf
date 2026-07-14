/*
 * =========================================================
 *  Wisecondor / gx-cnv Workflow (trio version)
 *
 *  CNV detection stack (parallel execution per group):
 *    WC   – Wisecondor    (orig / fetus / mom)
 *    WCX  – WisecondorX   (orig / fetus / mom, gender-aware)
 *    gx-cnv – Hybrid dual-track (orig / fetus / mom when refs available)
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
include { GXCNV_PLOT }        from '../modules/gxcnv'
include { GXCNV2_PREDICT }    from '../modules/gxcnv2'
include { GXCNV2_PLOT }       from '../modules/gxcnv2'
include { GXCNV1_PREDICT }    from '../modules/gxcnv1'
include { GXCNV1_PLOT }       from '../modules/gxcnv1'

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
        run_gxcnv2       // val: boolean — enable gxcnv2 (WCX-core, reuses RUN_WCX BED outputs)
        run_gxcnv1       // val: boolean — enable gxcnv1 (WC-core, reuses RUN_WC NPZ outputs)

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
        ch_wc_result = Channel.empty()
        if ( params.run_wc ) {
            RUN_WC(
                ch_bam_by_group,
                ch_config,
                labcode,
                analysisdir
            )
            ch_wc_result = RUN_WC.out.wc_result
        }

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

        // ── gx-cnv (orig / fetus / mom, mirroring WCX) ──────────────────
        // Reference layout:
        //   {ref_dir}/labs/{labcode}/GXCNV/{female|male}/reference.npz          ← orig
        //   {ref_dir}/labs/{labcode}/GXCNV/{female|male}_fetus/reference.npz    ← fetus (optional)
        //   {ref_dir}/labs/{labcode}/GXCNV/{female|male}_mom/reference.npz      ← mom   (optional)
        // If a group-specific reference does not exist the process will print a
        // warning and exit 0 with empty output, so the downstream join is safe.
        ch_gxcnv_calls      = Channel.empty()
        ch_gxcnv_comparison = Channel.empty()

        if ( run_gxcnv ) {
            // Build per-group BAM channels exactly like WCX
            ch_bam_all_groups = bam_trio.flatMap { t ->
                [
                    tuple(t[0], 'orig',  t[1], t[2]),
                    tuple(t[0], 'fetus', t[3], t[4]),
                    tuple(t[0], 'mom',   t[5], t[6]),
                ]
            }

            // Step 1: BAM → NPZ for all groups
            GXCNV_CONVERT(
                ch_bam_all_groups.map { sid, grp, bam, bai -> tuple("${sid}_${grp}", bam, bai) },
                gxcnv_bin_size,
                file('NO_FILE')
            )

            // Step 1b: GC inject (use orig WIG for all groups — same genome binning)
            GXCNV_GC_INJECT(
                GXCNV_CONVERT.out.npz,
                wig_norm_orig
            )

            // Restore (sample_id, group) structure from the composite key
            ch_gc_npz_by_group = GXCNV_GC_INJECT.out.npz
                .map { composite_id, npz ->
                    def parts = composite_id.toString().tokenize('_')
                    def grp   = parts[-1]            // orig / fetus / mom
                    def sid   = parts[0..-2].join('_')
                    tuple(sid, grp, npz)
                }

            GXCNV_PREDICT(
                ch_gc_npz_by_group.map { sid, grp, npz -> tuple("${sid}_${grp}", npz) },
                gender_txt,
                labcode,
                gxcnv_thresh_z,
                gxcnv_thresh_p,
                Channel.value('NA'),
                analysisdir
            )

            // Emit calls keyed back to (sample_id, group) for downstream join
            ch_gxcnv_calls = GXCNV_PREDICT.out.calls_tsv
                .map { composite_id, calls ->
                    def parts = composite_id.toString().tokenize('_')
                    def grp   = parts[-1]
                    def sid   = parts[0..-2].join('_')
                    tuple(sid, grp, calls)
                }

            // Step 5: Plot — genome-wide + per-chromosome figures
            ch_plot_input = GXCNV_PREDICT.out.bins_tsv
                .join( GXCNV_PREDICT.out.calls_tsv )
            GXCNV_PLOT( ch_plot_input, analysisdir )

            // Cross-check against WCX orig track (internal validation)
            if ( params.run_wcx ) {
                ch_wcx_orig = ch_wcx_result
                    .filter { sid, grp, p -> grp == 'orig' }
                    .map    { sid, grp, p -> tuple(sid, p) }

                ch_gxcnv_orig_calls = GXCNV_PREDICT.out.calls_tsv
                    .map { composite_id, calls ->
                        def parts = composite_id.toString().tokenize('_')
                        def grp   = parts[-1]
                        def sid   = parts[0..-2].join('_')
                        tuple(sid, grp, calls)
                    }
                    .filter { sid, grp, calls -> grp == 'orig' }
                    .map    { sid, grp, calls -> tuple(sid, calls) }

                ch_gxcnv_orig_regions = GXCNV_PREDICT.out.regions_tsv
                    .map { composite_id, regions ->
                        def parts = composite_id.toString().tokenize('_')
                        def grp   = parts[-1]
                        def sid   = parts[0..-2].join('_')
                        tuple(sid, grp, regions)
                    }
                    .filter { sid, grp, regions -> grp == 'orig' }
                    .map    { sid, grp, regions -> tuple(sid, regions) }

                ch_compare_input = ch_gxcnv_orig_calls
                    .join( ch_gxcnv_orig_regions )
                    .join( ch_wcx_orig, remainder: true )
                    .map { sid, calls, regions, aber ->
                        def aber_file = (aber == null) ? file('NO_FILE') : aber
                        tuple( sid, calls, regions, aber_file )
                    }
                GXCNV_COMPARE( ch_compare_input, analysisdir )
                ch_gxcnv_comparison = GXCNV_COMPARE.out.comparison
            }
        }

        // ── gxcnv2 (WCX-core, reuses RUN_WCX BED outputs — WCX result = gxcnv2) ──────
        // Wired directly to RUN_WCX; requires run_wcx=true.
        ch_gxcnv2_calls = Channel.empty()

        if ( run_gxcnv2 ) {
            if ( !params.run_wcx ) {
                log.warn "[gxcnv2] Skipped: run_wcx=false (gxcnv2 requires WCX BED outputs)"
            } else {
                ch_gxcnv2_input = RUN_WCX.out.wcx_beds
                    .map { sid, grp, bins_bed, segs_bed, aber_bed ->
                        tuple("${sid}_${grp}", bins_bed, segs_bed, aber_bed)
                    }

                GXCNV2_PREDICT(ch_gxcnv2_input, analysisdir)

                ch_gxcnv2_calls = GXCNV2_PREDICT.out.calls_tsv
                    .map { composite_id, calls ->
                        def parts = composite_id.toString().tokenize('_')
                        def grp   = parts[-1]
                        def sid   = parts[0..-2].join('_')
                        tuple(sid, grp, calls)
                    }

                ch_plot2_input = GXCNV2_PREDICT.out.bins_tsv
                    .join( GXCNV2_PREDICT.out.calls_tsv )
                GXCNV2_PLOT( ch_plot2_input, analysisdir )
            }
        }

        // ── gxcnv1 (WC-core, reuses RUN_WC NPZ outputs — WC result = gxcnv1) ───────
        // Wired directly to RUN_WC; requires run_wc=true.
        ch_gxcnv1_calls = Channel.empty()

        if ( run_gxcnv1 ) {
            if ( !params.run_wc ) {
                log.warn "[gxcnv1] Skipped: run_wc=false (gxcnv1 requires WC NPZ outputs)"
            } else {
                ch_gxcnv1_input = RUN_WC.out.wc_sample_npz
                    .join(RUN_WC.out.wc_out_npz, by: [0, 1])
                    .map { sid, grp, s_npz, out_npz ->
                        tuple("${sid}_${grp}", s_npz, out_npz)
                    }

                GXCNV1_PREDICT(ch_gxcnv1_input, analysisdir)

                ch_gxcnv1_calls = GXCNV1_PREDICT.out.calls_tsv
                    .map { composite_id, calls ->
                        def parts = composite_id.toString().tokenize('_')
                        def grp   = parts[-1]
                        def sid   = parts[0..-2].join('_')
                        tuple(sid, grp, calls)
                    }

                ch_plot1_input = GXCNV1_PREDICT.out.bins_tsv
                    .join( GXCNV1_PREDICT.out.calls_tsv )
                GXCNV1_PLOT( ch_plot1_input, analysisdir )
            }
        }

    emit:
        // Each element: (sample, group, result_path) — WC ∪ WCX
        wc_result        = ch_wc_result.mix( ch_wcx_result )
        gxcnv_calls      = ch_gxcnv_calls
        gxcnv_comparison = ch_gxcnv_comparison
        gxcnv2_calls     = ch_gxcnv2_calls
        gxcnv1_calls     = ch_gxcnv1_calls
}
