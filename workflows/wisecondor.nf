/*
 * =========================================================
 *  Wisecondor / gx-cnv Workflow
 *
 *  CNV detection stack (parallel execution):
 *    WC   – Wisecondor (legacy, per group: orig/fetus/mom)
 *    WCX  – WisecondorX (per group: orig/fetus/mom)
 *    gx-cnv – Hybrid dual-track (Track A Z-score CBS + Track B PCA Mahalanobis)
 *
 *  Validation mode:
 *    When gx-cnv reference is provided, gx-cnv runs in parallel with WCX.
 *    GXCNV_COMPARE produces a per-sample concordance report.
 *    Once gx-cnv is validated, WCX can be removed by setting
 *    params.run_wcx = false in nextflow.config.
 *
 *  Output:
 *    wc_result   – WC/WCX results (existing downstream consumers unchanged)
 *    gxcnv_calls – gx-cnv HIGH_RISK calls TSV
 *    gxcnv_comparison – concordance report vs WCX
 * =========================================================
 */

include { RUN_WC }          from '../modules/wisecondor'
include { RUN_WCX }         from '../modules/wisecondor'
include { GXCNV_CONVERT }   from '../modules/gxcnv'
include { GXCNV_PREDICT }   from '../modules/gxcnv'
include { GXCNV_COMPARE }   from '../modules/gxcnv'

workflow WC_WORKFLOW {
    take:
        sample_name      // val: sample identifier
        ch_bam           // tuple: (sample_id, bam, bai)
        ch_ff_result     // path: fetal_fraction result (FF_FINAL from ensemble)
        ch_config        // path: pipeline_config.json
        labcode          // val: lab identifier
        analysisdir      // val: output directory root
        gxcnv_reference  // path: pre-built gx-cnv reference .npz (or NO_FILE to skip)
        gxcnv_bin_size   // val: bin size for gxcnv convert (default: 100000)
        gxcnv_thresh_z   // val: Track A Z-score threshold (default: -3.0)
        gxcnv_thresh_p   // val: Track B p-value threshold (default: 0.05)

    main:
        ch_groups = Channel.of('orig', 'fetus', 'mom')

        // ── WC (Wisecondor) ───────────────────────────────────────────────────
        RUN_WC(
            sample_name,
            ch_groups,
            ch_bam,
            ch_ff_result,
            ch_config,
            labcode,
            analysisdir
        )

        // ── WCX (WisecondorX) ─────────────────────────────────────────────────
        // Controlled by params.run_wcx (default: true).
        // Set params.run_wcx = false once gx-cnv is validated.
        ch_wcx_result = Channel.empty()

        if ( params.run_wcx ) {
            RUN_WCX(
                sample_name,
                ch_groups,
                ch_bam,
                ch_ff_result,
                ch_config,
                labcode,
                analysisdir
            )
            ch_wcx_result = RUN_WCX.out.wcx_result
        }

        // ── gx-cnv (parallel, validation mode) ───────────────────────────────
        // Only runs when a pre-built reference NPZ is provided.
        ch_gxcnv_calls      = Channel.empty()
        ch_gxcnv_comparison = Channel.empty()

        if ( gxcnv_reference.name != 'NO_FILE' ) {

            // Step 1: BAM → NPZ
            GXCNV_CONVERT(
                ch_bam,
                gxcnv_bin_size,
                file('NO_FILE')   // no blacklist by default; override via params
            )

            // Step 2: Extract FF_FINAL for passing to gxcnv predict
            // ff_result TSV has column FF_FINAL; extract as a value
            ch_ff_val = ch_ff_result
                .map { sid, tsv ->
                    def ff = "NA"
                    try {
                        tsv.withReader { r ->
                            def header = null
                            r.eachLine { line ->
                                if ( line.startsWith("#") ) return
                                if ( header == null ) { header = line.split("\\t"); return }
                                def cols = line.split("\\t")
                                def idx  = header.findIndexOf { it == "FF_FINAL" }
                                if ( idx >= 0 ) ff = cols[idx]
                            }
                        }
                    } catch (e) { /* keep NA */ }
                    tuple( sid, ff )
                }

            // Step 3: Predict CNVs
            ch_predict_input = GXCNV_CONVERT.out.npz
                .join( ch_ff_val, remainder: true )
                .map { sid, npz, ff -> tuple( sid, npz, ff ?: "NA" ) }

            GXCNV_PREDICT(
                ch_predict_input.map { sid, npz, ff -> tuple(sid, npz) },
                gxcnv_reference,
                gxcnv_thresh_z,
                gxcnv_thresh_p,
                ch_predict_input.map { sid, npz, ff -> ff }
            )

            ch_gxcnv_calls = GXCNV_PREDICT.out.calls_tsv

            // Step 4: Compare gx-cnv vs WCX (only when WCX is also running)
            if ( params.run_wcx ) {
                // WCX aberrations file: assume it is published under analysisdir
                // Construct expected path from sample_name and analysisdir
                ch_wcx_aberrations = ch_wcx_result
                    .map { sid, result_dir ->
                        def aber = file("${result_dir}/aberrations.bed")
                        tuple( sid, aber )
                    }

                ch_compare_input = GXCNV_PREDICT.out.calls_tsv
                    .join( GXCNV_PREDICT.out.regions_tsv )
                    .join( ch_wcx_aberrations, remainder: true )
                    .map { sid, calls, regions, aber ->
                        def aber_file = (aber == null) ? file('NO_FILE') : aber
                        tuple( sid, calls, regions, aber_file )
                    }

                GXCNV_COMPARE( ch_compare_input )
                ch_gxcnv_comparison = GXCNV_COMPARE.out.comparison
            }
        }

    emit:
        wc_result          = RUN_WC.out.wc_result.mix( ch_wcx_result )
        gxcnv_calls        = ch_gxcnv_calls
        gxcnv_comparison   = ch_gxcnv_comparison
}
