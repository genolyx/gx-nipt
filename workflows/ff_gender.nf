/*
 * =========================================================
 *  Fetal Fraction & Gender Detection Workflow
 *
 *  FF estimation stack:
 *    YFF1        – Y-chromosome coverage ratio (gd_1)
 *    YFF2        – Adjusted YFF via wig normalisation (gd_2)
 *    seqFF       – Sequence-based FF from GC-corrected coverage
 *    Fragment FF – Short/long fragment ratio (gd_4)
 *    gx-FF       – LightGBM + DNN ensemble (NEW)
 *
 *  Final FF = GXFF_ENSEMBLE(seqFF, gx-FF)
 *    - Low-FF (<5%): gx-FF only (better trained for low-FF regime)
 *    - Normal:       weighted ensemble (gx-FF 60%, seqFF 40%)
 *    - gx-FF QC fail: fall back to seqFF
 *
 *  Gender decision uses FF_FINAL from ensemble.
 * =========================================================
 */

include { CALCULATE_YFF }        from '../modules/ff_gender'
include { CALCULATE_YFF2 }       from '../modules/ff_gender'
include { CALCULATE_SEQFF }      from '../modules/ff_gender'
include { CALCULATE_FRAGMENT_FF} from '../modules/ff_gender'
include { GENDER_DECISION }      from '../modules/ff_gender'
include { GXFF_PREDICT }         from '../modules/gxff'
include { GXFF_ENSEMBLE }        from '../modules/gxff'

workflow FF_GENDER_WORKFLOW {
    take:
        sample_name   // val: sample identifier
        ch_bam        // tuple: (sample_id, bam, bai)
        ch_bincount   // path: 50kb bin count file from HMMcopy (or NO_FILE)
        ch_config     // path: pipeline_config.json
        labcode       // val: lab identifier
        analysisdir   // val: output directory root
        gxff_model    // path: trained gx-FF model .pkl (or NO_FILE to skip)

    main:
        // ── YFF1: Y-chromosome coverage ratio ────────────────────────────────
        CALCULATE_YFF(
            sample_name,
            ch_bam,
            ch_config,
            labcode,
            analysisdir
        )

        // ── YFF2: Adjusted YFF using wig normalisation ────────────────────────
        CALCULATE_YFF2(
            sample_name,
            ch_config,
            labcode,
            analysisdir
        )

        // ── seqFF: Sequence-based FF ──────────────────────────────────────────
        CALCULATE_SEQFF(
            sample_name,
            ch_bam,
            ch_config,
            labcode,
            analysisdir
        )

        // ── Fragment-size based FF ────────────────────────────────────────────
        CALCULATE_FRAGMENT_FF(
            sample_name,
            ch_bam,
            ch_config,
            analysisdir
        )

        // ── gx-FF: LightGBM + DNN ensemble ───────────────────────────────────
        // Only runs when a trained model file is provided.
        // Falls back gracefully to seqFF-only if model is absent.
        ch_gxff_tsv = Channel.empty()

        if ( gxff_model.name != 'NO_FILE' ) {
            GXFF_PREDICT(
                ch_bam,
                ch_bincount,
                gxff_model
            )
            ch_gxff_tsv = GXFF_PREDICT.out.gxff_tsv
        }

        // ── Ensemble: combine seqFF + gx-FF ──────────────────────────────────
        // Join seqFF and gx-FF outputs on sample_id.
        // If gx-FF was skipped, the join produces an empty channel and
        // GXFF_ENSEMBLE falls back to seqFF-only mode.
        ch_seqff_tsv = CALCULATE_SEQFF.out.seqff_txt

        ch_ensemble_input = ch_seqff_tsv
            .join( ch_gxff_tsv.ifEmpty { tuple( sample_name, file('NO_FILE') ) }, remainder: true )
            .map { sid, seqff, gxff ->
                def gxff_file = (gxff == null || gxff.name == 'NO_FILE') ? file('NO_FILE') : gxff
                tuple( sid, seqff, gxff_file )
            }

        GXFF_ENSEMBLE( ch_ensemble_input )

        // ── Final gender decision ─────────────────────────────────────────────
        // Passes FF_FINAL from the ensemble to the gender decision step.
        GENDER_DECISION(
            sample_name,
            CALCULATE_YFF.out.yff_txt,
            CALCULATE_YFF2.out.yff2_txt,
            GXFF_ENSEMBLE.out.ff_tsv,    // replaces raw seqFF_txt
            CALCULATE_FRAGMENT_FF.out.frag_ff_txt,
            ch_config,
            analysisdir
        )

    emit:
        ff_result    = GENDER_DECISION.out.ff_result
        gender_txt   = GENDER_DECISION.out.gender_txt
        ff_ensemble  = GXFF_ENSEMBLE.out.ff_tsv   // for downstream logging/reporting
}
