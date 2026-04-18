#!/usr/bin/env nextflow
/*
 * =========================================================
 *  gx-nipt: Genolyx NIPT Analysis Pipeline (Nextflow)
 * =========================================================
 *  Author  : Hyukjung Kwon
 *  Version : 1.1.0
 *  License : MIT
 * =========================================================
 *
 * Pipeline overview:
 *
 *   FASTQ Input
 *       │
 *       ├─► [QC]  FastQC ──────────────────────────────────────────────────────►  Output_QC
 *       │
 *       ├─► [ALIGN] Downsample → BWA-MEM2 → Sort → Dedup → Proper-paired BAM
 *       │               │
 *       │               ├─► [QC] Qualimap → QC Filter ──────────────────────────►  Output_QC
 *       │               │
 *       │               ├─► [FF/Gender] YFF + YFF2 + seqFF + Fragment FF
 *       │               │              + gx-FF (LightGBM+DNN ensemble) [NEW]
 *       │               │              → Ensemble FF → Gender Decision ──────────►  Output_FF
 *       │               │
 *       │               ├─► [HMMcopy] readCounter (50kb + 10mb) → HMMcopy.R ───►  Output_hmmcopy
 *       │               │       │
 *       │               │       ├─► [EZD] orig / fetus / mom ──────────────────►  Output_EZD
 *       │               │       └─► [PRIZM] orig / fetus / mom ────────────────►  Output_PRIZM
 *       │               │
 *       │               ├─► [WC/WCX] Wisecondor + WisecondorX ───────────────►  Output_WC / Output_WCX
 *       │               │   [gx-cnv] Hybrid dual-track (parallel) [NEW] ───────►  Output_gxcnv
 *       │               │   [Compare] gx-cnv vs WCX concordance report [NEW] ──►  Output_cnv_comparison
 *       │               │
 *       │               └─► [MD] Microdeletion detection ─────────────────────►  Output_MD
 *       │
 *       └─► [REPORT] JSON + HTML ────────────────────────────────────────────►  Output_Result
 *
 * New in v1.1.0:
 *   - gx-FF: LightGBM + DNN ensemble FF estimation (parallel with seqFF)
 *     Set --gxff_model /path/to/model.pkl to enable.
 *     Falls back to seqFF-only when model is absent.
 *
 *   - gx-cnv: Hybrid dual-track CNV detection (parallel with WisecondorX)
 *     Set --gxcnv_reference /path/to/reference.npz to enable.
 *     Produces concordance report vs WCX for validation.
 *     Once validated, disable WCX with --run_wcx false.
 */

nextflow.enable.dsl = 2

// ─────────────────────────────────────────────────────────
// Import sub-workflows
// ─────────────────────────────────────────────────────────
include { ALIGN_WORKFLOW }      from './workflows/align'
include { QC_WORKFLOW }         from './workflows/qc'
include { FF_GENDER_WORKFLOW }  from './workflows/ff_gender'
include { HMMCOPY_WORKFLOW }    from './workflows/hmmcopy'
include { EZD_WORKFLOW }        from './workflows/ezd'
include { PRIZM_WORKFLOW }      from './workflows/prizm'
include { WC_WORKFLOW }         from './workflows/wisecondor'
include { MD_WORKFLOW }         from './workflows/microdeletion'
include { REPORT_WORKFLOW }     from './workflows/report'

// ─────────────────────────────────────────────────────────
// Parameter defaults (overridden by nextflow.config or CLI)
// ─────────────────────────────────────────────────────────
params.sample_name  = null
params.fastq_r1     = null
params.fastq_r2     = null
params.labcode      = null
params.age          = null
params.from_bam     = null   // optional: start from proper_paired.bam

params.root_dir     = null   // host root directory (e.g. /data/nipt)
params.work_dir     = null   // batch work directory (e.g. 250430_01)

params.outdir       = null   // resolved at runtime
params.analysisdir  = null   // resolved at runtime

params.algorithm_only = false
params.force          = false

// ── gx-FF parameters ─────────────────────────────────────
// Set to a .pkl model path to enable gx-FF.
// Leave as null to use seqFF-only (legacy behaviour).
params.gxff_model   = null   // e.g. /data/nipt/models/gxff_v1.pkl

// ── gx-cnv parameters ────────────────────────────────────
// Set to a .npz reference path to enable gx-cnv in parallel with WCX.
// Leave as null to skip gx-cnv entirely.
params.gxcnv_reference = null   // e.g. /data/nipt/refs/gxcnv_ref.npz
params.gxcnv_bin_size  = 100000 // bin size for gxcnv convert
params.gxcnv_thresh_z  = -3.0   // Track A Z-score threshold
params.gxcnv_thresh_p  = 0.05   // Track B p-value threshold

// When true, WisecondorX runs alongside gx-cnv for concordance validation.
// Set to false once gx-cnv is validated and WCX is no longer needed.
params.run_wcx = true

// ─────────────────────────────────────────────────────────
// Validate required parameters
// ─────────────────────────────────────────────────────────
def validate_params() {
    def errors = []
    if (!params.sample_name) errors << "  --sample_name is required"
    if (!params.labcode)     errors << "  --labcode is required"
    if (!params.root_dir)    errors << "  --root_dir is required"
    if (!params.work_dir)    errors << "  --work_dir is required"

    if (!params.from_bam && !params.algorithm_only) {
        if (!params.fastq_r1) errors << "  --fastq_r1 is required (or use --from_bam)"
        if (!params.fastq_r2) errors << "  --fastq_r2 is required (or use --from_bam)"
        if (!params.age)      errors << "  --age is required"
    }

    if (errors) {
        log.error "Parameter validation failed:\n" + errors.join("\n")
        System.exit(1)
    }
}

// ─────────────────────────────────────────────────────────
// Main workflow
// ─────────────────────────────────────────────────────────
workflow {
    validate_params()

    def sample_name  = params.sample_name
    def labcode      = params.labcode
    def root_dir     = params.root_dir
    def work_dir     = params.work_dir

    // Resolve output/analysis directories
    def outdir      = params.outdir      ?: "${root_dir}/output/${work_dir}/${sample_name}"
    def analysisdir = params.analysisdir ?: "${root_dir}/analysis/${work_dir}/${sample_name}"

    // ── Channel creation ──────────────────────────────────
    ch_config = Channel.fromPath("${root_dir}/config/${labcode}/pipeline_config.json", checkIfExists: true)

    // FASTQ or BAM input
    if (params.from_bam) {
        ch_bam   = Channel.fromPath(params.from_bam, checkIfExists: true)
        ch_fastq = Channel.empty()
    } else {
        def fastq_dir = "${root_dir}/fastq/${work_dir}/${sample_name}"
        ch_fastq = Channel.of([
            file("${fastq_dir}/${params.fastq_r1}", checkIfExists: true),
            file("${fastq_dir}/${params.fastq_r2}", checkIfExists: true)
        ])
        ch_bam = Channel.empty()
    }

    // ── Optional tool paths ───────────────────────────────
    // gx-FF model: use NO_FILE sentinel when not provided
    def gxff_model_path = params.gxff_model
        ? file(params.gxff_model, checkIfExists: true)
        : file('NO_FILE')

    // gx-cnv reference: use NO_FILE sentinel when not provided
    def gxcnv_ref_path = params.gxcnv_reference
        ? file(params.gxcnv_reference, checkIfExists: true)
        : file('NO_FILE')

    // ── Alignment & BAM generation ────────────────────────
    if (!params.algorithm_only) {
        ALIGN_WORKFLOW(
            sample_name,
            ch_fastq,
            ch_bam,
            ch_config,
            labcode,
            analysisdir
        )
        ch_proper_bam = ALIGN_WORKFLOW.out.proper_bam

        // ── QC (parallel with alignment output) ───────────
        QC_WORKFLOW(
            sample_name,
            ch_fastq,
            ch_proper_bam,
            ch_config,
            analysisdir
        )
        ch_qc_passed = QC_WORKFLOW.out.qc_result
    } else {
        // algorithm_only: skip alignment, use existing BAM
        ch_proper_bam = Channel.fromPath(
            "${analysisdir}/${sample_name}.proper_paired.bam",
            checkIfExists: true
        )
        ch_qc_passed = Channel.value(true)
    }

    // ── HMMcopy (parallel with FF/Gender) ─────────────────
    // Run first so bincount is available for gx-FF
    HMMCOPY_WORKFLOW(
        sample_name,
        ch_proper_bam,
        labcode,
        analysisdir
    )
    ch_norm_50kb  = HMMCOPY_WORKFLOW.out.norm_50kb
    ch_norm_10mb  = HMMCOPY_WORKFLOW.out.norm_10mb
    ch_count_10mb = HMMCOPY_WORKFLOW.out.count_10mb

    // The 50kb bin count file (HMMcopy readCounter output) is passed to
    // gx-FF so it can skip re-counting and use the already-computed bins.
    ch_bincount_50kb = HMMCOPY_WORKFLOW.out.bincount_50kb  // path or NO_FILE channel

    // ── Fetal Fraction & Gender ───────────────────────────
    FF_GENDER_WORKFLOW(
        sample_name,
        ch_proper_bam,
        ch_bincount_50kb,
        ch_config,
        labcode,
        analysisdir,
        gxff_model_path
    )
    ch_ff_result = FF_GENDER_WORKFLOW.out.ff_result

    // Log which FF method was used
    FF_GENDER_WORKFLOW.out.ff_ensemble
        .subscribe { sid, tsv ->
            log.info "[FF] Sample ${sid}: ensemble FF result -> ${tsv}"
        }

    // ── EZD (depends on HMMcopy 50kb) ─────────────────────
    EZD_WORKFLOW(
        sample_name,
        ch_norm_50kb,
        ch_ff_result,
        ch_config,
        labcode,
        analysisdir
    )
    ch_ezd_result = EZD_WORKFLOW.out.ezd_result

    // ── PRIZM (depends on HMMcopy 10mb count) ─────────────
    PRIZM_WORKFLOW(
        sample_name,
        ch_count_10mb,
        ch_ff_result,
        ch_config,
        labcode,
        analysisdir
    )
    ch_prizm_result = PRIZM_WORKFLOW.out.prizm_result

    // ── Wisecondor + gx-cnv (parallel) ────────────────────
    WC_WORKFLOW(
        sample_name,
        ch_proper_bam,
        ch_ff_result,
        ch_config,
        labcode,
        analysisdir,
        gxcnv_ref_path,
        params.gxcnv_bin_size,
        params.gxcnv_thresh_z,
        params.gxcnv_thresh_p
    )
    ch_wc_result         = WC_WORKFLOW.out.wc_result
    ch_gxcnv_calls       = WC_WORKFLOW.out.gxcnv_calls
    ch_gxcnv_comparison  = WC_WORKFLOW.out.gxcnv_comparison

    // Log gx-cnv concordance when available
    ch_gxcnv_comparison
        .subscribe { sid, tsv ->
            log.info "[gx-cnv] Sample ${sid}: concordance report -> ${tsv}"
        }

    // ── Microdeletion ──────────────────────────────────────
    MD_WORKFLOW(
        sample_name,
        ch_wc_result,
        ch_config,
        labcode,
        analysisdir
    )
    ch_md_result = MD_WORKFLOW.out.md_result

    // ── Final Report (JSON + HTML) ─────────────────────────
    REPORT_WORKFLOW(
        sample_name,
        ch_ezd_result,
        ch_prizm_result,
        ch_ff_result,
        ch_md_result,
        ch_qc_passed,
        ch_config,
        labcode,
        analysisdir,
        outdir
    )
}

// ─────────────────────────────────────────────────────────
// Workflow completion handler
// ─────────────────────────────────────────────────────────
workflow.onComplete {
    def status   = workflow.success ? "SUCCESS" : "FAILED"
    def duration = workflow.duration
    def gxff_status  = params.gxff_model       ? "ENABLED (${params.gxff_model})"      : "DISABLED (seqFF only)"
    def gxcnv_status = params.gxcnv_reference  ? "ENABLED (${params.gxcnv_reference})" : "DISABLED"
    def wcx_status   = params.run_wcx          ? "ENABLED" : "DISABLED"

    log.info """
    ╔══════════════════════════════════════════════════════╗
    ║  gx-nipt Pipeline Completed                         ║
    ║  Sample    : ${params.sample_name}
    ║  Status    : ${status}
    ║  Duration  : ${duration}
    ║  gx-FF     : ${gxff_status}
    ║  gx-cnv    : ${gxcnv_status}
    ║  WCX       : ${wcx_status}
    ║  Work Dir  : ${workflow.workDir}
    ╚══════════════════════════════════════════════════════╝
    """.stripIndent()
}

workflow.onError {
    log.error "Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}
