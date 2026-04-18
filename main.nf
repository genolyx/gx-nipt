#!/usr/bin/env nextflow
/*
 * =========================================================
 *  gx-nipt: Genolyx NIPT Analysis Pipeline (Nextflow)
 * =========================================================
 *  Author  : Hyukjung Kwon
 *  Version : 1.0.0
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
 *       │               ├─► [FF/Gender] YFF + YFF2 + SeqFF + Fragment FF ───────►  Output_FF
 *       │               │
 *       │               ├─► [HMMcopy] readCounter (50kb + 10mb) → HMMcopy.R ───►  Output_hmmcopy
 *       │               │       │
 *       │               │       ├─► [EZD] orig / fetus / mom ──────────────────►  Output_EZD
 *       │               │       └─► [PRIZM] orig / fetus / mom ────────────────►  Output_PRIZM
 *       │               │
 *       │               ├─► [WC/WCX] Wisecondor + WisecondorX ───────────────►  Output_WC / Output_WCX
 *       │               │
 *       │               └─► [MD] Microdeletion detection ─────────────────────►  Output_MD
 *       │
 *       └─► [REPORT] JSON + HTML ────────────────────────────────────────────►  Output_Result
 *
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
    // Config file channel
    ch_config = Channel.fromPath("${root_dir}/config/${labcode}/pipeline_config.json", checkIfExists: true)

    // FASTQ or BAM input
    if (params.from_bam) {
        ch_bam = Channel.fromPath(params.from_bam, checkIfExists: true)
        ch_fastq = Channel.empty()
    } else {
        def fastq_dir = "${root_dir}/fastq/${work_dir}/${sample_name}"
        ch_fastq = Channel.of([
            file("${fastq_dir}/${params.fastq_r1}", checkIfExists: true),
            file("${fastq_dir}/${params.fastq_r2}", checkIfExists: true)
        ])
        ch_bam = Channel.empty()
    }

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

    // ── Fetal Fraction & Gender (parallel) ────────────────
    FF_GENDER_WORKFLOW(
        sample_name,
        ch_proper_bam,
        ch_config,
        labcode,
        analysisdir
    )
    ch_ff_result = FF_GENDER_WORKFLOW.out.ff_result

    // ── HMMcopy (parallel) ────────────────────────────────
    HMMCOPY_WORKFLOW(
        sample_name,
        ch_proper_bam,
        labcode,
        analysisdir
    )
    ch_norm_50kb = HMMCOPY_WORKFLOW.out.norm_50kb
    ch_norm_10mb = HMMCOPY_WORKFLOW.out.norm_10mb
    ch_count_10mb = HMMCOPY_WORKFLOW.out.count_10mb

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

    // ── Wisecondor (parallel with EZD/PRIZM) ──────────────
    WC_WORKFLOW(
        sample_name,
        ch_proper_bam,
        ch_ff_result,
        ch_config,
        labcode,
        analysisdir
    )
    ch_wc_result = WC_WORKFLOW.out.wc_result

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
    def status = workflow.success ? "SUCCESS" : "FAILED"
    def duration = workflow.duration
    log.info """
    ╔══════════════════════════════════════════════════════╗
    ║  gx-nipt Pipeline Completed                         ║
    ║  Sample  : ${params.sample_name}
    ║  Status  : ${status}
    ║  Duration: ${duration}
    ║  Work Dir: ${workflow.workDir}
    ╚══════════════════════════════════════════════════════╝
    """.stripIndent()
}

workflow.onError {
    log.error "Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}
