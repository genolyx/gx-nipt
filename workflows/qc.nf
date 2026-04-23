/*
 * =========================================================
 *  QC Workflow
 *  FastQC → Qualimap → QC Filter
 * =========================================================
 */

include { FASTQC }           from '../modules/fastqc'
include { QUALIMAP }         from '../modules/qualimap'
include { QC_FILTER }        from '../modules/qualimap'

workflow QC_WORKFLOW {
    take:
        sample_name   // string
        ch_fastq      // channel: [r1, r2]
        ch_bam        // channel: proper_paired.bam
        ch_config     // channel: pipeline_config.json
        analysisdir   // string

    main:
        // FastQC runs in parallel with BAM generation
        if (!params.algorithm_only && !params.from_bam) {
            FASTQC(
                sample_name,
                ch_fastq,
                analysisdir
            )
        }

        // Qualimap on final proper_paired BAM
        QUALIMAP(
            sample_name,
            ch_bam,
            ch_config,
            analysisdir
        )

        // QC filter: check thresholds and write QC pass/fail
        QC_FILTER(
            sample_name,
            QUALIMAP.out.qc_txt,
            ch_config,
            analysisdir
        )

    emit:
        qc_result = QC_FILTER.out.qc_result
        qc_txt    = QUALIMAP.out.qc_txt
}
