/*
 * =========================================================
 *  Alignment Workflow
 *  FASTQ → Downsample → BWA-MEM2 → Sort → Dedup → Proper-paired BAM
 * =========================================================
 */

include { DOWNSAMPLE_FASTQ }         from '../modules/downsample'
include { BWA_ALIGN }                from '../modules/bwa'
include { SAMTOOLS_SORT_INDEX }      from '../modules/samtools'
include { PICARD_MARKDUP }           from '../modules/picard'
include { SAMTOOLS_FILTER_UNIQUE }   from '../modules/samtools'
include { SAMTOOLS_PROPER_PAIRED }   from '../modules/samtools'

workflow ALIGN_WORKFLOW {
    take:
        sample_name   // string
        ch_fastq      // channel: [r1, r2] or empty
        ch_bam        // channel: existing proper_paired.bam or empty
        ch_config     // channel: pipeline_config.json
        labcode       // string
        analysisdir   // string

    main:
        // ── Case 1: Start from existing proper_paired.bam ─
        if (params.from_bam) {
            ch_proper_bam = ch_bam
        } else {
            // ── Case 2: Full alignment pipeline ───────────
            // Step 1: Downsample FASTQ (if over max_fq_size)
            DOWNSAMPLE_FASTQ(
                sample_name,
                ch_fastq,
                ch_config
            )
            ch_ds_fastq = DOWNSAMPLE_FASTQ.out.fastq

            // Step 2: BWA-MEM2 alignment
            BWA_ALIGN(
                sample_name,
                ch_ds_fastq,
                ch_config,
                labcode,
                analysisdir
            )
            ch_raw_bam = BWA_ALIGN.out.bam

            // Step 3: Sort BAM
            SAMTOOLS_SORT_INDEX(
                sample_name,
                ch_raw_bam,
                'sorted',
                analysisdir
            )
            ch_sorted_bam = SAMTOOLS_SORT_INDEX.out.bam

            // Step 4: Picard MarkDuplicates
            PICARD_MARKDUP(
                sample_name,
                ch_sorted_bam,
                ch_config,
                analysisdir
            )
            ch_dedup_bam = PICARD_MARKDUP.out.bam

            // Step 5: Filter unique reads
            SAMTOOLS_FILTER_UNIQUE(
                sample_name,
                ch_dedup_bam,
                analysisdir
            )
            ch_unique_bam = SAMTOOLS_FILTER_UNIQUE.out.bam

            // Step 6: Filter proper-paired reads
            SAMTOOLS_PROPER_PAIRED(
                sample_name,
                ch_unique_bam,
                analysisdir
            )
            ch_proper_bam = SAMTOOLS_PROPER_PAIRED.out.bam
        }

    emit:
        proper_bam = ch_proper_bam
}
