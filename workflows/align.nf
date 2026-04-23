/*
 * =========================================================
 *  Alignment Workflow
 *  FASTQ → Downsample → BWA-MEM2 → Sort → Dedup → Proper-paired BAM
 *
 *  SSD Scratch Strategy (Strategy B):
 *    params.use_ssd = true  → 모든 BAM을 SSD에서 생성 후
 *                             proper_paired.bam만 HDD로 이동
 *    params.use_ssd = false → 기존 방식 (모두 HDD)
 * =========================================================
 */

include { DOWNSAMPLE_FASTQ }         from '../modules/downsample'
include { BWA_ALIGN }                from '../modules/bwa'
include { SAMTOOLS_SORT_INDEX }      from '../modules/samtools'
include { PICARD_MARKDUP }           from '../modules/picard'
include { SAMTOOLS_FILTER_UNIQUE }   from '../modules/samtools'
include { SAMTOOLS_PROPER_PAIRED }   from '../modules/samtools'
include { SAMTOOLS_SPLIT_FETUS_MOM } from '../modules/samtools'
include { SCRATCH_SETUP;
          SCRATCH_MOVE_FINAL }       from '../modules/scratch'

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
            ch_proper_bam     = ch_bam
            // When re-entering from an existing BAM, the index may or may
            // not be alongside. Try to pick it up; SPLIT process below
            // depends on .bai being present.
            ch_proper_bai = ch_bam.map { b ->
                def idx = file("${b}.bai")
                idx.exists() ? idx : null
            }
            ch_use_ssd_result = false

        } else {
            // ── Case 2: Full alignment pipeline ───────────

            // ── SSD Scratch 초기화 (use_ssd=true 시) ──────
            def scratch_dir = params.use_ssd
                ? (params.scratch_dir ?: '/tmp/nipt_scratch')
                : null

            if (params.use_ssd) {
                log.info "[SSD] Strategy B enabled. scratch_dir = ${scratch_dir}"
                SCRATCH_SETUP(sample_name, scratch_dir)
                ch_sample_scratch = SCRATCH_SETUP.out.sample_scratch_dir
            } else {
                log.info "[SSD] Strategy B disabled. Using default HDD work directory."
                ch_sample_scratch = Channel.value(null)
            }

            // ── Step 1: Downsample FASTQ ───────────────────
            DOWNSAMPLE_FASTQ(
                sample_name,
                ch_fastq,
                ch_config
            )
            ch_ds_fastq = DOWNSAMPLE_FASTQ.out.fastq

            // ── Step 2: BWA-MEM2 alignment ─────────────────
            // use_ssd=true 시 Nextflow workDir이 SSD를 가리키므로
            // 별도 경로 지정 없이 자동으로 SSD에 기록됨
            BWA_ALIGN(
                sample_name,
                ch_ds_fastq,
                ch_config,
                labcode,
                analysisdir
            )
            ch_raw_bam = BWA_ALIGN.out.bam

            // ── Step 3: Sort BAM ───────────────────────────
            SAMTOOLS_SORT_INDEX(
                sample_name,
                ch_raw_bam,
                'sorted',
                analysisdir
            )
            ch_sorted_bam = SAMTOOLS_SORT_INDEX.out.bam

            // ── Step 4: Picard MarkDuplicates ──────────────
            PICARD_MARKDUP(
                sample_name,
                ch_sorted_bam,
                ch_config,
                analysisdir
            )
            ch_dedup_bam = PICARD_MARKDUP.out.bam

            // ── Step 5: Filter unique reads ────────────────
            SAMTOOLS_FILTER_UNIQUE(
                sample_name,
                ch_dedup_bam,
                analysisdir
            )
            ch_unique_bam = SAMTOOLS_FILTER_UNIQUE.out.bam

            // ── Step 6: Filter proper-paired reads ─────────
            SAMTOOLS_PROPER_PAIRED(
                sample_name,
                ch_unique_bam,
                analysisdir
            )
            ch_proper_bam_ssd = SAMTOOLS_PROPER_PAIRED.out.bam
            ch_proper_bai_ssd = SAMTOOLS_PROPER_PAIRED.out.bai

            // ── Step 7: SSD → HDD 이동 (use_ssd=true 시) ──
            // proper_paired.bam만 HDD analysis 디렉토리로 복사 후
            // SSD의 임시 BAM 전체 삭제
            if (params.use_ssd) {
                SCRATCH_MOVE_FINAL(
                    sample_name,
                    ch_proper_bam_ssd,
                    ch_proper_bai_ssd,
                    analysisdir,
                    ch_sample_scratch
                )
                ch_proper_bam = SCRATCH_MOVE_FINAL.out.bam
                ch_proper_bai = SCRATCH_MOVE_FINAL.out.bai
            } else {
                ch_proper_bam = ch_proper_bam_ssd
                ch_proper_bai = ch_proper_bai_ssd
            }
        }

        // ── Step 8: TLEN-based split (orig / fetus / mom) ─────
        // proper_paired.bam → of_orig.bam + of_fetus.bam + of_mom.bam
        // Downstream HMMcopy / PRIZM / WC(X) expect all three.
        SAMTOOLS_SPLIT_FETUS_MOM(
            sample_name,
            ch_proper_bam,
            ch_proper_bai,
            analysisdir
        )

    emit:
        proper_bam = ch_proper_bam
        bam_trio   = SAMTOOLS_SPLIT_FETUS_MOM.out.trio
}
