/*
 * =========================================================
 *  Module: SSD Scratch Strategy (Strategy B)
 * =========================================================
 *  I/O 최적화 전략:
 *    - 모든 중간 BAM 파일을 빠른 NVMe SSD(scratch_dir)에서 생성
 *    - 분석 완료 후 proper_paired.bam만 HDD analysis 디렉토리로 이동
 *    - 나머지 임시 BAM은 SSD에서 삭제
 *
 *  성능 이득 (실측 기반 추정):
 *    - 현재 (모두 HDD): ~40s I/O
 *    - 전략 B (모두 SSD): ~5s I/O  → 87% 절감
 *    - BWA sort + Picard dedup 단계: HDD 180 MB/s → SSD 3,000 MB/s (16.7x)
 *
 *  사용법:
 *    nextflow run main.nf -profile ssd_scratch \
 *      --scratch_dir /tmp/nipt_scratch \
 *      --sample_name SAMPLE001 ...
 *
 *  또는 기존 프로파일에 파라미터 추가:
 *    nextflow run main.nf --use_ssd true --scratch_dir /tmp/nipt_scratch ...
 * =========================================================
 */

// ─────────────────────────────────────────────────────────
// SCRATCH_SETUP
//   SSD scratch 디렉토리 초기화 및 샘플별 작업 디렉토리 생성
// ─────────────────────────────────────────────────────────
process SCRATCH_SETUP {
    tag "${sample_name}"
    label 'process_low'

    input:
        val  sample_name
        val  scratch_dir

    output:
        val  "${scratch_dir}/${sample_name}", emit: sample_scratch_dir

    script:
        """
        set -euo pipefail

        SAMPLE_SCRATCH="${scratch_dir}/${sample_name}"
        mkdir -p "\${SAMPLE_SCRATCH}"

        # 이전 실행 잔여 파일 정리 (재실행 안전성)
        rm -f "\${SAMPLE_SCRATCH}"/*.bam "\${SAMPLE_SCRATCH}"/*.bai 2>/dev/null || true

        echo "[SCRATCH] Initialized: \${SAMPLE_SCRATCH}"
        echo "\${SAMPLE_SCRATCH}"
        """
}

// ─────────────────────────────────────────────────────────
// SCRATCH_MOVE_FINAL
//   SSD에서 생성된 proper_paired.bam을 HDD analysis 디렉토리로 이동
//   나머지 임시 BAM은 SSD에서 삭제
// ─────────────────────────────────────────────────────────
process SCRATCH_MOVE_FINAL {
    tag "${sample_name}"
    label 'process_low'

    input:
        val  sample_name
        path proper_paired_bam
        path proper_paired_bai
        val  analysisdir
        val  sample_scratch_dir

    output:
        path "${sample_name}.proper_paired.bam",     emit: bam
        path "${sample_name}.proper_paired.bam.bai", emit: bai

    script:
        """
        set -euo pipefail

        DEST="${analysisdir}"
        mkdir -p "\${DEST}"

        echo "[SCRATCH] Moving proper_paired.bam to HDD analysis directory..."
        echo "[SCRATCH]   From: ${proper_paired_bam}"
        echo "[SCRATCH]   To:   \${DEST}/${sample_name}.proper_paired.bam"

        # mv는 같은 파일시스템 내에서는 메타데이터 변경만 발생 (즉각적)
        # SSD → HDD 간 이동이므로 cp + rm으로 처리
        cp -p ${proper_paired_bam} "\${DEST}/${sample_name}.proper_paired.bam"
        cp -p ${proper_paired_bai} "\${DEST}/${sample_name}.proper_paired.bam.bai"

        echo "[SCRATCH] Move complete."

        # 임시 BAM 파일 정리 (SSD 공간 해제)
        echo "[SCRATCH] Cleaning up temporary BAMs from SSD..."
        rm -f "${sample_scratch_dir}"/*.bam \
              "${sample_scratch_dir}"/*.bai \
              "${sample_scratch_dir}"/*.bam.bai 2>/dev/null || true

        # 작업 디렉토리 자체는 유지 (로그 파일 보존)
        echo "[SCRATCH] Cleanup complete. SSD space released."

        # 출력 파일을 현재 작업 디렉토리에 심볼릭 링크로 노출
        ln -sf "\${DEST}/${sample_name}.proper_paired.bam" ./${sample_name}.proper_paired.bam
        ln -sf "\${DEST}/${sample_name}.proper_paired.bam.bai" ./${sample_name}.proper_paired.bam.bai
        """
}

// ─────────────────────────────────────────────────────────
// SCRATCH_CLEANUP_ON_FAILURE
//   파이프라인 실패 시 SSD 잔여 파일 정리
//   Nextflow workflow.onError 핸들러에서 호출
// ─────────────────────────────────────────────────────────
process SCRATCH_CLEANUP_ON_FAILURE {
    tag "${sample_name}"
    label 'process_low'
    errorStrategy 'ignore'   // 정리 실패가 전체 파이프라인을 막지 않도록

    input:
        val  sample_name
        val  sample_scratch_dir

    output:
        stdout

    script:
        """
        echo "[SCRATCH] Cleaning up after failure: ${sample_scratch_dir}"
        rm -rf "${sample_scratch_dir}" 2>/dev/null || true
        echo "[SCRATCH] Failure cleanup complete."
        """
}
