/*
 * =========================================================
 *  Module: Samtools operations
 * =========================================================
 */

process SAMTOOLS_SORT_INDEX {
    tag "${sample_name}"
    label 'process_medium'
    label 'nipt_docker'

    input:
        val  sample_name
        path bam
        val  suffix        // 'sorted', 'unique', etc.
        val  analysisdir

    output:
        path "${sample_name}.${suffix}.bam",     emit: bam
        path "${sample_name}.${suffix}.bam.bai", emit: bai

    script:
        def threads = task.cpus
        def mem     = task.memory.toGiga().intValue()
        """
        set -euo pipefail

        samtools sort \\
            -@ ${threads} \\
            -m ${mem}G \\
            -o ${sample_name}.${suffix}.bam \\
            ${bam}

        samtools index \\
            -@ ${threads} \\
            ${sample_name}.${suffix}.bam

        echo "[SAMTOOLS] Sort+Index complete: ${sample_name}.${suffix}.bam"
        """
}

process SAMTOOLS_FILTER_UNIQUE {
    tag "${sample_name}"
    label 'process_medium'
    label 'nipt_docker'

    input:
        val  sample_name
        path bam
        val  analysisdir

    output:
        path "${sample_name}.unique.bam",     emit: bam
        path "${sample_name}.unique.bam.bai", emit: bai

    script:
        def threads = task.cpus
        """
        set -euo pipefail

        # Filter: unique reads only (MQ >= 1, no secondary/supplementary)
        samtools view \\
            -@ ${threads} \\
            -b -q 1 \\
            -F 0x100 -F 0x800 \\
            ${bam} \\
        | samtools sort -@ ${threads} -o ${sample_name}.unique.bam -

        samtools index -@ ${threads} ${sample_name}.unique.bam

        echo "[SAMTOOLS] Unique filter complete: ${sample_name}.unique.bam"
        """
}

process SAMTOOLS_PROPER_PAIRED {
    tag "${sample_name}"
    label 'process_medium'
    label 'nipt_docker'

    input:
        val  sample_name
        path bam
        val  analysisdir

    output:
        path "${sample_name}.proper_paired.bam",     emit: bam
        path "${sample_name}.proper_paired.bam.bai", emit: bai

    script:
        def threads = task.cpus
        """
        set -euo pipefail

        # Filter: properly paired reads only
        samtools view \\
            -@ ${threads} \\
            -b -f 0x2 \\
            ${bam} \\
        | samtools sort -@ ${threads} -o ${sample_name}.proper_paired.bam -

        samtools index -@ ${threads} ${sample_name}.proper_paired.bam

        echo "[SAMTOOLS] Proper-paired filter complete: ${sample_name}.proper_paired.bam"
        """
}

/*
 * TLEN-based split: proper_paired.bam → of_orig / of_fetus / of_mom
 *
 *   of_orig    : proper_paired.bam을 그대로 재사용 (alias)
 *   of_fetus   : abs(tlen) <  160   (short cfDNA – 태아 유래 우세)
 *   of_mom     : abs(tlen) >  184   (long cfDNA  – 모체 유래 우세)
 *
 * ken-nipt의 awk -f {fetus,mom}.awk 로직과 동등. samtools 1.13+의
 * expression filter (-e)로 구현하여 awk 의존을 제거.
 */
process SAMTOOLS_SPLIT_FETUS_MOM {
    tag "${sample_name}"
    label 'process_medium'
    label 'nipt_docker'

    publishDir "${analysisdir}/${sample_name}", mode: 'copy', overwrite: true,
               pattern: "${sample_name}.of_*.bam*"

    input:
        val  sample_name
        path bam        // proper_paired.bam
        path bai        // proper_paired.bam.bai
        val  analysisdir

    output:
        tuple val(sample_name),
              path("${sample_name}.of_orig.bam"),   path("${sample_name}.of_orig.bam.bai"),
              path("${sample_name}.of_fetus.bam"),  path("${sample_name}.of_fetus.bam.bai"),
              path("${sample_name}.of_mom.bam"),    path("${sample_name}.of_mom.bam.bai"),
              emit: trio

    script:
        def threads   = task.cpus
        def fetus_max = 160       // same threshold as ken-nipt/bin/scripts/fetus.awk
        def mom_min   = 184       // same threshold as ken-nipt/bin/scripts/mom.awk
        """
        set -euo pipefail

        # 0) of_orig: proper_paired를 그대로 rename
        cp ${bam}                       ${sample_name}.of_orig.bam
        cp ${bai}                       ${sample_name}.of_orig.bam.bai

        # 1) of_fetus: |TLEN| < ${fetus_max}
        # awk col-9 (TLEN) squared < SIZE^2 — same as ken-nipt/bin/scripts/fetus.awk
        # samtools 1.13 does not support abs() in -e expressions; use awk instead.
        samtools view -h -@ ${threads} ${bam} | \\
            awk 'BEGIN{S2=${fetus_max}*${fetus_max}} /^@/{print;next} \$9*\$9<S2{print}' | \\
            samtools view -b -@ ${threads} -o ${sample_name}.of_fetus.bam
        samtools index -@ ${threads} ${sample_name}.of_fetus.bam

        # 2) of_mom: |TLEN| > ${mom_min}
        # awk col-9 (TLEN) squared > SIZE^2 — same as ken-nipt/bin/scripts/mom.awk
        samtools view -h -@ ${threads} ${bam} | \\
            awk 'BEGIN{S2=${mom_min}*${mom_min}} /^@/{print;next} \$9*\$9>S2{print}' | \\
            samtools view -b -@ ${threads} -o ${sample_name}.of_mom.bam
        samtools index -@ ${threads} ${sample_name}.of_mom.bam

        # Sanity counts (non-fatal; useful for trace log)
        for g in orig fetus mom; do
            n=\$(samtools view -c -@ ${threads} ${sample_name}.of_\${g}.bam)
            echo "[SPLIT] ${sample_name}.of_\${g}.bam reads=\${n}"
        done
        """
}
