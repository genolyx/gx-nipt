/*
 * =========================================================
 *  Module: HMMcopy
 *  readCounter → HMMcopy.R normalization
 * =========================================================
 */

process READCOUNTER_50KB {
    tag "${sample_name}"
    label 'process_low'
    label 'nipt_docker'

    publishDir "${analysisdir}/${sample_name}/Output_hmmcopy", mode: 'copy', overwrite: true

    input:
        val  sample_name
        path bam
        val  analysisdir

    output:
        path "${sample_name}.proper_paired.50kb.wig", emit: wig

    script:
        def gc_wig   = "/data/refs/hmmcopy/hg19.50kb.gc.wig"
        def map_wig  = "/data/refs/hmmcopy/hg19.50kb.map.wig"
        """
        set -euo pipefail

        readCounter \\
            --window 50000 \\
            --quality 20 \\
            --chromosome "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY" \\
            ${bam} \\
            > ${sample_name}.proper_paired.50kb.wig

        echo "[READCOUNTER 50kb] Complete for ${sample_name}"
        """
}

process READCOUNTER_10MB {
    tag "${sample_name}"
    label 'process_low'
    label 'nipt_docker'

    publishDir "${analysisdir}/${sample_name}/Output_hmmcopy", mode: 'copy', overwrite: true

    input:
        val  sample_name
        path bam
        val  analysisdir

    output:
        path "${sample_name}.proper_paired.10mb.wig", emit: wig

    script:
        """
        set -euo pipefail

        readCounter \\
            --window 10000000 \\
            --quality 20 \\
            --chromosome "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY" \\
            ${bam} \\
            > ${sample_name}.proper_paired.10mb.wig

        echo "[READCOUNTER 10mb] Complete for ${sample_name}"
        """
}

process HMMCOPY_R_50KB {
    tag "${sample_name}"
    label 'process_low'
    label 'nipt_docker'

    publishDir "${analysisdir}/${sample_name}/Output_hmmcopy", mode: 'copy', overwrite: true

    input:
        val  sample_name
        path wig
        val  analysisdir

    output:
        path "${sample_name}.proper_paired.50kb.wig.Normalization.txt", emit: norm_txt

    script:
        def gc_wig  = "/data/refs/hmmcopy/hg19.50kb.gc.wig"
        def map_wig = "/data/refs/hmmcopy/hg19.50kb.map.wig"
        def r_script = "/opt/gx-nipt/bin/scripts/HMMcopy.R"
        """
        set -euo pipefail

        Rscript ${r_script} \\
            ${wig} \\
            ${gc_wig} \\
            ${map_wig} \\
            ${sample_name}.proper_paired.50kb.wig.Normalization.txt \\
            50000

        echo "[HMMCOPY 50kb] Complete for ${sample_name}"
        """
}

process HMMCOPY_R_10MB {
    tag "${sample_name}"
    label 'process_low'
    label 'nipt_docker'

    publishDir "${analysisdir}/${sample_name}/Output_hmmcopy", mode: 'copy', overwrite: true

    input:
        val  sample_name
        path wig
        val  analysisdir

    output:
        path "${sample_name}.proper_paired.10mb.wig.Normalization.txt", emit: norm_txt

    script:
        def gc_wig  = "/data/refs/hmmcopy/hg19.10mb.gc.wig"
        def map_wig = "/data/refs/hmmcopy/hg19.10mb.map.wig"
        def r_script = "/opt/gx-nipt/bin/scripts/HMMcopy.R"
        """
        set -euo pipefail

        Rscript ${r_script} \\
            ${wig} \\
            ${gc_wig} \\
            ${map_wig} \\
            ${sample_name}.proper_paired.10mb.wig.Normalization.txt \\
            10000000

        echo "[HMMCOPY 10mb] Complete for ${sample_name}"
        """
}

process COUNT_10MB {
    tag "${sample_name}"
    label 'process_low'
    label 'nipt_docker'

    publishDir "${analysisdir}/${sample_name}/Output_hmmcopy", mode: 'copy', overwrite: true

    input:
        val  sample_name
        path wig
        val  analysisdir

    output:
        path "${sample_name}.proper_paired.10mb.count.txt", emit: count_txt

    script:
        """
        set -euo pipefail

        python3 /opt/gx-nipt/bin/scripts/utils/wig_to_count.py \\
            --wig ${wig} \\
            --output ${sample_name}.proper_paired.10mb.count.txt

        echo "[COUNT 10mb] Complete for ${sample_name}"
        """
}
