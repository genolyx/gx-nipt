# hg38 Migration Analysis

## Overview

The current GX-NIPT pipeline uses **hg19 (GRCh37)** as the reference genome. This document analyzes the benefits, risks, and migration path for transitioning to **hg38 (GRCh38)**.

---

## Benefits of hg38

| Category | hg19 | hg38 | Benefit |
|---|---|---|---|
| **Assembly completeness** | ~3.0 Gb | ~3.1 Gb | 5% more sequence, fewer gaps |
| **Centromere/telomere** | Unresolved (N-masked) | T2T-CHM13 patches available | Improved read mapping in pericentromeric regions |
| **ALT contigs** | Not included | Included (decoy sequences) | Reduced multi-mapping artifacts |
| **SNP/variant density** | dbSNP b137 | dbSNP b155+ | More accurate duplication marking |
| **Gene annotation** | GENCODE v19 | GENCODE v44 | More complete transcript models |
| **Tool compatibility** | Legacy tools | All modern tools | BWA-MEM2, GATK4, DeepVariant |

---

## Impact on NIPT Analysis

### Alignment Quality

hg38 includes **ALT contigs** and **decoy sequences** that capture reads which would otherwise multi-map in hg19. For NIPT cfDNA analysis:

- Reads from **pericentromeric regions** (chr13, chr18, chr21) map more accurately in hg38.
- The **chrY PAR1/PAR2** regions are better defined in hg38, which directly affects Y-chromosome FF calculation.
- **Segmental duplications** are better annotated, reducing false duplication marking.

### Fetal Fraction (YFF) Impact

The Y-chromosome pseudo-autosomal regions (PAR1, PAR2) overlap with chrX in hg19, causing systematic over-counting of chrY reads in females. hg38 separates PAR regions more cleanly:

```
hg19: chrY PAR1 = chrY:10,001-2,699,520 (overlaps chrX PAR1)
hg38: chrY PAR1 = chrY:10,001-2,781,479 (masked in chrY, counted only in chrX)
```

**Expected effect:** More accurate YFF in both male and female samples; reduced false-positive male calls.

### EZD / PRIZM Reference Compatibility

**Critical:** EZD and PRIZM references are **genome-version specific**. Migrating to hg38 requires:

1. Re-aligning all Normal and Positive training samples to hg38.
2. Rebuilding all EZD, PRIZM, and WCX references from scratch.
3. Recalculating all Threshold boundaries (UAR_min, UAR_max, Z_min, Z_max).

This is a **one-time migration cost** that cannot be avoided.

### Microdeletion (MD) Detection

hg38 provides more accurate coordinates for clinically relevant microdeletion regions:

| Region | hg19 coordinates | hg38 coordinates | Accuracy |
|---|---|---|---|
| 22q11.2 | chr22:18,900,000-21,500,000 | chr22:18,912,231-21,465,659 | Improved |
| 15q11-q13 | chr15:22,765,628-28,560,000 | chr15:22,835,778-28,649,634 | Improved |
| 5p15 (Cri-du-chat) | chr5:1-11,800,000 | chr5:1-11,600,000 | Similar |
| 1p36 | chr1:1-5,000,000 | chr1:1-5,000,000 | Similar |

---

## Migration Risks

| Risk | Severity | Mitigation |
|---|---|---|
| Reference rebuild cost | High | Plan 3–6 month parallel validation period |
| Tool version incompatibility | Medium | Pin tool versions in Docker image |
| BED file coordinate mismatch | High | Use liftOver or rebuild BED files from scratch |
| Existing result comparison | Medium | Run parallel hg19/hg38 on same samples for validation |
| chrY PAR masking behavior | Medium | Validate YFF on known male/female samples |

---

## Migration Recommendation

### Short-term (Current)

Continue with **hg19** for production. The existing reference data and validated thresholds represent significant investment.

### Medium-term (6–12 months)

Begin **parallel validation** with hg38:

1. Re-align 50+ Normal samples to hg38.
2. Build hg38 EZD/PRIZM/WCX references.
3. Run both pipelines on the same samples and compare results.
4. Validate FF accuracy improvement on known male samples.

### Long-term (12–24 months)

**Migrate to hg38** once:

- Parallel validation confirms equivalent or better sensitivity/specificity.
- All lab-specific reference files are rebuilt and validated.
- The Admin UI supports genome version tagging per reference set.

---

## Implementation Notes

### BWA-MEM2 Index

```bash
# Build hg38 index (requires ~60GB RAM)
bwa-mem2 index hg38.fa

# Recommended: use pre-built index from NCBI
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bwa_index.tar.gz
```

### Chromosome Naming

hg38 uses `chr1`-style naming (same as hg19 UCSC convention). No changes needed for chromosome name handling.

### BED File Liftover

```bash
# Liftover existing hg19 BED files to hg38
liftOver input_hg19.bed hg19ToHg38.over.chain output_hg38.bed unmapped.bed
```

---

## Conclusion

Migrating to hg38 is **recommended for new deployments** and provides measurable improvements in alignment accuracy, YFF estimation, and microdeletion detection. However, it requires a complete reference rebuild and parallel validation period. The GX-NIPT pipeline architecture (Nextflow + Docker) is designed to support both genome versions simultaneously via the `--genome` parameter.
