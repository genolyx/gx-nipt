# GX-NIPT Pipeline

**GX-NIPT** is a Nextflow-based Non-Invasive Prenatal Testing (NIPT) analysis pipeline, refactored and improved from the original `ken-nipt` pipeline. It provides a modular, reproducible, and scalable workflow for cfDNA-based chromosomal aneuploidy detection.

---

## Key Improvements over ken-nipt

| Category | Improvement |
|---|---|
| **Architecture** | Migrated from monolithic shell + Docker to Nextflow DSL2 modular pipeline |
| **Reproducibility** | Each process runs in an isolated container with pinned tool versions |
| **Parallelism** | Alignment, QC, and analysis steps run in parallel across samples |
| **Logging** | Structured per-process logs; Nextflow Tower / HTML report integration |
| **Bug Fixes** | Fixed pipeline completion marker race condition, f-string errors, exception handling gaps |
| **EZD** | Improved threshold calculation with Youden's J index; added ROC-based boundary validation |
| **PRIZM** | Added outlier-robust SD estimation (MAD-based) |
| **Fetal Fraction** | Combined Y-chromosome ratio (YFF) + fragment size-based FF with weighted ensemble |
| **Gender Detection** | Two-stage classification with confidence score; handles low-FF edge cases |
| **SCA** | Improved XXY/XYY/XO/XXX detection with FF-corrected linear model |
| **Reference** | Structured reference generation workflow with quality validation |
| **hg38** | Analysis and migration guide for hg38 transition (see `docs/hg38_migration.md`) |

---

## Pipeline Overview

```
FASTQ Input
    │
    ├─► [ALIGN] BWA-MEM2 → Samtools sort/index
    │
    ├─► [QC] FastQC + Qualimap → QC filter
    │
    ├─► [DOWNSAMPLE] Picard DownsampleSam (if needed)
    │
    ├─► [FF/GENDER] Y-chromosome FF + Fragment FF → Ensemble FF + Gender
    │
    ├─► [HMMCOPY] readCounter → HMMcopy R → Copy number segments
    │
    ├─► [EZD] UAR/Z-score calculation → Threshold decision
    │
    ├─► [PRIZM] Mean/SD-based Z-score → Aneuploidy call
    │
    ├─► [WISECONDOR] WisecondorX predict → Aneuploidy call
    │
    ├─► [MICRODELETION] Coverage-based MD detection
    │
    └─► [REPORT] JSON output → HTML review page → Portal upload
```

---

## Quick Start

### Prerequisites

- Nextflow >= 23.10
- Docker or Singularity
- Reference genome (hg19 or hg38)
- Lab-specific reference files (EZD, PRIZM, WCX npz)

### Installation

```bash
git clone https://github.com/genolyx/gx-nipt.git
cd gx-nipt
```

### Run

```bash
# Single sample
nextflow run main.nf \
  --labcode cordlife \
  --sample_id SAMPLE001 \
  --fastq_r1 /path/to/R1.fastq.gz \
  --fastq_r2 /path/to/R2.fastq.gz \
  --outdir /path/to/output \
  -profile docker

# Batch mode (CSV samplesheet)
nextflow run main.nf \
  --labcode cordlife \
  --samplesheet samples.csv \
  --outdir /path/to/output \
  -profile docker
```

### Samplesheet Format

```csv
sample_id,fastq_r1,fastq_r2
SAMPLE001,/data/SAMPLE001_R1.fastq.gz,/data/SAMPLE001_R2.fastq.gz
SAMPLE002,/data/SAMPLE002_R1.fastq.gz,/data/SAMPLE002_R2.fastq.gz
```

---

## Output Directory Structure

The output directory mirrors the structure expected by the Portal daemon:

```
{outdir}/{sample_id}/
├── bam/
│   ├── {sample_id}.sorted.bam
│   └── {sample_id}.sorted.bam.bai
├── qc/
│   ├── fastqc/
│   ├── qualimap/
│   └── qc_summary.json
├── ff_gender/
│   ├── ff_result.json
│   └── gender_result.json
├── hmmcopy/
│   ├── {sample_id}_readcounts.wig
│   └── {sample_id}_segments.csv
├── ezd/
│   └── ezd_result.json
├── prizm/
│   └── prizm_result.json
├── wisecondor/
│   └── wcx_result.json
├── microdeletion/
│   └── md_result.json
├── report/
│   ├── {sample_id}_result.json
│   └── {sample_id}_review.html
└── logs/
    ├── pipeline.log
    └── {step}_*.log
```

---

## Configuration

### Lab-specific Config

Lab configurations are stored in `conf/labs/{labcode}/pipeline_config.json`:

```json
{
  "qc_mapping_rate_min": 0.7,
  "qc_duplication_rate_max": 0.3,
  "qc_total_reads_min": 5000000,
  "ff_threshold_min": 2.0,
  "ff_threshold_low_risk": 4.0,
  "wc_binsize": 50000,
  "wcx_binsize": 200000,
  "ezd_qc_cutoff": 0.01,
  "prizm_qc_cutoff": 0.01
}
```

### Reference Files

gx-nipt keeps two kinds of reference data, separated by size:

**1. Lightweight lab config (tracked in the repo, copied per-run):**

```
conf/labs/{labcode}/
└── pipeline_config.json   # thresholds, QC cutoffs, module switches
```

`run_nipt.sh` seeds `${ROOT_DIR}/config/{labcode}/pipeline_config.json`
from this file on first run.

**2. Heavy reference data (NOT in the repo — provisioned on the host):**

`params.ref_dir` (default `/data/reference`) is bind-mounted read-only
into every container at the same path. Populate it with:

```
/data/reference/
├── genomes/hg19/hg19.fa      (+ .fai .0123 .amb .ann .bwt.2bit.64 .pac)
├── hmmcopy/hg19.{50kb,10mb}.{gc,map}.wig
├── models/seqff_model.pkl
└── labs/{labcode}/
    ├── WC/{group}/{group}_200k.npz
    ├── WCX/{group}/{M,F}_200k.npz
    ├── EZD/{group}/…
    ├── PRIZM/{group}/…
    └── bed/…
```

Override the root with `--ref_dir /path/to/reference` on the Nextflow
command line, or `--ref-dir …` / `NIPT_REF_DIR=…` when going through
`bin/run_nipt.sh` or `gx-daemon`.

See [docs/reference_generation.md](docs/reference_generation.md) for
how to build the lab-specific bundles.

---

## Reference Generation

See `bin/scripts/utils/reference/create_reference.py` and the [Reference Generation Guide](docs/reference_generation.md).

### Quick Reference Build

```bash
# Build EZD reference
python bin/scripts/utils/reference/create_reference.py \
  --labcode cordlife \
  --group orig \
  --sample_list samples_normal.txt \
  --algo EZD

# Build WCX reference
python bin/scripts/utils/reference/create_reference.py \
  --labcode cordlife \
  --group orig \
  --sample_list samples_normal.txt \
  --algo WCX \
  --gender M
```

---

## Algorithm Details

### EZD (Empirical Z-score Detection)

EZD computes per-chromosome UAR (Unique Autosomal Ratio) and Z-scores against a Normal reference distribution. Threshold boundaries (UAR_min, UAR_max, Z_min, Z_max) are derived using Youden's J index on ROC curves computed from Normal and Positive training samples.

**Key parameters:** `UAR_min`, `UAR_max`, `Z_min`, `Z_max`, `qc_cutoff`

### PRIZM

PRIZM uses a Mean/SD-based Z-score model with outlier-robust estimation (MAD). It provides a complementary aneuploidy call to EZD.

**Key parameters:** `qc_cutoff`, `zscore_cutoff`

### Fetal Fraction (FF)

FF is estimated using a weighted ensemble of:
1. **YFF** — Y-chromosome read ratio (males only)
2. **Fragment FF** — cfDNA fragment size distribution (fetal vs. maternal)

**Improvement:** The ensemble weight is dynamically adjusted based on gender confidence and read depth.

### Gender Detection

Two-stage classification:
1. **Stage 1** — Y-chromosome ratio threshold (GD1)
2. **Stage 2** — Fragment size profile validation (GD2)

Confidence score is reported alongside M/F classification.

### SCA (Sex Chromosome Aneuploidy)

FF-corrected linear model for XXY, XYY, XO, XXX detection. Each SCA type has independent slope/intercept/threshold parameters configurable via the Admin UI.

---

## Web Admin UI

The **GX-NIPT Admin** web application provides a browser-based interface for:

- Managing lab-specific Reference files (EZD, PRIZM, WCX)
- Editing Threshold parameters (UAR_min, UAR_max, Z_min, Z_max)
- Running the Threshold auto-calculation workflow (Youden's J)
- Managing Normal/Positive sample lists
- Running Reference generation jobs
- Configuring pipeline_config.json parameters per labcode
- Configuring SCA detection parameters (XXY, XYY, XO, XXX)
- Viewing audit logs of all configuration changes

See the [gx-nipt-admin repository](https://github.com/genolyx/gx-nipt-admin) for setup instructions.

---

## hg38 Migration

See [docs/hg38_migration.md](docs/hg38_migration.md) for a detailed analysis of the benefits and considerations for migrating from hg19 to hg38.

---

## Development

### Running Tests

```bash
# Python unit tests
python -m pytest tests/

# Nextflow pipeline test
nextflow run main.nf -profile test,docker
```

### Adding a New Lab

1. Create `conf/labs/{new_labcode}/pipeline_config.json`
2. Generate reference files using `create_reference.py`
3. Register the labcode in the Admin UI

---

## License

Proprietary — Genolyx Inc. All rights reserved.
