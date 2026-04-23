# Reference Generation Guide

## Overview

GX-NIPT requires lab-specific reference files for each algorithm (EZD, PRIZM, WCX, WC). References must be regenerated periodically as the sample pool grows to maintain accurate threshold boundaries.

---

## Host Layout and Docker Bind-Mount

Reference data is **not** baked into the `gx-nipt` image. The pipeline
reads everything from a single host directory controlled by
`params.ref_dir` (default `/data/reference`, overridable via
`--ref-dir` on `run_nipt.sh` or `--ref_dir` to Nextflow). That directory
is bind-mounted read-only into every container at the same path by the
Docker profile in `nextflow.config`.

Populate the host with this tree before running the pipeline:

```
/data/reference/
├── genomes/hg19/hg19.fa          (+ .fai .0123 .amb .ann .bwt.2bit.64 .pac)
├── hmmcopy/hg19.50kb.gc.wig      hg19.50kb.map.wig
│          hg19.10mb.gc.wig       hg19.10mb.map.wig
├── models/seqff_model.pkl
└── labs/<labcode>/
    ├── WC/<group>/<group>_200k.npz
    ├── WCX/<group>/{M,F}_200k.npz
    ├── EZD/<group>/…
    ├── PRIZM/<group>/…
    └── bed/…
```

`<labcode>` maps to the `--labcode` argument (`CORDLIFE`, `UCL`,
`VN`, …). `<group>` is the reference grouping used by the lab
(typically `orig`, `fetus`, `mom`, etc.; see each lab's
`conf/labs/<lab>/pipeline_config.json`).

> When invoked from `gx-daemon`, the daemon forwards `NIPT_REF_DIR` via
> `--ref-dir`, and the wrapper fails fast if the required genome/hmmcopy
> files or the `labs/<labcode>/` subtree are missing.

### Deferred: gx-FF model and gx-cnv reference

The new algorithms introduced in v1.1.0 ship **without** their own
reference artefacts:

| Algorithm | File type | Runtime flag | Status |
|-----------|-----------|--------------|--------|
| gx-FF     | `.pkl` LightGBM + DNN ensemble | `--gxff_model` / `--gxff-model` | deferred — trained in `genolyx/gx-FF` repo |
| gx-cnv    | `.npz` hybrid reference panel  | `--gxcnv_reference` / `--gxcnv-reference` (alias: `--gxcnv-model`, `--gxcnv_model`) | deferred — built in `genolyx/gx-cnv` repo |

These are **never required** to run the pipeline: `main.nf` substitutes
a `NO_FILE` sentinel when the flag is absent, `GXFF_ENSEMBLE` falls
back to seqFF-only FF estimation, and the gx-cnv pathway short-circuits
so WisecondorX stays on the critical path. The wrapper's reference
preflight does not include either file either.

Once the companion repos publish artefacts, drop them under
`${ref_dir}/models/gxff_model.pkl` and
`${ref_dir}/gxcnv/<panel>.npz` (or wherever you prefer) and pass the
path through `--gxff_model` / `--gxcnv_reference`, or set
`NIPT_GXFF_MODEL` / `NIPT_GXCNV_REFERENCE` (or the `..._MODEL` alias)
in `gx-daemon/.env` so every subsequent order picks them up
automatically.

Ensemble weighting (when gx-FF is enabled), as implemented in
`workflows/ff_gender.nf`:

- FF < 5 % — use gx-FF alone.
- FF ≥ 5 % — blend `0.6 * gx-FF + 0.4 * seqFF`.

---

## When to Regenerate References

| Trigger | Action |
|---|---|
| New lab onboarding | Generate initial references with minimum 30 Normal samples |
| Quarterly accumulation | Regenerate when Normal sample count increases by >20% |
| Threshold drift detected | Regenerate if false positive rate increases |
| Algorithm parameter change | Always regenerate after changing qc_cutoff or binsize |
| Genome version change | Full rebuild required (hg19 → hg38) |

---

## Minimum Sample Requirements

| Algorithm | Minimum Normal | Recommended | Positive (for EZD threshold) |
|---|---|---|---|
| EZD | 30 | 100+ | 10+ per chromosome |
| PRIZM | 30 | 100+ | N/A |
| WCX | 30 | 100+ | N/A |
| WC | 30 | 100+ | N/A |

---

## EZD Reference Generation

### Step 1: Prepare Sample List

Create a text file with one BAM path per line (Normal samples only):

```
/data/samples/SAMPLE001/bam/SAMPLE001.sorted.bam
/data/samples/SAMPLE002/bam/SAMPLE002.sorted.bam
...
```

### Step 2: Run Reference Builder

```bash
python bin/scripts/utils/reference/create_reference.py \
  --labcode cordlife \
  --group orig \
  --sample_list normal_samples.txt \
  --algo EZD \
  --output_dir conf/labs/cordlife/ezd/orig/
```

### Step 3: Calculate Thresholds (Youden's J)

After building the reference, calculate optimal thresholds using Normal + Positive samples:

```bash
python bin/scripts/utils/reference/create_reference.py \
  --labcode cordlife \
  --group orig \
  --normal_list normal_samples.txt \
  --positive_list positive_samples.txt \
  --algo EZD_THRESHOLD \
  --output_dir conf/labs/cordlife/ezd/orig/
```

This computes ROC curves for each chromosome and applies Youden's J index to determine:
- `UAR_min`, `UAR_max` — UAR-based boundaries
- `Z_min`, `Z_max` — Z-score-based boundaries

### Step 4: Validate

Review the generated thresholds in the Admin UI:
1. Navigate to **EZD Reference** → select labcode and group
2. Check that threshold values are reasonable (Z_max typically 2.5–4.0)
3. Review the chromosome-level distribution charts

---

## PRIZM Reference Generation

```bash
python bin/scripts/utils/reference/create_reference.py \
  --labcode cordlife \
  --group orig \
  --sample_list normal_samples.txt \
  --algo PRIZM \
  --output_dir conf/labs/cordlife/prizm/
```

PRIZM reference contains per-bin Mean and SD values. Outlier samples are automatically excluded using MAD-based filtering.

---

## WCX / WC Reference Generation

WCX and WC references must be built **separately for M and F** samples:

```bash
# Male reference
python bin/scripts/utils/reference/create_reference.py \
  --labcode cordlife \
  --group orig \
  --sample_list normal_male_samples.txt \
  --algo WCX \
  --gender M \
  --binsize 200000 \
  --output_dir conf/labs/cordlife/wcx/M/

# Female reference
python bin/scripts/utils/reference/create_reference.py \
  --labcode cordlife \
  --group orig \
  --sample_list normal_female_samples.txt \
  --algo WCX \
  --gender F \
  --binsize 200000 \
  --output_dir conf/labs/cordlife/wcx/F/
```

The output is a `.npz` file compatible with WisecondorX.

---

## Using the Admin UI for Reference Management

The **GX-NIPT Admin** web application provides a guided wizard for reference generation:

1. Navigate to **Reference Wizard**
2. Select Normal and/or Positive samples from the sample list
3. Choose algorithms (EZD, PRIZM, WCX, WC) and parameters
4. Review the job summary and submit
5. Monitor progress in **Audit Logs**

The wizard submits a job record to the database. The actual computation must be triggered on the analysis server using the job parameters.

---

## Quality Checks After Reference Generation

After generating new references, run the following validations:

### EZD Threshold Sanity Check

```python
import json

with open('conf/labs/cordlife/ezd/orig/ezd_reference.json') as f:
    ref = json.load(f)

for chr_name, data in ref['chromosomes'].items():
    z_range = data['Z_max'] - data['Z_min']
    assert 1.0 < z_range < 8.0, f"{chr_name}: Z range {z_range} is suspicious"
    print(f"{chr_name}: Z_min={data['Z_min']:.3f}, Z_max={data['Z_max']:.3f}, n={data['n_samples']}")
```

### WCX Reference Quality

```bash
python bin/scripts/utils/wisecondor_quality_check.py \
  --reference conf/labs/cordlife/wcx/M/reference.npz \
  --output_dir /tmp/wcx_qc/
```

---

## Troubleshooting

| Problem | Cause | Solution |
|---|---|---|
| Z_max < 2.0 for chr21 | Too few Positive samples | Add more T21 Positive samples |
| High false positive rate | Reference built from mixed samples | Remove Positive samples from Normal list |
| WCX reference fails | Insufficient samples | Minimum 30 samples per gender |
| Threshold drift over time | Sample population shift | Regenerate reference with recent samples |
