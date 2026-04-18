# Reference Generation Guide

## Overview

GX-NIPT requires lab-specific reference files for each algorithm (EZD, PRIZM, WCX, WC). References must be regenerated periodically as the sample pool grows to maintain accurate threshold boundaries.

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
