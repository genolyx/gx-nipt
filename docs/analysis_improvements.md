# Analysis Improvements: ken-nipt → gx-nipt

## 1. Bug Fixes

### 1.1 Pipeline Completion Marker Race Condition

**File:** `src/run_nipt.sh` (ken-nipt)

**Bug:** The `PIPELINE_MARKER` file check was performed immediately after `docker run`, but the Docker container writes the marker file asynchronously. In rare cases, the pipeline was reported as complete before the marker was written.

**Fix (gx-nipt):** Nextflow's process completion tracking replaces the marker file mechanism entirely. Each process reports its exit status directly to the Nextflow engine.

---

### 1.2 EZD Threshold Edge Case — Empty Positive Sample Set

**File:** `bin/scripts/utils/ezd_get_min_threshold.py` (ken-nipt)

**Bug:** When `positive_samples` was an empty list, the ROC calculation attempted `np.array([]) / np.array([])`, producing `NaN` thresholds that were silently written to the reference file.

**Fix (gx-nipt):**

```python
# Before (ken-nipt)
fpr, tpr, thresholds = roc_curve(y_true, y_score)

# After (gx-nipt)
if len(positive_samples) == 0:
    logger.warning(f"No positive samples for {chr_name}. Using default threshold.")
    return DEFAULT_THRESHOLD
fpr, tpr, thresholds = roc_curve(y_true, y_score)
if len(thresholds) < 2:
    logger.warning(f"Insufficient ROC points for {chr_name}. Using default threshold.")
    return DEFAULT_THRESHOLD
```

---

### 1.3 Gender Detector — Low-FF False Female Classification

**File:** `bin/scripts/modules/gender_detector.py` (ken-nipt)

**Bug:** Samples with FF < 2% had very low Y-chromosome read counts, causing male samples to be misclassified as female. The single-threshold approach (`y_ratio > GD1_THRESHOLD`) was insufficient for low-FF samples.

**Fix (gx-nipt — `ff_gender_improved.py`):**

```python
# Two-stage classification with FF-adjusted threshold
def classify_gender(y_ratio: float, ff: float, gd1: float, gd2: float) -> dict:
    # Stage 1: Standard threshold
    if y_ratio >= gd2:
        return {"gender": "M", "confidence": "high", "method": "stage2"}
    
    # Stage 2: FF-adjusted boundary for low-FF samples
    ff_adjusted_threshold = gd1 * max(0.5, min(1.0, ff / 4.0))
    if y_ratio >= ff_adjusted_threshold:
        return {"gender": "M", "confidence": "low", "method": "ff_adjusted"}
    
    return {"gender": "F", "confidence": "high" if y_ratio < gd1 * 0.3 else "medium", "method": "stage1"}
```

---

### 1.4 PRIZM Runner — SD Underestimation

**File:** `bin/scripts/modules/prizm_runner.py` (ken-nipt)

**Bug:** PRIZM used `np.std()` (population SD) instead of `np.std(ddof=1)` (sample SD) for reference building. With small sample sizes (N < 50), this underestimated SD by ~1%, leading to inflated Z-scores.

**Fix (gx-nipt):**

```python
# Before
sd = np.std(values)

# After — use sample SD and MAD-based outlier removal
from scipy.stats import median_abs_deviation
mad = median_abs_deviation(values, scale='normal')
filtered = values[np.abs(values - np.median(values)) < 3 * mad]
sd = np.std(filtered, ddof=1)
```

---

### 1.5 HMMcopy readCounter — Unmapped Read Inclusion

**File:** `src/run_nipt.sh` (ken-nipt)

**Bug:** The `readCounter` command was run without the `-q 20` quality filter, including low-quality and unmapped reads in the bin counts. This inflated read counts in low-complexity regions.

**Fix (gx-nipt — `modules/hmmcopy.nf`):**

```nextflow
script:
"""
readCounter \
    --window ${params.hmmcopy_binsize} \
    --quality 20 \
    --chromosome "${chroms}" \
    ${bam} > ${prefix}.wig
"""
```

---

## 2. Algorithm Improvements

### 2.1 EZD Threshold Calculation — Youden's J Index

**Previous approach:** Fixed percentile-based thresholds (e.g., mean ± 3SD of Normal distribution).

**Improved approach:** ROC-based threshold optimization using Youden's J index:

```
J = Sensitivity + Specificity - 1
Optimal threshold = argmax(J)
```

This maximizes the sum of sensitivity and specificity simultaneously, producing more balanced boundaries that adapt to the actual Normal/Positive sample distributions.

**Implementation:** `bin/scripts/utils/reference/create_reference.py` — `calculate_youden_threshold()`

---

### 2.2 Fetal Fraction — Weighted Ensemble

**Previous approach:** YFF (Y-chromosome ratio) only for males; no FF for females.

**Improved approach:** Weighted ensemble of YFF + Fragment-size FF:

```python
def ensemble_ff(yff: float, fragment_ff: float, gender: str, y_ratio: float) -> float:
    if gender == "M" and y_ratio > GD2_THRESHOLD:
        # High-confidence male: weight YFF more heavily
        w_yff = 0.7
        w_frag = 0.3
    elif gender == "M":
        # Low-confidence male: equal weight
        w_yff = 0.5
        w_frag = 0.5
    else:
        # Female: fragment FF only
        return fragment_ff
    
    return w_yff * yff + w_frag * fragment_ff
```

**Benefit:** More stable FF estimates, especially for samples near the gender boundary.

---

### 2.3 SCA Detection — FF Correction

**Previous approach:** Fixed slope/intercept per SCA type without FF correction.

**Improved approach:** FF-corrected linear model:

```
SCA_score = (observed_ratio - (slope * FF + intercept)) / SD_reference
```

This accounts for the fact that SCA signal strength scales with fetal fraction, reducing false positives in low-FF samples.

---

## 3. Workflow Improvements

### 3.1 Parallelism

| Step | ken-nipt | gx-nipt |
|---|---|---|
| Multi-sample | Sequential (one at a time) | Parallel (Nextflow executor) |
| Alignment | Single thread | Multi-thread (bwa-mem2 -t) |
| QC | After alignment only | Parallel with alignment |
| EZD/PRIZM/WCX | Sequential | Parallel |

### 3.2 Logging

| Feature | ken-nipt | gx-nipt |
|---|---|---|
| Log format | Unstructured text | Structured JSON per step |
| Log location | Single `pipeline.log` | Per-process logs in `logs/` |
| Error tracing | Manual grep | Nextflow HTML report |
| Progress tracking | `PIPELINE_MARKER` file | Nextflow process status |

### 3.3 Resource Management

Nextflow profiles define resource limits per process:

```groovy
process {
    withName: 'BWA_MEM' {
        cpus   = 8
        memory = '16 GB'
        time   = '4h'
    }
    withName: 'HMMCOPY_READCOUNTER' {
        cpus   = 2
        memory = '4 GB'
        time   = '1h'
    }
}
```

---

## 4. Output Compatibility

The output directory structure is **fully compatible** with the existing Portal daemon. The `report/` subdirectory contains the same `{sample_id}_result.json` format expected by the daemon's file watcher.

Key fields preserved:
- `result.json` — top-level aneuploidy calls (T21, T18, T13, SCA)
- `ezd_result.json` — per-chromosome UAR/Z-score values
- `prizm_result.json` — per-chromosome PRIZM Z-scores
- `wcx_result.json` — WisecondorX segment calls
- `ff_result.json` — FF and gender classification
- `qc_summary.json` — QC metrics (mapping rate, duplication rate, total reads)
