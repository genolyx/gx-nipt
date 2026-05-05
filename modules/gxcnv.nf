/**
 * gx-cnv Module
 *
 * Hybrid dual-track CNV detection engine for sWGS NIPT data.
 * Runs in parallel with WisecondorX. Results are compared and
 * both are included in the final report until gx-cnv is validated.
 *
 * Algorithm:
 *   Track A: Whole-genome Z-score normalisation + CBS segmentation (WisecondorX-inspired)
 *   Track B: Regional PCA denoising + Laplace-smoothed Mahalanobis distance (BinDel-inspired)
 *   Decision: AND-gate — HIGH_RISK only when both tracks independently exceed thresholds
 *
 * Commands:
 *   gxcnv convert  BAM → NPZ (GC-corrected bin counts)
 *   gxcnv newref   NPZ files → reference panel
 *   gxcnv predict  sample NPZ + reference → CNV calls
 */

// ─────────────────────────────────────────────────────────────────────────────
// Step 1: Convert BAM → NPZ
// ─────────────────────────────────────────────────────────────────────────────
process GXCNV_CONVERT {
    tag "${sample_id}"
    label 'process_medium'
    label 'nipt_docker'

    input:
    tuple val(sample_id), path(bam), path(bai)
    val   bin_size
    path  blacklist_bed

    output:
    tuple val(sample_id), path("${sample_id}.gxcnv.npz"), emit: npz
    path  "${sample_id}.gxcnv_convert.log",               emit: log

    script:
    def blacklist_arg = blacklist_bed.name != 'NO_FILE' ? "--blacklist ${blacklist_bed}" : ""
    """
    set -euo pipefail

    gxcnv convert \\
        ${bam} \\
        ${sample_id}.gxcnv.npz \\
        --bin-size ${bin_size} \\
        --min-mapq 1 \\
        ${blacklist_arg} \\
        2>&1 | tee ${sample_id}.gxcnv_convert.log

    if [ ! -s "${sample_id}.gxcnv.npz" ]; then
        echo "ERROR: gxcnv convert produced empty NPZ for ${sample_id}" >&2
        exit 1
    fi
    """

    stub:
    """
    python3 -c "import numpy as np; np.savez_compressed('${sample_id}.gxcnv.npz', bins=np.zeros(30000))"
    touch ${sample_id}.gxcnv_convert.log
    """
}


// ─────────────────────────────────────────────────────────────────────────────
// Step 1b: Inject GC fractions from HMMcopy 50kb WIG → NPZ
//
// gxcnv convert does not compute GC fractions from BAM input (bins[:,3] = NaN),
// so GC correction is skipped.  The reference panel was built with GC-corrected
// NPZ files; using uncorrected sample NPZ causes systematic Z-score artefacts.
//
// This process patches the NPZ in-place:
//   1. Load HMMcopy 50kb wig normalization file.
//   2. Aggregate adjacent 50 kb bins → 100 kb GC fractions (mean of valid bins).
//   3. Write GC fractions into bins[:,3].
//   4. Re-run polynomial GC correction (degree 3).
//   5. Re-normalise corrected counts to reads-per-valid-bin equivalent.
// ─────────────────────────────────────────────────────────────────────────────
process GXCNV_GC_INJECT {
    tag "${sample_id}"
    label 'process_low'
    label 'nipt_docker'

    input:
    tuple val(sample_id), path(sample_npz)
    path  wig_norm   // HMMcopy 50 kb normalization file (orig group)

    output:
    tuple val(sample_id), path("${sample_id}.gxcnv.gc.npz"), emit: npz
    path  "${sample_id}.gxcnv_gc_inject.log",                emit: log

    script:
    """
    set -euo pipefail

    python3 - <<'PYEOF'
import sys, logging
import numpy as np
import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [gc_inject] %(levelname)s %(message)s",
    handlers=[
        logging.FileHandler("${sample_id}.gxcnv_gc_inject.log"),
        logging.StreamHandler(sys.stdout),
    ],
)
logger = logging.getLogger(__name__)

# ── Load GC map from HMMcopy 50kb wig ────────────────────────────────────────
wig_path = "${wig_norm}"
try:
    df = pd.read_csv(
        wig_path, sep="\\t",
        names=["chr","start","end","reads","gc","map","valid",
               "ideal","cor.gc","cor.map","copy"],
        header=0, dtype={"chr": str},
    )
except Exception as e:
    logger.error("Cannot read wig %s: %s", wig_path, e)
    sys.exit(1)

df["gc"] = pd.to_numeric(df["gc"], errors="coerce")
df = df[df["gc"] > 0].copy()
df["start_100k"] = (df["start"].astype(int) // 100_000) * 100_000
gc_map = df.groupby(["chr", "start_100k"])["gc"].mean().reset_index()
gc_lookup = {(row.chr, int(row.start_100k)): float(row.gc)
             for row in gc_map.itertuples(index=False)}
logger.info("GC map loaded: %d 100kb entries", len(gc_lookup))

# ── Load NPZ ──────────────────────────────────────────────────────────────────
d      = dict(np.load("${sample_npz}", allow_pickle=True))
bins   = d["bins"]        # (N, 4): chrom_idx, start, end, gc
chroms = list(d["chroms"])
counts = d["counts"].astype(float)
mask   = d["mask"].astype(bool)

# ── Inject GC fractions ───────────────────────────────────────────────────────
gc_fractions = np.full(len(bins), np.nan, dtype=np.float32)
hit = 0
for i, (chrom_idx, start, end, _) in enumerate(bins):
    key = (chroms[int(chrom_idx)], int(start))
    if key in gc_lookup:
        gc_fractions[i] = gc_lookup[key]
        hit += 1

frac_hit = hit / max(len(bins), 1)
logger.info("GC hit rate: %d/%d (%.1f%%)", hit, len(bins), 100 * frac_hit)

if frac_hit < 0.3:
    logger.warning("Low GC hit rate — skipping GC correction; using raw counts")
    corrected = np.full(len(counts), np.nan)
    valid = mask & (counts > 0)
    corrected[valid] = counts[valid].astype(float)
    s = np.nansum(corrected[mask])
    if s > 0:
        corrected /= s / mask.sum()
    d["corrected"] = corrected
    bins[:, 3] = gc_fractions
    d["bins"] = bins
else:
    bins[:, 3] = gc_fractions
    d["bins"] = bins

    # ── Polynomial GC correction (degree 3) ──────────────────────────────────
    poly_degree = 3
    corrected   = np.full(len(counts), np.nan)
    valid = mask & (counts > 0) & np.isfinite(gc_fractions)

    if valid.sum() >= poly_degree + 1:
        gc_bins    = np.linspace(0, 1, 101)
        gc_idx     = np.clip(np.digitize(gc_fractions[valid], gc_bins) - 1, 0, 99)
        gc_centers = (gc_bins[:-1] + gc_bins[1:]) / 2
        gc_medians = np.full(100, np.nan)
        for i in range(100):
            in_bin = valid.nonzero()[0][gc_idx == i]
            if len(in_bin) >= 3:
                gc_medians[i] = np.median(counts[in_bin])
        finite = np.isfinite(gc_medians)
        if finite.sum() >= poly_degree + 1:
            coeffs    = np.polyfit(gc_centers[finite], gc_medians[finite], poly_degree)
            predicted = np.polyval(coeffs, gc_fractions[valid])
            predicted = np.where(predicted <= 0, np.nan, predicted)
            global_med = np.nanmedian(gc_medians[finite])
            adj = counts[valid].astype(float) / predicted * global_med
            corrected[valid] = np.where(np.isfinite(adj), adj, np.nan)
            logger.info("GC correction applied (poly degree %d)", poly_degree)
        else:
            corrected[valid] = counts[valid].astype(float)
            logger.warning("Not enough GC bins for poly fit — using raw counts")
    else:
        corrected[valid] = counts[valid].astype(float)
        logger.warning("Not enough valid bins — using raw counts")

    # Re-normalise
    s = np.nansum(corrected[mask])
    if s > 0:
        corrected = corrected / s * mask.sum()
    d["corrected"] = corrected

np.savez_compressed("${sample_id}.gxcnv.gc", **d)
logger.info("GC-injected NPZ saved: ${sample_id}.gxcnv.gc.npz")
PYEOF
    """

    stub:
    """
    python3 -c "import numpy as np; np.savez_compressed('${sample_id}.gxcnv.gc', bins=np.zeros((30000,4)), counts=np.zeros(30000), mask=np.ones(30000,dtype=bool), chroms=['chr1'], corrected=np.zeros(30000))"
    touch ${sample_id}.gxcnv_gc_inject.log
    """
}


// ─────────────────────────────────────────────────────────────────────────────
// Step 2: Build reference panel (run once, not per-sample)
// ─────────────────────────────────────────────────────────────────────────────
process GXCNV_NEWREF {
    tag "gxcnv_newref"
    label 'process_high'
    label 'nipt_docker'

    input:
    path npz_files   // collection of normal sample NPZ files
    val  output_name
    val  pca_variance

    output:
    path "${output_name}.npz", emit: reference
    path "${output_name}.newref.log", emit: log

    script:
    """
    set -euo pipefail

    gxcnv newref \\
        ${npz_files.join(' ')} \\
        -o ${output_name}.npz \\
        --pca-variance ${pca_variance} \\
        2>&1 | tee ${output_name}.newref.log

    if [ ! -s "${output_name}.npz" ]; then
        echo "ERROR: gxcnv newref produced empty reference" >&2
        exit 1
    fi
    """

    stub:
    """
    python3 -c "import numpy as np; np.savez_compressed('${output_name}.npz', ref=np.zeros((50,30000)))"
    touch ${output_name}.newref.log
    """
}


// ─────────────────────────────────────────────────────────────────────────────
// Step 3: Predict CNVs (per-sample)
// ─────────────────────────────────────────────────────────────────────────────
process GXCNV_PREDICT {
    tag "${sample_id}"
    label 'process_medium'
    label 'nipt_docker'

    publishDir "${analysisdir}/${sample_id}/gxcnv", mode: 'copy', pattern: "*.tsv", overwrite: true
    publishDir "${analysisdir}/${sample_id}/gxcnv", mode: 'copy', pattern: "*.txt", overwrite: true

    input:
    tuple val(sample_id), path(sample_npz)
    path  gender_txt      // <sample>.gender.txt — parsed to select female/male reference
    val   labcode
    val   thresh_z
    val   thresh_p
    val   fetal_fraction   // from FF ensemble step; "NA" if unavailable
    val   analysisdir

    output:
    tuple val(sample_id), path("${sample_id}_bins.tsv"),     emit: bins_tsv
    tuple val(sample_id), path("${sample_id}_segments.tsv"), emit: segments_tsv
    tuple val(sample_id), path("${sample_id}_calls.tsv"),    emit: calls_tsv
    tuple val(sample_id), path("${sample_id}_regions.tsv"),  emit: regions_tsv
    tuple val(sample_id), path("${sample_id}_qcmetrics.tsv"),emit: qcmetrics_tsv
    tuple val(sample_id), path("${sample_id}_sex.txt"),      emit: sex_txt
    path  "${sample_id}.gxcnv_predict.log",                  emit: log

    script:
    def ff_arg = fetal_fraction != 'NA' ? "--fetal-fraction ${fetal_fraction}" : ""
    """
    set -euo pipefail

    # ── Resolve gender-specific reference (mirrors WisecondorX logic) ──────────
    GENDER=\$(awk -F'[\\t ]+' 'tolower(\$1) == "final_gender" {
        g = toupper(\$2)
        if (g == "XY" || g == "MALE"   || g == "M") { print "male";   exit }
        if (g == "XX" || g == "FEMALE" || g == "F") { print "female"; exit }
    }' ${gender_txt})

    if [ -z "\${GENDER}" ]; then
        echo "[GXCNV] Could not parse final_gender from ${gender_txt}; defaulting to female." >&2
        GENDER="female"
    fi
    echo "[GXCNV] Detected gender: \${GENDER}"

    REF_NPZ="${params.ref_dir}/labs/${labcode}/GXCNV/\${GENDER}/reference.npz"

    if [ ! -f "\${REF_NPZ}" ]; then
        echo "ERROR: gx-cnv reference not found: \${REF_NPZ}" >&2
        exit 1
    fi

    # ── Run prediction ──────────────────────────────────────────────────────────
    gxcnv predict \\
        ${sample_npz} \\
        "\${REF_NPZ}" \\
        -o ${sample_id} \\
        --thresh-z ${thresh_z} \\
        --thresh-p ${thresh_p} \\
        ${ff_arg} \\
        2>&1 | tee ${sample_id}.gxcnv_predict.log

    # Validate required outputs
    for f in ${sample_id}_bins.tsv ${sample_id}_calls.tsv ${sample_id}_regions.tsv; do
        if [ ! -f "\$f" ]; then
            echo "ERROR: gxcnv predict missing output: \$f" >&2
            exit 1
        fi
    done
    """

    stub:
    """
    for suffix in bins segments calls regions qcmetrics; do
        echo -e "##gxcnv_version=0.1.0\\n#chrom\tstart\tend\tdual_call" > ${sample_id}_\${suffix}.tsv
    done
    echo -e "##gxcnv_version=0.1.0\\n#metric\tvalue\\npredicted_sex\tF" > ${sample_id}_sex.txt
    touch ${sample_id}.gxcnv_predict.log
    """
}


// ─────────────────────────────────────────────────────────────────────────────
// Step 4: Compare gx-cnv vs WisecondorX results
//         Produces a concordance summary for validation purposes.
// ─────────────────────────────────────────────────────────────────────────────
process GXCNV_COMPARE {
    tag "${sample_id}"
    label 'process_low'
    label 'nipt_docker'

    publishDir "${analysisdir}/${sample_id}/gxcnv", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id),
          path(gxcnv_calls),
          path(gxcnv_regions),
          path(wcx_aberrations)   // WisecondorX aberrations.bed or equivalent
    val   analysisdir

    output:
    tuple val(sample_id), path("${sample_id}.cnv_comparison.tsv"), emit: comparison
    path  "${sample_id}.cnv_comparison.log",                       emit: log

    script:
    """
    set -euo pipefail

    python3 - <<'PYEOF'
import sys
import csv
import logging

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [GXCNV_COMPARE] %(levelname)s %(message)s",
    handlers=[
        logging.FileHandler("${sample_id}.cnv_comparison.log"),
        logging.StreamHandler(sys.stdout),
    ],
)
logger = logging.getLogger(__name__)

sample_id = "${sample_id}"

# ── Parse gx-cnv HIGH_RISK calls ─────────────────────────────────────────────
gxcnv_calls = set()
try:
    with open("${gxcnv_calls}") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\\t")
            if len(parts) >= 4:
                region_name = parts[3] if len(parts) > 3 else f"{parts[0]}:{parts[1]}-{parts[2]}"
                gxcnv_calls.add(region_name)
    logger.info("gx-cnv HIGH_RISK calls: %s", gxcnv_calls)
except Exception as e:
    logger.warning("Could not parse gx-cnv calls: %s", e)

# ── Parse gx-cnv all regions ─────────────────────────────────────────────────
gxcnv_regions = {}
try:
    with open("${gxcnv_regions}") as f:
        reader = None
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                cols = line.lstrip("#").strip().split("\\t")
                reader = cols
                continue
            if reader is None:
                continue
            parts = line.strip().split("\\t")
            row = dict(zip(reader, parts))
            rname = row.get("region_name", "")
            if rname:
                gxcnv_regions[rname] = {
                    "gxcnv_dual_call": row.get("dual_call", ""),
                    "gxcnv_track_a_z": row.get("track_a_mean_z", "NA"),
                    "gxcnv_track_b_p": row.get("track_b_pvalue", "NA"),
                    "gxcnv_risk_pct":  row.get("risk_pct", "NA"),
                }
except Exception as e:
    logger.warning("Could not parse gx-cnv regions: %s", e)

# ── Parse WisecondorX aberrations ─────────────────────────────────────────────
wcx_calls = set()
try:
    with open("${wcx_aberrations}") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\\t")
            if len(parts) >= 4:
                wcx_calls.add(parts[3])
    logger.info("WisecondorX calls: %s", wcx_calls)
except Exception as e:
    logger.warning("Could not parse WisecondorX aberrations: %s", e)

# ── Concordance analysis ──────────────────────────────────────────────────────
all_regions = set(gxcnv_regions.keys()) | gxcnv_calls | wcx_calls

rows = []
for region in sorted(all_regions):
    gxcnv_call = "HIGH_RISK" if region in gxcnv_calls else "LOW_RISK"
    wcx_call   = "ABERRANT"  if region in wcx_calls   else "NORMAL"
    concordant = (gxcnv_call == "HIGH_RISK") == (wcx_call == "ABERRANT")
    status     = "CONCORDANT" if concordant else "DISCORDANT"

    region_info = gxcnv_regions.get(region, {})
    rows.append({
        "sample_id":       sample_id,
        "region":          region,
        "gxcnv_call":      gxcnv_call,
        "wcx_call":        wcx_call,
        "concordance":     status,
        "gxcnv_track_a_z": region_info.get("gxcnv_track_a_z", "NA"),
        "gxcnv_track_b_p": region_info.get("gxcnv_track_b_p", "NA"),
        "gxcnv_risk_pct":  region_info.get("gxcnv_risk_pct",  "NA"),
    })

# Summary stats
n_total       = len(rows)
n_concordant  = sum(1 for r in rows if r["concordance"] == "CONCORDANT")
n_gxcnv_only  = sum(1 for r in rows if r["gxcnv_call"] == "HIGH_RISK" and r["wcx_call"] == "NORMAL")
n_wcx_only    = sum(1 for r in rows if r["gxcnv_call"] == "LOW_RISK"  and r["wcx_call"] == "ABERRANT")

logger.info(
    "Concordance summary: total=%d concordant=%d gxcnv_only=%d wcx_only=%d",
    n_total, n_concordant, n_gxcnv_only, n_wcx_only,
)

# ── Write output ──────────────────────────────────────────────────────────────
fieldnames = ["sample_id", "region", "gxcnv_call", "wcx_call", "concordance",
              "gxcnv_track_a_z", "gxcnv_track_b_p", "gxcnv_risk_pct"]
with open("${sample_id}.cnv_comparison.tsv", "w", newline="") as out:
    out.write(f"## sample={sample_id}\\n")
    out.write(f"## n_regions={n_total}\\n")
    out.write(f"## n_concordant={n_concordant}\\n")
    out.write(f"## n_gxcnv_only={n_gxcnv_only}\\n")
    out.write(f"## n_wcx_only={n_wcx_only}\\n")
    writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\\t")
    writer.writeheader()
    writer.writerows(rows)

logger.info("Comparison written to ${sample_id}.cnv_comparison.tsv")
PYEOF
    """

    stub:
    """
    echo -e "## sample=${sample_id}\\n#sample_id\tregion\tgxcnv_call\twcx_call\tconcordance" > ${sample_id}.cnv_comparison.tsv
    touch ${sample_id}.cnv_comparison.log
    """
}
