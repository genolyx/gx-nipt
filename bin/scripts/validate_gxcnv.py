#!/usr/bin/env python3
"""
Batch validation: gx-cnv vs WisecondorX on existing GNCI samples.

For each sample it:
  1. gxcnv convert  BAM → NPZ
  2. GC inject      WIG → NPZ (polynomial GC correction)
  3. gxcnv predict  NPZ + reference → calls / figures
  4. Compare        gx-cnv calls vs WCX aberrations.bed

Writes:
  {out_dir}/validation_summary.tsv   – per-sample concordance metrics
  {out_dir}/validation_calls.tsv     – all HIGH_RISK gx-cnv calls
  {out_dir}/figures/                 – symlinks to per-sample PNG files
"""

import argparse
import csv
import glob
import logging
import multiprocessing as mp
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [validate_gxcnv] %(levelname)s %(message)s",
    stream=sys.stderr,
)
logger = logging.getLogger(__name__)

DOCKER_IMAGE = "gx-nipt:latest"
BIN_SIZE     = 100000
THRESH_Z     = -3.0
THRESH_P     = 0.05
REF_BASE     = "/home/ken/gx-nipt/refs/labs/cordlife/GXCNV"


# ── Helper: run shell command ──────────────────────────────────────────────────
def run(cmd, log_path=None, check=True):
    logger.debug("CMD: %s", " ".join(cmd))
    with open(log_path or os.devnull, "w") as flog:
        result = subprocess.run(cmd, stdout=flog, stderr=subprocess.STDOUT)
    if check and result.returncode != 0:
        raise RuntimeError(f"Command failed (rc={result.returncode}): {' '.join(cmd)}")
    return result.returncode


# ── Detect gender from gender.txt ─────────────────────────────────────────────
def detect_gender(gender_path):
    try:
        with open(gender_path) as f:
            for line in f:
                parts = line.strip().split()
                if parts and parts[0].lower() == "final_gender" and len(parts) >= 2:
                    g = parts[1].upper()
                    if g in ("XY", "MALE", "M"):
                        return "male"
                    if g in ("XX", "FEMALE", "F"):
                        return "female"
    except Exception:
        pass
    return "female"


# ── Parse WCX aberrations (coordinate-based) ──────────────────────────────────
def parse_wcx_calls(aber_path):
    """
    Returns list of (chrom, start, end, type) tuples.
    aberrations.bed: chr  start  end  ratio  zscore  type
    Header line starts with 'chr\t' (column name), skip it.
    """
    calls = []
    if not os.path.isfile(aber_path):
        return calls
    with open(aber_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("chr\t"):
                continue
            parts = line.split("\t")
            if len(parts) >= 3:
                try:
                    chrom = f"chr{parts[0]}" if not parts[0].startswith("chr") else parts[0]
                    start = int(parts[1])
                    end   = int(parts[2])
                    atype = parts[5] if len(parts) >= 6 else "unknown"
                    calls.append((chrom, start, end, atype))
                except ValueError:
                    continue
    return calls


# ── Parse gx-cnv regions TSV ──────────────────────────────────────────────────
def parse_gxcnv_calls(calls_path, regions_path):
    """
    Returns:
      high_risk: set of region_names called HIGH_RISK by gxcnv
      regions:   dict region_name -> row dict (includes chrom, start, end, scores)
    """
    high_risk = set()
    regions   = {}
    if os.path.isfile(calls_path):
        header = None
        with open(calls_path) as f:
            for line in f:
                if line.startswith("##"):
                    continue
                if line.startswith("#"):
                    header = line.lstrip("#").strip().split("\t")
                    continue
                if header is None:
                    continue
                parts = line.strip().split("\t")
                row   = dict(zip(header, parts))
                rname = row.get("region_name", "")
                if rname:
                    high_risk.add(rname)
    if os.path.isfile(regions_path):
        header = None
        with open(regions_path) as f:
            for line in f:
                if line.startswith("##"):
                    continue
                if line.startswith("#"):
                    header = line.lstrip("#").strip().split("\t")
                    continue
                if header is None:
                    continue
                parts = line.strip().split("\t")
                row   = dict(zip(header, parts))
                rname = row.get("region_name", "")
                if rname:
                    regions[rname] = row
    return high_risk, regions


# ── Overlap-based concordance: gxcnv region vs WCX aberrations ────────────────
def compare_by_overlap(gx_high_risk, gx_regions, wcx_calls):
    """
    For each gxcnv-monitored region, check whether any WCX aberration overlaps it.
    Returns per-region concordance list.
    """
    rows = []
    for rname, rinfo in gx_regions.items():
        try:
            r_chrom = rinfo.get("chrom", "")
            r_start = int(rinfo.get("start", 0))
            r_end   = int(rinfo.get("end",   0))
        except ValueError:
            continue

        # Check WCX overlap
        wcx_overlap = [
            c for c in wcx_calls
            if c[0] == r_chrom and c[2] > r_start and c[1] < r_end
        ]
        gxcnv_call = "HIGH_RISK" if rname in gx_high_risk else "LOW_RISK"
        wcx_call   = "ABERRANT"  if wcx_overlap            else "NORMAL"
        concordant = (gxcnv_call == "HIGH_RISK") == (wcx_call == "ABERRANT")

        rows.append({
            "region":      rname,
            "chrom":       r_chrom,
            "gxcnv_call":  gxcnv_call,
            "wcx_call":    wcx_call,
            "concordance": "CONCORDANT" if concordant else "DISCORDANT",
            "track_a_z":   rinfo.get("track_a_mean_z", "NA"),
            "track_b_p":   rinfo.get("track_b_pvalue",  "NA"),
            "risk_pct":    rinfo.get("risk_pct",          "NA"),
            "wcx_details": "; ".join(f"{c[0]}:{c[1]}-{c[2]}({c[3]})" for c in wcx_overlap),
        })
    return rows


# ── GC injection (Python, no Docker needed) ───────────────────────────────────
GC_INJECT_SCRIPT = r"""
import sys, numpy as np, pandas as pd

wig_path, npz_in, npz_out = sys.argv[1], sys.argv[2], sys.argv[3]

df = pd.read_csv(wig_path, sep='\t',
    names=['chr','start','end','reads','gc','map','valid','ideal','cor.gc','cor.map','copy'],
    header=0, dtype={'chr': str})
df['gc'] = pd.to_numeric(df['gc'], errors='coerce')
df = df[df['gc'] > 0]
df['s100k'] = (df['start'].astype(int) // 100_000) * 100_000
gc_map = df.groupby(['chr','s100k'])['gc'].mean()
gc_lookup = {(r[0], int(r[1])): r[2] for r in gc_map.reset_index().itertuples(index=False)}

d = dict(np.load(npz_in, allow_pickle=True))
bins   = d['bins']
chroms = list(d['chroms'])
counts = d['counts'].astype(float)
mask   = d['mask'].astype(bool)

gc_f = np.full(len(bins), np.nan, dtype=np.float32)
for i, (ci, s, e, _) in enumerate(bins):
    k = (chroms[int(ci)], int(s))
    if k in gc_lookup:
        gc_f[i] = gc_lookup[k]
bins[:, 3] = gc_f
d['bins'] = bins

poly_degree = 3
corrected   = np.full(len(counts), np.nan)
valid = mask & (counts > 0) & np.isfinite(gc_f)
if valid.sum() >= poly_degree + 1:
    gc_bins    = np.linspace(0, 1, 101)
    gc_idx     = np.clip(np.digitize(gc_f[valid], gc_bins) - 1, 0, 99)
    gc_centers = (gc_bins[:-1] + gc_bins[1:]) / 2
    gc_medians = np.full(100, np.nan)
    for i in range(100):
        ib = valid.nonzero()[0][gc_idx == i]
        if len(ib) >= 3:
            gc_medians[i] = np.median(counts[ib])
    finite = np.isfinite(gc_medians)
    if finite.sum() >= poly_degree + 1:
        coeffs    = np.polyfit(gc_centers[finite], gc_medians[finite], poly_degree)
        predicted = np.polyval(coeffs, gc_f[valid])
        predicted = np.where(predicted <= 0, np.nan, predicted)
        gm = np.nanmedian(gc_medians[finite])
        adj = counts[valid] / predicted * gm
        corrected[valid] = np.where(np.isfinite(adj), adj, np.nan)
    else:
        corrected[valid] = counts[valid]
else:
    corrected[valid] = counts[valid]
s = np.nansum(corrected[mask])
if s > 0:
    corrected = corrected / s * mask.sum()
d['corrected'] = corrected
np.savez_compressed(npz_out, **d)
print(f'[gc_inject] saved {npz_out}.npz')
"""


# ── Per-sample worker ─────────────────────────────────────────────────────────
def process_sample(args):
    (sample, bam, aber, wig, gender_path, tmp_dir, ref_base, n) = args
    sample_dir = Path(tmp_dir) / sample
    sample_dir.mkdir(parents=True, exist_ok=True)

    gender  = detect_gender(gender_path)
    ref_npz = f"{ref_base}/{gender}/reference.npz"
    if not os.path.isfile(ref_npz):
        return {"sample": sample, "error": f"Reference not found: {ref_npz}"}

    npz_raw = str(sample_dir / f"{sample}.gxcnv.npz")
    npz_gc  = str(sample_dir / f"{sample}.gxcnv.gc")
    prefix  = str(sample_dir / sample)
    log_cvt = str(sample_dir / "convert.log")
    log_gc  = str(sample_dir / "gc_inject.log")
    log_prd = str(sample_dir / "predict.log")

    try:
        # Step 1: convert BAM → NPZ (via Docker)
        run([
            "docker", "run", "--rm",
            "-v", f"{os.path.dirname(bam)}:{os.path.dirname(bam)}:ro",
            "-v", f"{str(sample_dir)}:{str(sample_dir)}:rw",
            DOCKER_IMAGE,
            "gxcnv", "convert", bam, npz_raw,
            "--bin-size", str(BIN_SIZE), "--min-mapq", "1",
        ], log_path=log_cvt)

        # Step 2: GC inject (local Python — faster than Docker overhead)
        gc_script = str(sample_dir / "gc_inject.py")
        with open(gc_script, "w") as f:
            f.write(GC_INJECT_SCRIPT)
        run(["python3", gc_script, wig, npz_raw, npz_gc], log_path=log_gc)

        # Step 3: predict (via Docker)
        run([
            "docker", "run", "--rm",
            "-v", f"{str(sample_dir)}:{str(sample_dir)}:rw",
            "-v", f"{ref_base}:{ref_base}:ro",
            DOCKER_IMAGE,
            "gxcnv", "predict",
            f"{npz_gc}.npz", ref_npz,
            "-o", prefix,
            "--thresh-z", str(THRESH_Z),
            "--thresh-p", str(THRESH_P),
        ], log_path=log_prd)

        # Step 4: compare (coordinate-based overlap)
        calls_path   = f"{prefix}_calls.tsv"
        regions_path = f"{prefix}_regions.tsv"
        gx_high_risk, gx_regions = parse_gxcnv_calls(calls_path, regions_path)
        wcx_calls = parse_wcx_calls(aber)

        region_rows  = compare_by_overlap(gx_high_risk, gx_regions, wcx_calls)
        n_total      = len(region_rows)
        n_concordant = sum(1 for r in region_rows if r["concordance"] == "CONCORDANT")
        n_gxcnv_only = sum(1 for r in region_rows if r["gxcnv_call"] == "HIGH_RISK" and r["wcx_call"] == "NORMAL")
        n_wcx_only   = sum(1 for r in region_rows if r["gxcnv_call"] == "LOW_RISK"  and r["wcx_call"] == "ABERRANT")
        concordance  = round(n_concordant / max(n_total, 1) * 100, 1)

        call_details = []
        for r in region_rows:
            if r["gxcnv_call"] == "HIGH_RISK" or r["wcx_call"] == "ABERRANT":
                call_details.append({
                    "sample":      sample,
                    "region":      r["region"],
                    "gxcnv":       r["gxcnv_call"],
                    "wcx":         r["wcx_call"],
                    "concordance": r["concordance"],
                    "track_a_z":   r["track_a_z"],
                    "track_b_p":   r["track_b_p"],
                    "risk_pct":    r["risk_pct"],
                    "wcx_details": r["wcx_details"],
                })

        return {
            "sample":      sample,
            "gender":      gender,
            "n_total":     n_total,
            "n_concordant": n_concordant,
            "n_gxcnv_only": n_gxcnv_only,
            "n_wcx_only":   n_wcx_only,
            "concordance":  concordance,
            "calls":        call_details,
            "figures": {
                "genome":  f"{prefix}_genome.png",
                "regions": f"{prefix}_regions.png",
                "qc":      f"{prefix}_qc.png",
            },
            "error": None,
        }

    except Exception as e:
        return {"sample": sample, "error": str(e)}


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    ap = argparse.ArgumentParser(description="Batch gx-cnv vs WCX validation")
    ap.add_argument("--analysis-dir",  default="/home/ken/ken-nipt/analysis")
    ap.add_argument("--ref-base",      default=REF_BASE)
    ap.add_argument("--out-dir",       default="/home/ken/gx-nipt/refs/gxcnv/validation")
    ap.add_argument("--max-samples",   type=int, default=200,
                    help="Max samples to validate (0=all)")
    ap.add_argument("--workers",       type=int, default=16)
    ap.add_argument("--tmp-dir",       default=None,
                    help="Temp dir for NPZ/figures (default: auto)")
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = out_dir / "figures"
    fig_dir.mkdir(exist_ok=True)

    # Collect samples
    logger.info("Scanning for samples in %s ...", args.analysis_dir)
    tasks = []
    for bam in sorted(glob.glob(f"{args.analysis_dir}/**/GNCI*.proper_paired.bam",
                                recursive=True)):
        parts   = bam.replace(args.analysis_dir + "/", "").split("/")
        work_dir, sample = parts[0], parts[1]
        adir    = f"{args.analysis_dir}/{work_dir}/{sample}"
        aber    = f"{adir}/Output_WCX/orig/{sample}.wcx.orig_aberrations.bed"
        wig     = f"{adir}/Output_hmmcopy/{sample}.of_orig.50kb.wig.Normalization.txt"
        gender  = f"{adir}/Output_FF/{sample}.gender.txt"
        if os.path.isfile(aber) and os.path.isfile(wig) and os.path.isfile(gender):
            tasks.append((sample, bam, aber, wig, gender))

    if args.max_samples > 0:
        tasks = tasks[:args.max_samples]
    logger.info("Processing %d samples with %d workers", len(tasks), args.workers)

    tmp_root = args.tmp_dir or tempfile.mkdtemp(prefix="gxcnv_val_")
    logger.info("Temp dir: %s", tmp_root)

    worker_args = [
        (s, b, a, w, g, tmp_root, args.ref_base, i)
        for i, (s, b, a, w, g) in enumerate(tasks)
    ]

    results = []
    with mp.Pool(processes=args.workers) as pool:
        for i, r in enumerate(pool.imap_unordered(process_sample, worker_args), 1):
            results.append(r)
            status = r.get("error") or f"concordance={r.get('concordance')}%"
            logger.info("[%d/%d] %s: %s", i, len(tasks), r["sample"], status)

    # Write summary TSV
    summary_path = out_dir / "validation_summary.tsv"
    with open(summary_path, "w", newline="") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=[
            "sample","gender","n_total","n_concordant","n_gxcnv_only",
            "n_wcx_only","concordance","error",
        ])
        w.writeheader()
        for r in sorted(results, key=lambda x: x["sample"]):
            w.writerow({k: r.get(k,"") for k in w.fieldnames})

    # Write calls TSV
    calls_path = out_dir / "validation_calls.tsv"
    all_calls  = [c for r in results for c in r.get("calls", [])]
    if all_calls:
        with open(calls_path, "w", newline="") as f:
            w = csv.DictWriter(f, delimiter="\t",
                               fieldnames=["sample","region","gxcnv","wcx",
                                           "concordance","track_a_z","track_b_p",
                                           "risk_pct","wcx_details"])
            w.writeheader()
            w.writerows(all_calls)

    # Symlink figures
    for r in results:
        if r.get("error"):
            continue
        for key, src in (r.get("figures") or {}).items():
            if src and os.path.isfile(src):
                dst = fig_dir / f"{r['sample']}_{key}.png"
                if not dst.exists():
                    dst.symlink_to(src)

    # Print summary statistics
    ok      = [r for r in results if not r.get("error")]
    errors  = [r for r in results if r.get("error")]
    concordances = [r["concordance"] for r in ok]
    avg_conc = sum(concordances) / max(len(concordances), 1)
    total_calls  = sum(r["n_gxcnv_only"] + r["n_wcx_only"] for r in ok)
    n_wcx_only   = sum(r["n_wcx_only"] for r in ok)
    n_gxcnv_only = sum(r["n_gxcnv_only"] for r in ok)

    logger.info("=" * 60)
    logger.info("Validation complete: %d samples processed, %d errors",
                len(results), len(errors))
    logger.info("Average concordance: %.1f%%", avg_conc)
    logger.info("WCX-only calls (gxcnv missed): %d", n_wcx_only)
    logger.info("gxcnv-only calls (WCX missed): %d", n_gxcnv_only)
    logger.info("Summary: %s", summary_path)
    logger.info("Calls:   %s", calls_path)
    logger.info("Figures: %s", fig_dir)

    if errors:
        logger.warning("Failed samples:")
        for r in errors:
            logger.warning("  %s: %s", r["sample"], r["error"])


if __name__ == "__main__":
    main()
