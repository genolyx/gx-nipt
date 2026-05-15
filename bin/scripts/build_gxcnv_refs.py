#!/usr/bin/env python3
"""
Build gx-cnv reference panels for fetus and mom groups.

For each labcode/gender/group combination:
  1. Collect existing GC-corrected NPZ files (from ken-nipt analysis dirs)
  2. gxcnv newref  NPZ files → reference.npz

Reference layout:
  {ref_dir}/labs/{labcode}/GXCNV/{gender}/reference.npz          ← orig (already built)
  {ref_dir}/labs/{labcode}/GXCNV/{gender}_fetus/reference.npz    ← fetus
  {ref_dir}/labs/{labcode}/GXCNV/{gender}_mom/reference.npz      ← mom

NPZ files are created on the fly using gxcnv convert + GC inject if not already present.
"""

import argparse
import glob
import logging
import os
import subprocess
import sys
import tempfile
import multiprocessing as mp
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [build_gxcnv_refs] %(levelname)s %(message)s",
    stream=sys.stderr,
)
logger = logging.getLogger(__name__)

DOCKER_IMAGE = "gx-nipt:latest"
BIN_SIZE     = 100000


# ── GC injection script (same as validate_gxcnv.py) ──────────────────────────
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
bins = d['bins']; chroms = list(d['chroms'])
counts = d['counts'].astype(float); mask = d['mask'].astype(bool)
gc_f = np.full(len(bins), np.nan, dtype=np.float32)
for i, (ci, s, e, _) in enumerate(bins):
    k = (chroms[int(ci)], int(s))
    if k in gc_lookup: gc_f[i] = gc_lookup[k]
bins[:, 3] = gc_f; d['bins'] = bins
poly_degree = 3; corrected = np.full(len(counts), np.nan)
valid = mask & (counts > 0) & np.isfinite(gc_f)
if valid.sum() >= poly_degree + 1:
    gc_bins = np.linspace(0, 1, 101); gc_idx = np.clip(np.digitize(gc_f[valid], gc_bins)-1,0,99)
    gc_centers = (gc_bins[:-1]+gc_bins[1:])/2; gc_medians = np.full(100, np.nan)
    for i in range(100):
        ib = valid.nonzero()[0][gc_idx==i]
        if len(ib)>=3: gc_medians[i] = np.median(counts[ib])
    finite = np.isfinite(gc_medians)
    if finite.sum() >= poly_degree+1:
        coeffs = np.polyfit(gc_centers[finite], gc_medians[finite], poly_degree)
        predicted = np.polyval(coeffs, gc_f[valid])
        predicted = np.where(predicted<=0, np.nan, predicted)
        gm = np.nanmedian(gc_medians[finite])
        adj = counts[valid]/predicted*gm
        corrected[valid] = np.where(np.isfinite(adj), adj, np.nan)
    else: corrected[valid] = counts[valid]
else: corrected[valid] = counts[valid]
s = np.nansum(corrected[mask])
if s > 0: corrected = corrected/s*mask.sum()
d['corrected'] = corrected
np.savez_compressed(npz_out, **d)
"""


def detect_gender(gender_path):
    try:
        with open(gender_path) as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 3 and parts[0] in ("gd_2", "final_gender"):
                    g = parts[2].upper()
                    return "male" if g in ("XY", "MALE", "M") else "female"
    except Exception:
        pass
    return "female"


def convert_and_gc(bam, wig, sample_id, out_dir):
    """Convert BAM → GC-injected NPZ. Returns path or None on failure."""
    npz_raw = f"{out_dir}/{sample_id}.gxcnv.npz"
    npz_gc  = f"{out_dir}/{sample_id}.gxcnv.gc"
    log_c   = f"{out_dir}/{sample_id}.convert.log"
    log_g   = f"{out_dir}/{sample_id}.gc.log"

    # Convert
    rc = subprocess.run([
        "docker", "run", "--rm",
        "-v", f"{os.path.dirname(bam)}:{os.path.dirname(bam)}:ro",
        "-v", f"{out_dir}:{out_dir}:rw",
        DOCKER_IMAGE,
        "gxcnv", "convert", bam, npz_raw,
        "--bin-size", str(BIN_SIZE), "--min-mapq", "1",
    ], stdout=open(log_c, "w"), stderr=subprocess.STDOUT).returncode
    if rc != 0:
        return None

    # GC inject
    script_path = f"{out_dir}/{sample_id}_gc_inject.py"
    with open(script_path, "w") as f:
        f.write(GC_INJECT_SCRIPT)
    rc = subprocess.run(
        ["python3", script_path, wig, npz_raw, npz_gc],
        stdout=open(log_g, "w"), stderr=subprocess.STDOUT
    ).returncode
    if rc != 0:
        return None

    return f"{npz_gc}.npz"


def _worker_convert(args):
    sample, bam, wig, out_dir = args
    result = convert_and_gc(bam, wig, sample, out_dir)
    return sample, result


def build_reference(npz_list, ref_path, pca_variance=0.95):
    """Run gxcnv newref on a list of NPZ files."""
    npz_dir  = os.path.dirname(npz_list[0])
    ref_name = os.path.splitext(ref_path)[0]
    log_path = ref_name + ".newref.log"
    os.makedirs(os.path.dirname(ref_path), exist_ok=True)

    rc = subprocess.run([
        "docker", "run", "--rm",
        "-v", f"{npz_dir}:{npz_dir}:rw",
        "-v", f"{os.path.dirname(ref_path)}:{os.path.dirname(ref_path)}:rw",
        DOCKER_IMAGE,
        "gxcnv", "newref",
        *npz_list,
        "-o", ref_name,
        "--pca-variance", str(pca_variance),
    ], stdout=open(log_path, "w"), stderr=subprocess.STDOUT).returncode

    if rc != 0:
        logger.error("gxcnv newref failed (rc=%d): %s", rc, log_path)
        return False
    logger.info("Reference built: %s  (n=%d samples)", ref_path, len(npz_list))
    return True


def main():
    ap = argparse.ArgumentParser(description="Build gxcnv references for fetus/mom groups")
    ap.add_argument("--labcode",      default="cordlife")
    ap.add_argument("--ref-dir",      default="/home/ken/gx-nipt/refs")
    ap.add_argument("--analysis-dir", default="/home/ken/ken-nipt/analysis")
    ap.add_argument("--groups",       nargs="+", default=["fetus", "mom"],
                    help="Groups to build (default: fetus mom)")
    ap.add_argument("--min-samples",  type=int, default=50,
                    help="Min samples per gender per group for reference")
    ap.add_argument("--max-samples",  type=int, default=150)
    ap.add_argument("--workers",      type=int, default=20)
    ap.add_argument("--pca-variance", type=float, default=0.95)
    ap.add_argument("--tmp-dir",      default=None)
    args = ap.parse_args()

    tmp_root = args.tmp_dir or tempfile.mkdtemp(prefix="gxcnv_ref_")
    logger.info("Temp dir: %s", tmp_root)

    # BAM suffix per group
    bam_suffix = {
        "orig":  ".proper_paired.bam",
        "fetus": ".of_fetus.bam",
        "mom":   ".of_mom.bam",
    }

    for group in args.groups:
        suffix = bam_suffix.get(group)
        if not suffix:
            logger.warning("Unknown group %s, skipping", group)
            continue

        logger.info("=== Building references for group: %s ===", group)
        tmp_grp = os.path.join(tmp_root, group)
        os.makedirs(tmp_grp, exist_ok=True)

        # Collect samples with all required files
        by_gender = {"male": [], "female": []}
        for bam in sorted(glob.glob(
            f"{args.analysis_dir}/**/GNCI*{suffix}", recursive=True
        )):
            parts    = bam.replace(args.analysis_dir + "/", "").split("/")
            work_dir, sample = parts[0], parts[1]
            adir     = f"{args.analysis_dir}/{work_dir}/{sample}"
            wig      = f"{adir}/Output_hmmcopy/{sample}.of_orig.50kb.wig.Normalization.txt"
            gender_f = f"{adir}/Output_FF/{sample}.gender.txt"
            if not (os.path.isfile(bam) and os.path.isfile(wig) and os.path.isfile(gender_f)):
                continue
            gender = detect_gender(gender_f)
            by_gender[gender].append((sample, bam, wig))

        logger.info("  Found: male=%d  female=%d",
                    len(by_gender["male"]), len(by_gender["female"]))

        for gender, samples in by_gender.items():
            ref_path = (f"{args.ref_dir}/labs/{args.labcode}/GXCNV"
                        f"/{gender}_{group}/reference.npz")

            if os.path.isfile(ref_path):
                logger.info("  [%s/%s] Reference already exists: %s", gender, group, ref_path)
                continue

            sel = samples[:args.max_samples]
            if len(sel) < args.min_samples:
                logger.warning("  [%s/%s] Only %d samples (need %d) — skipping",
                                gender, group, len(sel), args.min_samples)
                continue

            logger.info("  [%s/%s] Converting %d BAMs ...", gender, group, len(sel))
            tmp_gender = os.path.join(tmp_grp, gender)
            os.makedirs(tmp_gender, exist_ok=True)

            tasks = [(s, b, w, tmp_gender) for s, b, w in sel]
            npz_files = []
            with mp.Pool(args.workers) as pool:
                for i, (sid, npz) in enumerate(
                    pool.imap_unordered(_worker_convert, tasks), 1
                ):
                    if npz:
                        npz_files.append(npz)
                        logger.info("  [%d/%d] %s → OK", i, len(tasks), sid)
                    else:
                        logger.warning("  [%d/%d] %s → FAILED", i, len(tasks), sid)

            logger.info("  [%s/%s] %d/%d NPZ ready — building reference ...",
                        gender, group, len(npz_files), len(sel))

            if len(npz_files) < args.min_samples:
                logger.warning("  Not enough NPZ files (%d < %d) — skipping reference build",
                               len(npz_files), args.min_samples)
                continue

            build_reference(npz_files, ref_path, args.pca_variance)

    logger.info("Done. Reference layout:")
    for f in sorted(glob.glob(
        f"{args.ref_dir}/labs/{args.labcode}/GXCNV/**/reference.npz", recursive=True
    )):
        logger.info("  %s", f)


if __name__ == "__main__":
    main()
