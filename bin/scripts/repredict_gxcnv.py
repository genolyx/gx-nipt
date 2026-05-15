#!/usr/bin/env python3
"""
Re-run gxcnv predict on existing GC-injected NPZ files with different thresholds.
Much faster than full pipeline since convert+GC-inject steps are skipped.
"""
import argparse, csv, glob, os, subprocess, sys, multiprocessing as mp

DOCKER_IMAGE = "gx-nipt:latest"
REF_BASE     = "/home/ken/gx-nipt/refs/labs/cordlife/GXCNV"


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


def parse_wcx_calls(aber_path):
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
                    calls.append((chrom, int(parts[1]), int(parts[2]),
                                  parts[5] if len(parts) >= 6 else "unknown"))
                except ValueError:
                    pass
    return calls


def parse_gxcnv_regions(calls_path, regions_path):
    high_risk = set()
    regions   = {}
    for path, is_calls in [(calls_path, True), (regions_path, False)]:
        if not os.path.isfile(path):
            continue
        header = None
        with open(path) as f:
            for line in f:
                if line.startswith("##"):
                    continue
                if line.startswith("#"):
                    header = line.lstrip("#").strip().split("\t")
                    continue
                if header is None:
                    continue
                row = dict(zip(header, line.strip().split("\t")))
                rname = row.get("region_name", "")
                if rname:
                    if is_calls:
                        high_risk.add(rname)
                    else:
                        regions[rname] = row
    return high_risk, regions


def compare(gx_high_risk, gx_regions, wcx_calls):
    rows = []
    for rname, ri in gx_regions.items():
        try:
            rc = ri.get("chrom", "")
            rs = int(ri.get("start", 0))
            re = int(ri.get("end",   0))
        except ValueError:
            continue
        overlap  = [c for c in wcx_calls if c[0] == rc and c[2] > rs and c[1] < re]
        gc = "HIGH_RISK" if rname in gx_high_risk else "LOW_RISK"
        wc = "ABERRANT"  if overlap               else "NORMAL"
        rows.append({
            "region": rname, "gxcnv": gc, "wcx": wc,
            "concordance": "CONCORDANT" if (gc == "HIGH_RISK") == (wc == "ABERRANT") else "DISCORDANT",
            "track_a_z": ri.get("track_a_mean_z", "NA"),
            "track_b_p": ri.get("track_b_pvalue",  "NA"),
            "risk_pct":  ri.get("risk_pct", "NA"),
        })
    return rows


def worker(args):
    sample, npz_gc, aber, gender_path, ref_base, thresh_z, thresh_p, out_prefix = args
    gender  = detect_gender(gender_path)
    ref_npz = f"{ref_base}/{gender}/reference.npz"
    npz_dir = os.path.dirname(npz_gc)
    log     = f"{out_prefix}.predict.log"
    try:
        rc = subprocess.run([
            "docker", "run", "--rm",
            "-v", f"{npz_dir}:{npz_dir}:rw",
            "-v", f"{ref_base}:{ref_base}:ro",
            DOCKER_IMAGE,
            "gxcnv", "predict", npz_gc, ref_npz,
            "-o", out_prefix,
            "--thresh-z", str(thresh_z),
            "--thresh-p", str(thresh_p),
        ], stdout=open(log, "w"), stderr=subprocess.STDOUT).returncode

        if rc != 0:
            return {"sample": sample, "error": f"predict rc={rc}"}

        gx_hr, gx_reg = parse_gxcnv_regions(f"{out_prefix}_calls.tsv",
                                              f"{out_prefix}_regions.tsv")
        wcx = parse_wcx_calls(aber)
        rows = compare(gx_hr, gx_reg, wcx)
        n_total      = len(rows)
        n_concordant = sum(1 for r in rows if r["concordance"] == "CONCORDANT")
        n_gx_only    = sum(1 for r in rows if r["gxcnv"] == "HIGH_RISK" and r["wcx"] == "NORMAL")
        n_wcx_only   = sum(1 for r in rows if r["gxcnv"] == "LOW_RISK"  and r["wcx"] == "ABERRANT")
        conc_pct     = round(n_concordant / max(n_total, 1) * 100, 1)

        discordant = [dict(sample=sample, **r) for r in rows if r["concordance"] == "DISCORDANT"]
        calls      = [dict(sample=sample, **r) for r in rows if r["gxcnv"] == "HIGH_RISK"]
        return {
            "sample": sample, "gender": gender,
            "n_total": n_total, "n_concordant": n_concordant,
            "n_gxcnv_only": n_gx_only, "n_wcx_only": n_wcx_only,
            "concordance": conc_pct,
            "discordant": discordant, "calls": calls,
            "error": None,
        }
    except Exception as e:
        return {"sample": sample, "error": str(e)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tmp-dir",    default="/tmp/gxcnv_val_zl2i4t8t")
    ap.add_argument("--ref-base",   default=REF_BASE)
    ap.add_argument("--analysis-dir", default="/home/ken/ken-nipt/analysis")
    ap.add_argument("--out-dir",    default="/home/ken/gx-nipt/refs/gxcnv/validation_z25")
    ap.add_argument("--thresh-z",   type=float, default=-2.5)
    ap.add_argument("--thresh-p",   type=float, default=0.05)
    ap.add_argument("--workers",    type=int,   default=20)
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    import glob as _glob, logging
    logging.basicConfig(level=logging.INFO, stream=sys.stderr,
                        format="%(asctime)s %(levelname)s %(message)s")
    log = logging.getLogger()

    # Collect tasks from existing tmp dir
    tasks = []
    for sample_dir in sorted(_glob.glob(f"{args.tmp_dir}/GNCI*")):
        sample   = os.path.basename(sample_dir)
        npz_gc   = f"{sample_dir}/{sample}.gxcnv.gc.npz"
        if not os.path.isfile(npz_gc):
            continue
        aber = _glob.glob(f"{args.analysis_dir}/*/{sample}/Output_WCX/orig/{sample}.wcx.orig_aberrations.bed")
        gend = _glob.glob(f"{args.analysis_dir}/*/{sample}/Output_FF/{sample}.gender.txt")
        if not aber or not gend:
            continue
        out_prefix = f"{sample_dir}/{sample}_z25"
        tasks.append((sample, npz_gc, aber[0], gend[0], args.ref_base,
                      args.thresh_z, args.thresh_p, out_prefix))

    log.info("Re-predicting %d samples with thresh_z=%.1f ...", len(tasks), args.thresh_z)

    results = []
    with mp.Pool(args.workers) as pool:
        for i, r in enumerate(pool.imap_unordered(worker, tasks), 1):
            results.append(r)
            status = r.get("error") or f"concordance={r.get('concordance')}%"
            log.info("[%d/%d] %s: %s", i, len(tasks), r["sample"], status)

    # Summary
    ok = [r for r in results if not r.get("error")]
    avg_conc = sum(r["concordance"] for r in ok) / max(len(ok), 1)
    n_wcx_only   = sum(r["n_wcx_only"]   for r in ok)
    n_gxcnv_only = sum(r["n_gxcnv_only"] for r in ok)
    n_perfect    = sum(1 for r in ok if r["concordance"] == 100.0)

    log.info("="*60)
    log.info("thresh_z=%.1f  samples=%d  avg_concordance=%.1f%%",
             args.thresh_z, len(ok), avg_conc)
    log.info("100%% concordant: %d/%d", n_perfect, len(ok))
    log.info("WCX-only (missed by gxcnv): %d", n_wcx_only)
    log.info("gxcnv-only (not in WCX):    %d", n_gxcnv_only)

    # Write summary
    with open(f"{args.out_dir}/summary.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, delimiter="\t",
            fieldnames=["sample","gender","n_total","n_concordant",
                        "n_gxcnv_only","n_wcx_only","concordance","error"])
        w.writeheader()
        for r in sorted(results, key=lambda x: x["sample"]):
            w.writerow({k: r.get(k,"") for k in w.fieldnames})

    # Write discordant
    all_disc = [d for r in ok for d in r.get("discordant", [])]
    if all_disc:
        with open(f"{args.out_dir}/discordant.tsv", "w", newline="") as f:
            w = csv.DictWriter(f, delimiter="\t",
                fieldnames=["sample","region","gxcnv","wcx","concordance",
                            "track_a_z","track_b_p","risk_pct"])
            w.writeheader()
            w.writerows(all_disc)
        log.info("Discordant cases (%d):", len(all_disc))
        for d in all_disc:
            log.info("  %s | %s | gxcnv=%s wcx=%s z=%s p=%s risk=%s%%",
                     d["sample"], d["region"], d["gxcnv"], d["wcx"],
                     d["track_a_z"], d["track_b_p"], d["risk_pct"])

    # Write HIGH_RISK calls
    all_calls = [c for r in ok for c in r.get("calls", [])]
    if all_calls:
        with open(f"{args.out_dir}/calls.tsv", "w", newline="") as f:
            w = csv.DictWriter(f, delimiter="\t",
                fieldnames=["sample","region","gxcnv","wcx","concordance",
                            "track_a_z","track_b_p","risk_pct"])
            w.writeheader()
            w.writerows(all_calls)
        log.info("gxcnv HIGH_RISK calls: %d", len(all_calls))
        for c in all_calls:
            log.info("  %s | %s | wcx=%s z=%s p=%s risk=%s%%",
                     c["sample"], c["region"], c["wcx"],
                     c["track_a_z"], c["track_b_p"], c["risk_pct"])


if __name__ == "__main__":
    main()
