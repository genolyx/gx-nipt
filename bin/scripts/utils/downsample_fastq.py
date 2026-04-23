#!/usr/bin/env python3
"""
downsample_fastq.py
-------------------

Nextflow-facing helper that conditionally downsamples a paired FASTQ set.

* If the total read count (R1+R2) exceeds ``QC.max_fq_size`` in the lab
  config, the pair is subsampled to ``QC.downsample_size`` reads per mate
  with ``seqtk sample -s`` (seeded, reproducible across runs).
* Otherwise the input files are re-published verbatim under the expected
  output names.

CLI
~~~
    --r1      <path>                  forward read file (.fastq[.gz])
    --r2      <path>                  reverse read file (.fastq[.gz])
    --config  <pipeline_config.json>  lab pipeline configuration
    --sample  <sample_id>
    --outdir  <dir>                   output directory (created if missing)

Outputs (always gzipped):
    <outdir>/<sample>_R1.fastq.gz
    <outdir>/<sample>_R2.fastq.gz
    <outdir>/<sample>.fastq_check.txt
"""

import argparse
import gzip
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path


def _read_count(fastq_path: Path) -> int:
    """Count reads in a (possibly gzipped) FASTQ file."""
    opener = gzip.open if str(fastq_path).endswith(".gz") else open
    lines = 0
    with opener(str(fastq_path), "rt") as fh:
        for _ in fh:
            lines += 1
    return lines // 4


def _run(cmd: str) -> None:
    print(f"[downsample] $ {cmd}", flush=True)
    subprocess.check_call(cmd, shell=True, executable="/bin/bash")


def _seqtk_sample(src: Path, dst: Path, n: int, seed: int = 100) -> None:
    """Sample n reads with seqtk (reproducible) and gzip the result."""
    _run(f"seqtk sample -s {seed} {src} {n} | gzip -c > {dst}")


def _passthrough(src: Path, dst: Path) -> None:
    """Re-publish src as gzipped dst without modifying read count."""
    if str(src).endswith(".gz"):
        shutil.copyfile(src, dst)
    else:
        _run(f"gzip -c {src} > {dst}")


def main() -> int:
    ap = argparse.ArgumentParser(description="Conditional FASTQ downsampler")
    ap.add_argument("--r1", required=True, type=Path)
    ap.add_argument("--r2", required=True, type=Path)
    ap.add_argument("--config", required=True, type=Path)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--outdir", required=True, type=Path)
    args = ap.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    with open(args.config) as fh:
        config = json.load(fh)
    qc_cfg = config.get("QC", {})
    max_fq_size = int(qc_cfg.get("max_fq_size", 0))
    ds_size = int(qc_cfg.get("downsample_size", 0))

    r1_reads = _read_count(args.r1)
    r2_reads = _read_count(args.r2)
    total = r1_reads + r2_reads
    print(f"[downsample] {args.sample}: r1={r1_reads} r2={r2_reads} total={total} "
          f"max={max_fq_size} target={ds_size}")

    out_r1 = args.outdir / f"{args.sample}_R1.fastq.gz"
    out_r2 = args.outdir / f"{args.sample}_R2.fastq.gz"
    check_txt = args.outdir / f"{args.sample}.fastq_check.txt"

    if total == 0:
        with open(check_txt, "w") as fh:
            fh.write(f"0. FastQ size ({total}) : FAIL - empty input\n")
        print("[downsample] ERROR: empty FASTQ(s)", file=sys.stderr)
        return 1

    if max_fq_size > 0 and total > max_fq_size and ds_size > 0:
        # Downsample both mates with the same seed so pairs stay aligned.
        _seqtk_sample(args.r1, out_r1, ds_size)
        _seqtk_sample(args.r2, out_r2, ds_size)
        status = f"0. FastQ size ({total}) has been downsampled to {ds_size * 2}\n"
    else:
        _passthrough(args.r1, out_r1)
        _passthrough(args.r2, out_r2)
        status = f"0. FastQ size ({total}) : PASS\n"

    with open(check_txt, "w") as fh:
        fh.write(status)
    print(f"[downsample] {status.strip()}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
