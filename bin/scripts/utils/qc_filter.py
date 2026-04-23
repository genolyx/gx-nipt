#!/usr/bin/env python3
"""
qc_filter.py
------------

Apply per-lab QC thresholds (from ``pipeline_config.json`` > ``QC``)
against the parsed Qualimap summary emitted by parse_qualimap.py and
produce a single PASS/FAIL verdict file.

CLI
~~~
    --qc-txt  <sample>.qc.txt           parse_qualimap.py output
    --config  <pipeline_config.json>    lab pipeline configuration
    --sample  <sample_id>
    --output  <sample>.qc.filter.txt

Output: tab-separated, one check per line, plus a final ``overall`` row:

    <metric>    <value>    <threshold>    <comparison>    <status>
    ...
    overall     <n_fail>   fail_count     <=0             PASS|FAIL

``status`` is ``PASS`` / ``FAIL`` / ``SKIP`` (when either side is NA).
The script always exits 0 so downstream Nextflow processes can inspect
the file directly; exit code is reserved for unexpected errors.
"""

import argparse
import json
import sys
from pathlib import Path


def _num(value: str):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _check(name: str, actual, threshold, comparator: str):
    """Return (status, actual_str, threshold_str, comparator)."""
    if actual is None or threshold is None:
        return ("SKIP", "NA" if actual is None else str(actual),
                "NA" if threshold is None else str(threshold), comparator)
    if comparator == ">=":
        ok = actual >= threshold
    elif comparator == "<=":
        ok = actual <= threshold
    elif comparator == "between":
        low, high = threshold
        ok = low <= actual <= high
        return ("PASS" if ok else "FAIL",
                f"{actual}", f"{low}..{high}", "between")
    else:
        raise ValueError(f"Unknown comparator {comparator}")
    return ("PASS" if ok else "FAIL", f"{actual}", f"{threshold}", comparator)


def main() -> int:
    ap = argparse.ArgumentParser(description="Apply QC thresholds from lab config")
    ap.add_argument("--qc-txt", required=True, type=Path)
    ap.add_argument("--config", required=True, type=Path)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--output", required=True, type=Path)
    args = ap.parse_args()

    if not args.qc_txt.exists():
        print(f"[qc_filter] ERROR: missing QC input: {args.qc_txt}",
              file=sys.stderr)
        return 1

    metrics = {}
    for line in args.qc_txt.read_text().splitlines():
        if not line.strip() or "\t" not in line:
            continue
        key, value = line.split("\t", 1)
        metrics[key.strip()] = value.strip()

    with open(args.config) as fh:
        qc_cfg = json.load(fh).get("QC", {})

    checks = [
        ("number_of_reads",
         _num(metrics.get("number_of_reads")),
         _num(qc_cfg.get("number_of_reads")),
         ">="),
        ("number_of_mapped_reads",
         _num(metrics.get("number_of_mapped_reads")),
         _num(qc_cfg.get("number_of_mapped_reads")),
         ">="),
        ("mapping_rate",
         _num(metrics.get("mapping_rate")),
         _num(qc_cfg.get("mapping_rate")),
         ">="),
        ("duplication_rate",
         _num(metrics.get("duplication_rate")),
         _num(qc_cfg.get("duplication_rate")),
         "<="),
        ("gc_content",
         _num(metrics.get("gc_percentage")),
         (_num(qc_cfg.get("GC_content_min")), _num(qc_cfg.get("GC_content_max")))
             if qc_cfg.get("GC_content_min") is not None
                and qc_cfg.get("GC_content_max") is not None
             else None,
         "between"),
    ]

    fail_count = 0
    with open(args.output, "w") as fh:
        fh.write("sample_id\t{}\n".format(args.sample))
        for name, actual, threshold, comparator in checks:
            status, actual_s, threshold_s, cmp_s = _check(
                name, actual, threshold, comparator
            )
            if status == "FAIL":
                fail_count += 1
            fh.write(f"{name}\t{actual_s}\t{threshold_s}\t{cmp_s}\t{status}\n")
        overall = "PASS" if fail_count == 0 else "FAIL"
        fh.write(f"overall\t{fail_count}\tfail_count\t<=0\t{overall}\n")

    print(f"[qc_filter] {args.sample}: overall={overall} fail_count={fail_count}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
