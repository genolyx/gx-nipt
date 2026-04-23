#!/usr/bin/env python3
"""
parse_qualimap.py
-----------------

Parse Qualimap's ``genome_results.txt`` into a flat, key=value QC summary
that downstream steps (qc_filter.py, generate_json_output.py) can consume
without re-implementing the Qualimap parser logic.

CLI
~~~
    --input   <genome_results.txt>
    --sample  <sample_id>
    --output  <sample>.qc.txt

Output format (one key per line, tab-separated):
    sample_id               <sample>
    number_of_reads         <int>
    number_of_mapped_reads  <int>
    mapping_rate            <float, percent>
    duplication_rate        <float, percent>
    mean_coverage           <float>
    gc_percentage           <float>

Missing fields fall back to ``NA`` rather than aborting, so the QC filter
step can decide how to treat partial QC reports.
"""

import argparse
import re
import sys
from pathlib import Path


_PATTERNS = {
    "number_of_reads":        re.compile(r"number of reads\s*=\s*([0-9,\.]+)"),
    "number_of_mapped_reads": re.compile(r"number of mapped reads\s*=\s*([0-9,\.]+)"),
    "mapping_rate":           re.compile(r"number of mapped reads\s*=\s*[0-9,\.]+\s*\(([0-9\.]+)%\)"),
    "duplication_rate":       re.compile(r"duplication rate\s*=\s*([0-9\.]+)%?"),
    "mean_coverage":          re.compile(r"mean coverageData\s*=\s*([0-9\.Xx]+)"),
    "gc_percentage":          re.compile(r"GC percentage\s*=\s*([0-9\.]+)%?"),
}


def _clean_int(raw: str) -> str:
    return raw.replace(",", "").split(".")[0]


def _clean_float(raw: str) -> str:
    # Qualimap reports coverage as "12.34X" — strip trailing markers.
    return re.sub(r"[^0-9\.]", "", raw)


def parse(input_path: Path) -> dict:
    text = input_path.read_text()
    parsed = {}
    for key, pattern in _PATTERNS.items():
        match = pattern.search(text)
        if not match:
            parsed[key] = "NA"
            continue
        value = match.group(1)
        if key in ("number_of_reads", "number_of_mapped_reads"):
            parsed[key] = _clean_int(value)
        else:
            parsed[key] = _clean_float(value)
    return parsed


def main() -> int:
    ap = argparse.ArgumentParser(description="Parse Qualimap genome_results.txt")
    ap.add_argument("--input",  required=True, type=Path)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--output", required=True, type=Path)
    args = ap.parse_args()

    if not args.input.exists():
        print(f"[parse_qualimap] ERROR: input not found: {args.input}",
              file=sys.stderr)
        return 1

    metrics = parse(args.input)
    metrics["sample_id"] = args.sample

    ordered_keys = [
        "sample_id",
        "number_of_reads",
        "number_of_mapped_reads",
        "mapping_rate",
        "duplication_rate",
        "mean_coverage",
        "gc_percentage",
    ]
    with open(args.output, "w") as fh:
        for key in ordered_keys:
            fh.write(f"{key}\t{metrics.get(key, 'NA')}\n")

    print(f"[parse_qualimap] Wrote {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
