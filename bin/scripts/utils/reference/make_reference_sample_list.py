#!/usr/bin/env python3
import argparse
import glob
import os
import re
import json
from typing import Optional, Dict, List

# =========================
# QC parsing (TAB-based)
# =========================

WANT_QC_KEYS = {
    "number_of_reads",
    "number_of_mapped_reads",
    "mapping_rate",
    "duplication_rate",
    "mean_mapping_quality",
    "mean_coverageData",
    "GC_content",
}

# Sample folder name pattern will be set based on --prefix argument


def strip_bom(s: str) -> str:
    return s.replace("\ufeff", "").replace("\u200b", "")


def parse_qc_txt(qc_path: str) -> Dict[str, Optional[str]]:
    out = {k: None for k in WANT_QC_KEYS}

    if not os.path.exists(qc_path):
        return out

    with open(qc_path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = strip_bom(line).strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split("\t", 1)
            if len(parts) != 2:
                continue

            key = parts[0].strip().rstrip(":")
            val = parts[1].strip()

            if key in out:
                out[key] = val

    return out


def to_int(s: Optional[str]) -> Optional[int]:
    if s is None:
        return None
    try:
        return int(float(s.replace(",", "")))
    except Exception:
        return None


def to_float(s: Optional[str]) -> Optional[float]:
    if s is None:
        return None
    x = s.replace(",", "").replace("%", "")
    x = re.sub(r"[Xx]$", "", x)
    try:
        return float(x)
    except Exception:
        return None


# =========================
# FF / Gender parsing
# =========================

def parse_gender_gd2(path: str) -> str:
    if not os.path.exists(path):
        return ""
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = strip_bom(line).strip()
            if not line or line.lower().startswith("value"):
                continue
            parts = re.split(r"\s+", line)
            if len(parts) >= 3 and parts[0] == "gd_2":
                return parts[2]
    return ""


def parse_ff(path: str) -> Dict[str, str]:
    out = {"Fragment_FF": "", "YFF_2": "", "SeqFF": "", "M-SeqFF": ""}
    if not os.path.exists(path):
        return out
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = strip_bom(line).strip()
            if not line or line.lower() == "value":
                continue
            parts = re.split(r"\s+", line)
            if len(parts) >= 2 and parts[0] in out:
                out[parts[0]] = parts[1]
    return out


# =========================
# Report JSON parsing
# =========================

def parse_report_json(report_path: str) -> Dict[str, str]:
    """
    Extract Result, MDResult, Disease from <sample_id>_report.json
    """
    out = {"Result": "", "MDResult": "", "Disease": ""}
    if not os.path.exists(report_path):
        return out
    try:
        with open(report_path, "r", encoding="utf-8") as f:
            j = json.load(f)
        out["Result"] = str(j.get("Result", ""))
        out["MDResult"] = str(j.get("MDResult", ""))
        
        # Extract items from trisomy_result
        trisomy_result = j.get("trisomy_result", [])
        if trisomy_result and isinstance(trisomy_result, list):
            items = [item.get("item", "") for item in trisomy_result if isinstance(item, dict) and "item" in item]
            out["Disease"] = ", ".join(items) if items else ""
        
    except Exception:
        pass
    return out


def fmt(x: Optional[float], nd: int) -> str:
    return "" if x is None else f"{x:.{nd}f}"


# =========================
# Main
# =========================

def main():
    ap = argparse.ArgumentParser(
        description="Generate reference_sample_list.tsv (QC + FF + report Result/MDResult)."
    )
    ap.add_argument("--dirs", nargs="+", required=True,
                    help="analysis dirs, e.g. /home/ken/ken-nipt/analysis/2507")
    ap.add_argument("--out", default="reference_sample_list.tsv")
    ap.add_argument("--output-root", default="/home/ken/ken-nipt/output",
                    help="root directory containing report jsons")
    ap.add_argument("--prefix", default="GNMF",
                    help="sample ID prefix (default: GNMF)")
    args = ap.parse_args()
    
    # Create sample directory regex pattern based on prefix
    SAMPLE_DIR_RE = re.compile(rf"^{re.escape(args.prefix)}\d{{8}}$")

    header = [
        "month",
        "sample_id",
        "sample_dir",
        "number_of_reads",
        "number_of_mapped_reads",
        "mapping_rate(%)",
        "duplication_rate(%)",
        "mean_mapping_quality",
        "mean_coverageData(X)",
        "GC_content(%)",
        "fetal_gender(gd_2)",
        "Fragment_FF",
        "YFF_2",
        "SeqFF",
        "M-SeqFF",
        "Result",
        "Disease",
        "MDResult",
    ]

    rows: List[List[str]] = []

    for base in args.dirs:
        base = os.path.abspath(base)
        if not os.path.isdir(base):
            continue

        month = os.path.basename(base.rstrip("/"))

        sample_dirs = sorted(
            d for d in glob.glob(os.path.join(base, f"{args.prefix}*"))
            if os.path.isdir(d) and SAMPLE_DIR_RE.match(os.path.basename(d))
        )

        for sample_dir in sample_dirs:
            sample_id = os.path.basename(sample_dir)

            # QC
            qc_path = os.path.join(sample_dir, "Output_QC", f"{sample_id}.qc.txt")
            qc = parse_qc_txt(qc_path)

            reads = to_int(qc.get("number_of_reads"))
            mapped = to_int(qc.get("number_of_mapped_reads"))
            maprate = to_float(qc.get("mapping_rate"))
            duprate = to_float(qc.get("duplication_rate"))
            mapq = to_float(qc.get("mean_mapping_quality"))
            cov = to_float(qc.get("mean_coverageData"))
            gc = to_float(qc.get("GC_content"))

            # FF / gender
            ff_dir = os.path.join(sample_dir, "Output_FF")
            gender = parse_gender_gd2(os.path.join(ff_dir, f"{sample_id}.gender.txt"))
            ff = parse_ff(os.path.join(ff_dir, f"{sample_id}.fetal_fraction.txt"))

            # Report JSON
            report_path = os.path.join(
                args.output_root, month, sample_id, f"{sample_id}_report.json"
            )
            report = parse_report_json(report_path)

            rows.append([
                month,
                sample_id,
                sample_dir,
                "" if reads is None else str(reads),
                "" if mapped is None else str(mapped),
                fmt(maprate, 2),
                fmt(duprate, 2),
                fmt(mapq, 4),
                fmt(cov, 4),
                fmt(gc, 2),
                gender,
                ff["Fragment_FF"],
                ff["YFF_2"],
                ff["SeqFF"],
                ff["M-SeqFF"],
                report["Result"],
                report["Disease"],
                report["MDResult"],
            ])

    out = os.path.abspath(args.out)
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)

    with open(out, "w", encoding="utf-8") as w:
        w.write("\t".join(header) + "\n")
        for r in rows:
            w.write("\t".join(r) + "\n")

    print(f"[OK] wrote {out} (n={len(rows)})")


if __name__ == "__main__":
    main()
