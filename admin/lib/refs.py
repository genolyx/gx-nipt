"""Reference file inspection helpers."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from typing import Optional


def npz_info(path: Path) -> dict:
    """Return metadata about a WC/WCX .npz reference file."""
    info = {"path": str(path), "exists": path.exists(), "size_mb": None,
            "keys": [], "has_trained_cutoff": False, "has_is_nipt": False,
            "bins_per_chr_key": None, "shape_info": {}}
    if not path.exists():
        return info
    info["size_mb"] = round(path.stat().st_size / 1e6, 2)
    try:
        data = np.load(str(path), allow_pickle=True)
        info["keys"] = sorted(data.files)
        info["has_trained_cutoff"] = "trained_cutoff" in data.files
        info["has_is_nipt"] = "is_nipt" in data.files
        # detect bins_per_chr key variant
        for k in data.files:
            if "bins_per_chr" in k:
                info["bins_per_chr_key"] = k
                break
        for k in ["sample_ids", "reference", "bins_per_chr", "bins_per_chrA"]:
            if k in data.files:
                arr = data[k]
                info["shape_info"][k] = str(arr.shape)
    except Exception as e:
        info["error"] = str(e)
    return info


def ezd_thresholds(path: Path) -> Optional[pd.DataFrame]:
    """Load EZD thresholds TSV."""
    if not path.exists():
        return None
    try:
        return pd.read_csv(path, sep="\t")
    except Exception:
        return None


def prizm_csv_info(group_dir: Path) -> dict:
    """Return size and existence info for PRIZM CSV files."""
    expected = [
        "total_mean.csv", "total_sd.csv",
        "male_mean.csv",  "male_sd.csv",
        "female_mean.csv","female_sd.csv",
        "total_10mb_mean.csv", "total_10mb_sd.csv",
    ]
    result = {}
    for fname in expected:
        p = group_dir / fname
        result[fname] = {"exists": p.exists(),
                         "size_kb": round(p.stat().st_size / 1024, 1) if p.exists() else None}
    return result


def ref_status_table(ref_dir: Path, labs: list, groups: list) -> pd.DataFrame:
    """Build a summary DataFrame: lab × type × group → exists (bool)."""
    rows = []
    for lab in labs:
        for rtype in ["EZD", "PRIZM", "WC", "WCX"]:
            for group in groups:
                lab_path = ref_dir / "labs" / lab
                if rtype == "EZD":
                    exists = (lab_path / "EZD" / group / f"{group}_thresholds_new.tsv").exists()
                elif rtype == "PRIZM":
                    exists = (lab_path / "PRIZM" / group / "total_mean.csv").exists()
                elif rtype == "WC":
                    fname = {"orig": "orig_200k_proper_paired.npz",
                             "fetus": "fetus_200k_of.npz",
                             "mom": "mom_200k_of.npz"}[group]
                    exists = (lab_path / "WC" / fname).exists()
                elif rtype == "WCX":
                    fname = {"orig": "orig_200k_proper_paired.npz",
                             "fetus": "fetus_200k_of.npz",
                             "mom": "mom_200k_of.npz"}[group]
                    exists = (lab_path / "WCX" / fname).exists()
                rows.append({"Lab": lab, "Type": rtype, "Group": group, "Exists": exists})
    return pd.DataFrame(rows)


def sca_config(path: Path) -> Optional[dict]:
    """Load EZD sca_config.json."""
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text())
    except Exception:
        return None
