"""
NPZ Diagnostics — inspect WC/WCX reference files and patch missing keys.
"""
import sys, subprocess
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parents[1]))

import streamlit as st
import numpy as np
import pandas as pd

from lib.paths import REPO_ROOT, REF_DIR, PATCH_NPZ_SCRIPT
from lib.refs import npz_info
from lib.runner import python_cmd

st.set_page_config(page_title="NPZ Diagnostics | gx-nipt Admin", page_icon="🔍", layout="wide")
st.title("🔍 NPZ Diagnostics")
st.caption("Inspect and patch WC / WCX `.npz` reference files.")

# ─── Scan directory selector ───────────────────────────────────────────────────
st.subheader("Scan Directory")
scan_dir_str = st.text_input(
    "Directory to scan for .npz files",
    value=str(REPO_ROOT),
    help="Recursively searches for all .npz files under this path.",
)
scan_dir = Path(scan_dir_str)

if not scan_dir.exists():
    st.error(f"Directory not found: `{scan_dir}`")
    st.stop()

all_npz = sorted(p for p in scan_dir.rglob("*.npz") if ".tmp" not in p.name)

if not all_npz:
    st.warning(f"No .npz files found under `{scan_dir}`")
    st.stop()

# ─── Summary table ────────────────────────────────────────────────────────────
st.subheader("All NPZ Files")

rows = []
for p in all_npz:
    info = npz_info(p)
    try:
        rel = str(p.relative_to(scan_dir))
    except ValueError:
        rel = str(p)
    rows.append({
        "Path":             rel,
        "Size (MB)":        info["size_mb"],
        "trained_cutoff":   "✅" if info["has_trained_cutoff"] else "❌",
        "is_nipt":          "✅" if info["has_is_nipt"] else "❌",
        "bins_per_chr key": info["bins_per_chr_key"] or "—",
        "Keys (n)":         len(info["keys"]),
        "Error":            info.get("error", ""),
    })

df_all = pd.DataFrame(rows)

needs_patch = df_all[
    (df_all["trained_cutoff"] == "❌") | (df_all["is_nipt"] == "❌")
]

c1, c2, c3 = st.columns(3)
c1.metric("Total NPZ files", len(df_all))
c2.metric("Needs patch", len(needs_patch))
c3.metric("OK", len(df_all) - len(needs_patch))

st.dataframe(df_all, use_container_width=True, hide_index=True)

st.divider()

# ─── Detail inspector ─────────────────────────────────────────────────────────
st.subheader("Inspect Single File")

rel_paths = [r["Path"] for r in rows]
npz_choice = st.selectbox("Select NPZ", rel_paths)
chosen_path = scan_dir / npz_choice

if chosen_path.exists():
    info = npz_info(chosen_path)
    col1, col2, col3, col4 = st.columns(4)
    col1.metric("Size (MB)", info["size_mb"])
    col2.metric("trained_cutoff", "✅" if info["has_trained_cutoff"] else "❌")
    col3.metric("is_nipt",       "✅" if info["has_is_nipt"] else "❌")
    col4.metric("bins_per_chr key", info["bins_per_chr_key"] or "—")

    st.write("**All keys:**")
    st.code(", ".join(info["keys"]))

    if info.get("shape_info"):
        st.write("**Array shapes:**")
        st.json(info["shape_info"])

    if info.get("error"):
        st.error(info["error"])

    # Show sample_ids if available
    try:
        data = np.load(str(chosen_path), allow_pickle=True)
        if "sample_ids" in data.files:
            ids = data["sample_ids"]
            st.write(f"**Samples ({len(ids)}):**")
            st.dataframe(pd.DataFrame({"sample_id": ids}), height=200, use_container_width=True)
    except Exception as e:
        st.error(str(e))

st.divider()

# ─── Batch Patch ──────────────────────────────────────────────────────────────
st.subheader("Patch NPZ Files")
st.caption(
    "Adds `trained_cutoff = 0.003` and `is_nipt = True` to files missing these keys. "
    f"Script: `{PATCH_NPZ_SCRIPT}`"
)

if not PATCH_NPZ_SCRIPT.exists():
    st.warning(f"Patch script not found: `{PATCH_NPZ_SCRIPT}`")
else:
    patch_targets = st.multiselect(
        "Files to patch (pre-selected: those missing keys)",
        options=rel_paths,
        default=[r["Path"] for _, r in needs_patch.iterrows()],
    )
    tc_val  = st.number_input("trained_cutoff value", value=0.003, format="%.4f", step=0.001)
    dry_run = st.checkbox("Dry-run (show command only)")

    if st.button("🩹 Patch Selected Files", type="primary"):
        if not patch_targets:
            st.warning("No files selected.")
        else:
            log = ""
            log_area = st.empty()
            for rel in patch_targets:
                full = str(scan_dir / rel)
                cmd = [python_cmd(), str(PATCH_NPZ_SCRIPT), full,
                       "--trained-cutoff", str(tc_val)]
                if dry_run:
                    st.code(" ".join(cmd), language="bash")
                    continue
                result = subprocess.run(cmd, capture_output=True, text=True)
                status = "✅" if result.returncode == 0 else "❌"
                log += f"{status} {rel}\n{result.stdout}{result.stderr}\n"
                log_area.code(log, language="bash")

            if not dry_run:
                st.success("Patching complete. Refresh the table above to confirm.")
                st.rerun()
