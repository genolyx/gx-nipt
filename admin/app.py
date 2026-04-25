"""
gx-nipt Admin — Home
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

import streamlit as st
from lib.paths import REF_DIR, CONFIG_DIR, OUTPUT_DIR, LABS, REF_TYPES, GROUPS
from lib.refs import ref_status_table

st.set_page_config(
    page_title="gx-nipt Admin",
    page_icon="🧬",
    layout="wide",
)

st.title("🧬 gx-nipt Admin")
st.caption("Developer interface for managing references, lab configs, and pipeline runs.")

# ── Summary cards ─────────────────────────────────────────────────────────────
col1, col2, col3, col4 = st.columns(4)
with col1:
    st.metric("Labs", len(LABS))
with col2:
    completed = list(OUTPUT_DIR.glob("*/*/*.completed")) if OUTPUT_DIR.exists() else []
    st.metric("Completed Runs", len(completed))
with col3:
    npz_count = len(list(REF_DIR.glob("**/*.npz"))) if REF_DIR.exists() else 0
    st.metric("NPZ References", npz_count)
with col4:
    ezd_count = len(list(REF_DIR.glob("**/EZD/**/*_thresholds_new.tsv"))) if REF_DIR.exists() else 0
    st.metric("EZD Threshold Files", ezd_count)

st.divider()

# ── Reference status matrix ───────────────────────────────────────────────────
st.subheader("Reference Status Matrix")

if not REF_DIR.exists():
    st.warning(f"`{REF_DIR}` does not exist. Check `params.ref_dir` in nextflow.config.")
elif not LABS:
    st.warning(f"No lab configs found in `{CONFIG_DIR}`.")
else:
    df = ref_status_table(REF_DIR, LABS, GROUPS)
    # Pivot: rows = Lab+Group, cols = Type
    pivot = df.pivot_table(index=["Lab", "Group"], columns="Type", values="Exists", aggfunc="first")
    pivot = pivot.reindex(columns=REF_TYPES)

    def style_bool(val):
        if val is True:
            return "background-color:#d4edda; color:#155724"
        elif val is False:
            return "background-color:#f8d7da; color:#721c24"
        return ""

    st.dataframe(pivot.style.applymap(style_bool).format(
        {t: lambda v: "✅" if v else "❌" for t in REF_TYPES}
    ), use_container_width=True)

st.divider()

# ── Quick navigation ──────────────────────────────────────────────────────────
st.subheader("Quick Navigation")
cols = st.columns(5)
pages = [
    ("🏥 Lab Config",         "Edit pipeline_config.json per lab"),
    ("🧬 Reference Status",   "View reference files and EZD thresholds"),
    ("🔨 Reference Builder",  "Generate EZD / PRIZM / WC / WCX references"),
    ("🔍 NPZ Diagnostics",    "Inspect and patch NPZ reference files"),
    ("🩺 Analysis Monitor",   "Browse completed pipeline results"),
]
for col, (title, desc) in zip(cols, pages):
    with col:
        st.markdown(f"**{title}**")
        st.caption(desc)
