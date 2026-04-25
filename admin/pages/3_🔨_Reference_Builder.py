"""
Reference Builder — GUI front-end for create_reference.py
"""
import sys, subprocess, threading, queue, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parents[1]))

import streamlit as st
from lib.paths import LABS, REF_DIR, GROUPS, REF_SCRIPTS_DIR
from lib.runner import python_cmd

st.set_page_config(page_title="Reference Builder | gx-nipt Admin", page_icon="🔨", layout="wide")
st.title("🔨 Reference Builder")
st.caption("Generate EZD / PRIZM / WC / WCX reference files via `create_reference.py`.")

CREATE_REF = REF_SCRIPTS_DIR / "create_reference.py"
MAKE_LIST  = REF_SCRIPTS_DIR / "make_reference_sample_list.py"

if not CREATE_REF.exists():
    st.error(f"Script not found: `{CREATE_REF}`")
    st.stop()

# ── Step 1: Sample list ───────────────────────────────────────────────────────
st.subheader("Step 1 — Build Sample List")
with st.expander("Generate sample list from analysis directories", expanded=True):
    analysis_dirs = st.text_area(
        "Analysis directories (one per line)",
        placeholder="/home/ken/gx-nipt/analysis/2603\n/home/ken/gx-nipt/analysis/2602",
        height=100,
    )
    list_out = st.text_input("Output sample list path",
                             value=str(REF_DIR / "sample_list.txt"))
    prefix   = st.text_input("Prefix filter (optional)", value="")

    if st.button("📋 Generate Sample List"):
        dirs = [d.strip() for d in analysis_dirs.strip().splitlines() if d.strip()]
        if not dirs:
            st.warning("Enter at least one analysis directory.")
        elif not MAKE_LIST.exists():
            st.error(f"Script not found: `{MAKE_LIST}`")
        else:
            cmd = [python_cmd(), str(MAKE_LIST), "--out", list_out, "--dirs"] + dirs
            if prefix:
                cmd += ["--prefix", prefix]
            with st.spinner("Generating sample list…"):
                result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode == 0:
                st.success(f"Sample list written to `{list_out}`")
                st.code(result.stdout or "(no output)")
            else:
                st.error("Failed")
                st.code(result.stdout + result.stderr)

st.divider()

# ── Step 2: Reference generation ──────────────────────────────────────────────
st.subheader("Step 2 — Generate Reference")

col_left, col_right = st.columns(2)

with col_left:
    lab         = st.selectbox("Lab", sorted(LABS))
    sample_list = st.text_input("Sample list path",
                                value=str(REF_DIR / "sample_list.txt"))
    output_dir  = st.text_input("Output directory",
                                value=str(REF_DIR / "labs" / lab))
    ref_source  = st.text_input("Reference source (analysis root)",
                                value="/home/ken/gx-nipt/analysis")

with col_right:
    ref_types = st.multiselect("Algorithms", ["ezd", "prizm", "wc", "wcx"],
                               default=["ezd", "prizm"])
    groups    = st.multiselect("Groups", GROUPS, default=["orig"])
    preview   = st.checkbox("Preview only (dry-run)", value=False)

st.markdown("**Quality Filters**")
qcol1, qcol2, qcol3, qcol4 = st.columns(4)
min_seqff  = qcol1.number_input("Min seqFF",  value=2.0, step=0.5)
max_seqff  = qcol2.number_input("Max seqFF",  value=20.0, step=1.0)
min_map    = qcol3.number_input("Min mapping rate (%)", value=85.0, step=1.0)
max_dup    = qcol4.number_input("Max dup rate (%)",     value=40.0, step=1.0)

debug = st.checkbox("Debug output")

# Build command preview
cmd_parts = [
    python_cmd(), str(CREATE_REF),
    "--sample-list",        sample_list,
    "--labcode",            lab,
    "--ref-type",           *ref_types,
    "--groups",             *groups,
    "--output-dir",         output_dir,
    "--reference-source",   ref_source,
    "--min-seqff",          str(min_seqff),
    "--max-seqff",          str(max_seqff),
    "--min-mapping-rate",   str(min_map),
    "--max-duplication-rate", str(max_dup),
]
if preview:
    cmd_parts.append("--preview-only")
if debug:
    cmd_parts.append("--debug")

with st.expander("📋 Command preview"):
    st.code(" \\\n  ".join(cmd_parts), language="bash")

if st.button("🚀 Run Reference Generation", type="primary"):
    if not Path(sample_list).exists():
        st.warning(f"Sample list not found: `{sample_list}`")
    elif not ref_types:
        st.warning("Select at least one algorithm.")
    elif not groups:
        st.warning("Select at least one group.")
    else:
        log_area = st.empty()
        log_text = ""

        with st.spinner("Running…"):
            proc = subprocess.Popen(
                cmd_parts,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
            )
            for line in proc.stdout:
                log_text += line
                log_area.code(log_text[-6000:], language="bash")  # show last 6k chars
            proc.wait()

        if proc.returncode == 0:
            st.success("✅ Reference generation complete!")
        else:
            st.error(f"❌ Process exited with code {proc.returncode}")

        st.code(log_text[-4000:], language="bash")
