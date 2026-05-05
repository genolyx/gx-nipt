"""
Reference Builder — GUI front-end for reference generation scripts.

Tabs:
  1. EZD / PRIZM / WCX  — legacy create_reference.py
  2. gx-cnv             — build_gxcnv_reference.py  (BAM→NPZ + GC inject + newref)
  3. gx-FF              — build_gxff_model.py        (sample select + train)
"""
import sys, subprocess, threading, queue, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parents[1]))

import streamlit as st
from lib.paths import (
    LABS, REF_DIR, GROUPS, REF_SCRIPTS_DIR,
    BUILD_GXCNV_SCRIPT, BUILD_GXFF_SCRIPT, SCRIPTS_DIR,
)
from lib.runner import python_cmd

st.set_page_config(page_title="Reference Builder | gx-nipt Admin", page_icon="🔨", layout="wide")
st.title("🔨 Reference Builder")

tab_legacy, tab_gxcnv, tab_gxff = st.tabs(["📦 EZD / PRIZM / WCX", "🧬 gx-cnv", "🔬 gx-FF"])


# ─────────────────────────────────────────────────────────────────────────────
# Shared helper: streaming log runner
# ─────────────────────────────────────────────────────────────────────────────

def run_with_log(cmd: list, log_area):
    log_text = ""
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    for line in proc.stdout:
        log_text += line
        log_area.code(log_text[-8000:], language="bash")
    proc.wait()
    return proc.returncode, log_text


# ═════════════════════════════════════════════════════════════════════════════
# TAB 1 — Legacy (EZD / PRIZM / WCX)
# ═════════════════════════════════════════════════════════════════════════════

with tab_legacy:
    st.caption("Generate EZD / PRIZM / WC / WCX reference files via `create_reference.py`.")

    CREATE_REF = REF_SCRIPTS_DIR / "create_reference.py"
    MAKE_LIST  = REF_SCRIPTS_DIR / "make_reference_sample_list.py"

    if not CREATE_REF.exists():
        st.error(f"Script not found: `{CREATE_REF}`")
    else:
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
        st.subheader("Step 2 — Generate Reference")

        col_left, col_right = st.columns(2)
        with col_left:
            lab         = st.selectbox("Lab", sorted(LABS), key="leg_lab")
            sample_list = st.text_input("Sample list path",
                                        value=str(REF_DIR / "sample_list.txt"), key="leg_sl")
            output_dir  = st.text_input("Output directory",
                                        value=str(REF_DIR / "labs" / lab), key="leg_od")
            ref_source  = st.text_input("Reference source (analysis root)",
                                        value="/home/ken/gx-nipt/analysis", key="leg_rs")
        with col_right:
            ref_types = st.multiselect("Algorithms", ["ezd", "prizm", "wc", "wcx"],
                                       default=["ezd", "prizm"], key="leg_rt")
            groups    = st.multiselect("Groups", GROUPS, default=["orig"], key="leg_grp")
            preview   = st.checkbox("Preview only (dry-run)", value=False, key="leg_prev")

        st.markdown("**Quality Filters**")
        qc1, qc2, qc3, qc4 = st.columns(4)
        min_seqff = qc1.number_input("Min seqFF",  value=2.0,  step=0.5, key="leg_q1")
        max_seqff = qc2.number_input("Max seqFF",  value=20.0, step=1.0, key="leg_q2")
        min_map   = qc3.number_input("Min mapping rate (%)", value=85.0, step=1.0, key="leg_q3")
        max_dup   = qc4.number_input("Max dup rate (%)",     value=40.0, step=1.0, key="leg_q4")

        cmd_parts = [
            python_cmd(), str(CREATE_REF),
            "--sample-list", sample_list, "--labcode", lab,
            "--ref-type", *ref_types, "--groups", *groups,
            "--output-dir", output_dir, "--reference-source", ref_source,
            "--min-seqff", str(min_seqff), "--max-seqff", str(max_seqff),
            "--min-mapping-rate", str(min_map), "--max-duplication-rate", str(max_dup),
        ]
        if preview:
            cmd_parts.append("--preview-only")

        with st.expander("📋 Command preview"):
            st.code(" \\\n  ".join(cmd_parts), language="bash")

        if st.button("🚀 Run Reference Generation", type="primary", key="leg_run"):
            if not Path(sample_list).exists():
                st.warning(f"Sample list not found: `{sample_list}`")
            elif not ref_types or not groups:
                st.warning("Select algorithm and group.")
            else:
                log_area = st.empty()
                with st.spinner("Running…"):
                    rc, log_text = run_with_log(cmd_parts, log_area)
                st.success("✅ Done!") if rc == 0 else st.error(f"❌ Exit code {rc}")
                st.code(log_text[-4000:], language="bash")


# ═════════════════════════════════════════════════════════════════════════════
# TAB 2 — gx-cnv Reference Builder
# ═════════════════════════════════════════════════════════════════════════════

with tab_gxcnv:
    st.caption(
        "Build a gx-cnv reference panel (BAM→NPZ + HMMcopy GC injection + `gxcnv newref`).\n\n"
        "Requires `gxcnv` installed and sample data TSV with `sample_dir`, `fetal_gender(gd_2)`, QC columns."
    )

    if not BUILD_GXCNV_SCRIPT.exists():
        st.error(f"Script not found: `{BUILD_GXCNV_SCRIPT}`")
    else:
        col1, col2 = st.columns(2)
        with col1:
            cnv_lab       = st.selectbox("Lab", sorted(LABS), key="cnv_lab")
            cnv_sex       = st.selectbox("Fetal sex", ["female", "male"], key="cnv_sex")
            cnv_tsv       = st.text_input(
                "Sample data TSV",
                value="/home/ken/ken-nipt/bin/scripts/utils/reference/reference_sample_list_from_json_v3.tsv",
                key="cnv_tsv",
            )
            cnv_out       = st.text_input(
                "Output directory",
                value=str(REF_DIR / "labs" / cnv_lab / "GXCNV" / cnv_sex),
                key="cnv_out",
            )
            cnv_npz_cache = st.text_input(
                "NPZ cache directory",
                value=str(REF_DIR / "gxcnv" / f"npz_{cnv_sex}"),
                key="cnv_npz",
            )
        with col2:
            cnv_n       = st.number_input("Target # samples", value=100, min_value=20, step=10, key="cnv_n")
            cnv_threads = st.number_input("Parallel threads (BAM→NPZ)", value=8, min_value=1, key="cnv_thr")
            cnv_pca     = st.slider("PCA variance threshold", 0.80, 0.99, 0.95, step=0.01, key="cnv_pca")
            cnv_maxdup  = st.number_input("Max dup rate (%)", value=15.0, step=1.0, key="cnv_dup")
            cnv_minmap  = st.number_input("Min mapping rate (%)", value=85.0, step=1.0, key="cnv_map")

        cnv_cmd = [
            python_cmd(), str(BUILD_GXCNV_SCRIPT),
            "--sample-tsv",  cnv_tsv,
            "--labcode",     cnv_lab,
            "--sex",         cnv_sex,
            "--out-dir",     cnv_out,
            "--npz-cache",   cnv_npz_cache,
            "--n-samples",   str(cnv_n),
            "--threads",     str(cnv_threads),
            "--pca-variance",str(cnv_pca),
            "--max-dup",     str(cnv_maxdup),
            "--min-map",     str(cnv_minmap),
        ]

        with st.expander("📋 Command preview"):
            st.code(" \\\n  ".join(cnv_cmd), language="bash")

        if st.button("🚀 Build gx-cnv Reference", type="primary", key="cnv_run"):
            if not Path(cnv_tsv).exists():
                st.warning(f"TSV not found: `{cnv_tsv}`")
            else:
                log_area = st.empty()
                st.info(f"Building **{cnv_sex}** reference for **{cnv_lab}** — this may take 30–60 min.")
                with st.spinner("Running pipeline…"):
                    rc, log_text = run_with_log(cnv_cmd, log_area)
                if rc == 0:
                    ref_file = Path(cnv_out) / "reference.npz"
                    st.success(f"✅ Reference panel saved: `{ref_file}`")
                else:
                    st.error(f"❌ Process exited with code {rc}")
                st.code(log_text[-4000:], language="bash")


# ═════════════════════════════════════════════════════════════════════════════
# TAB 3 — gx-FF Model Builder
# ═════════════════════════════════════════════════════════════════════════════

with tab_gxff:
    st.caption(
        "Train a gx-FF ensemble model (LightGBM + DNN) from male-fetus NIPT samples.\n\n"
        "Uses HMMcopy wig (coverage features) + BAM (fragment features).\n"
        "Requires the patched `pipeline.py` with `BAM_PATH` column support."
    )

    if not BUILD_GXFF_SCRIPT.exists():
        st.error(f"Script not found: `{BUILD_GXFF_SCRIPT}`")
    else:
        col1, col2 = st.columns(2)
        with col1:
            ff_tsv     = st.text_input(
                "Sample data TSV",
                value="/home/ken/ken-nipt/bin/scripts/utils/reference/reference_sample_list_from_json_v3.tsv",
                key="ff_tsv",
            )
            ff_out     = st.text_input(
                "Output directory (model)",
                value=str(REF_DIR / "gxff" / "model"),
                key="ff_out",
            )
            ff_n       = st.number_input("Target # training samples", value=170, min_value=30, step=10, key="ff_n")
            ff_min_ff  = st.number_input("Min YFF_2 (%) for male samples", value=4.0, step=0.5, key="ff_minff")
        with col2:
            ff_features = st.multiselect("Features", ["coverage", "fragment", "nucleosome"],
                                         default=["coverage", "fragment"], key="ff_feat")
            ff_threads  = st.number_input("Parallel threads", value=8, min_value=1, key="ff_thr")
            ff_cv       = st.number_input("CV folds", value=5, min_value=2, key="ff_cv")
            ff_augment  = st.checkbox("Augment low-FF samples", value=True, key="ff_aug")
            ff_lfw      = st.number_input("Low-FF sample weight", value=3.0, step=0.5, key="ff_lfw")

        ff_cmd = [
            python_cmd(), str(BUILD_GXFF_SCRIPT),
            "--sample-tsv",    ff_tsv,
            "--out-dir",       ff_out,
            "--n-samples",     str(ff_n),
            "--min-ff",        str(ff_min_ff),
            "--features",      *ff_features,
            "--threads",       str(ff_threads),
            "--cv-folds",      str(ff_cv),
            "--low-ff-weight", str(ff_lfw),
        ]
        if not ff_augment:
            ff_cmd.append("--no-augment")

        with st.expander("📋 Command preview"):
            st.code(" \\\n  ".join(ff_cmd), language="bash")

        st.warning(
            "⚠️ Training with coverage + fragment features takes **~3–4 hours** for 170 samples. "
            "Consider running in a terminal for long jobs.",
            icon="⏳",
        )

        col_btn1, col_btn2 = st.columns(2)
        config_only = col_btn1.button("📋 Generate Config Only", key="ff_cfg_only")
        run_train   = col_btn2.button("🚀 Train Model", type="primary", key="ff_run")

        if config_only:
            cmd = ff_cmd + ["--config-only"]
            with st.spinner("Generating config…"):
                result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode == 0:
                cfg_path = Path(ff_out) / "training_config.tsv"
                st.success(f"Training config saved: `{cfg_path}`")
                if cfg_path.exists():
                    import pandas as _pd
                    df_cfg = _pd.read_csv(cfg_path, sep="\t")
                    st.dataframe(df_cfg.head(10), use_container_width=True)
            else:
                st.error("Failed")
                st.code(result.stdout + result.stderr)

        if run_train:
            if not Path(ff_tsv).exists():
                st.warning(f"TSV not found: `{ff_tsv}`")
            elif not ff_features:
                st.warning("Select at least one feature type.")
            else:
                log_area = st.empty()
                st.info("Training started — output streams below.")
                with st.spinner("Training… (this will take a while)"):
                    rc, log_text = run_with_log(ff_cmd, log_area)
                if rc == 0:
                    metrics_path = Path(ff_out) / "cv_metrics.tsv"
                    st.success(f"✅ Model saved: `{Path(ff_out) / 'gxff_model.pkl'}`")
                    if metrics_path.exists():
                        import pandas as _pd
                        m = _pd.read_csv(metrics_path, sep="\t")
                        st.dataframe(m, use_container_width=True)
                else:
                    st.error(f"❌ Process exited with code {rc}")
                st.code(log_text[-4000:], language="bash")
