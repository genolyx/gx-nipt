"""
Reference Status Dashboard — per-lab file existence, EZD thresholds, PRIZM info.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parents[1]))

import streamlit as st
import pandas as pd
import plotly.graph_objects as go

from lib.paths import REF_DIR, LABS, GROUPS, lab_ref_dir
from lib.refs import ref_status_table, ezd_thresholds, prizm_csv_info, npz_info

st.set_page_config(page_title="Reference Status | gx-nipt Admin", page_icon="🧬", layout="wide")
st.title("🧬 Reference Status")

if not LABS:
    st.error("No labs found.")
    st.stop()

lab = st.selectbox("Lab", sorted(LABS))
ref_base = lab_ref_dir(lab)

# ─── Top-level existence matrix ───────────────────────────────────────────────
st.subheader("File Existence Matrix")
df = ref_status_table(REF_DIR, [lab], GROUPS)
pivot = df.pivot_table(index="Group", columns="Type", values="Exists", aggfunc="first")
pivot = pivot.reindex(columns=["EZD", "PRIZM", "WC", "WCX"])

st.dataframe(
    pivot.style.format(
        {t: (lambda v: "✅" if v else "❌") for t in ["EZD","PRIZM","WC","WCX"]}
    ),
    use_container_width=True,
)

st.divider()

# ─── Detail tabs ──────────────────────────────────────────────────────────────
tab_ezd, tab_prizm, tab_wc, tab_wcx = st.tabs(["EZD", "PRIZM", "WC", "WCX"])

# ── EZD ──────────────────────────────────────────────────────────────────────
with tab_ezd:
    st.subheader("EZD Thresholds")
    group = st.selectbox("Group", GROUPS, key="ezd_group")
    tsv_path = ref_base / "EZD" / group / f"{group}_thresholds_new.tsv"
    df_ezd = ezd_thresholds(tsv_path)
    if df_ezd is None:
        st.warning(f"Not found: `{tsv_path}`")
    else:
        st.caption(f"`{tsv_path}` — {len(df_ezd)} rows")
        st.dataframe(df_ezd, use_container_width=True)

        # Chart: Z_min / Z_max per chromosome
        if {"chr", "Z_min", "Z_max"}.issubset(df_ezd.columns):
            fig = go.Figure()
            fig.add_trace(go.Bar(name="Z_min", x=df_ezd["chr"], y=df_ezd["Z_min"],
                                 marker_color="#5b9bd5"))
            fig.add_trace(go.Bar(name="Z_max", x=df_ezd["chr"], y=df_ezd["Z_max"],
                                 marker_color="#ed7d31"))
            fig.update_layout(
                barmode="group", title="Z_min / Z_max per Chromosome",
                xaxis_title="Chromosome", yaxis_title="Z-score",
                height=380,
                paper_bgcolor="rgba(0,0,0,0)",
                plot_bgcolor="rgba(0,0,0,0)",
                font_color="#fafafa",
            )
            fig.update_xaxes(gridcolor="rgba(255,255,255,0.15)", zerolinecolor="rgba(255,255,255,0.3)")
            fig.update_yaxes(gridcolor="rgba(255,255,255,0.15)", zerolinecolor="rgba(255,255,255,0.3)")
            st.plotly_chart(fig, use_container_width=True)

        # sca_config
        sca_path = ref_base / "EZD" / group / "sca_config.json"
        if sca_path.exists():
            import json
            with st.expander("sca_config.json"):
                st.json(json.loads(sca_path.read_text()))

# ── PRIZM ─────────────────────────────────────────────────────────────────────
with tab_prizm:
    st.subheader("PRIZM Reference CSVs")
    group = st.selectbox("Group", GROUPS, key="prizm_group")
    prizm_dir = ref_base / "PRIZM" / group
    info = prizm_csv_info(prizm_dir)
    rows = [{"File": k, "Exists": "✅" if v["exists"] else "❌",
             "Size (KB)": v["size_kb"] or "—"} for k, v in info.items()]
    st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)

    # Preview one CSV
    csv_files = sorted(prizm_dir.glob("*.csv")) if prizm_dir.exists() else []
    if csv_files:
        chosen = st.selectbox("Preview CSV", [f.name for f in csv_files], key="prizm_csv")
        try:
            df_p = pd.read_csv(prizm_dir / chosen, header=None)
            st.caption(f"Shape: {df_p.shape}")
            st.dataframe(df_p.head(10), use_container_width=True)
        except Exception as e:
            st.error(str(e))

# ── WC ────────────────────────────────────────────────────────────────────────
with tab_wc:
    st.subheader("WC (Wisecondor) NPZ References")
    wc_dir = ref_base / "WC"
    npz_names = {"orig": "orig_200k_proper_paired.npz",
                 "fetus": "fetus_200k_of.npz",
                 "mom": "mom_200k_of.npz"}
    for g, fname in npz_names.items():
        path = wc_dir / fname
        info = npz_info(path)
        with st.expander(f"{g} — `{fname}` {'✅' if info['exists'] else '❌'}"):
            if not info["exists"]:
                st.warning("File not found.")
            else:
                c1, c2, c3 = st.columns(3)
                c1.metric("Size (MB)", info["size_mb"])
                c2.metric("trained_cutoff", "✅" if info["has_trained_cutoff"] else "❌")
                c3.metric("is_nipt",       "✅" if info["has_is_nipt"] else "❌")
                st.write("**Keys:**", ", ".join(info["keys"]))
                if info.get("error"):
                    st.error(info["error"])

# ── WCX ───────────────────────────────────────────────────────────────────────
with tab_wcx:
    st.subheader("WCX (WisecondorX) NPZ References")
    wcx_dir = ref_base / "WCX"
    wcx_files = {
        "orig M":    "orig_M_200k_proper_paired.npz",
        "orig F":    "orig_F_200k_proper_paired.npz",
        "fetus M":   "fetus_M_200k_of.npz",
        "fetus F":   "fetus_F_200k_of.npz",
        "mom":       "mom_200k_of.npz",
    }
    for label, fname in wcx_files.items():
        path = wcx_dir / fname
        info = npz_info(path)
        with st.expander(f"{label} — `{fname}` {'✅' if info['exists'] else '❌'}"):
            if not info["exists"]:
                st.warning("File not found.")
            else:
                c1, c2, c3 = st.columns(3)
                c1.metric("Size (MB)", info["size_mb"])
                c2.metric("trained_cutoff", "✅" if info["has_trained_cutoff"] else "❌")
                c3.metric("is_nipt",       "✅" if info["has_is_nipt"] else "❌")
                st.write("**Keys:**", ", ".join(info["keys"]))
                if info.get("shape_info"):
                    st.write("**Shapes:**", info["shape_info"])
                if info.get("error"):
                    st.error(info["error"])
