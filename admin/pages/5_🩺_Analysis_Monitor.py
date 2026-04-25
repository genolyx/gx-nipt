"""
Analysis Monitor — browse completed pipeline results.
"""
import sys, json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parents[1]))

import streamlit as st
import pandas as pd
import plotly.express as px

from lib.paths import OUTPUT_DIR

st.set_page_config(page_title="Analysis Monitor | gx-nipt Admin", page_icon="🩺", layout="wide")
st.title("🩺 Analysis Monitor")
st.caption("Browse completed gx-nipt pipeline results.")

if not OUTPUT_DIR.exists():
    st.warning(f"Output directory not found: `{OUTPUT_DIR}`")
    st.stop()

# ─── Collect all completed runs ───────────────────────────────────────────────
completed_markers = sorted(OUTPUT_DIR.glob("*/*/*.completed"))

if not completed_markers:
    st.info("No completed runs found yet.")
    st.stop()

rows = []
for marker in completed_markers:
    sample_id  = marker.stem          # GNCI26030001
    work_dir   = marker.parent.name   # GNCI26030001
    batch_dir  = marker.parent.parent.name  # 2603
    result_json = marker.parent / "Output_Result" / f"{sample_id}.result.json"
    completed_at = marker.read_text().strip()

    row = {
        "Batch": batch_dir,
        "Sample": sample_id,
        "Completed": completed_at,
        "Result JSON": result_json.exists(),
        "HTML Report": (marker.parent / "Output_Result" / f"{sample_id}.review.html").exists(),
    }

    # Parse result.json for key fields
    if result_json.exists():
        try:
            data  = json.loads(result_json.read_text())
            nipt  = data.get("NIPT", {})
            fr    = nipt.get("final_results", {})
            row.update({
                "Gender":          fr.get("fetal_gender", "—"),
                "FF (YFF)":        fr.get("fetal_fraction_yff", "—"),
                "FF (seqFF)":      fr.get("fetal_fraction_seqff", "—"),
                "Sample Bias":     fr.get("sample_bias", "—"),
                "Trisomy Result":  str(fr.get("final_trisomy_result") or "—"),
                "MD Result":       str(fr.get("final_md_result") or "—"),
            })
        except Exception:
            pass

    rows.append(row)

df = pd.DataFrame(rows).sort_values(["Batch", "Sample"], ascending=[False, True])

# ─── Filters ──────────────────────────────────────────────────────────────────
col_f1, col_f2 = st.columns(2)
batches = ["All"] + sorted(df["Batch"].unique(), reverse=True)
sel_batch  = col_f1.selectbox("Batch", batches)
search     = col_f2.text_input("Search sample ID")

filtered = df.copy()
if sel_batch != "All":
    filtered = filtered[filtered["Batch"] == sel_batch]
if search:
    filtered = filtered[filtered["Sample"].str.contains(search, case=False)]

# ─── Summary metrics ──────────────────────────────────────────────────────────
mc1, mc2, mc3, mc4 = st.columns(4)
mc1.metric("Total Runs", len(filtered))
mc2.metric("With JSON",  int(filtered["Result JSON"].sum()))
mc3.metric("With HTML",  int(filtered["HTML Report"].sum()))
if "Trisomy Result" in filtered.columns:
    abnormal = filtered[filtered["Trisomy Result"].str.contains("Trisomy|Detected", na=False)]
    mc4.metric("Trisomy Flags", len(abnormal))

# ─── Main table ───────────────────────────────────────────────────────────────
display_cols = ["Batch","Sample","Completed","Gender","FF (seqFF)","Sample Bias",
                "Trisomy Result","MD Result","Result JSON","HTML Report"]
display_cols = [c for c in display_cols if c in filtered.columns]

st.dataframe(
    filtered[display_cols],
    use_container_width=True,
    hide_index=True,
)

st.divider()

# ─── Detail viewer ────────────────────────────────────────────────────────────
st.subheader("Sample Detail")
if filtered.empty:
    st.info("No samples to show.")
else:
    sample_choice = st.selectbox(
        "Select sample",
        filtered["Sample"].tolist(),
        format_func=lambda s: f"{filtered[filtered['Sample']==s]['Batch'].values[0]} / {s}",
    )
    sel_row  = filtered[filtered["Sample"] == sample_choice].iloc[0]
    batch    = sel_row["Batch"]
    out_path = OUTPUT_DIR / batch / sample_choice

    tab_json, tab_qc, tab_files = st.tabs(["📊 Result JSON", "📈 QC", "📁 Files"])

    with tab_json:
        result_json = out_path / "Output_Result" / f"{sample_choice}.result.json"
        if result_json.exists():
            try:
                data = json.loads(result_json.read_text())
                nipt = data.get("NIPT", {})
                fr   = nipt.get("final_results", {})

                c1, c2, c3 = st.columns(3)
                c1.metric("Gender",     fr.get("fetal_gender","—"))
                c2.metric("FF (seqFF)", fr.get("fetal_fraction_seqff","—"))
                c3.metric("Sample Bias",fr.get("sample_bias","—"))

                st.markdown("**Trisomy Results**")
                tris = nipt.get("trisomy_results", [])
                if tris:
                    st.dataframe(pd.DataFrame(tris), use_container_width=True, hide_index=True)
                else:
                    st.info("No trisomy results in JSON.")

                with st.expander("Full JSON"):
                    st.json(data)
            except Exception as e:
                st.error(str(e))
        else:
            st.warning("result.json not found.")

    with tab_qc:
        qc_txt = out_path / "Output_QC" / f"{sample_choice}.qc.filter.txt"
        if qc_txt.exists():
            try:
                df_qc = pd.read_csv(qc_txt, sep="\t", header=None,
                                    names=["metric","value","threshold","op","status"])
                st.dataframe(df_qc, use_container_width=True, hide_index=True)
            except Exception:
                st.code(qc_txt.read_text())
        else:
            st.warning("QC filter file not found.")

    with tab_files:
        if out_path.exists():
            all_files = sorted(out_path.rglob("*"))
            file_rows = [{"Path": str(f.relative_to(out_path)),
                          "Size (KB)": round(f.stat().st_size/1024,1)}
                         for f in all_files if f.is_file()]
            st.dataframe(pd.DataFrame(file_rows), use_container_width=True,
                         hide_index=True, height=400)
        else:
            st.warning(f"Output directory not found: `{out_path}`")

# ─── FF distribution chart (if multiple samples) ──────────────────────────────
if "FF (seqFF)" in filtered.columns and len(filtered) > 1:
    st.divider()
    st.subheader("Fetal Fraction Distribution")

    df_ff = filtered[["Sample", "FF (seqFF)", "FF (YFF)"]].copy() if "FF (YFF)" in filtered.columns \
            else filtered[["Sample", "FF (seqFF)"]].copy()
    for col in ["FF (seqFF)", "FF (YFF)"]:
        if col in df_ff.columns:
            df_ff[col] = pd.to_numeric(df_ff[col], errors="coerce")

    ff_vals = df_ff["FF (seqFF)"].dropna()

    _chart_layout = dict(
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        font_color="#fafafa",
    )
    _axis_style = dict(gridcolor="rgba(255,255,255,0.15)", zerolinecolor="rgba(255,255,255,0.3)")

    if not ff_vals.empty:
        chart_col, table_col = st.columns([3, 2])

        with chart_col:
            if len(ff_vals) < 10:
                # Few samples: bar chart per sample is clearer than histogram
                fig = px.bar(
                    df_ff.dropna(subset=["FF (seqFF)"]),
                    x="Sample", y="FF (seqFF)",
                    title="seqFF per Sample",
                    text_auto=".2f",
                )
                fig.update_traces(textposition="outside")
            else:
                fig = px.histogram(ff_vals, nbins=20, labels={"value": "FF (seqFF)"},
                                   title="seqFF Distribution")
            fig.update_layout(**_chart_layout)
            fig.update_xaxes(**_axis_style)
            fig.update_yaxes(**_axis_style)
            st.plotly_chart(fig, use_container_width=True)

        with table_col:
            st.markdown("**Raw FF values**")
            st.dataframe(df_ff, use_container_width=True, hide_index=True)
