"""
Lab Config Editor — view and edit pipeline_config.json per labcode.
"""
import sys, json, copy
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parents[1]))

import streamlit as st
from lib.paths import LABS, lab_config_path, CONFIG_DIR

st.set_page_config(page_title="Lab Config | gx-nipt Admin", page_icon="🏥", layout="wide")
st.title("🏥 Lab Config Editor")
st.caption("View and edit `pipeline_config.json` for each labcode.")

if not LABS:
    st.error(f"No lab configs found in `{CONFIG_DIR}`.")
    st.stop()

# ── Lab selector ──────────────────────────────────────────────────────────────
lab = st.selectbox("Select Lab", sorted(LABS))
cfg_path = lab_config_path(lab)

if not cfg_path.exists():
    st.error(f"Config not found: `{cfg_path}`")
    st.stop()

config = json.loads(cfg_path.read_text())

# Sections hidden from tabs (still preserved when saving)
HIDDEN_SECTIONS = {"WCFF"}
# Sections rendered with 3-column layout
THREE_COL_SECTIONS = {"QC", "FF_Gender_Config", "WCX"}

# ── Tabs per section ──────────────────────────────────────────────────────────
visible_sections = [s for s in config.keys() if s not in HIDDEN_SECTIONS]
tabs = st.tabs(visible_sections)

edited = copy.deepcopy(config)


def _parse_number(raw: str, original):
    """Parse a text input back to int or float, matching original type."""
    try:
        parsed = float(raw)
        if isinstance(original, int) and parsed == int(parsed):
            return int(parsed)
        return parsed
    except ValueError:
        return None


for tab, section in zip(tabs, visible_sections):
    with tab:
        data = config[section]
        if not isinstance(data, dict):
            st.write(data)
            continue

        st.markdown(f"#### {section}")
        new_vals = {}

        # ── MD_Target sections: show bed + regions 1-4 clearly ─────────────────
        if section.startswith("MD_Target"):
            if "bed" in data:
                new_vals["bed"] = st.text_input(
                    "bed", value=data["bed"], key=f"{lab}_{section}_bed"
                )

            region_keys = sorted(k for k in data if k.startswith("region"))
            if region_keys:
                region_cols = st.columns(len(region_keys))
                for col, rkey in zip(region_cols, region_keys):
                    with col:
                        st.markdown(f"**{rkey}**")
                        rval = data[rkey]
                        new_rval = {}
                        for subk, subv in rval.items():
                            display = str(int(subv)) if float(subv) == int(float(subv)) else str(subv)
                            raw = st.text_input(subk, value=display,
                                                key=f"{lab}_{section}_{rkey}_{subk}")
                            parsed = _parse_number(raw, subv)
                            if parsed is None:
                                st.error(f"Invalid: {raw}")
                                new_rval[subk] = subv
                            else:
                                new_rval[subk] = parsed
                        new_vals[rkey] = new_rval

        # ── All other sections ──────────────────────────────────────────────────
        else:
            ncols = 3 if section in THREE_COL_SECTIONS else 2
            cols = st.columns(ncols)
            for i, (key, val) in enumerate(data.items()):
                with cols[i % ncols]:
                    if isinstance(val, bool):
                        new_vals[key] = st.checkbox(
                            key, value=val, key=f"{lab}_{section}_{key}"
                        )
                    elif isinstance(val, (int, float)):
                        raw = st.text_input(
                            key, value=str(val), key=f"{lab}_{section}_{key}"
                        )
                        parsed = _parse_number(raw, val)
                        if parsed is None:
                            st.error(f"Invalid number for {key}")
                            new_vals[key] = val
                        else:
                            new_vals[key] = parsed
                    elif isinstance(val, str):
                        new_vals[key] = st.text_input(
                            key, value=val, key=f"{lab}_{section}_{key}"
                        )
                    elif isinstance(val, dict):
                        st.markdown(f"**{key}**")
                        st.json(val, expanded=False)
                        new_vals[key] = val
                    else:
                        st.write(f"{key}: {val}")
                        new_vals[key] = val

        edited[section] = new_vals

# ── Diff & Save ───────────────────────────────────────────────────────────────
st.divider()
col_diff, col_save = st.columns([3, 1])

with col_diff:
    if st.button("🔍 Preview Changes"):
        original_str = json.dumps(config, indent=2)
        edited_str   = json.dumps(edited, indent=2)
        if original_str == edited_str:
            st.info("No changes detected.")
        else:
            st.code(edited_str, language="json")

with col_save:
    if st.button("💾 Save Config", type="primary"):
        cfg_path.write_text(json.dumps(edited, indent=4))
        st.success(f"Saved → `{cfg_path}`")
        st.rerun()

st.divider()

# ── Raw JSON viewer ───────────────────────────────────────────────────────────
with st.expander("📄 Raw JSON"):
    st.json(config)

# ── Add new lab ───────────────────────────────────────────────────────────────
with st.expander("➕ Add New Lab (copy from existing)"):
    src_lab  = st.selectbox("Copy from", sorted(LABS), key="new_lab_src")
    new_name = st.text_input("New lab code (e.g. newlab)")
    if st.button("Create"):
        if not new_name:
            st.warning("Enter a lab code.")
        elif new_name in LABS:
            st.warning(f"Lab `{new_name}` already exists.")
        else:
            new_path = CONFIG_DIR / new_name / "pipeline_config.json"
            new_path.parent.mkdir(parents=True, exist_ok=True)
            src_cfg = json.loads(lab_config_path(src_lab).read_text())
            new_path.write_text(json.dumps(src_cfg, indent=4))
            st.success(f"Created `{new_path}`. Reload to edit.")
