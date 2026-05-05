"""
Centralised path definitions for the gx-nipt Admin UI.
All pages import from here so changing one variable fixes everything.
"""
from pathlib import Path

# ── Repo root (two levels up from admin/lib/) ────────────────────────────────
REPO_ROOT  = Path(__file__).resolve().parents[2]

# ── Runtime directories (outside git, on host) ───────────────────────────────
REF_DIR    = REPO_ROOT / "refs"
CONFIG_DIR = REPO_ROOT / "config"
OUTPUT_DIR = REPO_ROOT / "output"
LOG_DIR    = REPO_ROOT / "log"

# ── Scripts ──────────────────────────────────────────────────────────────────
SCRIPTS_DIR          = REPO_ROOT / "bin" / "scripts"
REF_SCRIPTS_DIR      = SCRIPTS_DIR / "utils" / "reference"
PATCH_NPZ_SCRIPT     = SCRIPTS_DIR / "utils" / "patch_npz.py"
RUN_NIPT_SH          = REPO_ROOT / "src" / "run_nipt.sh"
BUILD_GXCNV_SCRIPT   = SCRIPTS_DIR / "build_gxcnv_reference.py"
BUILD_GXFF_SCRIPT    = SCRIPTS_DIR / "build_gxff_model.py"
GXFF_PATCH_DIR       = REPO_ROOT / "bin" / "gxff_patch"

# ── Lab helpers ───────────────────────────────────────────────────────────────
LABS = [d.name for d in CONFIG_DIR.iterdir() if d.is_dir()] if CONFIG_DIR.exists() else []
GROUPS   = ["orig", "fetus", "mom"]
REF_TYPES = ["EZD", "PRIZM", "WC", "WCX"]

def lab_config_path(labcode: str) -> Path:
    return CONFIG_DIR / labcode / "pipeline_config.json"

def lab_ref_dir(labcode: str) -> Path:
    return REF_DIR / "labs" / labcode
