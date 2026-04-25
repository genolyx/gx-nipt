#!/usr/bin/env bash
# gx-nipt Admin UI launcher
# Usage: bash admin/start.sh [port]
set -euo pipefail

PORT=${1:-8501}
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

# Install deps if needed
pip install -q -r "${SCRIPT_DIR}/requirements.txt"

echo "Starting gx-nipt Admin on http://0.0.0.0:${PORT}"
streamlit run "${SCRIPT_DIR}/app.py" \
    --server.port "${PORT}" \
    --server.address 0.0.0.0 \
    --server.headless true \
    --browser.gatherUsageStats false
