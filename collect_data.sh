#!/usr/bin/env bash
# collect_data.sh — Set up Python environment and run extract_roi_values.py in batch mode.
#
# Usage:
#   bash collect_data.sh <project_dir> [output.csv]
#
# Arguments:
#   project_dir   Directory containing subject folders (required)
#   output.csv    Output filename (optional, default: cat12_hypothalamus_data.csv)

set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ---------------------------------------------------------------------------
# Arguments
# ---------------------------------------------------------------------------
if [ $# -lt 1 ]; then
    echo "Usage: bash collect_data.sh <project_dir> [output.csv]"
    exit 1
fi

PROJECT_DIR="$1"
OUTPUT_CSV="${2:-cat12_hypothalamus_data.csv}"

if [ ! -d "$PROJECT_DIR" ]; then
    echo "Error: project directory not found: $PROJECT_DIR"
    exit 1
fi

# ---------------------------------------------------------------------------
# Check Python 3
# ---------------------------------------------------------------------------
if ! command -v python3 &>/dev/null; then
    echo "Error: python3 not found. Install Python 3.8+ and re-run."
    exit 1
fi
echo "[ok] Python 3: $(python3 --version)"

# ---------------------------------------------------------------------------
# Create / activate virtual environment
# ---------------------------------------------------------------------------
VENV_DIR="$REPO_DIR/data_collection_scripts/venv"
if [ ! -d "$VENV_DIR" ]; then
    echo "[install] Creating virtual environment at $VENV_DIR ..."
    python3 -m venv "$VENV_DIR"
    echo "[ok] Virtual environment created."
else
    echo "[skip] Virtual environment already exists."
fi

# shellcheck disable=SC1091
source "$VENV_DIR/bin/activate"

# ---------------------------------------------------------------------------
# Install dependencies
# ---------------------------------------------------------------------------
echo "[install] Installing Python dependencies..."
pip install --quiet -r "$REPO_DIR/data_collection_scripts/requirements.txt"
echo "[ok] Dependencies installed."

# ---------------------------------------------------------------------------
# Warn if hypothalamus atlas is missing
# ---------------------------------------------------------------------------
ATLAS_PATH="$(grep 'HYPOTHALAMUS_ATLAS_PATH' "$REPO_DIR/data_collection_scripts/extract_roi_values.py" \
    | head -1 | sed 's/.*= *"\(.*\)".*/\1/')"

if [ ! -f "$ATLAS_PATH" ]; then
    echo ""
    echo "WARNING: Hypothalamus atlas not found at:"
    echo "  $ATLAS_PATH"
    echo "Run 'bash install.sh' first to download libraries and patch paths."
    echo ""
fi

# ---------------------------------------------------------------------------
# Run extraction
# ---------------------------------------------------------------------------
echo ""
echo "Running extract_roi_values.py on: $PROJECT_DIR"
echo "Output CSV: $OUTPUT_CSV"
echo ""
python3 "$REPO_DIR/data_collection_scripts/extract_roi_values.py" \
    "$PROJECT_DIR" --batch -o "$OUTPUT_CSV"

echo ""
echo "=== Done ==="
echo "Results written to: $(pwd)/$OUTPUT_CSV"
