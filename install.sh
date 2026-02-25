#!/usr/bin/env bash
# install.sh — Download SPM12 + CAT12, copy custom ROIs, and patch hardcoded paths.
# Run once from the repo root before submitting any jobs.

set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "=== RunCat12 install ==="
echo "Repo root: $REPO_DIR"
echo ""

# ---------------------------------------------------------------------------
# 1. Create libs/ parent directory
# ---------------------------------------------------------------------------
mkdir -p "$REPO_DIR/libs/spm12"

# ---------------------------------------------------------------------------
# 2. Download and extract SPM12 r7487
# ---------------------------------------------------------------------------
SPM_DIR="$REPO_DIR/libs/spm12/spm12"
if [ -d "$SPM_DIR" ] && [ -n "$(ls -A "$SPM_DIR")" ]; then
    echo "[skip] SPM12 already present at $SPM_DIR"
else
    echo "[install] Downloading SPM12 r7487..."
    TMP_SPM=$(mktemp -d)
    curl -fL "https://github.com/spm/spm12/archive/refs/tags/r7487.zip" -o "$TMP_SPM/spm12.zip"
    echo "[install] Extracting SPM12..."
    unzip -q "$TMP_SPM/spm12.zip" -d "$TMP_SPM"
    rm "$TMP_SPM/spm12.zip"
    # GitHub source archives extract to a single subdirectory (e.g. spm12-r7487/)
    EXTRACTED_SPM=$(find "$TMP_SPM" -maxdepth 1 -mindepth 1 -type d | head -1)
    mv "$EXTRACTED_SPM" "$SPM_DIR"
    rm -rf "$TMP_SPM"
    echo "[ok] SPM12 installed at $SPM_DIR"
fi

# Ensure toolbox dir exists now that SPM12 is in place
mkdir -p "$SPM_DIR/toolbox"

# ---------------------------------------------------------------------------
# 3. Download and extract CAT12 12.8.2
# ---------------------------------------------------------------------------
CAT_DIR="$SPM_DIR/toolbox/cat12"
if [ -d "$CAT_DIR" ] && [ -n "$(ls -A "$CAT_DIR")" ]; then
    echo "[skip] CAT12 already present at $CAT_DIR"
else
    echo "[install] Downloading CAT12 12.8.2..."
    TMP_CAT=$(mktemp -d)
    curl -fL "https://github.com/ChristianGaser/cat12/releases/download/12.8.2/cat12.8.2.zip" -o "$TMP_CAT/cat12.zip"
    echo "[install] Extracting CAT12..."
    unzip -q "$TMP_CAT/cat12.zip" -d "$TMP_CAT/extracted"
    rm "$TMP_CAT/cat12.zip"
    # Release zips may have a single subdirectory or files at the root
    TOP_DIRS=$(find "$TMP_CAT/extracted" -maxdepth 1 -mindepth 1 -type d | wc -l)
    TOP_FILES=$(find "$TMP_CAT/extracted" -maxdepth 1 -mindepth 1 -type f | wc -l)
    if [ "$TOP_DIRS" -eq 1 ] && [ "$TOP_FILES" -eq 0 ]; then
        EXTRACTED_CAT=$(find "$TMP_CAT/extracted" -maxdepth 1 -mindepth 1 -type d | head -1)
        mv "$EXTRACTED_CAT" "$CAT_DIR"
    else
        mv "$TMP_CAT/extracted" "$CAT_DIR"
    fi
    rm -rf "$TMP_CAT"
    echo "[ok] CAT12 installed at $CAT_DIR"
fi

# ---------------------------------------------------------------------------
# 4. Copy custom ROIs into the CAT12 templates directory
# ---------------------------------------------------------------------------
TEMPLATES_DIR="$CAT_DIR/templates_MNI152NLin2009cAsym"
mkdir -p "$TEMPLATES_DIR"

echo "[install] Copying custom ROIs into $TEMPLATES_DIR ..."
cp "$REPO_DIR/resampled_template_ROIs/rhypothalamusAtlas_template_v2.nii" \
   "$TEMPLATES_DIR/hypothalamusAtlas.nii"
cp "$REPO_DIR/resampled_template_ROIs/rJHU.nii" \
   "$TEMPLATES_DIR/JHU.nii"
echo "[ok] ROIs copied."

# ---------------------------------------------------------------------------
# 5. Patch run_new_cat_normseg.m
# ---------------------------------------------------------------------------
MATLAB_FILE="$REPO_DIR/run_new_cat_normseg.m"
OLD_PREFIX="/isilon/datalake/riipl/original/DEMONco/Hellcat-12.9/"

if grep -qF "$OLD_PREFIX" "$MATLAB_FILE"; then
    cp "$MATLAB_FILE" "${MATLAB_FILE}.bak"
    sed -i "s|${OLD_PREFIX}|${REPO_DIR}/|g" "$MATLAB_FILE"
    echo "[ok] Patched $MATLAB_FILE (backup: ${MATLAB_FILE}.bak)"
else
    echo "[skip] $MATLAB_FILE — old prefix not found (already patched?)"
fi

# ---------------------------------------------------------------------------
# 6. Patch create_batch_files.sh
# ---------------------------------------------------------------------------
BATCH_FILE="$REPO_DIR/create_batch_files.sh"

cp "$BATCH_FILE" "${BATCH_FILE}.bak"

# Patch the cd line — matches both the original placeholder and any previously-patched line
if grep -qE "cd .* # (UPDATE THIS|Updated by install\.sh)" "$BATCH_FILE"; then
    sed -i "s|cd .* # .*|cd ${REPO_DIR} # Updated by install.sh|" "$BATCH_FILE"
    echo "[ok] Patched cd path in $BATCH_FILE"
else
    echo "[warn] $BATCH_FILE — cd line not found in expected format; skipping cd patch"
fi

# Patch the matlab command to include SPM12 addpath
if grep -q "addpath.*run_new_cat_normseg" "$BATCH_FILE"; then
    # Already has addpath — update the path in place
    sed -i "s|addpath('[^']*'); run_new_cat_normseg|addpath('${SPM_DIR}'); run_new_cat_normseg|g" "$BATCH_FILE"
    echo "[ok] Updated SPM12 addpath in matlab command in $BATCH_FILE"
elif grep -q 'matlab -r "run_new_cat_normseg' "$BATCH_FILE"; then
    # No addpath yet — insert it before run_new_cat_normseg
    sed -i "s|matlab -r \"run_new_cat_normseg|matlab -r \"addpath('${SPM_DIR}'); run_new_cat_normseg|g" "$BATCH_FILE"
    echo "[ok] Added SPM12 addpath to matlab command in $BATCH_FILE"
else
    echo "[warn] $BATCH_FILE — matlab command not found in expected format; skipping addpath patch"
fi
echo "       (backup: ${BATCH_FILE}.bak)"

# ---------------------------------------------------------------------------
# 7. Patch data_collection_scripts/extract_roi_values.py
# ---------------------------------------------------------------------------
PY_FILE="$REPO_DIR/data_collection_scripts/extract_roi_values.py"
NEW_ATLAS_PATH="${TEMPLATES_DIR}/hypothalamusAtlas.nii"

if grep -q 'HYPOTHALAMUS_ATLAS_PATH' "$PY_FILE"; then
    cp "$PY_FILE" "${PY_FILE}.bak"
    sed -i "s|HYPOTHALAMUS_ATLAS_PATH = \".*\"|HYPOTHALAMUS_ATLAS_PATH = \"${NEW_ATLAS_PATH}\"|" "$PY_FILE"
    echo "[ok] Patched $PY_FILE (backup: ${PY_FILE}.bak)"
else
    echo "[skip] $PY_FILE — HYPOTHALAMUS_ATLAS_PATH not found"
fi

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
echo ""
echo "=== Done ==="
echo ""
echo "Installed:"
echo "  SPM12 : $SPM_DIR"
echo "  CAT12 : $CAT_DIR"
echo "  Atlas : $TEMPLATES_DIR/hypothalamusAtlas.nii"
echo "  JHU   : $TEMPLATES_DIR/JHU.nii"
echo ""
echo "Patched files:"
echo "  $MATLAB_FILE"
echo "  $BATCH_FILE"
echo "  $PY_FILE"
echo ""
echo "ACTION REQUIRED:"
echo "  Open create_batch_files.sh and update the SLURM account:"
echo "    #SBATCH -A ansir-users  →  #SBATCH -A <your-hpc-account>"
echo ""
echo "Verify the patches look correct before submitting any jobs."
