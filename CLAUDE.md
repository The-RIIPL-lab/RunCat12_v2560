# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a neuroimaging post-processing pipeline for RIIPL data. It runs CAT12 v2560 segmentation and normalization on T1w MRI data, co-registers and normalizes diffusion (DTI, NODDI) and perfusion (ASL) images, exports atlases to native space, and collects ROI metrics via Python.

**Target system:** DEMON HPC (SLURM), Matlab 2022b+, CAT12 v2560, SPM12.

## Required Customizations Before Use

Two files must be updated with local paths/account:

- `create_batch_files.sh`: Update `#SBATCH -A ansir-users` (HPC account) and `cd /isilon/datalake/...` (working directory)
- `run_new_cat_normseg.m`: Update `addpath(...)` calls pointing to SPM12/CAT12, `label_map_path`, `template_path`, and the TPM path inside `run_cat12_segmentation()`
- `data_collection_scripts/extract_roi_values.py`: Update `HYPOTHALAMUS_ATLAS_PATH`

## Workflow

### Step 1: Create SLURM jobs

```bash
for x in `cat ./subject_lists/<subject_list>.txt`; do
    if [ -d "/<project_folder>/${x}" ] && [ ! -d "/<project_folder>/${x}/nifti/cat12_v2560" ]; then
        bash ./create_batch_files.sh /<project_folder>/${x}
    fi
done
```

Jobs are written to `jobs/`, logs go to `logs/`.

### Step 2: Submit jobs (run inside a tmux session)

```bash
ls jobs/ | grep sbatch | xargs -I {} -P 10 sbatch -W jobs/{}
```

### Step 3: Install Python dependencies and collect data

```bash
pip install -r ./data_collection_scripts/requirements.txt

# Single subject
python ./data_collection_scripts/extract_roi_values.py <subject_dir>

# Batch (parent directory containing multiple subject folders)
python ./data_collection_scripts/extract_roi_values.py <parent_dir> --batch -o output.csv
```

## Architecture

### MATLAB Pipeline (`run_new_cat_normseg.m`)

Entry point called from SLURM jobs: `matlab -r "run_new_cat_normseg('<subject_dir>'); exit;"`

The function orchestrates 5 steps for each subject:

1. **CAT12 segmentation** – runs `spm_jobman` to segment/normalize the T1w image, producing deformation fields (`y_*`, `iy_*`) and tissue maps (`p1*`, `p2*`, `p3*`) in `nifti/cat12_v2560/mri/`
2. **DTI co-registration & normalization** – co-registers DTI metrics to native T1, moves `r*` files to `mri/DTI_native/`, then warps them to MNI space (`wr*`) in `mri/DTI/`
3. **NODDI co-registration & normalization** – same pattern as DTI, output in `mri/NODDI_native/` and `mri/NODDI/`
4. **ASL co-registration & normalization** – handles multiple ASL runs (`ASL_native_run1/`, etc.) and normalized outputs in `mri/ASL/`
5. **Atlas-to-native-space export** – applies inverse deformation field to warp template atlases (neuromorphometrics, AAL3, JHU, hypothalamusAtlas, etc.) into native space with `w` prefix in `nifti/cat12_v2560/`

Helper file-discovery functions: `getASLFiles.m`, `getDTIFiles.m`, `getNODDIFiles.m`

**Subject cohort prefixes** (used to find T1w images):
- MARVEL: `M*-tfl3d*ns.nii`
- ADRC: `3*-tfl3d116ns.nii`
- SWITCH: `S*-tfl3d*ns.nii`

DTI files are identified by `*_ECC_*.nii` naming; NODDI by `*-NODDI_*.nii`; ASL by `M*_M0_masked.nii` + `*_CBF*.nii`/`*_CVR*.nii`. All three filter out already-processed files (starting with lowercase `r`, `w`, `m`).

**Output file prefix conventions:**
- `r*` – co-registered to native T1 space
- `wr*` – normalized to MNI template space
- `w*` (atlases) – warped from template to native space using inverse deformation
- `p1*`, `p2*`, `p3*` – GM, WM, CSF tissue probability maps

QC PNG images are saved to `nifti/cat12_v2560/QC_registration/`.

### Python Data Collection (`data_collection_scripts/`)

- **`extract_roi_values.py`** – main script; expects `<subject>/nifti/cat12_v2560/` to exist. Extracts: brain mask volume, hypothalamus ROI 1 & 2 volumes (from `whypothalamusAtlas.nii`), GM/WM/CSF volumes (thresholded at 0.5), and native/normalized DTI+NODDI mean values per hypothalamus ROI. Outputs a CSV.
- **`brain_data_exploration.py`** – visualization/analysis of the extracted CSV (cross-sectional; averages longitudinal subjects).

### Template ROIs (`resampled_template_ROIs/`)

Contains pre-resampled atlas files in MNI space used for normalized-space diffusion analysis, including the hypothalamus atlas (`rhypothalamusAtlas_template_v2.nii`).
