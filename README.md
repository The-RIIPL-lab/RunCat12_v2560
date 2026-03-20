# Post Processing Script for RIIPL Data

This project is a post-processing pipeline for existing RIIPL data using MATLAB and CAT12 r2560. It creates native-space labels and volume measurements, optional surface maps, and collects diffusion metrics from hypothalamus ROIs.

## Requirements

*   MATLAB 2022b+
*   CAT12 r2560 (downloaded automatically by `install.sh`)
*   SPM12 r7487 (downloaded automatically by `install.sh`)
*   SLURM HPC system (designed for the DEMON HPC)
*   Python 3.8+ (for data collection scripts)

## Setup

Run the install script once after cloning to download SPM12 and CAT12 r2560, copy custom ROI atlases, and patch all hardcoded paths:

```bash
bash install.sh
```

Then open `create_batch_files.sh` and update the one line that cannot be auto-detected:

```
#SBATCH -A ansir-users   →   #SBATCH -A <your-hpc-account>
```

## Workflow

### STEP 1: Create Jobs

Run the following loop to generate one `.sbatch` file per subject into `jobs/`:

```bash
for x in $(cat ./subject_lists/baseline_ADRC_50_subjects.txt); do
    if [ -d "/<your project folder>/${x}" ] && [ ! -d "/<your project folder>/${x}/nifti/cat12_v2560" ]; then
        bash ./create_batch_files.sh "/<your project folder>/${x}"
    fi
done
```

`create_batch_files.sh` also accepts an optional second argument to write outputs to a separate directory:

```bash
bash ./create_batch_files.sh /path/to/SUBJECT_ID /path/to/output_root
```

### STEP 2: Run Jobs

Run in a tmux session from the `jobs/` directory (10 jobs at a time):

```bash
cd jobs && ls | grep sbatch | xargs -I {} -P 10 sbatch -W {}
```

Logs are written to `logs/` (one `.out` and `.err` per subject).

### STEP 3: Collect Data

Install Python dependencies and run the extraction script in one step:

```bash
bash collect_data.sh /path/to/project/ [output.csv]
```

Or run the Python script directly:

```bash
# Single subject
python ./data_collection_scripts/extract_roi_values.py /path/to/SUBJECT_ID

# Batch (entire project directory)
python ./data_collection_scripts/extract_roi_values.py --batch /path/to/project/ -o results.csv --verbose
```

QC overlay images (hypothalamus ROIs on native T1):

```bash
python ./data_collection_scripts/hypothalamus_roi_qc.py --batch /path/to/project/
```

## Processing Options

`run_new_cat_normseg` accepts optional name-value flags after the two positional arguments. All flags default to values that match the original behaviour.

| Flag | Default | Description |
|---|---|---|
| `'hypothalamus'` | `true` | Include the hypothalamus atlas in the native-space warp and atlas list |
| `'surface'` | `false` | Enable CAT12 surface mapping and surface measures |
| `'qc_images'` | `false` | Generate SPM registration QC PNG images in `QC_registration/` |

When `'surface'` is `true`, a `surface_aal3_stats.csv` file is written to the CAT12 output folder containing per-ROI cortical thickness, gyrification, and sulcal depth for all AAL3 regions (left and right hemisphere).

To pass flags via a SLURM batch file, edit the MATLAB command in `create_batch_files.sh`:

```bash
# Enable surface mapping
matlab -r "addpath('<spm_dir>'); run_new_cat_normseg('$1', '$output_root', 'surface', true); exit;"

# Disable hypothalamus atlas, enable QC images
matlab -r "addpath('<spm_dir>'); run_new_cat_normseg('$1', '$output_root', 'hypothalamus', false, 'qc_images', true); exit;"
```

## Output Structure

```
SUBJECT_ID/
└── nifti/
    └── cat12_v2560/
        ├── w*.nii                  ← native-space atlas labelmaps (w prefix)
        ├── surface_aal3_stats.csv  ← surface ROI stats (if 'surface' = true)
        ├── QC_registration/        ← registration QC images (if 'qc_images' = true)
        └── mri/
            ├── p1*.nii             ← GM probability map
            ├── p2*.nii             ← WM probability map
            ├── p3*.nii             ← CSF probability map
            ├── m*.nii              ← bias-corrected native T1
            ├── DTI_native/         ← coregistered DTI (r* prefix)
            ├── DTI/                ← normalised DTI (wr* prefix)
            ├── NODDI_native/       ← coregistered NODDI (r* prefix)
            └── NODDI/              ← normalised NODDI (wr* prefix)
```

## `extract_roi_values.py` Options

| Flag | Default | Description |
|---|---|---|
| `--batch` | off | Iterate over all subdirectories of `input_path` |
| `-o / --output` | `cat12_hypothalamus_data.csv` | Output CSV filename |
| `--hyp-atlas` | auto (repo-relative) | Override path to the hypothalamus atlas template |
| `--verbose` | off | Print per-subject status during batch runs |

**Output columns:**

| Column group | Description |
|---|---|
| `subject_id`, `processing_status` | Subject name and outcome |
| `brain_mask_volume_mm3` | Whole-brain mask volume |
| `hypothalamus_roi1_volume_mm3`, `hypothalamus_roi2_volume_mm3` | Hypothalamus left/right volumes |
| `gm_volume_mm3`, `wm_volume_mm3`, `csf_volume_mm3` | Tissue volumes |
| `native_DTI_{FA,MD,L1,L2,L3}_hyp_roi{1,2}_mean` | Mean DTI metrics in native-space hypothalamus |
| `native_NODDI_{ICVF,ISOVF,OD}_hyp_roi{1,2}_mean` | Mean NODDI metrics in native-space hypothalamus |
| `normalized_DTI_{FA,MD,L1,L2,L3}_hyp_roi{1,2}_mean` | Mean DTI metrics in template-space hypothalamus |
| `normalized_NODDI_{ICVF,ISOVF,OD}_hyp_roi{1,2}_mean` | Mean NODDI metrics in template-space hypothalamus |

## Acknowledgements

*   [CAT12](https://neuro-jena.github.io/cat/) (r2560)
*   [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/) (r7487)
*   MATLAB
