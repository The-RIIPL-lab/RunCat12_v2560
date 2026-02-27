It looks like the `README.md` file wasn’t found in the current directory, so I couldn’t update it directly.

Here’s the updated version of your README with the additional details you requested. You can copy and paste this into your file:

***

# Post Processing Script for RIIPL Data

This project is a post-processing script designed to be run on existing RIIPL processed data using Matlab and Cat12 v2560. The intent is to use updated software to create native space labels and volume measurements, surface maps, and processing reports.

## Requirements

*   Matlab 2022b+
*   CAT12 v2560
*   SPM12
*   SLURM HPC system (designed for the DEMON HPC)

## Setup

1.  Ensure you have the required software installed:
    *   Matlab 2022b or later
    *   CAT12 v2560
    *   SPM12

2.  This script is intended to be run on a SLURM HPC system. Make sure you have access to such a system, preferably the DEMON HPC.

## Customization

The `create_batch_files.sh` script will require some customization:

*   Update the script path to point to your specific setup.
*   Add your HPC account name.

## Creating and Submitting Batch Jobs

### Quick setup

Run the install script once to download SPM12 + CAT12 and patch all hardcoded paths:

```bash
bash install.sh
```

After that, open `create_batch_files.sh` and update the one line that cannot be auto-detected:

```
#SBATCH -A ansir-users   →   #SBATCH -A <your-hpc-account>
```

### Using create_batch_files.sh

The script takes a single subject directory as its argument and writes a
ready-to-submit SLURM `.sbatch` file into a `jobs/` folder:

```bash
bash ./create_batch_files.sh /path/to/project/SUBJECT_ID
```

To process a list of subjects in bulk (skipping any already processed):

```bash
for x in $(cat ./subject_lists/my_subjects.txt); do
    if [ -d "/project/${x}" ] && [ ! -d "/project/${x}/nifti/cat12_v2560" ]; then
        bash ./create_batch_files.sh "/project/${x}"
    fi
done
```

This creates one `.sbatch` file per subject in `jobs/`.
To submit all jobs (10 at a time) in a tmux session:

```bash
cd jobs && ls | grep sbatch | xargs -I {} -P 10 sbatch -W {}
```

Logs are written to the `logs/` directory (one `.out` and `.err` per subject).

## Usage

1.  Customize the `create_batch_files.sh` script as described above.
2.  Submit the job to the SLURM HPC system.

## Output

The script will generate:

*   Native space labels and volume measurements
*   Surface maps
*   Processing reports

## Acknowledgements

This project uses:

*   Matlab
*   CAT12
*   SPM12

***

## Detailed Steps to Run the Code

The software works like this:

1.  You create your jobs
2.  You run the jobs in the DEMON HPC system
3.  You use Python scripts to collect your data

### STEP 0: Get the Software

*   Clone the repository:
    ```bash
    git clone https://github.com/The-RIIPL-lab/RunCat12_v2560
    ```

*   Update the code for your location:
    *   Update the `create_batch_files` script to match your account number.
    *   Update your working directory.

### STEP 1: Create Your Jobs

Run the following loop to create jobs:

```bash
for x in `cat ./subject_lists/baseline_ADRC_50_subjects.txt`; do 
    if [ -d "/<your project folder>/${x}" ] && [ ! -d "/<your project folder>/${x}/nifti/cat12_v2560" ]; then 
        bash ./create_batch_files.sh /<your project folder>/${x}
    fi
done
```

> This will create jobs in the `job/` directory.

### STEP 2: Run Your Jobs

(I recommend running these in a tmux session)

```bash
ls | grep sbatch | xargs -I {} -P 10 sbatch -W {}
```

> This will run 10 jobs at once until the list is complete.

### STEP 3: Collect Data Using Python Scripts

Install the required Python packages:

```bash
pip install -r ./data_collection_scripts/requirements.txt
```

#### `extract_roi_values.py`

Extracts volumes and diffusion metrics from CAT12 outputs and writes a single CSV file.

**Single subject:**

```bash
python ./data_collection_scripts/extract_roi_values.py /path/to/SUBJECT_ID
```

**Batch (entire project directory):**

```bash
python ./data_collection_scripts/extract_roi_values.py --batch /path/to/project/
```

The script automatically discovers the `cat12_v*` output folder inside each
subject's `nifti/` directory, so no path configuration is needed.

**All options:**

| Flag | Default | Description |
|---|---|---|
| `--batch` | off | Iterate over all subdirectories of `input_path` |
| `-o / --output` | `cat12_hypothalamus_data.csv` | Output CSV filename |
| `--hyp-atlas` | auto (repo-relative) | Override path to the hypothalamus atlas template |
| `--verbose` | off | Print per-subject status during batch runs |

**Expected subject directory structure:**

```
SUBJECT_ID/
└── nifti/
    └── cat12_v<build>/          ← discovered automatically
        ├── wbrainmask_T1.nii
        ├── whypothalamusAtlas.nii
        └── mri/
            ├── p1*.nii          ← GM probability map
            ├── p2*.nii          ← WM probability map
            ├── p3*.nii          ← CSF probability map
            ├── DTI_native/      ← coregistered DTI (r* prefix)
            ├── DTI/             ← normalised DTI (wr* prefix)
            ├── NODDI_native/    ← coregistered NODDI (r* prefix)
            └── NODDI/           ← normalised NODDI (wr* prefix)
```

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