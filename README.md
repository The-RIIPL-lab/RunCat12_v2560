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

*   Create the Python environment:
    ```bash
    pip install -r ./data_collection_scripts/requirements.txt
    ```

*   Update the `HYPOTHALAMUS_ATLAS_PATH` and `JHU_ATLAS_PATH` in the `extract_roi_values.py` script.

*   Run the script:
    ```bash
    python ./data_collection_scripts/extract_roi_values.py <your project folder>
    ```