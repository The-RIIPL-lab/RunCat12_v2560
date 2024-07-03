
# Post Processing Script for RIIPL Data

This project is a post-processing script designed to be run on existing RIIPL processed data using Matlab and Cat12 v2560. The intent is to use updated software to create native space labels and volume measurements, surface maps, and processing reports.

## Requirements

- Matlab 2022b+
- CAT12 v2560
- SPM12
- SLURM HPC system (designed for the DEMON HPC)

## Setup

1. Ensure you have the required software installed:
   - Matlab 2022b or later
   - CAT12 v2560
   - SPM12

2. This script is intended to be run on a SLURM HPC system. Make sure you have access to such a system, preferably the DEMON HPC.

## Customization

The `create_batch_files.sh` script will require some customization:
- Update the script path to point to your specific setup.
- Add your HPC account name.

## Usage

1. Customize the `create_batch_files.sh` script as described above.
2. Submit the job to the SLURM HPC system.

## Output

The script will generate:
- Native space labels and volume measurements
- Surface maps
- Processing reports

## Acknowledgements

This project uses:
- [Matlab](https://www.mathworks.com/products/matlab.html)
- [CAT12](http://www.neuro.uni-jena.de/cat/)
- [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/)
