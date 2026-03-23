#!/bin/bash

# Usage: create_batch_files.sh <subject_dir> [output_root] [--hypothalamus true|false] [--surface true|false] [--qc-images true|false]

# Check if the path is valid
if [ ! -d "$1" ]; then
    echo "Invalid path: $1"
    exit 1
fi

# Check for the existence of "jobs" and "logs" directories, create them if they do not exist
[ ! -d "jobs" ] && mkdir jobs
[ ! -d "logs" ] && mkdir logs

# Get the base name of the input directory
base_name=$(basename $1)

# Optional second positional argument: separate output root directory
# If it starts with '--' it is a flag, not an output root
if [[ "${2:-}" == --* ]] || [[ -z "${2:-}" ]]; then
    output_root=""
    shift 1
else
    output_root="$2"
    shift 2
fi

# Parse optional flags
matlab_flags=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --hypothalamus)
            matlab_flags="${matlab_flags}, 'hypothalamus', ${2}"
            shift 2
            ;;
        --surface)
            matlab_flags="${matlab_flags}, 'surface', ${2}"
            shift 2
            ;;
        --qc-images)
            matlab_flags="${matlab_flags}, 'qc_images', ${2}"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Create a sbatch file in the jobs folder
cat << EOF > jobs/CAT_${base_name}.sbatch
#!/bin/tcsh
#SBATCH --job-name=CAT_${base_name:0:9}
#SBATCH --output=`pwd`/logs/${base_name}.out
#SBATCH --error=`pwd`/logs/${base_name}.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH -p defq
#SBATCH -A ansir-users  # UPDATE THIS!
#SBATCH -W
module load matlab
cd /home/richard/RunCat12 # Updated by install.sh
matlab -r "addpath('/home/richard/RunCat12/libs/spm12/spm12'); run_new_cat_normseg('$1', '$output_root'${matlab_flags}); exit;"
wait
EOF
