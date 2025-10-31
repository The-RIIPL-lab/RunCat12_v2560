#!/bin/bash

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

# Create a sbatch file in the jobs folder
cat << EOF > jobs/CAT_${base_name}.sbatch
#!/bin/tcsh
#SBATCH --job-name=CAT_${base_name:0:9}
#SBATCH --output=`pwd`/logs/${base_name}.out
#SBATCH --error=`pwd`/logs/${base_name}.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=20:00
#SBATCH -p defq
#SBATCH -A ansir-users
#SBATCH -W
module load matlab
cd /isilon/datalake/riipl/original/DEMONco/Hellcat-12.9/
matlab -r "run_new_cat_normseg('$1'); exit;"
wait
EOF
