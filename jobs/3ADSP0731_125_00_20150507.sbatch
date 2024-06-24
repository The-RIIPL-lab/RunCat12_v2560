#!/bin/tcsh
#SBATCH --job-name=3ADSP0731_125_00_20150507
#SBATCH --output=/isilon/datalake/riipl/scratch/ADRC/Hellcat-12.9/logs/3ADSP0731_125_00_20150507.out
#SBATCH --error=/isilon/datalake/riipl/scratch/ADRC/Hellcat-12.9/logs/3ADSP0731_125_00_20150507.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=20:00
#SBATCH -p defq
#SBATCH -A ansir-users
module load matlab
matlab -r "run_new_cat_normseg('/isilon/datalake/riipl/original/ADPPR/3ADSP0731_125_00_20150507'); exit;"
