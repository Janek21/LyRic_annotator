#!/bin/bash

#SBATCH --output=logs/run/%x_%j.out
#SBATCH --error=logs/run/%x_%j.err
#SBATCH --qos=normal
#SBATCH --job-name=lyric

#SBATCH --mem=24G
#SBATCH --time=300

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4

#record start
start_time=$(date +%s)
echo ">STARTING at $(date)"

module load CMake/3.29.3-GCCcore-12.3.0 
module load Python/3.13.5-GCCcore-14.3.0

source ~/bin/snakemake/bin/activate
cpus="${SLURM_CPUS_PER_TASK:-2}"

#run from the repo root, pointing snakemake at the species clone
species_name="$1"
species_dir="$(realpath "$species_name")"
echo "Running LyRic for $species_name"

snakemake --snakefile "$species_dir/workflow/Snakefile" --directory "$species_dir" --configfile "$species_dir/config/default.yaml" --unlock
snakemake --snakefile "$species_dir/workflow/Snakefile" --directory "$species_dir" --configfile "$species_dir/config/default.yaml" --cores $cpus --keep-going

#record memory usage
cgroup_dir=$(awk -F: '{print $NF}' /proc/self/cgroup)
peak_mem=`cat /sys/fs/cgroup$cgroup_dir/memory.peak`
peak_mem_mb=$(awk "BEGIN {printf \"%.2f\", $peak_mem / 1048576}") #transfer to mb
echo ">Peak memory was $peak_mem_mb MegaBytes"

#record end
elapsed_time=$(( $(date +%s) - start_time ))
echo "It takes $((elapsed_time / 60 )) minutes"
echo ">ENDING at $(date)"

