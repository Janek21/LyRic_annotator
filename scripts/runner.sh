#!/bin/bash

#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --qos=normal
#SBATCH --job-name=lyric

#SBATCH --mem=42G
#SBATCH --time=500

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8

#record start
echo ">STARTING at $(date)"

module load CMake/3.29.3-GCCcore-12.3.0 
module load Python/3.13.5-GCCcore-14.3.0

source ~/bin/snakemake/bin/activate

snakemake --unlock
snakemake --cores all --configfile config/default.yaml --keep-going

#record memory usage
cgroup_dir=$(awk -F: '{print $NF}' /proc/self/cgroup)
peak_mem=`cat /sys/fs/cgroup$cgroup_dir/memory.peak`
peak_mem_mb=$(awk "BEGIN {printf \"%.2f\", $peak_mem / 1048576}") #transfer to mb
echo ">Peak memory was $peak_mem_mb MegaBytes"

#record end
echo ">ENDING at $(date)"

