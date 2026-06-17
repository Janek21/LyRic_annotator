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

#use the CPUs and memory granted by the SLURM allocation as Snakemake's budget;
#fall back to small defaults when run outside SLURM (e.g. a local/interactive run)
cpus="${SLURM_CPUS_PER_TASK:-2}"
if [ -n "${SLURM_MEM_PER_NODE:-}" ]; then
	mem_mb="$SLURM_MEM_PER_NODE"
elif [ -n "${SLURM_MEM_PER_CPU:-}" ]; then
	mem_mb=$(( SLURM_MEM_PER_CPU * cpus ))
else
	#must stay >= the largest per-rule mem_mb reservation in the profile
	#(longReadMapping reserves 5000), or Snakemake refuses to schedule it
	mem_mb=8000
fi

#run from the repo root, pointing snakemake at the species clone
species_name="$1"
species_dir="$(realpath "$species_name")"
echo "Running LyRic for $species_name"

#drop the whole species if any download is incomplete:
#0-byte *.fastq.gz (no dw)
#*.fastq.gz.part / *.fastq.gz.uniq (killed mid-write)
if find "$species_name/data/fastq" -maxdepth 1 -type f \( -name '*.fastq.gz' -size 0 -o -name '*.fastq.gz.part' -o -name '*.fastq.gz.uniq' \) | grep -q .; then
	echo "Incomplete download in $species_name/data/fastq (empty .fastq.gz or leftover .part/.uniq); removing $species_name."
	rm -rf "$species_name"
	exit 1
fi

snakemake --snakefile "$species_dir/workflow/Snakefile" --directory "$species_dir" --configfile "$species_dir/config/default.yaml" --unlock
snakemake --snakefile "$species_dir/workflow/Snakefile" --directory "$species_dir" --configfile "$species_dir/config/default.yaml" --cores "$cpus" --resources mem_mb="$mem_mb" --keep-going

#record memory usage
cgroup_dir=$(awk -F: '{print $NF}' /proc/self/cgroup)
peak_mem=`cat /sys/fs/cgroup$cgroup_dir/memory.peak`
peak_mem_mb=$(awk "BEGIN {printf \"%.2f\", $peak_mem / 1048576}") #transfer to mb
echo ">Peak memory was $peak_mem_mb MegaBytes"

#record end
elapsed_time=$(( $(date +%s) - start_time ))
echo "It takes $((elapsed_time / 60 )) minutes"
echo ">ENDING at $(date)"

