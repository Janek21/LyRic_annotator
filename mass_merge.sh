#!/bin/bash
#SBATCH --job-name=lyric_mass_merge
#SBATCH --qos=normal
#SBATCH --time=150
#SBATCH --mem=16G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --output=logs/eval/merge/%x.%A_%a.out
#SBATCH --error=logs/eval/merge/%x.%A_%a.err
#
# Mass merge + evaluate every species already set up in the repo root.
# Assumes the earlier stages (download, run, evaluation.sh) have finished:
# a species is "ready" when evaluation.sh produced output/files/merged_<sp>_ann.gff.
# merge_evaluation.sh itself skips species that lack a reference annotation.
#
# Usage (from the repo root):
#   bash mass_merge.sh [busco_db]
# The first call discovers the ready species and submits a throttled SLURM array;
# each array task then runs merge_evaluation.sh on one species.

busco_db="${1:-/no_backup/rg/references/busco_downloads}"
species_list="${2:-}"

#--- bootstrap: not yet inside the array -> discover species and submit the array ---
if [ -z "${SLURM_ARRAY_TASK_ID:-}" ]; then
	mkdir -p logs/eval/merge
	species_list="logs/eval/merge/mass_merge_species.txt"
	: > "$species_list"

	shopt -s nullglob
	for d in */; do
		d="${d%/}"
		#sp derived exactly as evaluation.sh names the merged file (plain 2nd token)
		sp=$(echo "$d" | cut -f2 -d"_")
		[ -s "$d/output/files/merged_${sp}_ann.gff" ] && echo "$d" >> "$species_list"
	done
	shopt -u nullglob

	count=$(wc -l < "$species_list")
	if [ "$count" -eq 0 ]; then
		echo "No species ready to merge (no */output/files/merged_*_ann.gff found)."
		exit 1
	fi
	echo "Found $count species ready to merge:"
	cat "$species_list"

	#throttle to 10 concurrent tasks to avoid hammering NCBI entrez (get_busco_db / genetic code)
	script_path="$(realpath "$0")"
	sbatch --array=0-$((count - 1))%10 "$script_path" "$busco_db" "$species_list"
	exit 0
fi

#--- array task: merge + evaluate the species on this line ---
start_time=$(date +%s)
echo ">STARTING at $(date)"

#let merge_evaluation.sh / busco pick up the task allocation
export SLURM_CPUS_PER_TASK

species_name=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$species_list")
if [ -z "$species_name" ]; then
	echo "ERROR: no species at index $SLURM_ARRAY_TASK_ID in $species_list"
	exit 1
fi

echo "Merging + evaluating $species_name"
bash scripts/merge_evaluation.sh "$species_name" "$busco_db"

elapsed_time=$(( $(date +%s) - start_time ))
echo "It takes $((elapsed_time / 60)) minutes"
echo ">ENDING at $(date)"
