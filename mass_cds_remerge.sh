#!/bin/bash
#SBATCH --job-name=lyric_mass_cds_remerge
#SBATCH --qos=normal
#SBATCH --time=480
#SBATCH --mem=16G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --output=logs/eval/merge/%x.%A_%a.out
#SBATCH --error=logs/eval/merge/%x.%A_%a.err
#
# Backfill CDS prediction + remerge for every species already run by LyRic.
#
# Existing species were evaluated before CDS integration, so they only have the
# exon-only prediction (output/files/merged_<sp>_ann.gff) and no CDS-augmented
# annotation. This driver re-runs the canonical evaluation chain on each one, in
# order, within a single array task:
#   1. evaluation.sh        - rebuilds the merge, runs the TD2 ORF prediction, lifts
#                             the ORFs to genome CDS (CDSmerged_<sp>_ann.gff), reruns
#                             BUSCO and rewrites the gene/transcript model counts.
#                             Relocates the CDS-augmented annotation as the new
#                             summary/pred/<sp>_<taxon>_pred.gff.
#   2. merge_evaluation.sh  - merges the (now CDS-augmented) LyRic annotation with the
#                             reference, reruns BUSCO + model counting on the merge, and
#                             relocates summary/merge/pred/<sp>_<taxon>_mergedRef.gff.
#                             Self-skips species that have no own reference annotation.
#
# "Ready" = LyRic (snakemake) produced read models, i.e. the species working dir still
# has output/mappings/mergedReads/*.gff | ont_*.gz (same test as check_runs_status.sh's
# has_run and the input evaluation.sh re-merges from). Species whose working dir was
# deleted (only summary/ artifacts survive) cannot be re-evaluated and are skipped.
#
# Usage (from the repo root):
#   bash mass_cds_remerge.sh [busco_db]
# The first call discovers the ready species and submits a throttled SLURM array;
# each array task then runs evaluation.sh + merge_evaluation.sh on one species.

busco_db="${1:-/no_backup/rg/references/busco_downloads}"
species_list="${2:-}"

#--- bootstrap: not yet inside the array -> discover species and submit the array ---
if [ -z "${SLURM_ARRAY_TASK_ID:-}" ]; then
	mkdir -p logs/eval/merge
	species_list="logs/eval/merge/mass_cds_remerge_species.txt"
	: > "$species_list"

	shopt -s nullglob
	for d in */; do
		d="${d%/}"
		#ready = LyRic read models still present (evaluation.sh re-merges from these)
		lyric_out="$d/output/mappings/mergedReads"
		models=("$lyric_out"/*.gff "$lyric_out"/ont_*.gz)
		[ "${#models[@]}" -gt 0 ] && echo "$d" >> "$species_list"
	done
	shopt -u nullglob

	count=$(wc -l < "$species_list")
	if [ "$count" -eq 0 ]; then
		echo "No species ready to re-evaluate (no */output/mappings/mergedReads/*.gff found)."
		exit 1
	fi
	echo "Found $count species ready to re-evaluate (CDS) + remerge:"
	cat "$species_list"

	#throttle to 10 concurrent tasks to avoid hammering NCBI entrez (get_busco_db / genetic code)
	script_path="$(realpath "$0")"
	sbatch --array=0-$((count - 1))%10 "$script_path" "$busco_db" "$species_list"
	exit 0
fi

#--- array task: re-evaluate (CDS) then remerge the species on this line ---
start_time=$(date +%s)
echo ">STARTING at $(date)"

#let evaluation.sh / merge_evaluation.sh / busco pick up the task allocation
export SLURM_CPUS_PER_TASK

species_name=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$species_list")
if [ -z "$species_name" ]; then
	echo "ERROR: no species at index $SLURM_ARRAY_TASK_ID in $species_list"
	exit 1
fi

echo "=== Stage 1: evaluation.sh (CDS prediction + eval + counts) for $species_name ==="
bash evaluation.sh "$species_name" "$busco_db"

echo "=== Stage 2: merge_evaluation.sh (remerge + eval + counts) for $species_name ==="
bash scripts/merge_evaluation.sh "$species_name" "$busco_db"

elapsed_time=$(( $(date +%s) - start_time ))
echo "It takes $((elapsed_time / 60)) minutes"
echo ">ENDING at $(date)"
