#!/bin/bash
#SBATCH --job-name=lyric_mass_pred_merge
#SBATCH --qos=normal
#SBATCH --time=30
#SBATCH --mem=4G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/eval/merge/%x_%j.out
#SBATCH --error=logs/eval/merge/%x_%j.err
#
# Backfill the prediction relocation for already-annotated species, then merge them all.
#
# Two jobs, in order:
#   1. move + link  - for every species that evaluation.sh already annotated, relocate the
#                     prediction output/files/merged_<sp>_ann.gff into summary/pred/ and hardlink
#                     it back. New evaluation.sh runs already do this; this backfills the species
#                     that were annotated before that step existed. Idempotent (skips ones already
#                     relocated), so it is safe to re-run.
#   2. mass merge   - hand off to mass_merge.sh, which submits the throttled SLURM array that runs
#                     merge_evaluation.sh on each ready species (that step relocates its own
#                     mergedRef annotation into summary/merge/pred/).
#
# "Ready" = evaluation.sh produced output/files/merged_<sp>_ann.gff (same test as mass_merge.sh).
#
# Usage (from the repo root):
#   sbatch mass_pred_merge.sh [busco_db]

set -uo pipefail

busco_db="${1:-/no_backup/rg/references/busco_downloads}"

mkdir -p logs/eval/merge

pred_dir="summary/pred"
mkdir -p "$pred_dir"

#--- 1. backfill: relocate each species' prediction into summary/pred and link back ---
relocated=0
skipped=0
missing_taxon=0

shopt -s nullglob
for d in */; do
	d="${d%/}"
	#sp derived exactly as evaluation.sh names the merged file (plain 2nd token)
	sp=$(echo "$d" | cut -f2 -d"_")
	merged="$d/output/files/merged_${sp}_ann.gff"
	[ -s "$merged" ] || continue   #not annotated yet, nothing to relocate

	#taxon id = most repeated id in the taxon column (same derivation as evaluation.sh)
	taxonID=""
	[ -s "$d/srr_select.tsv" ] && taxonID=$(cut "$d/srr_select.tsv" -f4 | sort | uniq -c | sort -nr | awk '{print $2}' | head -n1)
	if [ -z "$taxonID" ]; then
		echo "  ! $d: cannot resolve taxonID (missing/empty srr_select.tsv); skipping relocation"
		missing_taxon=$((missing_taxon + 1))
		continue
	fi

	pred_dest="$pred_dir/${d}_${taxonID}_pred.gff"
	#already relocated (species file is the same inode as the summary copy) -> nothing to do
	if [ "$merged" -ef "$pred_dest" ]; then
		skipped=$((skipped + 1))
		continue
	fi

	rm -f "$pred_dest"             #refresh on reruns
	mv "$merged" "$pred_dest"      #relocate the prediction into the central summary tree
	ln "$pred_dest" "$merged"      #link it back so the original species location stays valid
	echo "  + $d -> $pred_dest"
	relocated=$((relocated + 1))
done
shopt -u nullglob

echo "Prediction relocation: $relocated moved, $skipped already done, $missing_taxon skipped (no taxonID)."

#--- 2. merge every ready species (mass_merge.sh submits the throttled SLURM array) ---
echo "Submitting mass merge..."
bash mass_merge.sh "$busco_db"
