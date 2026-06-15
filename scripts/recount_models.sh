#!/usr/bin/env bash
#SBATCH --job-name=lyric_recount
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=30
#SBATCH --output=logs/eval/%x_%j.out
#SBATCH --error=logs/eval/%x_%j.err
# Recompute gene/transcript model counts for every already-relocated annotation,
# without rerunning evaluation.sh / merge_evaluation.sh (no AGAT, BUSCO or TD2).
#
# It rewrites the count files that the old col3-only counter populated with 0 on
# tmerge GTF-style output. Two annotation flavours are handled by count_models:
#   - prediction          summary/pred/<sp>_<taxon>_pred.gff       -> summary/counts/
#   - merged-with-reference summary/merge/pred/<sp>_<taxon>_mergedRef.gff -> summary/merge/counts/
#
# Usage: bash scripts/recount_models.sh   (or sbatch scripts/recount_models.sh)

set -euo pipefail
echo ">STARTING at $(date)"

#count gene + transcript models in one gffread pass. --keep-genes normalises any
#input (tmerge GTF or AGAT GFF3) into gene + transcript records: real gene
#features are preserved (AGAT-clustered loci, not the per-transcript gene_id),
#and one gene + transcript is synthesised per id when the input has no such
#feature. Counting the normalised col3 feature types then works for every
#annotation flavour. Prints "<gene_count>\t<transcript_count>".
count_both() {  # $1=gff
	{ gffread "$1" --keep-genes -o - 2>/dev/null || true; } | awk -F'\t' '
		/^#/ { next }
		$3 ~ /^([A-Za-z_]*gene)$/                { g++; next }
		$3 ~ /^(transcript|mRNA|[A-Za-z_]*RNA)$/ { t++ }
		END { printf "%d\t%d\n", g, t }'
}

# recount_dir <gff-glob-dir> <gff-suffix> <counts-dir>
# For each <stem><suffix> GFF, write <stem>_gc.txt and <stem>_tc.txt into counts-dir.
recount_dir() {
	local pred_dir="$1" suffix="$2" counts_dir="$3"
	[ -d "$pred_dir" ] || { echo "  (no $pred_dir, skipping)"; return; }
	mkdir -p "$counts_dir"

	shopt -s nullglob
	local gffs=("$pred_dir"/*"$suffix")
	shopt -u nullglob
	if [ "${#gffs[@]}" -eq 0 ]; then
		echo "  (no *$suffix files in $pred_dir)"
		return
	fi

	local gff stem counts genes transcripts
	for gff in "${gffs[@]}"; do
		[ -s "$gff" ] || { echo "  SKIP empty $gff"; continue; }
		stem=$(basename "$gff" "$suffix")            #<species>_<taxon>
		counts=$(count_both "$gff")
		genes=${counts%%$'\t'*}
		transcripts=${counts##*$'\t'}
		echo "$genes" > "$counts_dir/${stem}_gc.txt"
		echo "$transcripts" > "$counts_dir/${stem}_tc.txt"
		printf '  %-55s genes=%-7s transcripts=%s\n' "$stem" "$genes" "$transcripts"
	done
}

echo "== Prediction annotations =="
recount_dir "summary/pred" "_pred.gff" "summary/counts"

echo "== Merged-with-reference annotations =="
recount_dir "summary/merge/pred" "_mergedRef.gff" "summary/merge/counts"

echo "Recount completed!"
echo ">ENDING at $(date)"
