#!/bin/bash
# One-off migration: bring existing busco_summary/ entries up to the current standard.
#   1. rename  busco_summary/<name>.json  ->  busco_summary/<name>_busco.json
#   2. (re)generate gene + transcript model counts from each species' merged annotation:
#        busco_summary/<name>_gc.txt  and  busco_summary/<name>_tc.txt
#
# <name> is the json basename, i.e. <species_name>_<taxonID>. The merged annotation is
# looked up at <species_name>/output/files/merged_<sp>_ann.gff (sp = 2nd field of the name),
# matching evaluation.sh. Run from the repo root. JSONs are always renamed; counts are
# skipped (with a warning) when the species' merged annotation is not present.
#
# Usage: bash scripts/backfill_busco_summary.sh [busco_summary_dir]

busco_summary_dir="${1:-busco_summary}"

shopt -s nullglob
for json in "$busco_summary_dir"/*.json; do
	base=$(basename "$json" .json)

	#skip entries already migrated
	if [[ "$base" == *_busco ]]; then
		continue
	fi

	#split the basename into species_name (+ taxon id, if the name already carries one)
	if [[ "$base" =~ ^(.+)_([0-9]+)$ ]]; then
		species_name="${BASH_REMATCH[1]}"
		taxonID="${BASH_REMATCH[2]}"
	else
		species_name="$base"
		taxonID=""
	fi

	#recover the taxon id from srr_select.tsv when the filename lacks one
	#(same lookup as evaluation.sh: most repeated id in the taxon column)
	if [ -z "$taxonID" ] && [ -f "$species_name/srr_select.tsv" ]; then
		taxonID=$(cut "$species_name/srr_select.tsv" -f4|sort|uniq -c|sort -nr|awk '{print $2}'|head -n1)
	fi

	#canonical stem: <species_name>_<taxonID> (drop the trailing _ if no taxon could be found)
	stem="${species_name}_${taxonID}"
	stem="${stem%_}"

	#1. rename to the <stem>_busco.json standard
	mv -v "$json" "$busco_summary_dir/${stem}_busco.json"

	#2. gene + transcript counts from the merged annotation
	sp=$(echo "$species_name" | cut -f2 -d"_")
	merged="$species_name/output/files/merged_${sp}_ann.gff"

	if [ -f "$merged" ]; then
		gene_count=$(cut -f3 "$merged" | grep -cxF "gene" || true)
		transcript_count=$(cut -f3 "$merged" | grep -cxE 'transcript|mRNA' || true)
		echo "$gene_count" > "$busco_summary_dir/${stem}_gc.txt"
		echo "$transcript_count" > "$busco_summary_dir/${stem}_tc.txt"
		echo "  $stem -> Gene models: $gene_count | Transcript models: $transcript_count"
	else
		echo "  $stem -> merged annotation not found at $merged; renamed json, skipped counts"
	fi
done
shopt -u nullglob

echo "Backfill complete."
