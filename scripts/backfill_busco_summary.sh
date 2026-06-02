#!/bin/bash
# One-off migration: bring existing BUSCO summary entries up to the current standard.
#   - taxon-specific (lineage) jsons  -> summary/busco_lineage/<stem>_Lbusco.json
#   - eukaryote jsons                 -> summary/busco_eukaryote/<stem>_Ebusco.json
#   - (re)generate gene + transcript counts -> summary/counts/<stem>_gc.txt / _tc.txt
#
# Sources migrated into summary/busco_lineage (treated as lineage results):
#   - legacy summary/busco/*.json (the single-lineage layout this repo used before)
#   - any raw summary/busco_lineage/*.json not yet suffixed
# summary/busco_eukaryote/*.json is migrated in place.
#
# <stem> is <species_name>_<taxonID>; the taxon id is recovered from
# <species_name>/srr_select.tsv when the existing filename lacks one. The merged annotation
# is looked up at <species_name>/output/files/merged_<sp>_ann.gff (sp = 2nd field of the name),
# matching evaluation.sh. Run from the repo root. JSONs are always renamed; counts are skipped
# (with a warning) when the merged annotation is absent, and computed at most once per species.
#
# Usage: bash scripts/backfill_busco_summary.sh

counts_dir="summary/counts"
mkdir -p "$counts_dir"

#species whose counts were already generated this run (avoid recomputing per lineage)
declare -A counted

#canonical_stem <basename> -> sets STEM=<species_name>_<taxonID> and SPECIES_NAME=<species_name>
canonical_stem() {
	local base="$1" species taxon
	#strip an existing _busco suffix left over from an earlier migration
	base="${base%_busco}"
	if [[ "$base" =~ ^(.+)_([0-9]+)$ ]]; then
		species="${BASH_REMATCH[1]}"
		taxon="${BASH_REMATCH[2]}"
	else
		species="$base"
		taxon=""
	fi
	#recover the taxon id from srr_select.tsv when the filename lacks one
	if [ -z "$taxon" ] && [ -f "$species/srr_select.tsv" ]; then
		taxon=$(cut "$species/srr_select.tsv" -f4|sort|uniq -c|sort -nr|awk '{print $2}'|head -n1)
	fi
	SPECIES_NAME="$species"
	STEM="${species}_${taxon}"
	STEM="${STEM%_}"   #drop trailing _ if no taxon could be found
}

generate_counts() {
	local species="$1" stem="$2"
	[ -n "${counted[$stem]:-}" ] && return
	counted[$stem]=1
	local sp merged gene transcript
	sp=$(echo "$species" | cut -f2 -d"_")
	merged="$species/output/files/merged_${sp}_ann.gff"
	if [ -f "$merged" ]; then
		gene=$(cut -f3 "$merged" | grep -cxF "gene" || true)
		transcript=$(cut -f3 "$merged" | grep -cxE 'transcript|mRNA' || true)
		echo "$gene" > "$counts_dir/${stem}_gc.txt"
		echo "$transcript" > "$counts_dir/${stem}_tc.txt"
		echo "  $stem -> Gene models: $gene | Transcript models: $transcript"
	else
		echo "  $stem -> merged annotation not found at $merged; skipped counts"
	fi
}

#migrate <src_dir> <dest_dir> <suffix>   (suffix is Lbusco or Ebusco)
migrate() {
	local src="$1" dest="$2" suffix="$3"
	[ -d "$src" ] || return
	mkdir -p "$dest"
	shopt -s nullglob
	for json in "$src"/*.json; do
		base=$(basename "$json" .json)
		#skip entries already at the target standard
		[[ "$base" == *_"$suffix" ]] && continue
		canonical_stem "$base"
		mv -v "$json" "$dest/${STEM}_${suffix}.json"
		generate_counts "$SPECIES_NAME" "$STEM"
	done
	shopt -u nullglob
}

#lineage: legacy summary/busco first, then any raw files already in summary/busco_lineage
migrate "summary/busco"          "summary/busco_lineage"   "Lbusco"
migrate "summary/busco_lineage"  "summary/busco_lineage"   "Lbusco"
#eukaryote: in place
migrate "summary/busco_eukaryote" "summary/busco_eukaryote" "Ebusco"

#remove the now-empty legacy dir if everything moved out
rmdir summary/busco 2>/dev/null && echo "Removed empty legacy summary/busco"

echo "Backfill complete."
