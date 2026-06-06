#!/bin/bash
# Step 0 normalization: keep ../data/species sources always compressed.
# Compresses any plain genome (*_genomic.fna) / annotation (*GC*.gff) source in
# place with pigz, skipping files that already have a .gz sibling. Safe to re-run.
#
# Usage: bash scripts/recompress_sources.sh [species_root]

species_root="${1:-../data/species}"

#collect plain sources that lack a .gz counterpart
to_compress=()
while IFS= read -r f; do
	[ -e "$f.gz" ] || to_compress+=("$f")
done < <(find "$species_root"/*/GC* -type f \( -name "*_genomic.fna" -o -name "*GC*.gff" \))

if [ "${#to_compress[@]}" -eq 0 ]; then
	echo "Nothing to compress; all sources already gzipped."
	exit 0
fi

echo "Compressing ${#to_compress[@]} source file(s) with pigz..."
#pigz replaces each file with file.gz; run nproc compressions in parallel
printf '%s\n' "${to_compress[@]}" | xargs -r -P "$(nproc)" pigz -v
echo "Done. Sources in $species_root are now compressed."
