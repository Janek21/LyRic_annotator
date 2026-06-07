#!/bin/bash
# Classify every species listed in dataspecie.txt by how far it got in LyRic.
# Buckets (mutually exclusive):
#   NO FOLDER   - species in the list with no working dir (never started)
#   NOT RUN     - folder exists but LyRic (snakemake) produced no read models yet
#   NOT EVAL    - LyRic ran but the annotation was not merged/evaluated
#   DONE        - LyRic ran and the merged annotation exists
#
# Folder names are matched the same way the pipeline creates them: the raw name
# from dataspecie.txt is sanitized with  sed -E 's/[^A-Za-z0-9._-]+/_/g'
# (see lyric_template.sh), so e.g. "Entamoeba_histolytica_HM-1:IMSS" -> "..._HM-1_IMSS".
#
# Usage (from the repo root):
#   bash check_against_dataspecie.sh [dataspecie.txt] [base_dir]
# Defaults:
#   dataspecie.txt = /home/jj/Desktop/Data_science/CRG/TFM2/projects/busco_references/dataspecie.txt
#   base_dir       = .  (where the species working dirs live)

shopt -s nullglob

datafile="${1:-/home/jj/Desktop/Data_science/CRG/TFM2/projects/busco_references/dataspecie.txt}"
base_dir="${2:-.}"

if [ ! -s "$datafile" ]; then
	echo "Species list not found or empty: $datafile" >&2
	exit 1
fi

#--- stage detection (mirrors check_lyric_status.sh) ---

has_run() {  #LyRic/snakemake produced per-read models
	local lyric_out="$1/output/mappings/mergedReads"
	local models=("$lyric_out"/*.gff "$lyric_out"/ont_*.gz)
	[ "${#models[@]}" -gt 0 ]
}

has_eval() {  #merged annotation produced by evaluation.sh
	local merged=("$1"/output/files/merged_*_ann.gff)
	[ "${#merged[@]}" -gt 0 ] && [ -s "${merged[0]}" ]
}

#--- classify ---

no_folder=()
not_run=()
not_eval=()
done_sp=()

while IFS= read -r raw || [ -n "$raw" ]; do
	raw="${raw%$'\r'}"                 # strip stray CR
	[ -z "${raw// }" ] && continue     # skip blank lines
	work_name=$(printf '%s' "$raw" | sed -E 's/[^A-Za-z0-9._-]+/_/g')
	sp_dir="$base_dir/$work_name"

	if [ ! -d "$sp_dir" ]; then
		no_folder+=("$raw")
	elif ! has_run "$sp_dir"; then
		not_run+=("$raw")
	elif ! has_eval "$sp_dir"; then
		not_eval+=("$raw")
	else
		done_sp+=("$raw")
	fi
done < "$datafile"

#--- report ---

print_group() {
	local title="$1"; shift
	printf '\n== %s (%d) ==\n' "$title" "$#"
	if [ "$#" -eq 0 ]; then
		echo "  (none)"
	else
		printf '  %s\n' "$@"
	fi
}

total=$(( ${#no_folder[@]} + ${#not_run[@]} + ${#not_eval[@]} + ${#done_sp[@]} ))
echo "Species list: $datafile"
echo "Base dir:     $base_dir"
echo "Total species: $total"

print_group "NO FOLDER (never started)"        "${no_folder[@]}"
print_group "FOLDER, LyRic NOT RUN"            "${not_run[@]}"
print_group "RUN but NOT EVALUATED"            "${not_eval[@]}"
print_group "DONE (run + evaluated)"           "${done_sp[@]}"
