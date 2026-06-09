#!/bin/bash
# Classify every species listed in dataspecie.txt by how far it got in LyRic.
# Buckets (mutually exclusive, in pipeline order):
#   NO FOLDER          - species in the list with no working dir (never started)
#   NOT RUN            - folder exists but LyRic (snakemake) produced no read models yet
#   NOT EVALUATED      - LyRic ran but the prediction was not merged/evaluated
#   EVAL, NOT MERGED   - prediction done, but not yet merged with the reference annotation
#   NO REFERENCE       - prediction done, merge skipped (no reference annotation available)
#   MERGED, NOT EVAL   - merged with the reference, but BUSCO on the merged annotation not done
#   MERGED + EVALUATED - merged with the reference and BUSCO-evaluated (fully complete)
#
# Stage artifacts (mirrors evaluation.sh / merge_evaluation.sh; merged/eval outputs are
# collected under summary/ so they survive deletion of the species working dir):
#   run             $sp_dir/output/mappings/mergedReads/*.gff | ont_*.gz
#   evaluated       $sp_dir/output/files/merged_*_ann.gff      (prediction; hardlinked into summary/pred/)
#   no reference    summary/merge/no_reference.txt             (species listed here)
#   merged          $sp_dir/output/files/mergedRef_*_ann.gff | summary/merge/pred/<sp>_*_mergedRef.gff
#   merged+eval     summary/merge/busco_lineage/<sp>_*_Lbusco.json
#
# Folder names are matched the same way the pipeline creates them: the raw name
# from dataspecie.txt is sanitized with  sed -E 's/[^A-Za-z0-9._-]+/_/g'
# (see lyric_template.sh), so e.g. "Entamoeba_histolytica_HM-1:IMSS" -> "..._HM-1_IMSS".
#
# Usage (from the repo root):
#   bash check_runs_status.sh [dataspecie.txt] [base_dir]
# Defaults:
#   dataspecie.txt = /home/jj/Desktop/Data_science/CRG/TFM2/projects/busco_references/dataspecie.txt
#   base_dir       = .  (where the species working dirs and summary/ live)

shopt -s nullglob

datafile="${1:-dataspecie.txt}"
base_dir="${2:-.}"
summary_dir="$base_dir/summary"

if [ ! -s "$datafile" ]; then
	echo "Species list not found or empty: $datafile" >&2
	exit 1
fi

#--- stage detection (mirrors evaluation.sh / merge_evaluation.sh) ---

has_run() {  #LyRic/snakemake produced per-read models
	local lyric_out="$1/output/mappings/mergedReads"
	local models=("$lyric_out"/*.gff "$lyric_out"/ont_*.gz)
	[ "${#models[@]}" -gt 0 ]
}

has_eval() {  #prediction annotation produced by evaluation.sh
	local merged=("$1"/output/files/merged_*_ann.gff)
	[ "${#merged[@]}" -gt 0 ] && [ -s "${merged[0]}" ]
}

is_no_ref() {  #merge legitimately skipped (no reference annotation); $1 = sanitized work name
	[ -f "$summary_dir/merge/no_reference.txt" ] && grep -qxF "$1" "$summary_dir/merge/no_reference.txt"
}

has_merged() {  #LyRic + reference merged by merge_evaluation.sh; $1 = sp_dir, $2 = work name
	local f=("$1"/output/files/mergedRef_*_ann.gff "$summary_dir"/merge/pred/"$2"_*_mergedRef.gff)
	[ "${#f[@]}" -gt 0 ] && [ -s "${f[0]}" ]
}

has_merged_eval() {  #BUSCO completed on the merged-with-reference annotation; $1 = work name
	local j=("$summary_dir"/merge/busco_lineage/"$1"_*_Lbusco.json)
	[ "${#j[@]}" -gt 0 ] && [ -s "${j[0]}" ]
}

#--- classify ---

no_folder=()
not_run=()
not_eval=()
not_merged=()
no_ref=()
merged_only=()
merged_eval=()

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
	elif is_no_ref "$work_name"; then
		no_ref+=("$raw")
	elif ! has_merged "$sp_dir" "$work_name"; then
		not_merged+=("$raw")
	elif ! has_merged_eval "$work_name"; then
		merged_only+=("$raw")
	else
		merged_eval+=("$raw")
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

total=$(( ${#no_folder[@]} + ${#not_run[@]} + ${#not_eval[@]} + ${#no_ref[@]} \
	+ ${#not_merged[@]} + ${#merged_only[@]} + ${#merged_eval[@]} ))
echo "Species list: $datafile"
echo "Base dir:     $base_dir"
echo "Total species: $total"

print_group "NO FOLDER (never started)"          "${no_folder[@]}"
print_group "FOLDER, LyRic NOT RUN"              "${not_run[@]}"
print_group "RUN but NOT EVALUATED"              "${not_eval[@]}"
print_group "EVALUATED but NOT MERGED"           "${not_merged[@]}"
print_group "NO REFERENCE (merge skipped)"       "${no_ref[@]}"
print_group "MERGED but NOT EVALUATED"           "${merged_only[@]}"
print_group "MERGED + EVALUATED (complete)"      "${merged_eval[@]}"
