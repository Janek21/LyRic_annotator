#!/usr/bin/env bash
#SBATCH --job-name=lyric_merged_eval
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=90
#SBATCH --output=logs/eval/%x_%j.out
#SBATCH --error=logs/eval/%x_%j.err
# Merge the LyRic annotation with the reference annotation and evaluate the result.
# Runs only for species that have a reference annotation (data/input/Annotation.gff).
# Evaluation = BUSCO protein completeness + gene/transcript counts (no gffcompare).
#
# Usage: bash scripts/merge_evaluation.sh <species_name> [busco_db]

echo ">STARTING at $(date)"

species_name="$1"
busco_db="${2:-/no_backup/rg/references/busco_downloads}"
data_root="${3:-../data/species}"   #where the per-species source assemblies/annotations live
cpus="${SLURM_CPUS_PER_TASK:-$(nproc)}"

sp=$(echo "$species_name"|cut -f2 -d"_")
echo "$sp"

#activate the shared conda env (agat, gffread, TD2, busco)
source $(conda info --base)/etc/profile.d/conda.sh
conda activate buscomania

tmp_files="$species_name/output/files"
mkdir -p "$tmp_files"

#reference annotation placed by lyric_template.sh in data/input/Annotation.gff. check if it can be used by lookig if gff exists on ../data/species
shopt -s nullglob
own_ref_src=("$data_root/${species_name}"*/GC*/"${species_name}"*GC*.gff.gz)
shopt -u nullglob

ref_gff="$species_name/data/input/Annotation.gff"
#full LyRic annotation produced by evaluation.sh
lyric_gff="$tmp_files/merged_${sp}_ann.gff"

#only species with their OWN reference annotation are merged; the rest are logged and skipped
merge_summary_dir="summary/merge"
mkdir -p "$merge_summary_dir"
if [ "${#own_ref_src[@]}" -eq 0 ] || [ ! -s "${own_ref_src[0]}" ]; then
	echo "No own reference annotation for $species_name under $data_root (only a fallback or none); skipping merge."
	#record once so reruns don't pile up duplicate lines
	grep -qxF "$species_name" "$merge_summary_dir/no_reference.txt" 2>/dev/null \
		|| echo "$species_name" >> "$merge_summary_dir/no_reference.txt"
	exit 0
fi
if [ ! -s "$ref_gff" ]; then
	echo "Own reference source exists but $ref_gff is missing; re-run lyric_template.sh setup. Aborting."
	exit 1
fi
if [ ! -s "$lyric_gff" ]; then
	echo "LyRic annotation $lyric_gff missing; run evaluation.sh first. Aborting."
	exit 1
fi

#per-task AGAT config so parallel jobs don't collide on agat_config.yaml
agat_cfg="$tmp_files/agat_merge_${species_name}_${SLURM_ARRAY_TASK_ID:-$$}.yaml"
agat config --expose --no-log --output "$agat_cfg" >/dev/null 2>&1
trap 'rm -f "$agat_cfg"' EXIT

#merge LyRic + reference into one non-redundant annotation
merged_ref="$tmp_files/mergedRef_${sp}_ann.gff"
agat_sp_merge_annotations.pl --gff "$lyric_gff" --gff "$ref_gff" --config "$agat_cfg" --out "$merged_ref"
echo "Merged LyRic + reference at $merged_ref"

#count gene and transcript models in the merged annotation
counts_dir="$merge_summary_dir/counts"
mkdir -p "$counts_dir"
#taxon id = most repeated id in the taxon column (same as evaluation.sh)
taxonID=$(cut "$species_name/srr_select.tsv" -f4|sort|uniq -c|sort -nr|awk '{print $2}'|head -n1)
gene_count=$(cut -f3 "$merged_ref" | grep -cxF "gene" || true)
transcript_count=$(cut -f3 "$merged_ref" | grep -cxE 'transcript|mRNA' || true)
echo "$gene_count" > "$counts_dir/${species_name}_${taxonID}_gc.txt"
echo "$transcript_count" > "$counts_dir/${species_name}_${taxonID}_tc.txt"
echo "      Gene models: $gene_count | Transcript models: $transcript_count"

#keep one isoform per gene before translating to proteins
agat_sp_keep_longest_isoform.pl --gff "$merged_ref" --config "$agat_cfg" --out "$tmp_files/longestRef_${sp}_ann.gff"
echo "Found longest isoforms."

#resolve the NCBI nuclear genetic code for this taxon (codon table for non-standard cases)
gcode=$(python3 scripts/get_genetic_code.py -e "ibdyjsayzcllkyvjkc@nespf.com" -k "${NCBI_API_KEY:-}" -t "$taxonID" 2>/dev/null)
if ! [[ "$gcode" =~ ^[0-9]+$ ]]; then
	echo ">Could not resolve genetic code for taxon $taxonID; defaulting to table 1."
	gcode=1
fi
echo "Translation table for $taxonID: $gcode"

#generate the transcriptome with gffread, using the uncompressed genome copied into the species dir
shortname=$(python3 scripts/LyRic_setup.py shortname -s "$species_name")
td_work="$tmp_files/transdecoder_merge_work"
mkdir -p "$td_work"
transcripts_abs="$(realpath "$tmp_files/transcriptsRef_$sp.fa")"
gffread "$tmp_files/longestRef_${sp}_ann.gff" -g "$species_name/data/fasta/$shortname.fa" -w "$transcripts_abs"

(cd "$td_work" && #move to folder for TD2 execution ONLY
	#find ORFs in transcripts
	TD2.LongOrfs -t "$transcripts_abs" -O . -G "$gcode"
	#select most probable ORFs to create proteins
	TD2.Predict -t "$transcripts_abs" -O . -G "$gcode"
)
prot_file="$tmp_files/protRef_$sp.fa"
mv "$td_work/transcriptsRef_$sp.fa.TD2.pep" "$prot_file"
echo "TransDecoder proteins in $prot_file"

##run busco on the merged proteome (taxon-specific + eukaryote lineages)
res_lineage="$species_name/output/busco_mergedRef_lineage"
res_euk="$species_name/output/busco_mergedRef_eukaryote"
odb_version="odb12"
busco_lineage_dir="$merge_summary_dir/busco_lineage"
busco_euk_dir="$merge_summary_dir/busco_eukaryote"

rm -rf "$res_lineage" "$res_euk"
mkdir -p "$res_lineage" "$res_euk" "$busco_lineage_dir" "$busco_euk_dir"

#taxon-specific (custom) lineage
busco_lineage=$(python3 scripts/get_busco_db.py -e "ibdyjsayzcllkyvjkc@nespf.com" -t "$taxonID" -b "$busco_db/file_versions.tsv" -v "$odb_version")
echo "BUSCO lineage for $taxonID is $busco_lineage"
euk_lineage="eukaryota_${odb_version}"

#1. taxon-specific lineage busco -> summary/merge/busco_lineage/<stem>_Lbusco.json
busco -m protein -i "$prot_file" --download_path "$busco_db" -l "$busco_lineage" -c "$cpus" -f --out_path "${species_name}/output" -o busco_mergedRef_lineage --tar
lineage_json="$busco_lineage_dir/${species_name}_${taxonID}_Lbusco.json"
mv "$res_lineage"/*.json "$lineage_json"

#2. eukaryote lineage busco -> summary/merge/busco_eukaryote/<stem>_Ebusco.json
busco -m protein -i "$prot_file" --download_path "$busco_db" -l "$euk_lineage" -c "$cpus" -f --out_path "${species_name}/output" -o busco_mergedRef_eukaryote --tar
euk_json="$busco_euk_dir/${species_name}_${taxonID}_Ebusco.json"
mv "$res_euk"/*.json "$euk_json"

#merged-with-reference annotation relocated to summary/, hardlinked back to the species location
#(canonical inode lives in summary/merge/pred so the species folder can be removed safely)
pred_dir="$merge_summary_dir/pred"
mkdir -p "$pred_dir"
pred_dest="$pred_dir/${species_name}_${taxonID}_mergedRef.gff"
rm -f "$pred_dest"                 #refresh on reruns
mv "$merged_ref" "$pred_dest"     #relocate the merged annotation into the central summary tree
ln "$pred_dest" "$merged_ref"     #link it back so the original species location stays valid
echo "Merged-with-reference annotation collected into $pred_dir/"

rm -rf agat_log_*
echo "Analysis completed!"
echo ">ENDING at $(date)"
