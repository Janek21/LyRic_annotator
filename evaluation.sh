#!/usr/bin/env bash
#SBATCH --job-name=lyric_eval
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=90
#SBATCH --output=logs/eval/%x_%j.out
#SBATCH --error=logs/eval/%x_%j.err

echo ">STARTING at $(date)"

species_name="$1"
busco_db="${2:-/no_backup/rg/references/busco_downloads}"

sp=$(echo "$species_name"|cut -f2 -d"_")
echo "$sp"

#cativate busco conda env
source $(conda info --base)/etc/profile.d/conda.sh
conda activate buscomania

#create storing folders and variables
tmp_files="$species_name/output/files"
lyric_out="$species_name/output/mappings/mergedReads"

rm -rf "$tmp_files"
mkdir -p "$tmp_files"

#per-task AGAT config so parallel jobs don't collide on agat_config.yaml
agat_cfg="$tmp_files/agat_${species_name}_${SLURM_ARRAY_TASK_ID:-$$}.yaml"
agat config --expose --no-log --output "$agat_cfg" >/dev/null 2>&1
trap 'rm -f "$agat_cfg"' EXIT

##decompressions
#decompress gffs
find "$lyric_out" -type f -name "ont_*.gz"|xargs -r -P $(nproc) unpigz -df  #$(nproc) unpigz -df #"$SLURM_CPUS_PER_TASK" unpigz -df
#genome stays gzipped in ../data/species; the uncompressed copy in data/fasta is used instead

#rename for long file names
#Removes the prefix and sufix(minNreads (N varies)) and replaces it with nothing ('')
rename "ont_HpreCap_0+_" "" "$lyric_out"/ont_HpreCap_0+_[DSE]RR*.gff
rename 's/\.HiSS\.tmerge\.min[0-9]+reads\.splicing_status-all\.endSupport-all//' "$lyric_out"/*.gff

#for f in "$lyric_out"/*.gff; do
#	new=$(echo "$f" | sed -E 's/\.HiSS\.tmerge\.min[0-9]+reads\.splicing_status-all\.endSupport-all//')
#	[ "$f" != "$new" ] && mv "$f" "$new"
#done

#detect number of files in folder(if 1 only, dont merge)
shopt -s nullglob
gff_files=("$lyric_out"/*.gff)
shopt -u nullglob
file_count=${#gff_files[@]}

echo "FC is $file_count"

if [ "$file_count" -eq 1 ]; then
	#if only 1 file, rename it for the rest of the pipeline
	cp "${gff_files[0]}" "$tmp_files/merged_${sp}_ann.gff"
	echo "Copied files at $tmp_files/merged_${sp}_ann.gff"
else
	#merge gffs
	agat_sp_merge_annotations.pl --gff "$lyric_out" --config "$agat_cfg" --out "$tmp_files/merged_${sp}_ann.gff"
	echo "Merged files at $tmp_files/merged_${sp}_ann.gff"
fi

#count gene and transcript models in the merged annotation
counts_dir="summary/counts"
mkdir -p "$counts_dir"
merged="$tmp_files/merged_${sp}_ann.gff"
#taxon id = most repeated id in the taxon column (drives the genetic code and BUSCO lineage)
taxonID=$(cut "$species_name/srr_select.tsv" -f4|sort|uniq -c|sort -nr|awk '{print $2}'|head -n1)
gene_count=$(cut -f3 "$merged" | grep -cxF "gene" || true)
transcript_count=$(cut -f3 "$merged" | grep -cxE 'transcript|mRNA' || true)
echo "$gene_count" > "$counts_dir/${species_name}_${taxonID}_gc.txt"
echo "$transcript_count" > "$counts_dir/${species_name}_${taxonID}_tc.txt"
echo "      Gene models: $gene_count | Transcript models: $transcript_count"

#get longest isoform
agat_sp_keep_longest_isoform.pl --gff "$tmp_files/merged_${sp}_ann.gff" --config "$agat_cfg" --out "$tmp_files/longest_${sp}_ann.gff"
echo "Found longest isoforms."

#resolve the NCBI nuclear genetic code for this taxon (codon table for on-standard cases)
echo "KEY: $NCBI_API_KEY"
gcode=$(python3 scripts/get_genetic_code.py -e "ibdyjsayzcllkyvjkc@nespf.com" -k "${NCBI_API_KEY:-}" -t "$taxonID" 2>/dev/null)
if ! [[ "$gcode" =~ ^[0-9]+$ ]]; then
	echo ">Could not resolve genetic code for taxon $taxonID; defaulting to table 1."
	gcode=1
fi
echo "Translation table for $taxonID: $gcode"

#transform to proteins (sequences with premature stops or frameshifts will be translated exactly as your in gff3+computationally better)

#ensure files are generated in particular folders(no naming clash)
td_work="$tmp_files/transdecoder_work"
mkdir -p "$td_work"
transcripts_abs="$(realpath "$tmp_files/transcripts_$sp.fa")"

#generate transcriptome with gffread, using the uncompressed genome copied into the species dir
shortname=$(python3 scripts/LyRic_setup.py shortname -s "$species_name")
gffread "$tmp_files/longest_${sp}_ann.gff" -g "$species_name/data/fasta/$shortname.fa" -w "$transcripts_abs"

(cd "$td_work" && #move to folder for TD2 execution ONLY
	#Find ORFs in transcripts
	TD2.LongOrfs -t "$transcripts_abs" -O . -G "$gcode"
	#Select most probable ORFs to create proteins
	TD2.Predict -t "$transcripts_abs" -O . -G "$gcode" #-O is output of ORFs
)

#move TD2 prot files to correct folders
mv "$td_work/transcripts_$sp.fa.TD2.pep" "$tmp_files/prot_$sp.fa"

echo "TransDecoder proteins in $tmp_files/prot_$sp.fa"

##run busco inline (joined from the former scripts/busco_evaluation.sh so the whole
##evaluation runs as one job; no extra sbatch dependency to chain)
cpus="${SLURM_CPUS_PER_TASK:-$(nproc)}"
res_lineage="$species_name/output/busco_res_lineage"
res_euk="$species_name/output/busco_res_eukaryote"
odb_version="odb12"
busco_lineage_dir="summary/busco_lineage"
busco_euk_dir="summary/busco_eukaryote"

rm -rf "$res_lineage" "$res_euk"
mkdir -p "$res_lineage" "$res_euk" "$busco_lineage_dir" "$busco_euk_dir"

#get the taxon-specific (custom) lineage (taxonID resolved above for the genetic code)
busco_lineage=$(python3 scripts/get_busco_db.py -e "ibdyjsayzcllkyvjkc@nespf.com" -t "$taxonID" -b "$busco_db/file_versions.tsv" -v "$odb_version")
echo "BUSCO lineage for $taxonID is $busco_lineage"
#fixed eukaryote lineage (run for every species alongside the custom one)
euk_lineage="eukaryota_${odb_version}"

#1. taxon-specific lineage busco -> summary/busco_lineage/<stem>_Lbusco.json
busco -m protein -i "$tmp_files/prot_$sp.fa" --download_path "$busco_db" -l "$busco_lineage" -c "$cpus" -f --out_path "${species_name}/output" -o busco_res_lineage --tar
lineage_json="$busco_lineage_dir/${species_name}_${taxonID}_Lbusco.json"
mv "$res_lineage"/*.json "$lineage_json"
ln -sfv "$(realpath "$lineage_json")" "$res_lineage/${species_name}_${taxonID}_Lbusco.json"
busco --plot "$busco_lineage_dir"

#2. eukaryote lineage busco -> summary/busco_eukaryote/<stem>_Ebusco.json
busco -m protein -i "$tmp_files/prot_$sp.fa" --download_path "$busco_db" -l "$euk_lineage" -c "$cpus" -f --out_path "${species_name}/output" -o busco_res_eukaryote --tar
euk_json="$busco_euk_dir/${species_name}_${taxonID}_Ebusco.json"
mv "$res_euk"/*.json "$euk_json"
ln -sfv "$(realpath "$euk_json")" "$res_euk/${species_name}_${taxonID}_Ebusco.json"
busco --plot "$busco_euk_dir"

#predicted (merged) annotation relocated to summary/, hardlinked back to the species location
#(canonical inode lives in summary/pred so the species folder can be removed safely)
pred_dir="summary/pred"
mkdir -p "$pred_dir"
pred_dest="$pred_dir/${species_name}_${taxonID}_pred.gff"
rm -f "$pred_dest"                 #refresh on reruns
mv "$merged" "$pred_dest"         #relocate the prediction into the central summary tree
ln "$pred_dest" "$merged"         #link it back so the original species location stays valid
echo "Predicted annotation collected into $pred_dir/"

rm -rf agat_log_*
echo "Analysis completed!"
echo ">ENDING at $(date)"
