#!/bin/bash

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
agat_cfg="$tmp_files/agat_${sp}_${SLURM_ARRAY_TASK_ID:-$$}.yaml"
agat config --expose --output "$agat_cfg" >/dev/null 2>&1
trap 'rm -f "$agat_cfg"' EXIT

##decompressions
#decompress gffs
find "$lyric_out" -type f -name "ont_*.gz"|xargs -r -P $(nproc) unpigz -df  #$(nproc) unpigz -df #"$SLURM_CPUS_PER_TASK" unpigz -df
#decompress fna if they are compressed still
find ../data/species/"$species_name"*/GC* -type f -name "GC*_genomic.fna.gz"|xargs -r -P $(nproc) unpigz -df  #$(nproc) unpigz -df #"$SLURM_CPUS_PER_TASK" unpigz -df

#rename for long file names
#Removes the prefix and sufix and replaces it with nothing ('')
rename "ont_HpreCap_0+_" "" "$lyric_out"/ont_HpreCap_0+_[DSE]RR*.gff
rename ".HiSS.tmerge.min2reads.splicing_status-all.endSupport-all" "" "$lyric_out"/*.gff

#detect number of files in folder(if 1 only ,dont merge)
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
#taxon id = most repeated id in the taxon column (same as busco_evaluation.sh)
taxonID=$(cut "$species_name/srr_select.tsv" -f4|sort|uniq -c|sort -nr|awk '{print $2}'|head -n1)
gene_count=$(cut -f3 "$merged" | grep -cxF "gene" || true)
transcript_count=$(cut -f3 "$merged" | grep -cxE 'transcript|mRNA' || true)
echo "$gene_count" > "$counts_dir/${species_name}_${taxonID}_gc.txt"
echo "$transcript_count" > "$counts_dir/${species_name}_${taxonID}_tc.txt"
echo "      Gene models: $gene_count | Transcript models: $transcript_count"

#get longest isoform
agat_sp_keep_longest_isoform.pl --gff "$tmp_files/merged_${sp}_ann.gff" --config "$agat_cfg" --out "$tmp_files/longest_${sp}_ann.gff"
echo "Found longest isoforms."

#resolve the NCBI nuclear genetic code for this taxon. gffread only extracts
#nucleotide transcripts (table-independent); the translation table matters at
#the ORF step, so it is passed to TransDecoder. Non-standard codes (e.g.
#ciliates like Paramecium/Tetrahymena use table 6, where TAA/TAG code for Gln,
#not stop) would otherwise be mistranslated under the default table 1.
gcode=$(python3 scripts/get_genetic_code.py -e "ibdyjsayzcllkyvjkc@nespf.com" -k "${NCBI_API_KEY:-}" -t "$taxonID" 2>/dev/null)
if ! [[ "$gcode" =~ ^[0-9]+$ ]]; then
	echo ">Could not resolve genetic code for taxon $taxonID; defaulting to table 1."
	gcode=1
fi
echo "Translation table for $taxonID: $gcode"

#i# transform to proteins (sequences with premature stops or frameshifts will be translated exactly as your in gff3+computationally better)
#generate transcriptome with gffread
gffread "$tmp_files/longest_${sp}_ann.gff" -g ../data/species/"$species_name"*/GC*/GC*.fna -w "$tmp_files/transcripts_$sp.fa"

#Find ORFs in transcripts
TD2.LongOrfs -t "$tmp_files/transcripts_$sp.fa" -O "$tmp_files/transdecoder_work" -G "$gcode"

#Select most probable ORFs to create proteins
TD2.Predict -t "$tmp_files/transcripts_$sp.fa" -O "$tmp_files/transdecoder_work" -G "$gcode" #-O is output of ORFs

#move TD2 files to correct folders(as prot and to log)
mv "./transcripts_$sp.fa.TD2.pep" "$tmp_files/prot_$sp.fa"
mv ./*.fa.TD2.* "$tmp_files/transdecoder_work"

echo "TransDecoder proteins in $tmp_files/prot_$sp.fa"

mkdir -p logs

##in Case of massExecution
#bash scripts/busco_evaluation.sh "$species_name"
#exit

##run busco
sbatch \
	--job-name="busco_${sp}" \
	--cpus-per-task=4 \
	--mem=16G \
	--output="logs/%x_%j.out" \
	--error="logs/%x_%j.err" \
	--time=60 \
	scripts/busco_evaluation.sh "$species_name" "$busco_db"

##run gffcompare
sbatch \
	--job-name="gffcmp_${sp}" \
	--cpus-per-task=2 \
	--mem=4G \
	--output="logs/%x_%j.out" \
	--error="logs/%x_%j.err" \
	--time=10 \
	scripts/gffcompare_evaluation.sh "$species_name"

rm -rf agat_log_*
echo "Analysis completed!"
