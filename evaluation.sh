#!/bin/bash

species_name="$1"
busco_db="${2:-/no_backup/rg/references/busco_downloads}"

sp=$(echo "$species_name"|cut -f2 -d"_")
echo "$sp"

#cativate busco conda env
source $(conda info --base)/etc/profile.d/conda.sh
conda activate buscomania

#create storing folders
mkdir -p "$species_name/output/files" 


tmp_files="$species_name/output/files"
lyric_out="$species_name/output/mappings/mergedReads"

##decompressions
#decompress gffs
find "$lyric_out" -type f -name "ont_*.gz"|xargs -r -P $(nproc) unpigz -df #$(nproc) unpigz -df #"$SLURM_CPUS_PER_TASK" unpigz -df
#decompress fna if they are compressed still
find ../data/species/"$species_name"*/GCA* -type f -name "GCA*_genomic.fna.gz"|xargs -r -P $(nproc) unpigz -df #$(nproc) unpigz -df #"$SLURM_CPUS_PER_TASK" unpigz -df

#rename for long file names
#Removes the prefix and sufix and replaces it with nothing ('')
rename.ul "ont_HpreCap_0+_" "" "$lyric_out"/ont_HpreCap_0+_[DSE]RR*.gff
rename.ul ".HiSS.tmerge.min2reads.splicing_status-all.endSupport-all" "" "$lyric_out"/*.gff

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
	agat_sp_merge_annotations.pl --gff "$lyric_out" --out "$tmp_files/merged_${sp}_ann.gff"
	echo "Merged files at $tmp_files/merged_${sp}_ann.gff"
fi

#get longest isoform
agat_sp_keep_longest_isoform.pl --gff "$tmp_files/merged_${sp}_ann.gff" --out "$tmp_files/longest_${sp}_ann.gff"
echo "Found longest isoforms."

#transform to transcripts
gffread "$tmp_files/longest_${sp}_ann.gff" -g ../data/species/"$species_name"*/GCA*/GCA*.fna -w "$tmp_files/trsc_$sp.fa"
echo "Transcript files at $tmp_files/trsc_$sp.fa"

mkdir -p logs

##run busco
sbatch \
	--job-name="busco_${sp}" \
	--cpus-per-task=4 \
	--mem=16G \
	--output="logs/%x_%j.out" \
	--error="logs/%x_%j.err" \
	--time=90 \
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


