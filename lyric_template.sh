#!/bin/bash

species_name="$1"
longread_protists_db="${2:-../data/longread_protists.tsv}"

sp=$(echo "$species_name"|cut -f2 -d"_")
sp_extra=$(echo "$species_name"|cut -f3 -d"_")
echo "$sp"

source $(conda info --base)/etc/profile.d/conda.sh
conda activate buscomania

#clone templane for non-human annotation
git clone -v https://github.com/Janek21/LyRic_nonhuman "$species_name"

#select sra for specie
srr_list="$species_name/srr_list.tsv"

#attempt maximum specificity (term 2+3 of name)
echo "Searching for $sp+$sp_extra"
search_res=$(grep -i "$sp" "$longread_protists_db" | grep -i "$sp_extra")

if [ -z "$search_res" ]; then
    echo "No match found for both terms. Falling back to: $sp"
    search_res=$(grep -i "$sp" "$longread_protists_db")
fi

echo "$search_res" > "$species_name/full_srr.tsv"
#select best SRRs
python3 scripts/SRA_selector.py -i "$species_name/full_srr.tsv" -o "$species_name/srr_select.tsv" -s "$srr_list" -e error_species.txt -t 15 -m 8

#if nothing survived the filtering (file missing or empty), clean up and abort this species
#[Species, SRA, size] rows are logged to error_species.txt by SRA_selector.py
srr_count=$(wc -l < "$species_name/srr_select.tsv" 2>/dev/null || echo 0)
if [ ! -s "$species_name/srr_select.tsv" ]; then
	echo "No SRA selected for $species_name; logged to error_species.txt. Removing $species_name and aborting."
	rm -rf "$species_name"
	exit 1
fi
echo "Selected SRA are $srr_count"


#Modify git for the current species and samples
shortname=$(python3 scripts/LyRic_setup.py shortname -s "$species_name")

#config.default.yaml, per the species name
python3 scripts/LyRic_setup.py config -s "$species_name" -o "$species_name/config/default.yaml"
#decompress the genome sequence into the working dir (sources in ../data/species stay gzipped)
python3 scripts/LyRic_setup.py file_transfer -s "$species_name" -i ../data/species/"$species_name"*/GC*/GC*_genomic.fna.gz -o "$species_name/data/fasta/$shortname.fa"
#the genome must be a non-empty FASTA, otherwise the pipeline later dies on indexing
genome_fa="$species_name/data/fasta/$shortname.fa"
if [ ! -s "$genome_fa" ] || [ "$(head -c1 "$genome_fa")" != ">" ]; then
	echo "Genome $shortname.fa missing or not a valid FASTA for $species_name; aborting."
	exit 1
fi
#copy the genome annotation (decompresses the gzipped source onto the plain Annotation.gff)
python3 scripts/LyRic_setup.py file_transfer -s "$species_name" -i ../data/species/"$species_name"*/GC*/"$species_name"*GC*.gff.gz -o "$species_name/data/input/Annotation.gff"
#if no annotation was produced, find the closest related species that has one
python3 scripts/annotation_fallback.py -s "$species_name" -d "$longread_protists_db" -r "../data/species" -o "$species_name/data/input/Annotation.gff"
#set up the sample annotations
python3 scripts/LyRic_setup.py annotate_config -s "$species_name" -i "$srr_list" -o "$species_name/data/sample_annotations.tsv"

#generate empty files
mkdir -p "$species_name/data/input"
touch "$species_name/data/input/fakeCAGE.bed" "$species_name/data/input/fakeDHS.bed"

mkdir -p "$species_name/data/fastq"
#cp scripts/srr_dw.sh $species_name

#SRR downloader
srr_count=$(wc -l < "$srr_list")
array_max=$((srr_count - 1))

dl_jobid=$(sbatch --parsable \
	--job-name="srr_download_${sp}" \
	--output="logs/dw/%x_%A.%a.out" \
	--error="logs/dw/%x_%A.%a.err" \
	--array=0-${array_max} \
	scripts/srr_dw.sh "$species_name")
echo "Download array submitted: job $dl_jobid"
#parsable line so lyric_prepare.sh can chain the next stages on this job
echo "DOWNLOAD_JOBID=$dl_jobid"

echo "LyRic is ready to execute"


