#!/bin/bash

species_name="$1"
longread_protists_db="${2:-../data/longread_protists.tsv}"

sp=$(echo "$species_name"|cut -f2 -d"_")
echo "$sp"

source $(conda info --base)/etc/profile.d/conda.sh
conda activate buscomania

#clone templane for non-human annotation
git clone -v https://github.com/Janek21/LyRic_nonhuman "$species_name"

#select sra for specie
srr_list="$species_name/srr_list.tsv"
grep "$sp" "$longread_protists_db" > "$species_name/full_srr.tsv"
#select best SRRs
python3 scripts/SRA_selector.py -i "$species_name/full_srr.tsv" -o "$species_name/srr_select.tsv" -s "$srr_list" -t 10 -m 6
echo "SRA are $(wc -l $species_name/srr_select.tsv) at $species_name/srr_select.tsv"

#Modify git for the current species and samples
shortname=$(python3 scripts/LyRic_setup.py shortname -s "$species_name")

#decompress files if any remain compressed
find ../data/species/"$species_name"*/GCA* -type f -name "*GCA*.gz"|xargs -r -P $(nproc) unpigz -df
echo "Decompressed genome and reference anotation files."

#config.default.yaml, per the species name
python3 scripts/LyRic_setup.py config -s "$species_name" -o "$species_name/config/default.yaml"
#copy and compress the genome sequence
python3 scripts/LyRic_setup.py file_transfer -s "$species_name" -i ../data/species/"$species_name"*/GCA*/GCA*_genomic.fna -o "$species_name/data/fasta/$shortname.fa.gz"
#copy the genome annotation
python3 scripts/LyRic_setup.py file_transfer -s "$species_name" -i ../data/species/"$species_name"*/GCA*/"$species_name"*GCA*.gff -o "$species_name/data/input/Annotation.gff"
#set up the sample annotations
python3 scripts/LyRic_setup.py annotate_config -s "$species_name" -i "$srr_list" -o "$species_name/data/sample_annotations.tsv"

#generate empty files
mkdir -p "$species_name/data/input"
touch "$species_name/data/input/fakeCAGE.bed" "$species_name/data/input/fakeDHS.bed"

cp scripts/runner.sh "$species_name/runner.sh"

mkdir -p "$species_name/data/fastq"
#cp scripts/srr_dw.sh $species_name

#SRR downloader
srr_count=$(wc -l < "$srr_list")
array_max=$((srr_count - 1))

sbatch \
	--job-name="srr_download_${sp}" \
	--output="logs/%x.%A_%a.out" \
	--error="logs/%x.%A_%a.err" \
	--array=0-${array_max} \
	scripts/srr_dw.sh "$species_name"

echo "LyRic is ready to execute"


