#!/bin/bash

raw_name="$1"
longread_protists_db="${2:-../data/longread_protists.tsv}"

#snakemake breaks on specieal characters(:) > sanitize folder name, but keep the raw name to match source
species_name=$(printf '%s' "$raw_name" | sed -E 's/[^A-Za-z0-9._-]+/_/g')

sp=$(echo "$raw_name"|cut -f2 -d"_")
sp_extra=$(echo "$raw_name"|cut -f3 -d"_")

#'sp.'/'cf.'/'aff.' are placeholders, not a real epithet (e.g. Schizochytrium_sp._CCTCC_M209059) fall back to the strain so the SRA is specific instead of matching every 'sp' substring.
case "$sp" in
	sp|sp.|cf|cf.|aff|aff.)
	sp="$sp_extra"
	sp_extra=$(echo "$raw_name"|cut -f4 -d"_")
	;;
esac
echo "$sp"

source $(conda info --base)/etc/profile.d/conda.sh
conda activate buscomania

#clone templane for non-human annotation
rm -rf "$species_name"
git clone -v https://github.com/Janek21/LyRic_nonhuman "$species_name"

#select sra for specie
srr_list="$species_name/srr_list.tsv"

#match the organism-name column (field 5) exactly, or as a "<name> <strain>" prefix (avoid gracilis vs neogracilis))
#name first (strain-specific), then (genus + species)
as_words=$(echo "$raw_name" | tr '_' ' ')
binom=$(echo "$as_words" | awk '{print $1, $2}')
echo "Searching for organism '$as_words' (binomial '$binom')"
search_res=$(awk -F'\t' -v q="$as_words" 'BEGIN{q=tolower(q)} {o=tolower($5)} o==q || index(o, q" ")==1' "$longread_protists_db")

if [ -z "$search_res" ]; then
    echo "No exact organism match. Falling back to binomial '$binom'."
    search_res=$(awk -F'\t' -v q="$binom" 'BEGIN{q=tolower(q)} {o=tolower($5)} o==q || index(o, q" ")==1' "$longread_protists_db")
fi

echo "$search_res" > "$species_name/full_srr.tsv"
#select best SRRs
python3 scripts/SRA_selector.py -i "$species_name/full_srr.tsv" -o "$species_name/srr_select.tsv" -s "$srr_list" -e error_species.txt -t 15 -m 18

#if no SRA survive filtering (file empty), clean up and abort species
#[Species, SRA, size] rows logged to error_species.txt by SRA_selector.py
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
python3 scripts/LyRic_setup.py file_transfer -s "$species_name" -i ../data/species/"$raw_name"*/GC*/GC*_genomic.fna* -o "$species_name/data/fasta/$shortname.fa"
#the genome must be a non-empty FASTA, otherwise the pipeline later dies on indexing
genome_fa="$species_name/data/fasta/$shortname.fa"
if [ ! -s "$genome_fa" ] || [ "$(head -c1 "$genome_fa")" != ">" ]; then
	echo "Genome $shortname.fa missing or not a valid FASTA for $species_name; aborting."
	exit 1
fi
#copy the genome annotation (decompresses the gzipped source onto the plain Annotation.gff)
python3 scripts/LyRic_setup.py file_transfer -s "$species_name" -i ../data/species/"$raw_name"*/GC*/"$raw_name"*GC*.gff* -o "$species_name/data/input/Annotation.gff"
#if no annotation was produced, find the closest related species that has one
python3 scripts/annotation_fallback.py -s "$raw_name" -d "$longread_protists_db" -r "../data/species" -o "$species_name/data/input/Annotation.gff"
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
	--job-name="srr_download" \
	--dependency=singleton \
	--output="logs/dw/srr_download_${sp}.%A_%a.out" \
	--error="logs/dw/srr_download_${sp}.%A_%a.err" \
	--array=0-${array_max}%5 \
	scripts/srr_dw.sh "$species_name")
echo "Download array submitted: job $dl_jobid"
#parsable line so lyric_prepare.sh can chain the next stages on this job
echo "DOWNLOAD_JOBID=$dl_jobid"

echo "LyRic is ready to execute"


