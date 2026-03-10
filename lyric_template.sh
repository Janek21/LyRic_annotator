#!/bin/bash

species_name="$1"
longread_protists_db="$2"

sp=$(echo $species_name|cut -f2 -d"_")
echo $sp

#clone templane for non-human annotation
git clone https://github.com/Janek21/LyRic_nonhuman $species_name

#select sra for specie
srr_list="$species_name/srr_list.tsv"
grep "$sp" "$longread_protists_db" > "$species_name/full_srr.tsv"
#select best SRRs
python3 scripts/SRA_selector.py -i "$species_name/full_srr.tsv" -o "$species_name/srr_select.tsv" -s "$species_name/srr_list.tsv" -t 2
echo "SRA are $(wc -l $species_name/srr_select.tsv) at $species_name/srr_select.tsv"

#Modify git for the current species and samples
shortname=$(python3 scripts/LyRic_setup.py shortname -s $species_name)
#config.default.yaml, per the species name
python3 scripts/LyRic_setup.py config -s $species_name -o $species_name/config/default.yaml
#copy and compress the genome sequence
python3 scripts/LyRic_setup.py file_transfer -s $species_name -i ../data/species/$species_name/raw_*_gn.fa -o $species_name/data/fasta/$shortname.fa.gz
#copy the genome annotation
python3 scripts/LyRic_setup.py file_transfer -s $species_name -i ../data/species/$species_name/raw_*_ann.gff -o $species_name/data/input/Annotation.gff
#set up the sample annotations
python3 scripts/LyRic_setup.py annotate_config -s $species_name -i $srr_list -o $species_name/data/sample_annotations.tsv

#generate empty files
mkdir -p $species_name/data/input
touch $species_name/data/input/fakeCAGE.bed $species_name/data/input/fakeDHS.bed

#sbatch srr_dw.sh $species_name

echo "LyRic is ready to execute"


