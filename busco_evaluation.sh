#!/bin/bash

species_name="$1"
#if no 2nd argument is given, it uses /no_backup...
busco_db="${2:-/no_backup/rg/references/busco_downloads}"

sp=$(echo $species_name|cut -f2 -d"_")
echo $sp

#cativate busco conda env
source $(conda info --base)/etc/profile.d/conda.sh
conda activate buscomania

#create storing folders
mkdir -p $species_name/output/busco_res
mkdir -p $species_name/output/files 
mkdir -p busco_summary

tmp_files="$species_name/output/files"
lyric_out="$species_name/output/mappings/mergedReads"

##decompressions
#decompress gffs
find "$lyric_out" -type f -name "ont_*.gz"|xargs -r -P "$SLURM_CPUS_PER_TASK" unpigz -df
#decompress fna if they are compressed still
find ../data/species/"$species_name"/GCA* -type f -name "GCA*_genomic.fna.gz"|xargs -r -P "$SLURM_CPUS_PER_TASK" unpigz -df

#rename for long file names
#Removes the prefix and sufix and replaces it with nothing ('')
rename "ont_HpreCap_0+_" "" $lyric_output/ont_HpreCap_0+_[SE]RR*.gff
rename ".HiSS.tmerge.min2reads.splicing_status-all.endSupport-all" "" $lyric_output/*.gff
#merge gffs
agat_sp_merge_annotations.pl --gff $lyric_out --out $tmp_files/merged_${sp}_ann.gff
echo "Merged files at $tmp_files/merged_$sp.gff"

#get longest isoform
agat_sp_keep_longest_isoform.pl --gff $tmp_files/merged_${sp}_ann.gff --out $tmp_files/longest_${sp}_ann.gff
echo "Found longest isoforms."

#transform to transcripts
gffread $tmp_files/longest_${sp}_ann.gff -g ../data/species/$species_name/GCA*/GCA*_genomic.fna -w $tmp_files/trsc_$sp.fa
echo "Transcript files at $tmp_files/trsc_$sp.fa"

##run busco

#get taxon id form SRR list(get most repeated id in taxon column for specie)
taxonID=$(cut $species_name/srr_select.tsv -f4|sort|uniq -c|sort -nr|awk '{print $2}'|head -n1)
#get lineage
busco_lineage=$(python3 scripts/get_busco_db.py -e ibdyjsayzcllkyvjkc@nespf.com -t $taxonID -b $busco_db/file_versions.tsv -v odb12)
echo "BUSCO lineage for $taxonID is $busco_lineage"

#Run busco
busco -m transcriptome -i $tmp_files/trsc_$sp.fa --download_path $busco_db -l $busco_lineage -c "$SLURM_CPUS_PER_TASK" -f --out_path $species_name/output/busco_res --tar

#summary for all
ln -vf $species_name/output/busco_res/*$sp*json busco_summary
busco --plot busco_summary



