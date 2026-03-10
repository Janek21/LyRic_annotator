#!/bin/bash

species_name="$1"

sp=$(echo $species_name|cut -f2 -d"_")
echo $sp

#create storing folders
mkdir -p $species_name/output/busco_res
mkdir -p $species_name/output/files 
mkdir -p busco_summary

tmp_files="$species_name/output/files"
lyric_out="$species_name/output/mappings/mergedReads"

#merge gffs #@.gz?
agat_sp_merge_annotations.pl --gff $lyric_out/ont_*.gff --out $tmp_files/merged_${sp}.gff
echo "Merged files at $tmp_files/merged_$sp.gff"

#get longest isoform
agat_sp_keep_longest_isoform.pl --gff $tmp_files/merged_$sp.gff --out $tmp_files/longest_${sp}.gff
echo "Found longest isoforms."

#transform to transcripts
gffread $tmp_files/longestFP_${sp}_ann.gff -g raw_Pvivax_gn.fa -w $tmp_files/trsc_$sp.fa
echo "Transcript files at $tmp_files/trsc_$sp.fa"

##run busco

#get taxon id form SRR list(get most repeated id in taxon column for specie)
taxonID=$(cut $species_name/srr_select.tsv -f4|sort|uniq -c|sort -nr|awk "{print $2}"|head -n1)
#get lineage
busco_lineage=$(python3 get_busco_db.py -e ibdyjsayzcllkyvjkc@nespf.com -t $taxonID -b ../../data/busco_downloads/file_versions.tsv -v odb12)
echo "BUSCO lineage for $taxonID is $busco_lineage"

#Run busco
busco -m transcriptome -i $tmp_files/trsc_$sp.fa --download_path ../../data/busco_downloads/ -l $busco_lineage -c $(($(nproc))) -f --out_path $species_name/output/busco_res -o $sp --tar

#summary for all
ln -vf $species_name/output/busco_res/*$sp*json busco_summary
busco --plots busco_summary



