#!/bin/bash

echo ">STARTING at $(date)"

species_name="$1"
#if no 2nd argument is given, it uses /no_backup...
busco_db="${2:-/no_backup/rg/references/busco_downloads}"
cpus="${SLURM_CPUS_PER_TASK:-$(nproc)}"

#species shortname
sp=$(echo "$species_name"|cut -f2 -d"_")

#cativate busco conda env
source $(conda info --base)/etc/profile.d/conda.sh
conda activate buscomania

#create storing folders and variables
tmp_files="$species_name/output/files"
res_lineage="$species_name/output/busco_res_lineage"
res_euk="$species_name/output/busco_res_eukaryote"

odb_version="odb12"
busco_lineage_dir="summary/busco_lineage"
busco_euk_dir="summary/busco_eukaryote"

rm -rf "$res_lineage" "$res_euk"
mkdir -p "$res_lineage" "$res_euk" "$busco_lineage_dir" "$busco_euk_dir"

##run busco

#get taxon id form SRR list(get most repeated id in taxon column for specie)
taxonID=$(cut "$species_name/srr_select.tsv" -f4|sort|uniq -c|sort -nr|awk '{print $2}'|head -n1)
echo "TAXON IS: $taxonID"
#get the taxon-specific (custom) lineage
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

#record memory usage
cgroup_dir=$(awk -F: '{print $NF}' /proc/self/cgroup)
peak_mem=`cat /sys/fs/cgroup$cgroup_dir/memory.peak`
peak_mem_mb=$(awk "BEGIN {printf \"%.2f\", $peak_mem / 1048576}") #transfer to mb
echo ">Peak memory was $peak_mem_mb MegaBytes"

#record end
echo ">ENDING at $(date)"
