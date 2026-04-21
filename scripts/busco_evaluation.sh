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
res_folder="$species_name/output/busco_res"

rm -rf "$res_folder"
mkdir -p "$res_folder"
mkdir -p busco_summary

##run busco

#get taxon id form SRR list(get most repeated id in taxon column for specie)
taxonID=$(cut "$species_name/srr_select.tsv" -f4|sort|uniq -c|sort -nr|awk '{print $2}'|head -n1)
echo "TAXON IS: $taxonID"
#get lineage
busco_lineage=$(python3 scripts/get_busco_db.py -e "ibdyjsayzcllkyvjkc@nespf.com" -t "$taxonID" -b "$busco_db/file_versions.tsv" -v odb12)
echo "BUSCO lineage for $taxonID is $busco_lineage"

#Run busco
busco -m protein -i "$tmp_files/prot_$sp.fa" --download_path "$busco_db" -l "$busco_lineage" -c "$cpus" -f --out_path "${species_name}/output" -o busco_res --tar

#summary for all
mv "$res_folder"/*json "$res_folder/${species_name}_${taxonID}.json"
ln -vf "$res_folder"/*json busco_summary
busco --plot busco_summary

#record memory usage
cgroup_dir=$(awk -F: '{print $NF}' /proc/self/cgroup)
peak_mem=`cat /sys/fs/cgroup$cgroup_dir/memory.peak`
peak_mem_mb=$(awk "BEGIN {printf \"%.2f\", $peak_mem / 1048576}") #transfer to mb
echo ">Peak memory was $peak_mem_mb MegaBytes"

#record end
echo ">ENDING at $(date)"
