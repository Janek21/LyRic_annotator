#!/bin/bash

echo ">STARTING at $(date)"

species_name="$1"

#activate busco conda env
source $(conda info --base)/etc/profile.d/conda.sh
conda activate buscomania

sp=$(echo "$species_name"|cut -f2 -d"_")
echo "$sp"

#create storing folders and variables
gffcmp_dir="$species_name/output/gffcmp"
log_dir="$species_name/output/compare_logs"

mkdir -p "$gffcmp_dir"
mkdir -p "$log_dir"
mkdir -p gffcmp_summary

#reference annotation
ref_gff=$(realpath "../data/species/$species_name"*/GCA*/*GCA*.gff)
#busco_evaluation cleaned annotation
pred_gff="$species_name/output/files/longest_${sp}_ann.gff"

echo "Running gffcompare for $species_name"
echo "Reference at: $ref_gff"
echo "Predicted at: $pred_gff"

#all will be computed at logs folder
prefix="$log_dir/${species_name}-Lycmp"

gffcompare -r "$ref_gff" "$pred_gff" -o "$prefix"

#move stats files to summary folders
shopt -s extglob
mv "$log_dir"/*.stats "$gffcmp_dir"/
ln -vf "$gffcmp_dir"/*.stats gffcmp_summary

echo "done_${species_name}"



#record memory usage
cgroup_dir=$(awk -F: '{print $NF}' /proc/self/cgroup)
peak_mem=`cat /sys/fs/cgroup$cgroup_dir/memory.peak`
peak_mem_mb=$(awk "BEGIN {printf \"%.2f\", $peak_mem / 1048576}") #transfer to mb
echo ">Peak memory was $peak_mem_mb MegaBytes"

#record end
echo ">ENDING at $(date)"
