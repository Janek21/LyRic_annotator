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
gffcmp_summary_dir="summary/gffcmp"

mkdir -p "$gffcmp_dir"
mkdir -p "$log_dir"
mkdir -p "$gffcmp_summary_dir"

#reference annotation (the uncompressed annotation already placed in the species working dir)
ref_gff="$species_name/data/input/Annotation.gff"
#busco_evaluation cleaned annotation
pred_gff="$species_name/output/files/longest_${sp}_ann.gff"

echo "Running gffcompare for $species_name"
echo "Reference at: $ref_gff"
echo "Predicted at: $pred_gff"

#fix gff if needed(replace ? in strand column for .)
questionPresence_ref=$(cut -f7 $ref_gff|grep -Fx "?"|wc -l)

if [ "$questionPresence_ref" -ne 0 ]; then
	echo "Replacing ? strand symbol in reference annotation"
	awk -F'\t' 'BEGIN{OFS="\t"} {$7=gensub(/\?/, ".", "g", $7); print}' "$ref_gff" > "$species_name/output/files/newRef_${species_name}.gff"
	#replace reference variable
	ref_gff="$species_name/output/files/newRef_${species_name}.gff"
fi

questionPresence_pred=$(cut -f7 $pred_gff|grep -Fx "?"|wc -l)

if [ "$questionPresence_pred" -ne 0 ]; then
	echo "Replacing ? strand in predicted annotation"
        awk -F'\t' 'BEGIN{OFS="\t"} {$7=gensub(/\?/, ".", "g", $7); print}' "$pred_gff" > "$species_name/output/files/newLongest_${sp}_ann.gff"
        #replace pred variable
	pred_gff="$species_name/output/files/newLongest_${sp}_ann.gff"
fi


#all will be computed at logs folder
prefix="$log_dir/${species_name}-Lycmp"

gffcompare -r "$ref_gff" "$pred_gff" -o "$prefix"

#move stats files to summary folders
shopt -s extglob
mv "$log_dir"/*.stats "$gffcmp_dir"/
ln -vf "$gffcmp_dir"/*.stats "$gffcmp_summary_dir"

echo "done_${species_name}"



#record memory usage
cgroup_dir=$(awk -F: '{print $NF}' /proc/self/cgroup)
peak_mem=`cat /sys/fs/cgroup$cgroup_dir/memory.peak`
peak_mem_mb=$(awk "BEGIN {printf \"%.2f\", $peak_mem / 1048576}") #transfer to mb
echo ">Peak memory was $peak_mem_mb MegaBytes"

#record end
echo ">ENDING at $(date)"
