#!/bin/bash

#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

#SBATCH --job-name=heavy_srrDW

#SBATCH --qos=normal
#SBATCH --time=650

#SBATCH --mem=6G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2

echo ">STARTING at $(date)"

module load SRA-Toolkit

species_name="$1"
cd $species_name

dw_link="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR100/051/ERR10032851/ERR10032851.fastq.gz"
SRRid="ERR10032851"

echo "--- Processing $SRRid ---"
	
#define things
out_dir="data/fastq"
SRR_path="$out_dir/$SRRid.fastq.gz"
result_file="$out_dir/ont_HpreCap_0+_$SRRid.fastq.gz"
	
	
if [ -f "$result_file" ]; then #check if file has alreacy been downloaded
	echo "Skipping $SRRid: $result_file already exists."
		
else
	#Sownload the file
	#wget -nc -P "${out_dir}/" "$dw_link"
		
	#check existence each time to avoid getting stuck in errors
	if [ -f "$SRR_path" ]; then
		echo "downloaded $SRRid"

		#rename the file to specific format
		mv "$SRR_path" "$result_file"
		echo "Complete: new file is $result_file"
	else
		echo "Error: $SRR_path not found. Download may have failed."
	fi
fi


# Record memory usage (at the end of all 4 downloads)
cgroup_dir=$(awk -F: '{print $NF}' /proc/self/cgroup)
# Check if the path exists to avoid errors on different cgroup versions
if [ -f "/sys/fs/cgroup$cgroup_dir/memory.peak" ]; then
	peak_mem=$(cat "/sys/fs/cgroup$cgroup_dir/memory.peak")
	peak_mem_mb=$(awk "BEGIN {printf \"%.2f\", $peak_mem / 1048576}")
	echo ">Peak memory was $peak_mem_mb MegaBytes"
fi

# Record end
echo ">ENDING at $(date)"
