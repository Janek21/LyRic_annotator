#!/bin/bash

#SBATCH --output=logs/%x.%A_%a.out
#SBATCH --error=logs/%x.%A_%a.err

#SBATCH --job-name=srr_download

#SBATCH --qos=normal
#SBATCH --time=120

#SBATCH --mem=16G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4

#SBATCH --array=0-1

#22 arrays of 4=can process up to 88 lines, we have 86(last task ahas 2 lines only


echo ">STARTING at $(date)"

module load SRA-Toolkit

species_name="$1"
cd $species_name

#How many lines to process per array task
LINES_PER_TASK=1

#Line range for this specific array ID
#If ID=0, START=1, END=4.
START_LINE=$(( SLURM_ARRAY_TASK_ID * LINES_PER_TASK + 1 ))
END_LINE=$(( (SLURM_ARRAY_TASK_ID + 1) * LINES_PER_TASK ))

echo "Processing lines $START_LINE through $END_LINE from $species_name/srr_list.tsv"

#Extract block of 4 SRR IDs and loop them
selectedSRRs=$(sed -n "${START_LINE},${END_LINE}p" srr_list.tsv)


for SRRid in $selectedSRRs; do
	echo "--- Processing $SRRid ---"
	
	#define things
	out_dir="data/fastq"
	SRR_path="$out_dir/$SRRid.fastq"
	result_file="$out_dir/ont_HpreCap_0+_$SRRid.fastq.gz"
	
	
	if [ -f "$result_file" ]; then #check if file has alreacy been downloaded
		echo "Skipping $SRRid: $result_file already exists."
		
	else
		#Sownload the file
		fasterq-dump $SRRid -O "${out_dir}/" -e "$SLURM_CPUS_PER_TASK"
		
		#check existence each time to avoid getting stuck in errors
		if [ -f "$SRR_path" ]; then
		echo "downloaded $SRRid"
		
		#zip file
		pigz -9 -p "$SLURM_CPUS_PER_TASK" "$SRR_path"
		
		#rename the file to specific format
		mv "$SRR_path.gz" "$result_file"
		echo "Complete: new file is $result_file"
		else
		echo "Error: $SRR_path not found. Download may have failed."
		fi
	fi
done

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
