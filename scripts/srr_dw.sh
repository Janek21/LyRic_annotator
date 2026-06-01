#!/bin/bash

#SBATCH --output=logs/%x.%A_%a.out
#SBATCH --error=logs/%x.%A_%a.err

#SBATCH --job-name=srr_download

#SBATCH --qos=normal
#SBATCH --time=180

#SBATCH --mem=4G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2

echo ">STARTING at $(date)"

set -euo pipefail

species_name="$1"
cd "$species_name"

out_dir="data/fastq"
mkdir -p "$out_dir"

# accession for this array task (line SLURM_ARRAY_TASK_ID + 1)
SRRid=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" srr_list.tsv)

if [ -z "$SRRid" ]; then
	echo "ERROR: no accession at index $SLURM_ARRAY_TASK_ID in $species_name/srr_list.tsv"
	exit 1
fi

result_file="$out_dir/ont_HpreCap_0+_$SRRid.fastq.gz"

if [ -f "$result_file" ]; then #check if file has already been downloaded
	echo "Skipping $SRRid: $result_file already exists."
else
	echo "--- Processing $SRRid ---"

	#Resolve the FASTQ FTP URL from ENA (returns paths without ftp:// prefix)
	FTP_URL=$(curl -sf "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${SRRid}&result=read_run&fields=fastq_ftp" \
		| tail -n +2 | cut -f2 | tr ';' '\n' | head -1)

	if [ -z "$FTP_URL" ]; then
		echo "Error: could not retrieve FTP URL for $SRRid from ENA."
		exit 1
	fi

	#Download to a temp file, move on success (avoids leaving a partial file that
	#the skip-check above would mistake for a finished download)
	tmp_file="$out_dir/$SRRid.fastq.gz.part"
	wget -q -O "$tmp_file" "ftp://${FTP_URL}"
	mv "$tmp_file" "$result_file"
	echo "Complete: new file is $result_file"
fi

# Record memory usage at the end
cgroup_dir=$(awk -F: '{print $NF}' /proc/self/cgroup)
# Check if the path exists to avoid errors on different cgroup versions
if [ -f "/sys/fs/cgroup$cgroup_dir/memory.peak" ]; then
	peak_mem=$(cat "/sys/fs/cgroup$cgroup_dir/memory.peak")
	peak_mem_mb=$(awk "BEGIN {printf \"%.2f\", $peak_mem / 1048576}")
	echo ">Peak memory was $peak_mem_mb MegaBytes"
fi

# Record end
echo ">ENDING at $(date)"
