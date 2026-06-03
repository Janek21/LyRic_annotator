#!/bin/bash
# Orchestrate the full LyRic pipeline for one species, chaining the stages with
# SLURM job dependencies so this script returns immediately:
#   1. lyric_template.sh  - set up the species dir + submit the ENA download array
#   2. runner.sh          - run the LyRic (snakemake) pipeline, after downloads finish
#   3. evaluation.sh      - merge + evaluate the annotation, after the pipeline finishes
#
# Usage: bash lyric_prepare.sh <Genus_species[_extra]> [longread_db] [busco_db]

set -euo pipefail

species_name="$1"
longread_db="${2:-../data/longread_protists.tsv}"
busco_db="${3:-/no_backup/rg/references/busco_downloads}"

sp=$(echo "$species_name" | cut -f2 -d"_")

mkdir -p logs

#1.Set up the species and submit the download
echo "=== Stage 1: lyric_template.sh ($species_name) ==="
#stream lyric_template.sh output live to the shell while capturing it to extract the job id
template_log=$(mktemp)
bash lyric_template.sh "$species_name" "$longread_db" 2>&1 | tee "$template_log"
dl_jobid=$(sed -n 's/^DOWNLOAD_JOBID=//p' "$template_log" | tail -1)
rm -f "$template_log"

#2. Run the pipeline once the downloads finish.(from species_name dir)
echo "=== Stage 2: runner.sh ==="
mkdir -p "$species_name/logs"
echo "runner.sh will start after download job $dl_jobid"
run_jobid=$(sbatch --parsable \
	--job-name="lyric_${sp}" \
	--chdir="$species_name" \
	--dependency=afterok:"$dl_jobid" \
	--qos=normal \
	--cpus-per-task=6 \
	--mem=36G \
	--time=500 \
	--output="logs/%x_%j.out" \
	--error="logs/%x_%j.err" \
	"$species_name/runner.sh")
echo "Pipeline submitted: job $run_jobid"

#3. Evaluate once the pipeline finishes (evaluation.sh runs from the repo root).
echo "=== Stage 3: evaluation.sh ==="
ev_jobid=$(sbatch --parsable \
	--job-name="eval_ly_${sp}" \
	--dependency=afterok:"$run_jobid" \
	--cpus-per-task=4 \
	--mem=12G \
	--time=90 \
	--output="logs/eval/%x_%j.out" \
	--error="logs/eval/%x_%j.err" \
	evaluation.sh "$species_name" "$busco_db")
echo "Evaluation submitted: job $ev_jobid (starts after job $run_jobid)"

echo "Pipeline chain for $species_name: download(${dl_jobid:-none}) -> run($run_jobid) -> eval($ev_jobid)"
