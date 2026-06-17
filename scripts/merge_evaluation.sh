#!/usr/bin/env bash
#SBATCH --job-name=lyric_merged_eval
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=90
#SBATCH --output=logs/eval/%x_%j.out
#SBATCH --error=logs/eval/%x_%j.err
# Merge the LyRic annotation with the reference annotation and evaluate the result.
# Runs only for species that have a reference annotation (data/input/Annotation.gff).
# Evaluation = BUSCO protein completeness + gene/transcript counts (no gffcompare).
#
# Usage: bash scripts/merge_evaluation.sh <species_name> [busco_db]

echo ">STARTING at $(date)"

species_name="$1"
busco_db="${2:-/no_backup/rg/references/busco_downloads}"
data_root="${3:-../data/species}"   #where the per-species source assemblies/annotations live
cpus="${SLURM_CPUS_PER_TASK:-$(nproc)}"

sp=$(echo "$species_name"|cut -f2 -d"_")
echo "$sp"

#count gene + transcript models in one gffread pass. --keep-genes normalises any
#input (tmerge GTF or AGAT GFF3) into gene + transcript records: real gene
#features are preserved (AGAT-clustered loci, not the per-transcript gene_id),
#and one gene + transcript is synthesised per id when the input has no such
#feature. Prints "<gene_count>\t<transcript_count>".
count_models() {  # $1=gff
	{ gffread "$1" --keep-genes -o - 2>/dev/null || true; } | awk -F'\t' '
		/^#/ { next }
		$3 ~ /^([A-Za-z_]*gene)$/                { g++; next }
		$3 ~ /^(transcript|mRNA|[A-Za-z_]*RNA)$/ { t++ }
		END { printf "%d\t%d\n", g, t }'
}

#activate the shared conda env (agat, gffread, TD2, busco)
source $(conda info --base)/etc/profile.d/conda.sh
conda activate buscomania

tmp_files="$species_name/output/files"
mkdir -p "$tmp_files"

#reference annotation placed by lyric_template.sh in data/input/Annotation.gff. check if it can be used by lookig if gff exists on ../data/species
shopt -s nullglob
own_ref_src=("$data_root/${species_name}"*/GC*/"${species_name}"*GC*.gff.gz)
shopt -u nullglob

ref_gff="$species_name/data/input/Annotation.gff"
#LyRic annotation produced by evaluation.sh: prefer the CDS-augmented annotation
#(exons+CDS); fall back to the exon-only merge for species evaluated before CDS
#integration so reruns of older results keep working.
lyric_gff="$tmp_files/CDSmerged_${sp}_ann.gff"
[ -s "$lyric_gff" ] || lyric_gff="$tmp_files/merged_${sp}_ann.gff"

#only species with their OWN reference annotation are merged; the rest are logged and skipped
merge_summary_dir="summary/merge"
mkdir -p "$merge_summary_dir"
if [ "${#own_ref_src[@]}" -eq 0 ] || [ ! -s "${own_ref_src[0]}" ]; then
	echo "No own reference annotation for $species_name under $data_root (only a fallback or none); skipping merge."
	#record once so reruns don't pile up duplicate lines
	grep -qxF "$species_name" "$merge_summary_dir/no_reference.txt" 2>/dev/null \
		|| echo "$species_name" >> "$merge_summary_dir/no_reference.txt"
	exit 0
fi
if [ ! -s "$ref_gff" ]; then
	echo "Own reference source exists but $ref_gff is missing; re-run lyric_template.sh setup. Aborting."
	exit 1
fi
if [ ! -s "$lyric_gff" ]; then
	echo "LyRic annotation $lyric_gff missing; run evaluation.sh first. Aborting."
	exit 1
fi

#per-task AGAT config so parallel jobs don't collide on agat_config.yaml
agat_cfg="$tmp_files/agat_merge_${species_name}_${SLURM_ARRAY_TASK_ID:-$$}.yaml"
agat config --expose --no-log --output "$agat_cfg" >/dev/null 2>&1
trap 'rm -f "$agat_cfg"' EXIT

#merge LyRic + reference into one non-redundant annotation
merged_ref="$tmp_files/mergedRef_${sp}_ann.gff"
agat_sp_merge_annotations.pl --gff "$lyric_gff" --gff "$ref_gff" --config "$agat_cfg" --out "$merged_ref"
echo "Merged LyRic + reference at $merged_ref"

#count gene and transcript models in the merged annotation
counts_dir="$merge_summary_dir/counts"
mkdir -p "$counts_dir"
#taxon id = most repeated id in the taxon column (same as evaluation.sh)
taxonID=$(cut "$species_name/srr_select.tsv" -f4|sort|uniq -c|sort -nr|awk '{print $2}'|head -n1)
counts=$(count_models "$merged_ref")
gene_count=${counts%%$'\t'*}
transcript_count=${counts##*$'\t'}
echo "      Gene models: $gene_count | Transcript models: $transcript_count"

#resolve the NCBI nuclear genetic code for this taxon (codon table for non-standard cases)
gcode=$(python3 scripts/get_genetic_code.py -e "ibdyjsayzcllkyvjkc@nespf.com" -k "${NCBI_API_KEY:-}" -t "$taxonID" 2>/dev/null)
if ! [[ "$gcode" =~ ^[0-9]+$ ]]; then
	echo ">Could not resolve genetic code for taxon $taxonID; defaulting to table 1."
	gcode=1
fi
echo "Translation table for $taxonID: $gcode"

#proteins for BUSCO via the same ORF inference as the prediction step
#(scripts/infer_cds.sh). The merged annotation already carries CDS - the CDS file
#is one of the merge inputs - so only the longest-isoform proteome is produced
#here, with no CDS splice.
shortname=$(python3 scripts/LyRic_setup.py shortname -s "$species_name")
genome_fa="$species_name/data/fasta/$shortname.fa"

#genome size = total assembly length (exact; sum of contig lengths, incl. N gaps)
fai="${genome_fa}.fai"
if [ -s "$fai" ]; then
	genome_size=$(cut -f2 "$fai" | awk '{s+=$1} END{print s+0}')
elif [[ "$genome_fa" == *.gz ]]; then
	genome_size=$(pigz -dcp "${SLURM_CPUS_PER_TASK:-$(nproc)}" "$genome_fa" \
		| awk '/^>/{next} {s+=length($0)} END{print s+0}')
else
	genome_size=$(awk '/^>/{next} {s+=length($0)} END{print s+0}' "$genome_fa")
fi
echo "      Genome size: ${genome_size} bp"

td_work="$tmp_files/transdecoder_merge_work"
prot_file="$tmp_files/protRef_$sp.fa"
bash scripts/infer_cds.sh "$merged_ref" "$genome_fa" "$gcode" "$td_work" "$prot_file"
echo "TransDecoder proteins (longest isoform per gene) in $prot_file"

# ── derived metrics + one consolidated per-species file ──────────────
# transcriptome_transcripts = records that entered ORF calling ($td_work/transcripts.fa);
# coding_transcripts        = unique source transcripts with >=1 TD2 ORF (full .pep,
#                             not the longest-isoform proteome), .pN suffix stripped.
tx_fa="$td_work/transcripts.fa"
pep_fa="$td_work/transcripts.fa.TD2.pep"
transcriptome_tx=$(grep -c '^>' "$tx_fa" 2>/dev/null || echo 0)
coding_tx=$( { grep '^>' "$pep_fa" 2>/dev/null || true; } \
	| awk '{sub(/^>/,"",$1); sub(/\.p[0-9]+$/,"",$1); print $1}' | sort -u | wc -l)

metrics_file="$counts_dir/${species_name}_${taxonID}_metrics.tsv"
awk -v sp="${species_name}_${taxonID}" \
	-v gc="$gene_count" -v tc="$transcript_count" -v gs="$genome_size" \
	-v cod="$coding_tx" -v ntx="$transcriptome_tx" 'BEGIN {
		OFS = "\t"
		gd  = (gs  > 0) ? sprintf("%.2f", gc / (gs / 1e6)) : "NA"
		td  = (gs  > 0) ? sprintf("%.2f", tc / (gs / 1e6)) : "NA"
		ipg = (gc  > 0) ? sprintf("%.3f", tc / gc)         : "NA"
		cf  = (ntx > 0) ? sprintf("%.4f", cod / ntx)       : "NA"
		print "species", "gene_count", "transcript_count", "genome_size_bp", \
		      "coding_transcripts", "transcriptome_transcripts", \
		      "gene_density_per_mb", "transcript_density_per_mb", \
		      "isoforms_per_gene", "coding_fraction"
		print sp, gc, tc, gs, cod, ntx, gd, td, ipg, cf
	}' > "$metrics_file"
read -r _ _ _ _ _ _ gd td ipg cf < <(tail -1 "$metrics_file")
echo "      Metrics: gene_density=${gd}/Mb transcript_density=${td}/Mb isoforms/gene=${ipg} coding_fraction=${cf}"
echo "      Consolidated metrics -> $metrics_file"

##run busco on the merged proteome (taxon-specific + eukaryote lineages)
res_lineage="$species_name/output/busco_mergedRef_lineage"
res_euk="$species_name/output/busco_mergedRef_eukaryote"
odb_version="odb12"
busco_lineage_dir="$merge_summary_dir/busco_lineage"
busco_euk_dir="$merge_summary_dir/busco_eukaryote"

rm -rf "$res_lineage" "$res_euk"
mkdir -p "$res_lineage" "$res_euk" "$busco_lineage_dir" "$busco_euk_dir"

#taxon-specific (custom) lineage
busco_lineage=$(python3 scripts/get_busco_db.py -e "ibdyjsayzcllkyvjkc@nespf.com" -t "$taxonID" -b "$busco_db/file_versions.tsv" -v "$odb_version")
echo "BUSCO lineage for $taxonID is $busco_lineage"
euk_lineage="eukaryota_${odb_version}"

#1. taxon-specific lineage busco -> summary/merge/busco_lineage/<stem>_Lbusco.json
busco -m protein -i "$prot_file" --download_path "$busco_db" -l "$busco_lineage" -c "$cpus" -f --out_path "${species_name}/output" -o busco_mergedRef_lineage --tar
lineage_json="$busco_lineage_dir/${species_name}_${taxonID}_Lbusco.json"
mv "$res_lineage"/*.json "$lineage_json"

#2. eukaryote lineage busco -> summary/merge/busco_eukaryote/<stem>_Ebusco.json
busco -m protein -i "$prot_file" --download_path "$busco_db" -l "$euk_lineage" -c "$cpus" -f --out_path "${species_name}/output" -o busco_mergedRef_eukaryote --tar
euk_json="$busco_euk_dir/${species_name}_${taxonID}_Ebusco.json"
mv "$res_euk"/*.json "$euk_json"

#merged-with-reference annotation relocated to summary/, hardlinked back to the species location
#(canonical inode lives in summary/merge/pred so the species folder can be removed safely)
pred_dir="$merge_summary_dir/pred"
mkdir -p "$pred_dir"
pred_dest="$pred_dir/${species_name}_${taxonID}_mergedRef.gff"
rm -f "$pred_dest"                 #refresh on reruns
mv "$merged_ref" "$pred_dest"     #relocate the merged annotation into the central summary tree
ln "$pred_dest" "$merged_ref"     #link it back so the original species location stays valid
echo "Merged-with-reference annotation collected into $pred_dir/"

rm -rf agat_log_*
echo "Analysis completed!"
echo ">ENDING at $(date)"
