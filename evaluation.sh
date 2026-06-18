#!/usr/bin/env bash
#SBATCH --job-name=lyric_eval
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=90
#SBATCH --output=logs/eval/%x_%j.out
#SBATCH --error=logs/eval/%x_%j.err

echo ">STARTING at $(date)"

species_name="$1"
busco_db="${2:-/no_backup/rg/references/busco_downloads}"

set -euo pipefail

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

#cativate busco conda env
source $(conda info --base)/etc/profile.d/conda.sh
conda activate buscomania

#create storing folders and variables
tmp_files="$species_name/output/files"
lyric_out="$species_name/output/mappings/mergedReads"

rm -rf "$tmp_files"
mkdir -p "$tmp_files"

#per-task AGAT config so parallel jobs don't collide on agat_config.yaml
agat_cfg="$tmp_files/agat_${species_name}_${SLURM_ARRAY_TASK_ID:-$$}.yaml"
agat config --expose --no-log --output "$agat_cfg" >/dev/null 2>&1
trap 'rm -f "$agat_cfg"' EXIT

##decompressions
#decompress gffs
find "$lyric_out" -type f -name "ont_*.gz"|xargs -r -P $(nproc) unpigz -df  #$(nproc) unpigz -df #"$SLURM_CPUS_PER_TASK" unpigz -df
#genome stays gzipped in ../data/species; the uncompressed copy in data/fasta is used instead

#rename for long file names
#Removes the prefix and sufix(minNreads (N varies))
for f in "$lyric_out"/ont_HpreCap_0+_[DSE]RR*.gff; do
    [ -e "$f" ] || continue
    b=$(basename "$f"); b=${b#ont_HpreCap_0+_}; b=${b%%.HiSS.*}
    mv -- "$f" "$lyric_out/$b.gff"
done

#detect number of files in folder(if 1 only, dont merge)
shopt -s nullglob
gff_files=("$lyric_out"/*.gff)
shopt -u nullglob
file_count=${#gff_files[@]}

echo "FC is $file_count"

if [ "$file_count" -eq 1 ]; then
	#normalise gff
	agat_convert_sp_gxf2gxf.pl -g "${gff_files[0]}" --config "$agat_cfg" -o "$tmp_files/merged_${sp}_ann.gff"
	echo "Normalized single file at $tmp_files/merged_${sp}_ann.gff"
else
	#merge gffs
	agat_sp_merge_annotations.pl --gff "$lyric_out" --config "$agat_cfg" --out "$tmp_files/merged_${sp}_ann.gff"
	echo "Merged files at $tmp_files/merged_${sp}_ann.gff"
fi

#count gene and transcript models in the merged annotation
counts_dir="summary/counts"
mkdir -p "$counts_dir"
merged="$tmp_files/merged_${sp}_ann.gff"
#taxon id = most repeated id in the taxon column (drives the genetic code and BUSCO lineage)
taxonID=$(cut "$species_name/srr_select.tsv" -f4|sort|uniq -c|sort -nr|awk '{print $2}'|head -n1)
counts=$(count_models "$merged")
gene_count=${counts%%$'\t'*}
transcript_count=${counts##*$'\t'}
echo "      Gene models: $gene_count | Transcript models: $transcript_count"

#resolve the NCBI nuclear genetic code for this taxon (codon table for on-standard cases)
echo "KEY: $NCBI_API_KEY"
gcode=$(python3 scripts/get_genetic_code.py -e "ibdyjsayzcllkyvjkc@nespf.com" -k "${NCBI_API_KEY:-}" -t "$taxonID" 2>/dev/null)
if ! [[ "$gcode" =~ ^[0-9]+$ ]]; then
	echo ">Could not resolve genetic code for taxon $taxonID; defaulting to table 1."
	gcode=1
fi
echo "Translation table for $taxonID: $gcode"

#predict ORFs with TD2 and (1) splice the resulting CDS onto the exon-only LyRic
#annotation (tmerge models exons only), (2) emit one protein per gene (longest
#isoform) for BUSCO. See scripts/infer_cds.sh.
td_work="$tmp_files/transdecoder_work"
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

cds_merged="$tmp_files/CDSmerged_${sp}_ann.gff"
prot_file="$tmp_files/prot_$sp.fa"
bash scripts/infer_cds.sh "$merged" "$genome_fa" "$gcode" "$td_work" "$prot_file" "$cds_merged"
echo "CDS-augmented annotation written to $cds_merged"
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

##run busco inline (joined from the former scripts/busco_evaluation.sh so the whole
##evaluation runs as one job; no extra sbatch dependency to chain)
cpus="${SLURM_CPUS_PER_TASK:-$(nproc)}"
res_lineage="$species_name/output/busco_res_lineage"
res_euk="$species_name/output/busco_res_eukaryote"
odb_version="odb12"
busco_lineage_dir="summary/busco_lineage"
busco_euk_dir="summary/busco_eukaryote"

rm -rf "$res_lineage" "$res_euk"
mkdir -p "$res_lineage" "$res_euk" "$busco_lineage_dir" "$busco_euk_dir"

#get the taxon-specific (custom) lineage (taxonID resolved above for the genetic code)
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

#predicted (CDS-augmented) annotation relocated to summary/, hardlinked back to the species location
#(canonical inode lives in summary/pred so the species folder can be removed safely).
#the exon-only $merged stays in place: mass_merge/mass_pred_merge gate on it and
#merge_evaluation falls back to it for species evaluated before CDS integration.
pred_dir="summary/pred"
mkdir -p "$pred_dir"
pred_dest="$pred_dir/${species_name}_${taxonID}_pred.gff"
rm -f "$pred_dest"                 #refresh on reruns
mv "$cds_merged" "$pred_dest"     #relocate the CDS-augmented prediction into the central summary tree
ln "$pred_dest" "$cds_merged"     #link it back so the original species location stays valid
echo "Predicted annotation collected into $pred_dir/"

rm -rf agat_log_*
echo "Analysis completed!"
echo ">ENDING at $(date)"
