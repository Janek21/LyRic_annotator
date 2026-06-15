#!/usr/bin/env bash
# Predict ORFs with TD2 and turn them into:
#   1. a longest-isoform protein FASTA for BUSCO (always)
#   2. a CDS-augmented copy of the input annotation (only if out_cds_gff is given)
#
# LyRic/tmerge annotations model exons only; TD2 predicts the ORF per transcript,
# which is lifted back to the genome to recover CDS features. The same ORF set is
# filtered to one (longest) isoform per gene so BUSCO sees a non-redundant proteome.
#
# Usage:
#   infer_cds.sh <annotation_gff> <genome_fa> <gcode> <workdir> <out_prot_fa> [out_cds_gff]
# Pass out_cds_gff only when the input is exon-only (no pre-existing CDS): the
# predicted CDS rows are concatenated onto it. Omit it to emit proteins alone.

set -euo pipefail

ann_gff="$1"
genome_fa="$2"
gcode="$3"
workdir="$4"
out_prot="$5"
out_cds="${6:-}"

td2_util="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/TD2/util"
mkdir -p "$workdir"

#normalised GTF (always carries gene_id/transcript_id) + transcriptome over the
#full annotation, so ids line up across the transcriptome, the TD2 ORFs and the
#alignment gff3.
norm_gtf="$workdir/normalized.gtf"
transcripts="$(realpath "$workdir/transcripts.fa")"
gffread "$ann_gff" -T -o "$norm_gtf"
gffread "$norm_gtf" -g "$genome_fa" -w "$transcripts"

#predict ORFs; TD2 writes transcripts.fa.TD2.{gff3,pep} into the workdir
( cd "$workdir" \
	&& TD2.LongOrfs -t "$transcripts" -O . -G "$gcode" \
	&& TD2.Predict  -t "$transcripts" -O . -G "$gcode" )
td2_gff3="$workdir/transcripts.fa.TD2.gff3"
td2_pep="$workdir/transcripts.fa.TD2.pep"

#longest isoform per gene, from exon spans in the normalised GTF (ids match the
#TD2 .pep records exactly)
longest_ids="$workdir/longest_ids.txt"
awk -F'\t' '
	$3=="exon" {
		t=""; g="";
		if (match($9, /transcript_id "[^"]+/)) t=substr($9, RSTART+15, RLENGTH-15)
		if (match($9, /gene_id "[^"]+/))       g=substr($9, RSTART+9,  RLENGTH-9)
		if (t=="") next
		if (g=="") g=t
		len[t]+=$5-$4+1; gene[t]=g
	}
	END { for (t in len) print gene[t]"\t"t"\t"len[t] }' "$norm_gtf" \
	| sort -k1,1 -k3,3nr | awk -F'\t' '!seen[$1]++ {print $2}' > "$longest_ids"

#keep only ORFs whose transcript is a longest isoform (strip the .pN ORF suffix)
awk 'NR==FNR { keep[$1]=1; next }
	/^>/ { tid=substr($1,2); sub(/\.p[0-9]+$/,"",tid); out=(tid in keep) }
	out' "$longest_ids" "$td2_pep" > "$out_prot"

#optional: lift the transcript-space ORFs onto the genome and splice the CDS rows
#onto the (exon-only) input annotation
if [ -n "$out_cds" ]; then
	aln_gff3="$workdir/alignment.gff3"
	cds_gff3="$workdir/cds.gff3"
	perl "$td2_util/gtf_to_alignment_gff3.pl" "$norm_gtf" > "$aln_gff3"
	perl "$td2_util/cdna_alignment_orf_to_genome_orf.pl" "$td2_gff3" "$aln_gff3" "$transcripts" > "$cds_gff3"
	{ grep -P "\tCDS\t" "$cds_gff3" || true; } | cat "$ann_gff" - > "$out_cds"
fi
