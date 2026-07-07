#!/bin/bash
# Build the concatenated dual-genome reference: prefix each parent's contigs with a
# parental name (chr1 -> B73_chr1 / Mo17_chr1), then concatenate both assemblies and
# both GTFs into one reference. GTF seqnames are prefixed to match; gene_id/transcript_id
# are left unchanged so paired-feature tables keep matching. Handles .gz inputs.

set -euo pipefail

# ":?" errors on unset OR empty
PARENT1_FA=${PARENT1_FA:?set PARENT1_FA}
PARENT2_FA=${PARENT2_FA:?set PARENT2_FA}
PARENT1_GTF=${PARENT1_GTF:?set PARENT1_GTF}
PARENT2_GTF=${PARENT2_GTF:?set PARENT2_GTF}
PARENT1_PREFIX=${PARENT1_PREFIX:?set PARENT1_PREFIX}
PARENT2_PREFIX=${PARENT2_PREFIX:?set PARENT2_PREFIX}
OUTDIR=${OUTDIR:-resources/references/${PARENT1_PREFIX}_${PARENT2_PREFIX}}

# identical prefixes would re-create the contig collision the design avoids
[ "$PARENT1_PREFIX" != "$PARENT2_PREFIX" ] || {
    echo "ERROR: PARENT1_PREFIX and PARENT2_PREFIX must differ (both '$PARENT1_PREFIX')." >&2
    exit 1
}

P1="${PARENT1_PREFIX}_"
P2="${PARENT2_PREFIX}_"
mkdir -p "$OUTDIR"
FA="$OUTDIR/concatenated.fa"
GTF="$OUTDIR/concatenated.gtf"
REPORT="$OUTDIR/reference_build_report.tsv"

read_any() { case "$1" in *.gz) zcat "$1" ;; *) cat "$1" ;; esac; }

# prefix + concatenate FASTA (header lines only)
echo "[$(date)] Prefixing and concatenating FASTA..."
read_any "$PARENT1_FA" | awk -v p="$P1" '/^>/{sub(/^>/, ">" p)} 1'  >  "$FA"
read_any "$PARENT2_FA" | awk -v p="$P2" '/^>/{sub(/^>/, ">" p)} 1'  >> "$FA"

# prefix + concatenate GTF (seqname = col 1)
echo "[$(date)] Prefixing and concatenating GTF..."
read_any "$PARENT1_GTF" | awk -v p="$P1" 'BEGIN{FS=OFS="\t"} /^#/{print;next} {$1=p$1} 1'  >  "$GTF"
read_any "$PARENT2_GTF" | awk -v p="$P2" 'BEGIN{FS=OFS="\t"} /^#/{print;next} {$1=p$1} 1'  >> "$GTF"

echo "[$(date)] Indexing..."
samtools faidx "$FA"

count_contigs() { awk -v p="$1" 'index($1,p)==1{n++} END{print n+0}' "$FA.fai"; }
sum_bp()        { awk -v p="$1" 'index($1,p)==1{s+=$2} END{print s+0}' "$FA.fai"; }

{
    echo -e "metric\tvalue"
    echo -e "parent1_prefix\t$PARENT1_PREFIX"
    echo -e "parent2_prefix\t$PARENT2_PREFIX"
    echo -e "parent1_contigs\t$(count_contigs "$P1")"
    echo -e "parent2_contigs\t$(count_contigs "$P2")"
    echo -e "parent1_bp\t$(sum_bp "$P1")"
    echo -e "parent2_bp\t$(sum_bp "$P2")"
    echo -e "total_bp\t$(awk '{s+=$2} END{print s+0}' "$FA.fai")"
} > "$REPORT"

echo "[$(date)] Done."
echo "  FASTA : $FA"
echo "  GTF   : $GTF"
echo "  Report: $REPORT"
