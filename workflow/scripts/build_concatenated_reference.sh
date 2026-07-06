#!/bin/bash
# build_concatenated_reference.sh
# Build the hybrid (concatenated) reference from two parental assemblies.
#
# This is the FIRST real pipeline module. It takes two parental assemblies that
# independently name their contigs (and therefore collide - both may have chr1),
# adds a per-parent prefix to every contig name, and concatenates them into one
# hybrid reference. After this step, P1 sequence lives on P1_* contigs and P2
# sequence on P2_* contigs, so reads can later be split by contig prefix.
#
# Inputs:
#   P1.fa  P2.fa     parental assemblies (bare contig names, may collide)
#   P1.gtf P2.gtf    parental annotations (bare contig names)
#
# Outputs:
#   results/reference/concatenated.fa
#   results/reference/concatenated.gtf
#   results/reference/reference_build_report.tsv
#
# Prefix convention: contig chr1 of parent 1 becomes P1_chr1, etc.
# Prefixes are parameters (default P1 / P2) so real runs can use B73 / Mo17.
#
# Requires: samtools, awk

set -euo pipefail

# ---- parameters -------------------------------------------------------------
P1_FA=${P1_FA:-P1.fa}
P2_FA=${P2_FA:-P2.fa}
P1_GTF=${P1_GTF:-P1.gtf}
P2_GTF=${P2_GTF:-P2.gtf}
P1_PREFIX=${P1_PREFIX:-P1}
P2_PREFIX=${P2_PREFIX:-P2}
SEP=${SEP:-_}                       # separator between prefix and contig name
OUTDIR=${OUTDIR:-results/reference}

for cmd in samtools awk; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "ERROR: required command not found in PATH: $cmd" >&2
        exit 1
    fi
done

if [ "$P1_PREFIX" = "$P2_PREFIX" ]; then
    echo "ERROR: P1_PREFIX and P2_PREFIX must differ (both are '$P1_PREFIX')." >&2
    echo "       Identical prefixes would re-create contig-name collisions." >&2
    exit 1
fi

for f in "$P1_FA" "$P2_FA" "$P1_GTF" "$P2_GTF"; do
    [ -f "$f" ] || { echo "ERROR: input not found: $f" >&2; exit 1; }
done

mkdir -p "$OUTDIR"
CAT_FA="$OUTDIR/concatenated.fa"
CAT_GTF="$OUTDIR/concatenated.gtf"
REPORT="$OUTDIR/reference_build_report.tsv"

# ---- prefix + concatenate FASTA --------------------------------------------
# Only rename header lines (lines starting with >). Sequence lines untouched.
echo "[$(date)] Prefixing and concatenating FASTA..."
awk -v p="${P1_PREFIX}${SEP}" '/^>/{sub(/^>/, ">" p)} {print}' "$P1_FA"  >  "$CAT_FA"
awk -v p="${P2_PREFIX}${SEP}" '/^>/{sub(/^>/, ">" p)} {print}' "$P2_FA"  >> "$CAT_FA"

# ---- prefix + concatenate GTF ----------------------------------------------
# Only column 1 (seqname) is prefixed. Comment lines (#...) are passed through.
echo "[$(date)] Prefixing and concatenating GTF..."
awk -v p="${P1_PREFIX}${SEP}" 'BEGIN{FS=OFS="\t"} /^#/{print; next} {$1=p $1; print}' "$P1_GTF"  >  "$CAT_GTF"
awk -v p="${P2_PREFIX}${SEP}" 'BEGIN{FS=OFS="\t"} /^#/{print; next} {$1=p $1; print}' "$P2_GTF"  >> "$CAT_GTF"

# ---- index ------------------------------------------------------------------
echo "[$(date)] Indexing concatenated FASTA..."
samtools faidx "$CAT_FA"

# ---- build report -----------------------------------------------------------
echo "[$(date)] Writing build report..."
P1_CONTIGS=$(awk -v p="${P1_PREFIX}${SEP}" 'BEGIN{FS="\t"} $1 ~ "^" p' "$CAT_FA.fai" | wc -l | tr -d ' ')
P2_CONTIGS=$(awk -v p="${P2_PREFIX}${SEP}" 'BEGIN{FS="\t"} $1 ~ "^" p' "$CAT_FA.fai" | wc -l | tr -d ' ')
P1_BP=$(awk -v p="${P1_PREFIX}${SEP}" 'BEGIN{FS="\t";s=0} $1 ~ "^" p {s+=$2} END{print s+0}' "$CAT_FA.fai")
P2_BP=$(awk -v p="${P2_PREFIX}${SEP}" 'BEGIN{FS="\t";s=0} $1 ~ "^" p {s+=$2} END{print s+0}' "$CAT_FA.fai")
TOTAL_CONTIGS=$((P1_CONTIGS + P2_CONTIGS))
TOTAL_BP=$((P1_BP + P2_BP))

{
    echo -e "# HyARQ reference build report"
    echo -e "# generated: $(date)"
    echo -e "metric\tvalue"
    echo -e "P1_prefix\t${P1_PREFIX}${SEP}"
    echo -e "P2_prefix\t${P2_PREFIX}${SEP}"
    echo -e "P1_contigs\t${P1_CONTIGS}"
    echo -e "P2_contigs\t${P2_CONTIGS}"
    echo -e "total_contigs\t${TOTAL_CONTIGS}"
    echo -e "P1_bp\t${P1_BP}"
    echo -e "P2_bp\t${P2_BP}"
    echo -e "total_bp\t${TOTAL_BP}"
    echo -e "#"
    echo -e "# per-contig detail"
    echo -e "contig\tlength\tsource_parent"
    awk -v p="${P1_PREFIX}${SEP}" 'BEGIN{FS="\t"} $1 ~ "^" p {print $1"\t"$2"\tP1"}' "$CAT_FA.fai"
    awk -v p="${P2_PREFIX}${SEP}" 'BEGIN{FS="\t"} $1 ~ "^" p {print $1"\t"$2"\tP2"}' "$CAT_FA.fai"
} > "$REPORT"

# ---- console summary --------------------------------------------------------
echo ""
echo "[$(date)] Reference build complete."
echo "  FASTA : $CAT_FA"
echo "  GTF   : $CAT_GTF"
echo "  Report: $REPORT"
echo ""
echo "Contigs in hybrid reference:"
cut -f1,2 "$CAT_FA.fai" | sed 's/^/  /'
