#!/bin/bash
# build_concatenated_reference.sh
# Build the hybrid (concatenated) reference from two parental assemblies.
#
# This is the FIRST real pipeline module. It takes two parental assemblies that
# independently name their contigs (and therefore collide - both may have chr1),
# adds a per-parent prefix to every contig name, and concatenates them into one
# hybrid reference. After this step, parent1 sequence lives on contigs prefixed
# with the user-provided parent1 name and parent2 sequence lives on contigs
# prefixed with the user-provided parent2 name.
#
# Inputs:
#   PARENT1_FA / PARENT2_FA   parental assemblies (bare contig names, may collide)
#   PARENT1_GTF / PARENT2_GTF parental annotations (bare contig names)
#   PARENT1_PREFIX / PARENT2_PREFIX parental names to write into output contigs
#
# Outputs:
#   results/reference/concatenated.fa
#   results/reference/concatenated.gtf
#   results/reference/reference_build_report.tsv
#
# Prefix convention: with PARENT1_PREFIX=B73, contig chr1 becomes B73_chr1.
# parent1/parent2 are internal roles only; output names use the user-provided
# parental prefixes.
# GTF seqnames are always prefixed to match FASTA contigs. GTF gene_id and
# transcript_id attributes are preserved by default so user-provided gene pair
# tables can keep matching the final annotation.
#
# Requires: samtools, awk, perl

set -euo pipefail

# ---- parameters -------------------------------------------------------------
REFERENCE_ID=${REFERENCE_ID:-}
PARENT1_FA=${PARENT1_FA:-${P1_FA:-P1.fa}}
PARENT2_FA=${PARENT2_FA:-${P2_FA:-P2.fa}}
PARENT1_GTF=${PARENT1_GTF:-${P1_GTF:-P1.gtf}}
PARENT2_GTF=${PARENT2_GTF:-${P2_GTF:-P2.gtf}}
PARENT1_PREFIX=${PARENT1_PREFIX:-${P1_PREFIX:-}}
PARENT2_PREFIX=${PARENT2_PREFIX:-${P2_PREFIX:-}}
SEP=${SEP:-_}                       # separator between prefix and contig name
PREFIX_GTF_ATTRIBUTES=${PREFIX_GTF_ATTRIBUTES:-false}
OUTDIR=${OUTDIR:-results/reference}

for cmd in samtools awk perl; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "ERROR: required command not found in PATH: $cmd" >&2
        exit 1
    fi
done

if [ -z "$PARENT1_PREFIX" ] || [ -z "$PARENT2_PREFIX" ]; then
    echo "ERROR: set PARENT1_PREFIX and PARENT2_PREFIX to real parental names." >&2
    echo "       Backward-compatible aliases P1_PREFIX and P2_PREFIX are also accepted." >&2
    echo "       Example: PARENT1_PREFIX=B73 PARENT2_PREFIX=Mo17" >&2
    exit 1
fi

if [ "$PARENT1_PREFIX" = "$PARENT2_PREFIX" ]; then
    echo "ERROR: PARENT1_PREFIX and PARENT2_PREFIX must differ (both are '$PARENT1_PREFIX')." >&2
    echo "       Identical prefixes would re-create contig-name collisions." >&2
    exit 1
fi

case "$PREFIX_GTF_ATTRIBUTES" in
    true|false) ;;
    *)
        echo "ERROR: PREFIX_GTF_ATTRIBUTES must be true or false. Got: $PREFIX_GTF_ATTRIBUTES" >&2
        exit 1
        ;;
esac

for f in "$PARENT1_FA" "$PARENT2_FA" "$PARENT1_GTF" "$PARENT2_GTF"; do
    [ -f "$f" ] || { echo "ERROR: input not found: $f" >&2; exit 1; }
done

PARENT1_FULL_PREFIX="${PARENT1_PREFIX}${SEP}"
PARENT2_FULL_PREFIX="${PARENT2_PREFIX}${SEP}"
GTF_ATTRIBUTE_MODE=preserved
if [ "$PREFIX_GTF_ATTRIBUTES" = "true" ]; then
    GTF_ATTRIBUTE_MODE=prefixed
fi

check_fasta_not_prefixed() {
    local fasta=$1
    local first_prefixed
    first_prefixed=$(awk -v p1="$PARENT1_FULL_PREFIX" -v p2="$PARENT2_FULL_PREFIX" '
        /^>/ {
            name = substr($1, 2)
            if (index(name, p1) == 1 || index(name, p2) == 1) {
                print name
                exit
            }
        }' "$fasta")
    if [ -n "$first_prefixed" ]; then
        echo "ERROR: input FASTA appears to be already parent-prefixed: $fasta" >&2
        echo "       First prefixed contig seen: $first_prefixed" >&2
        echo "       Refusing to add ${PARENT1_FULL_PREFIX}/${PARENT2_FULL_PREFIX} a second time." >&2
        exit 1
    fi
}

prefix_gtf() {
    local input_gtf=$1
    local parent_prefix=$2
    HYARQ_PARENT_PREFIX="$parent_prefix" \
    HYARQ_PREFIX_GTF_ATTRIBUTES="$PREFIX_GTF_ATTRIBUTES" \
    perl -e '
        use strict;
        use warnings;
        my $prefix = $ENV{"HYARQ_PARENT_PREFIX"};
        my $prefix_attrs = $ENV{"HYARQ_PREFIX_GTF_ATTRIBUTES"};
        while (my $line = <>) {
            chomp $line;
            if ($line =~ /^#/) {
                print $line, "\n";
                next;
            }
            my @fields = split(/\t/, $line, -1);
            if (@fields >= 9) {
                $fields[0] = $prefix . $fields[0];
                if ($prefix_attrs eq "true") {
                    $fields[8] =~ s/(gene_id\s+")([^"]*)(")/$1$prefix$2$3/g;
                    $fields[8] =~ s/(transcript_id\s+")([^"]*)(")/$1$prefix$2$3/g;
                }
                print join("\t", @fields), "\n";
            } else {
                print $line, "\n";
            }
        }
    ' "$input_gtf"
}

count_cross_parent_gene_id_duplicates() {
    awk -F'\t' -v p1="$PARENT1_FULL_PREFIX" -v p2="$PARENT2_FULL_PREFIX" '
        function gene_id(attrs,    start, rest, endq) {
            start = index(attrs, "gene_id \"")
            if (start == 0) return ""
            rest = substr(attrs, start + 9)
            endq = index(rest, "\"")
            if (endq == 0) return ""
            return substr(rest, 1, endq - 1)
        }
        /^#/ || NF < 9 { next }
        {
            gid = gene_id($9)
            if (gid == "") next
            if (index($1, p1) == 1) parent1_seen[gid] = 1
            else if (index($1, p2) == 1) parent2_seen[gid] = 1
        }
        END {
            n = 0
            for (gid in parent1_seen) {
                if (gid in parent2_seen) n++
            }
            print n + 0
        }' "$CAT_GTF"
}

check_fasta_not_prefixed "$PARENT1_FA"
check_fasta_not_prefixed "$PARENT2_FA"

mkdir -p "$OUTDIR"
CAT_FA="$OUTDIR/concatenated.fa"
CAT_GTF="$OUTDIR/concatenated.gtf"
REPORT="$OUTDIR/reference_build_report.tsv"

# ---- prefix + concatenate FASTA --------------------------------------------
# Only rename header lines (lines starting with >). Sequence lines untouched.
echo "[$(date)] Prefixing and concatenating FASTA..."
awk -v p="$PARENT1_FULL_PREFIX" '/^>/{sub(/^>/, ">" p)} {print}' "$PARENT1_FA"  >  "$CAT_FA"
awk -v p="$PARENT2_FULL_PREFIX" '/^>/{sub(/^>/, ">" p)} {print}' "$PARENT2_FA"  >> "$CAT_FA"

# ---- prefix + concatenate GTF ----------------------------------------------
# Column 1 (seqname) is always prefixed. gene_id/transcript_id attributes are
# preserved unless PREFIX_GTF_ATTRIBUTES=true.
echo "[$(date)] Prefixing and concatenating GTF..."
prefix_gtf "$PARENT1_GTF" "$PARENT1_FULL_PREFIX"  >  "$CAT_GTF"
prefix_gtf "$PARENT2_GTF" "$PARENT2_FULL_PREFIX"  >> "$CAT_GTF"

# ---- index ------------------------------------------------------------------
echo "[$(date)] Indexing concatenated FASTA..."
samtools faidx "$CAT_FA"

# ---- build report -----------------------------------------------------------
echo "[$(date)] Writing build report..."
PARENT1_CONTIGS=$(awk -v p="$PARENT1_FULL_PREFIX" 'BEGIN{FS="\t"} $1 ~ "^" p' "$CAT_FA.fai" | wc -l | tr -d ' ')
PARENT2_CONTIGS=$(awk -v p="$PARENT2_FULL_PREFIX" 'BEGIN{FS="\t"} $1 ~ "^" p' "$CAT_FA.fai" | wc -l | tr -d ' ')
PARENT1_BP=$(awk -v p="$PARENT1_FULL_PREFIX" 'BEGIN{FS="\t";s=0} $1 ~ "^" p {s+=$2} END{print s+0}' "$CAT_FA.fai")
PARENT2_BP=$(awk -v p="$PARENT2_FULL_PREFIX" 'BEGIN{FS="\t";s=0} $1 ~ "^" p {s+=$2} END{print s+0}' "$CAT_FA.fai")
TOTAL_CONTIGS=$((PARENT1_CONTIGS + PARENT2_CONTIGS))
TOTAL_BP=$((PARENT1_BP + PARENT2_BP))
DUPLICATE_GENE_IDS=$(count_cross_parent_gene_id_duplicates)

{
    echo -e "# HyARQ reference build report"
    echo -e "# generated: $(date)"
    echo -e "metric\tvalue"
    if [ -n "$REFERENCE_ID" ]; then
        echo -e "reference_id\t${REFERENCE_ID}"
    fi
    echo -e "parent1_fasta\t${PARENT1_FA}"
    echo -e "parent2_fasta\t${PARENT2_FA}"
    echo -e "parent1_gtf\t${PARENT1_GTF}"
    echo -e "parent2_gtf\t${PARENT2_GTF}"
    echo -e "parent1_prefix\t${PARENT1_PREFIX}"
    echo -e "parent2_prefix\t${PARENT2_PREFIX}"
    echo -e "gtf_attributes_prefixed\t${PREFIX_GTF_ATTRIBUTES}"
    echo -e "duplicate_gene_ids_detected\t${DUPLICATE_GENE_IDS}"
    echo -e "parent1_contigs\t${PARENT1_CONTIGS}"
    echo -e "parent2_contigs\t${PARENT2_CONTIGS}"
    echo -e "parent1_bp\t${PARENT1_BP}"
    echo -e "parent2_bp\t${PARENT2_BP}"
    echo -e "total_contigs\t${TOTAL_CONTIGS}"
    echo -e "total_bp\t${TOTAL_BP}"
    echo -e "#"
    echo -e "# per-contig detail"
    echo -e "contig\tlength\tsource_parent"
    awk -v p="$PARENT1_FULL_PREFIX" 'BEGIN{FS="\t"} $1 ~ "^" p {print $1"\t"$2"\tparent1"}' "$CAT_FA.fai"
    awk -v p="$PARENT2_FULL_PREFIX" 'BEGIN{FS="\t"} $1 ~ "^" p {print $1"\t"$2"\tparent2"}' "$CAT_FA.fai"
} > "$REPORT"

if [ "$PREFIX_GTF_ATTRIBUTES" = "false" ] && [ "$DUPLICATE_GENE_IDS" -gt 0 ]; then
    echo "ERROR: detected ${DUPLICATE_GENE_IDS} gene_id value(s) present in both parent annotations." >&2
    echo "       Duplicate gene IDs can cause featureCounts to merge parent-specific features." >&2
    echo "       Provide parent annotations with unique gene IDs, provide disambiguated GTFs," >&2
    echo "       or rerun with PREFIX_GTF_ATTRIBUTES=true." >&2
    echo "       Report written before exit: $REPORT" >&2
    exit 1
fi

# ---- console summary --------------------------------------------------------
echo ""
echo "[$(date)] Reference build complete."
echo "  FASTA : $CAT_FA"
echo "  GTF   : $CAT_GTF"
echo "  Report: $REPORT"
echo "  GTF attributes: $GTF_ATTRIBUTE_MODE"
echo ""
echo "Contigs in hybrid reference:"
cut -f1,2 "$CAT_FA.fai" | sed 's/^/  /'
