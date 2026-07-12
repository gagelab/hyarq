#!/bin/bash
# count_rna_genes.sh
# Count RNA-seq fragments over paired parent1/parent2 genes on the concatenated reference.
#
# Required inputs:
#   BAM          coordinate-sorted STAR BAM for one RNA library
#   GTF          concatenated annotation on prefixed parental contigs
#   GENE_PAIRS   TSV pairing parent1 and parent2 gene IDs
#   PARENT1_PREFIX  contig prefix for parent1, provided as an environment variable
#   PARENT2_PREFIX  contig prefix for parent2, provided as an environment variable
#
# Outputs (OUTDIR, default results/counts/rna/${LIBRARY_ID}):
#   gene_counts.tsv
#   filter_summary.tsv
#
# Read-pair assignment is based only on mapped contig prefixes. A pair is
# counted only when both primary mates map to the same haplotype.

set -euo pipefail

SAMTOOLS=${SAMTOOLS:-samtools}
FEATURECOUNTS=${FEATURECOUNTS:-featureCounts}

die() {
    echo "ERROR: $*" >&2
    exit 1
}

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}

require_cmd() {
    command -v "$1" >/dev/null 2>&1 || die "required command not found in PATH: $1"
}

for cmd in "$SAMTOOLS" "$FEATURECOUNTS" awk python3; do
    require_cmd "$cmd"
done

REP=${REP:-1}
SAMPLE=${SAMPLE:-rna}
LIBRARY_ID=${LIBRARY_ID:-${SAMPLE}_rep${REP}}
BAM=${BAM:-results/mapping/rna/${LIBRARY_ID}.Aligned.sortedByCoord.out.bam}
GTF=${GTF:-results/reference/concatenated.gtf}
GENE_PAIRS=${GENE_PAIRS:-gene_pairs.tsv}
OUTDIR=${OUTDIR:-results/counts/rna/${LIBRARY_ID}}
COUNT_TMPDIR=${COUNT_TMPDIR:-${RNA_TMPDIR:-rna_count_tmp/${LIBRARY_ID}}}
THREADS=${THREADS:-4}
MAPQ_MIN=${MAPQ_MIN:-255}
FEATURE_TYPE=${FEATURE_TYPE:-exon}
GROUP_ATTR=${GROUP_ATTR:-gene_id}
STRANDNESS=${STRANDNESS:-0}

[ -n "${PARENT1_PREFIX:-}" ] || die "required environment variable not set: PARENT1_PREFIX"
[ -n "${PARENT2_PREFIX:-}" ] || die "required environment variable not set: PARENT2_PREFIX"
PARENT1_CONTIG_PREFIX="${PARENT1_PREFIX}_"
PARENT2_CONTIG_PREFIX="${PARENT2_PREFIX}_"

log "RNA gene counting: ${LIBRARY_ID}"
echo "  BAM        = $BAM"
echo "  GTF        = $GTF"
echo "  GENE_PAIRS = $GENE_PAIRS"
echo "  PARENT1_PREFIX = $PARENT1_PREFIX"
echo "  PARENT2_PREFIX = $PARENT2_PREFIX"
echo "  OUTDIR     = $OUTDIR"

for f in "$BAM" "$GTF" "$GENE_PAIRS"; do
    [ -f "$f" ] || die "input not found: $f"
done

mkdir -p "$OUTDIR" "$COUNT_TMPDIR"

log "Checking GTF feature type and parent contig prefixes"
n_features=$(awk -F'\t' -v t="$FEATURE_TYPE" '$3==t' "$GTF" | wc -l | tr -d ' ')
[ "$n_features" -gt 0 ] || die "no '${FEATURE_TYPE}' features in $GTF"
n_parent1_features=$(awk -F'\t' -v p="$PARENT1_CONTIG_PREFIX" -v t="$FEATURE_TYPE" '$3==t && index($1,p)==1' "$GTF" | wc -l | tr -d ' ')
n_parent2_features=$(awk -F'\t' -v p="$PARENT2_CONTIG_PREFIX" -v t="$FEATURE_TYPE" '$3==t && index($1,p)==1' "$GTF" | wc -l | tr -d ' ')
[ "$n_parent1_features" -gt 0 ] || die "no '${FEATURE_TYPE}' features on parent1 contigs with prefix ${PARENT1_CONTIG_PREFIX} in $GTF"
[ "$n_parent2_features" -gt 0 ] || die "no '${FEATURE_TYPE}' features on parent2 contigs with prefix ${PARENT2_CONTIG_PREFIX} in $GTF"

log "Filtering BAM to unique, properly paired, primary alignments"
TOTAL_BAM_RECORDS=$("$SAMTOOLS" view -c "$BAM")
"$SAMTOOLS" view -b -q "$MAPQ_MIN" -f 0x2 -F 0x904 "$BAM" \
    > "$COUNT_TMPDIR/filtered.bam"
FILTERED_RECORDS=$("$SAMTOOLS" view -c "$COUNT_TMPDIR/filtered.bam")

log "Name-sorting filtered BAM for read-pair classification"
"$SAMTOOLS" sort -n -@ "$THREADS" -o "$COUNT_TMPDIR/filtered.namesorted.bam" \
    "$COUNT_TMPDIR/filtered.bam"

log "Classifying read pairs by same-haplotype concordance"
: > "$COUNT_TMPDIR/keep_names.txt"
"$SAMTOOLS" view "$COUNT_TMPDIR/filtered.namesorted.bam" \
    | awk -v parent1_prefix="$PARENT1_CONTIG_PREFIX" -v parent2_prefix="$PARENT2_CONTIG_PREFIX" \
          -v keep="$COUNT_TMPDIR/keep_names.txt" \
          -v counts="$COUNT_TMPDIR/filter_counts.tsv" '
    function haplo(contig) {
        if (index(contig, parent1_prefix) == 1) return "parent1"
        if (index(contig, parent2_prefix) == 1) return "parent2"
        return "NA"
    }
    function flush(name, n, h1, h2,    hap) {
        groups++
        if (n != 2) { nontwo++; return }
        if (h1 == "NA" || h2 == "NA") { unknown++; return }
        if (h1 != h2) { cross++; return }
        hap = h1
        same++
        if (hap == "parent1") parent1f++; else parent2f++
        print name >> keep
    }
    {
        name = $1
        contig = $3
        if (name != cur) {
            if (cur != "") flush(cur, cnt, ha, hb)
            cur = name
            cnt = 0
            ha = ""
            hb = ""
        }
        cnt++
        if (cnt == 1) ha = haplo(contig)
        else if (cnt == 2) hb = haplo(contig)
    }
    END {
        if (cur != "") flush(cur, cnt, ha, hb)
        print "read_groups_total\t" groups+0 > counts
        print "same_haplotype_pairs_kept\t" same+0 >> counts
        print "cross_haplotype_pairs_excluded\t" cross+0 >> counts
        print "unknown_haplotype_pairs_excluded\t" unknown+0 >> counts
        print "non_two_record_groups_excluded\t" nontwo+0 >> counts
        print "parent1_fragments\t" parent1f+0 >> counts
        print "parent2_fragments\t" parent2f+0 >> counts
    }'

[ -s "$COUNT_TMPDIR/keep_names.txt" ] || {
    echo "ERROR: no same-haplotype RNA read pairs remain after filtering." >&2
    echo "       Inspect $COUNT_TMPDIR/filter_counts.tsv." >&2
    exit 1
}

log "Extracting same-haplotype pairs"
"$SAMTOOLS" view -b -N "$COUNT_TMPDIR/keep_names.txt" "$COUNT_TMPDIR/filtered.namesorted.bam" \
    > "$COUNT_TMPDIR/concordant.namesorted.bam"
"$SAMTOOLS" sort -@ "$THREADS" -o "$COUNT_TMPDIR/concordant.sorted.bam" \
    "$COUNT_TMPDIR/concordant.namesorted.bam"

log "Splitting concordant BAM by parent prefix"
"$SAMTOOLS" view -h "$COUNT_TMPDIR/concordant.sorted.bam" \
    | awk -v p="$PARENT1_CONTIG_PREFIX" 'BEGIN { OFS="\t" } /^@/ { print; next } index($3, p) == 1 { print }' \
    | "$SAMTOOLS" view -b -o "$COUNT_TMPDIR/parent1.bam" -
"$SAMTOOLS" view -h "$COUNT_TMPDIR/concordant.sorted.bam" \
    | awk -v p="$PARENT2_CONTIG_PREFIX" 'BEGIN { OFS="\t" } /^@/ { print; next } index($3, p) == 1 { print }' \
    | "$SAMTOOLS" view -b -o "$COUNT_TMPDIR/parent2.bam" -

[ "$("$SAMTOOLS" view -c "$COUNT_TMPDIR/parent1.bam")" -gt 0 ] || die "no parent1 alignments after prefix split"
[ "$("$SAMTOOLS" view -c "$COUNT_TMPDIR/parent2.bam")" -gt 0 ] || die "no parent2 alignments after prefix split"

log "Counting genes per haplotype with featureCounts"
"$FEATURECOUNTS" -T "$THREADS" -p --countReadPairs -B -C -s "$STRANDNESS" \
    -t "$FEATURE_TYPE" -g "$GROUP_ATTR" \
    -a "$GTF" -o "$COUNT_TMPDIR/parent1_counts.txt" "$COUNT_TMPDIR/parent1.bam" \
    > "$COUNT_TMPDIR/featureCounts.parent1.log" 2>&1
"$FEATURECOUNTS" -T "$THREADS" -p --countReadPairs -B -C -s "$STRANDNESS" \
    -t "$FEATURE_TYPE" -g "$GROUP_ATTR" \
    -a "$GTF" -o "$COUNT_TMPDIR/parent2_counts.txt" "$COUNT_TMPDIR/parent2.bam" \
    > "$COUNT_TMPDIR/featureCounts.parent2.log" 2>&1
[ -s "$COUNT_TMPDIR/parent1_counts.txt" ] || die "featureCounts produced no parent1 output; see $COUNT_TMPDIR/featureCounts.parent1.log"
[ -s "$COUNT_TMPDIR/parent2_counts.txt" ] || die "featureCounts produced no parent2 output; see $COUNT_TMPDIR/featureCounts.parent2.log"

log "Pairing parent1/parent2 gene counts"
python3 - "$COUNT_TMPDIR/parent1_counts.txt" "$COUNT_TMPDIR/parent2_counts.txt" \
          "$GENE_PAIRS" "$OUTDIR/gene_counts.tsv" <<'PYEOF'
import csv
import sys

parent1_counts_path, parent2_counts_path, pairs_path, output_path = sys.argv[1:5]


def load_featurecounts(path):
    counts = {}
    with open(path) as source:
        for line_number, line in enumerate(source, start=1):
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if not fields or fields[0] == "Geneid":
                continue
            gene_id = fields[0]
            if gene_id in counts:
                raise SystemExit(f"ERROR: duplicate gene ID in {path}: {gene_id}")
            try:
                counts[gene_id] = int(fields[-1])
            except ValueError as error:
                raise SystemExit(
                    f"ERROR: non-integer featureCounts value in {path}:{line_number}"
                ) from error
    return counts


parent1_counts = load_featurecounts(parent1_counts_path)
parent2_counts = load_featurecounts(parent2_counts_path)

with open(pairs_path, newline="") as source:
    reader = csv.reader(source, delimiter="\t")
    try:
        header = next(reader)
    except StopIteration:
        raise SystemExit(f"ERROR: missing header in {pairs_path}")
    expected_header = ["parent1_gene_id", "parent2_gene_id"]
    expected_header_text = "\t".join(expected_header)
    if header != expected_header:
        raise SystemExit(
            f"ERROR: {pairs_path} header must be exactly: {expected_header_text}"
        )
    rows = list(reader)

with open(output_path, "w", newline="") as output:
    writer = csv.writer(output, delimiter="\t", lineterminator="\n")
    writer.writerow(("parent1_gene_id", "parent2_gene_id", "parent1_count", "parent2_count", "total_count"))
    for line_number, row in enumerate(rows, start=2):
        if len(row) != 2:
            raise SystemExit(f"ERROR: expected two columns in {pairs_path}:{line_number}")
        parent1_gene, parent2_gene = row
        if parent1_gene not in parent1_counts:
            raise SystemExit(f"ERROR: parent1 gene absent from featureCounts output: {parent1_gene}")
        if parent2_gene not in parent2_counts:
            raise SystemExit(f"ERROR: parent2 gene absent from featureCounts output: {parent2_gene}")
        parent1_count = parent1_counts[parent1_gene]
        parent2_count = parent2_counts[parent2_gene]
        writer.writerow((parent1_gene, parent2_gene, parent1_count, parent2_count, parent1_count + parent2_count))
PYEOF

{
    echo -e "metric\tvalue"
    echo -e "assay\trna"
    echo -e "sample\t${SAMPLE}"
    echo -e "rep\t${REP}"
    echo -e "library_id\t${LIBRARY_ID}"
    echo -e "bam\t${BAM}"
    echo -e "gtf\t${GTF}"
    echo -e "gene_pairs\t${GENE_PAIRS}"
    echo -e "mapq_min\t${MAPQ_MIN}"
    echo -e "strandness\t${STRANDNESS}"
    echo -e "feature_type\t${FEATURE_TYPE}"
    echo -e "group_attr\t${GROUP_ATTR}"
    echo -e "parent1_prefix\t${PARENT1_PREFIX}"
    echo -e "parent2_prefix\t${PARENT2_PREFIX}"
    echo -e "total_bam_records\t${TOTAL_BAM_RECORDS}"
    echo -e "kept_after_filter_records\t${FILTERED_RECORDS}"
    cat "$COUNT_TMPDIR/filter_counts.tsv"
} > "$OUTDIR/filter_summary.tsv"

show() {
    if command -v column >/dev/null 2>&1; then
        column -t "$1"
    else
        cat "$1"
    fi
}

echo ""
log "RNA gene counting complete: ${OUTDIR}"
echo ""
echo "--- filter_summary.tsv ---"
sed 's/^/  /' "$OUTDIR/filter_summary.tsv"
echo ""
echo "--- gene_counts.tsv ---"
show "$OUTDIR/gene_counts.tsv" | sed 's/^/  /'

rm -rf "$COUNT_TMPDIR"
