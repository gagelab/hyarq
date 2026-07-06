#!/bin/bash
# count_rna_genes.sh
# Count RNA-seq fragments over paired P1/P2 genes on the concatenated reference.
#
# Required inputs:
#   BAM          coordinate-sorted STAR BAM for one RNA library
#   GTF          concatenated annotation on prefixed P1/P2 contigs
#   GENE_PAIRS   TSV pairing P1 and P2 gene IDs
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
P1_PREFIX=${P1_PREFIX:-P1_}
P2_PREFIX=${P2_PREFIX:-P2_}
THREADS=${THREADS:-4}
MAPQ_MIN=${MAPQ_MIN:-255}
FEATURE_TYPE=${FEATURE_TYPE:-exon}
GROUP_ATTR=${GROUP_ATTR:-gene_id}
STRANDNESS=${STRANDNESS:-0}

log "RNA gene counting: ${LIBRARY_ID}"
echo "  BAM        = $BAM"
echo "  GTF        = $GTF"
echo "  GENE_PAIRS = $GENE_PAIRS"
echo "  OUTDIR     = $OUTDIR"

for f in "$BAM" "$GTF" "$GENE_PAIRS"; do
    [ -f "$f" ] || die "input not found: $f"
done

mkdir -p "$OUTDIR" "$COUNT_TMPDIR"

log "Checking GTF feature type and parent contig prefixes"
n_features=$(awk -F'\t' -v t="$FEATURE_TYPE" '$3==t' "$GTF" | wc -l | tr -d ' ')
[ "$n_features" -gt 0 ] || die "no '${FEATURE_TYPE}' features in $GTF"
n_p1_features=$(awk -F'\t' -v p="$P1_PREFIX" -v t="$FEATURE_TYPE" '$3==t && index($1,p)==1' "$GTF" | wc -l | tr -d ' ')
n_p2_features=$(awk -F'\t' -v p="$P2_PREFIX" -v t="$FEATURE_TYPE" '$3==t && index($1,p)==1' "$GTF" | wc -l | tr -d ' ')
[ "$n_p1_features" -gt 0 ] || die "no '${FEATURE_TYPE}' features on ${P1_PREFIX} contigs in $GTF"
[ "$n_p2_features" -gt 0 ] || die "no '${FEATURE_TYPE}' features on ${P2_PREFIX} contigs in $GTF"

log "Filtering BAM to unique, properly paired, primary alignments"
TOTAL_BAM_RECORDS=$("$SAMTOOLS" view -c "$BAM")
"$SAMTOOLS" view -b -q "$MAPQ_MIN" -f 0x2 -F 0x904 "$BAM" \
    > "$COUNT_TMPDIR/filtered.bam"
FILTERED_RECORDS=$("$SAMTOOLS" view -c "$COUNT_TMPDIR/filtered.bam")

log "Name-sorting filtered BAM for read-pair classification"
"$SAMTOOLS" sort -n -@ "$THREADS" -o "$COUNT_TMPDIR/filtered.namesorted.bam" \
    "$COUNT_TMPDIR/filtered.bam"

log "Classifying read pairs by same-haplotype concordance"
"$SAMTOOLS" view "$COUNT_TMPDIR/filtered.namesorted.bam" \
    | awk -v p1="$P1_PREFIX" -v p2="$P2_PREFIX" \
          -v keep="$COUNT_TMPDIR/keep_names.txt" \
          -v counts="$COUNT_TMPDIR/filter_counts.tsv" '
    function haplo(contig) {
        if (index(contig, p1) == 1) return "P1"
        if (index(contig, p2) == 1) return "P2"
        return "NA"
    }
    function flush(name, n, h1, h2,    hap) {
        groups++
        if (n != 2) { nontwo++; return }
        if (h1 == "NA" || h2 == "NA") { unknown++; return }
        if (h1 != h2) { cross++; return }
        hap = h1
        same++
        if (hap == "P1") p1f++; else p2f++
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
        print "P1_fragments\t" p1f+0 >> counts
        print "P2_fragments\t" p2f+0 >> counts
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
    | awk -v p="$P1_PREFIX" 'BEGIN { OFS="\t" } /^@/ { print; next } index($3, p) == 1 { print }' \
    | "$SAMTOOLS" view -b -o "$COUNT_TMPDIR/P1.bam" -
"$SAMTOOLS" view -h "$COUNT_TMPDIR/concordant.sorted.bam" \
    | awk -v p="$P2_PREFIX" 'BEGIN { OFS="\t" } /^@/ { print; next } index($3, p) == 1 { print }' \
    | "$SAMTOOLS" view -b -o "$COUNT_TMPDIR/P2.bam" -

[ "$("$SAMTOOLS" view -c "$COUNT_TMPDIR/P1.bam")" -gt 0 ] || die "no P1 alignments after prefix split"
[ "$("$SAMTOOLS" view -c "$COUNT_TMPDIR/P2.bam")" -gt 0 ] || die "no P2 alignments after prefix split"

log "Counting genes per haplotype with featureCounts"
"$FEATURECOUNTS" -T "$THREADS" -p --countReadPairs -B -C -s "$STRANDNESS" \
    -t "$FEATURE_TYPE" -g "$GROUP_ATTR" \
    -a "$GTF" -o "$COUNT_TMPDIR/P1_counts.txt" "$COUNT_TMPDIR/P1.bam" \
    > "$COUNT_TMPDIR/featureCounts.P1.log" 2>&1
"$FEATURECOUNTS" -T "$THREADS" -p --countReadPairs -B -C -s "$STRANDNESS" \
    -t "$FEATURE_TYPE" -g "$GROUP_ATTR" \
    -a "$GTF" -o "$COUNT_TMPDIR/P2_counts.txt" "$COUNT_TMPDIR/P2.bam" \
    > "$COUNT_TMPDIR/featureCounts.P2.log" 2>&1
[ -s "$COUNT_TMPDIR/P1_counts.txt" ] || die "featureCounts produced no P1 output; see $COUNT_TMPDIR/featureCounts.P1.log"
[ -s "$COUNT_TMPDIR/P2_counts.txt" ] || die "featureCounts produced no P2 output; see $COUNT_TMPDIR/featureCounts.P2.log"

log "Pairing P1/P2 gene counts"
python3 - "$COUNT_TMPDIR/P1_counts.txt" "$COUNT_TMPDIR/P2_counts.txt" \
          "$GENE_PAIRS" "$OUTDIR/gene_counts.tsv" <<'PYEOF'
import csv
import sys

p1_counts_path, p2_counts_path, pairs_path, output_path = sys.argv[1:5]


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


def choose_column(header, candidates, path):
    for candidate in candidates:
        if candidate in header:
            return candidate
    raise SystemExit(
        f"ERROR: missing required column in {path}; expected one of: {', '.join(candidates)}"
    )


p1_counts = load_featurecounts(p1_counts_path)
p2_counts = load_featurecounts(p2_counts_path)

with open(pairs_path, newline="") as source:
    reader = csv.DictReader(source, delimiter="\t")
    if reader.fieldnames is None:
        raise SystemExit(f"ERROR: missing header in {pairs_path}")
    id_col = choose_column(reader.fieldnames, ("gene_pair_id", "pair_id", "gene_id"), pairs_path)
    p1_col = choose_column(reader.fieldnames, ("P1_gene_id", "P1_gene", "P1_feature"), pairs_path)
    p2_col = choose_column(reader.fieldnames, ("P2_gene_id", "P2_gene", "P2_feature"), pairs_path)
    rows = list(reader)

seen_ids = set()
with open(output_path, "w", newline="") as output:
    writer = csv.writer(output, delimiter="\t", lineterminator="\n")
    writer.writerow(("gene_pair_id", "P1_gene_id", "P2_gene_id", "P1_count", "P2_count", "total_count"))
    for row in rows:
        pair_id = row[id_col]
        p1_gene = row[p1_col]
        p2_gene = row[p2_col]
        if not pair_id or pair_id in seen_ids:
            raise SystemExit(f"ERROR: blank or duplicate gene pair ID in {pairs_path}: {pair_id!r}")
        if p1_gene not in p1_counts:
            raise SystemExit(f"ERROR: P1 gene absent from featureCounts output: {p1_gene}")
        if p2_gene not in p2_counts:
            raise SystemExit(f"ERROR: P2 gene absent from featureCounts output: {p2_gene}")
        seen_ids.add(pair_id)
        p1_count = p1_counts[p1_gene]
        p2_count = p2_counts[p2_gene]
        writer.writerow((pair_id, p1_gene, p2_gene, p1_count, p2_count, p1_count + p2_count))
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
    echo -e "p1_prefix\t${P1_PREFIX}"
    echo -e "p2_prefix\t${P2_PREFIX}"
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
