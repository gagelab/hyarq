#!/bin/bash
# count_atac_regions.sh
# Count ATAC-seq fragments over paired P1/P2 regions on the concatenated
# reference.
#
# Required inputs:
#   BAM            coordinate-sorted STAR BAM for one ATAC library
#   REGIONS        BED with prefixed P1/P2 region coordinates and names in col 4
#   REGION_PAIRS   TSV pairing P1 and P2 region names
#
# Outputs (OUTDIR, default results/counts/atac/${LIBRARY_ID}):
#   region_counts.tsv
#   filter_summary.tsv
#
# Fragment assignment is based only on mapped contig prefixes. A fragment is
# counted only when both mates form one same-contig, same-haplotype interval.

set -euo pipefail

SAMTOOLS=${SAMTOOLS:-samtools}
BEDTOOLS=${BEDTOOLS:-bedtools}

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

for cmd in "$SAMTOOLS" "$BEDTOOLS" awk python3 sort; do
    require_cmd "$cmd"
done

REP=${REP:-1}
SAMPLE=${SAMPLE:-atac}
LIBRARY_ID=${LIBRARY_ID:-${SAMPLE}_rep${REP}}
BAM=${BAM:-results/mapping/atac/${LIBRARY_ID}.Aligned.sortedByCoord.out.bam}
REGIONS=${REGIONS:-${PEAKS:-atac_regions.bed}}
REGION_PAIRS=${REGION_PAIRS:-region_pairs.tsv}
OUTDIR=${OUTDIR:-results/counts/atac/${LIBRARY_ID}}
COUNT_TMPDIR=${COUNT_TMPDIR:-${ATAC_TMPDIR:-atac_count_tmp/${LIBRARY_ID}}}
P1_PREFIX=${P1_PREFIX:-P1_}
P2_PREFIX=${P2_PREFIX:-P2_}
THREADS=${THREADS:-4}
MAPQ_MIN=${MAPQ_MIN:-255}
MIN_OVERLAP_FRAC=${MIN_OVERLAP_FRAC:-}

log "ATAC region counting: ${LIBRARY_ID}"
echo "  BAM          = $BAM"
echo "  REGIONS      = $REGIONS"
echo "  REGION_PAIRS = $REGION_PAIRS"
echo "  OUTDIR       = $OUTDIR"

for f in "$BAM" "$REGIONS" "$REGION_PAIRS"; do
    [ -f "$f" ] || die "input not found: $f"
done

mkdir -p "$OUTDIR" "$COUNT_TMPDIR"

INTERSECT_F=()
if [ -n "$MIN_OVERLAP_FRAC" ]; then
    INTERSECT_F=(-f "$MIN_OVERLAP_FRAC")
fi

log "Filtering BAM to unique, properly paired, primary alignments"
TOTAL_BAM_RECORDS=$("$SAMTOOLS" view -c "$BAM")
"$SAMTOOLS" view -b -q "$MAPQ_MIN" -f 0x2 -F 0x904 "$BAM" \
    > "$COUNT_TMPDIR/filtered.bam"
FILTERED_RECORDS=$("$SAMTOOLS" view -c "$COUNT_TMPDIR/filtered.bam")

log "Name-sorting filtered BAM for BEDPE conversion"
"$SAMTOOLS" sort -n -@ "$THREADS" -o "$COUNT_TMPDIR/filtered.namesorted.bam" \
    "$COUNT_TMPDIR/filtered.bam"

log "Converting paired alignments to BEDPE"
"$BEDTOOLS" bamtobed -bedpe -i "$COUNT_TMPDIR/filtered.namesorted.bam" \
    > "$COUNT_TMPDIR/pairs.bedpe" 2> "$COUNT_TMPDIR/bamtobed.stderr.log"
[ -s "$COUNT_TMPDIR/pairs.bedpe" ] || {
    echo "ERROR: BEDPE conversion produced no records." >&2
    echo "       Check $COUNT_TMPDIR/bamtobed.stderr.log and the filtered BAM." >&2
    exit 1
}

log "Building same-haplotype fragment intervals"
python3 - "$COUNT_TMPDIR/pairs.bedpe" "$P1_PREFIX" "$P2_PREFIX" \
          "$COUNT_TMPDIR/fragments.bed" "$COUNT_TMPDIR/filter_counts.tsv" <<'PYEOF'
import sys

bedpe, p1_prefix, p2_prefix, fragments_path, counts_path = sys.argv[1:6]


def haplotype(contig):
    if contig.startswith(p1_prefix):
        return "P1"
    if contig.startswith(p2_prefix):
        return "P2"
    return None


counts = {
    "total_bedpe_pairs": 0,
    "same_haplotype_fragments_kept": 0,
    "cross_contig_or_unmapped_pairs_excluded": 0,
    "unknown_haplotype_pairs_excluded": 0,
    "cross_haplotype_pairs_excluded": 0,
    "P1_fragments": 0,
    "P2_fragments": 0,
}

with open(bedpe) as source, open(fragments_path, "w") as output:
    for line_number, line in enumerate(source, start=1):
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 7:
            continue
        chrom1, start1, end1, chrom2, start2, end2, name = fields[:7]
        counts["total_bedpe_pairs"] += 1

        if chrom1 != chrom2 or start1 == "-1" or start2 == "-1":
            counts["cross_contig_or_unmapped_pairs_excluded"] += 1
            continue

        hap1 = haplotype(chrom1)
        hap2 = haplotype(chrom2)
        if hap1 is None or hap2 is None:
            counts["unknown_haplotype_pairs_excluded"] += 1
            continue
        if hap1 != hap2:
            counts["cross_haplotype_pairs_excluded"] += 1
            continue

        try:
            start = min(int(start1), int(start2))
            end = max(int(end1), int(end2))
        except ValueError as error:
            raise SystemExit(
                f"ERROR: non-integer BEDPE coordinates at {bedpe}:{line_number}"
            ) from error
        if start < 0 or end <= start:
            counts["cross_contig_or_unmapped_pairs_excluded"] += 1
            continue

        output.write(f"{chrom1}\t{start}\t{end}\t{name}\t{hap1}\n")
        counts["same_haplotype_fragments_kept"] += 1
        counts[f"{hap1}_fragments"] += 1

with open(counts_path, "w") as output:
    for metric, value in counts.items():
        output.write(f"{metric}\t{value}\n")
PYEOF

[ -s "$COUNT_TMPDIR/fragments.bed" ] || {
    echo "ERROR: no same-haplotype fragments remain after filtering." >&2
    echo "       Inspect $COUNT_TMPDIR/filter_counts.tsv and $COUNT_TMPDIR/pairs.bedpe." >&2
    exit 1
}

awk -F'\t' '$5=="P1"' "$COUNT_TMPDIR/fragments.bed" \
    | sort -k1,1 -k2,2n > "$COUNT_TMPDIR/P1_fragments.bed"
awk -F'\t' '$5=="P2"' "$COUNT_TMPDIR/fragments.bed" \
    | sort -k1,1 -k2,2n > "$COUNT_TMPDIR/P2_fragments.bed"
awk -v prefix="$P1_PREFIX" 'BEGIN { FS=OFS="\t" } $1 ~ "^" prefix' "$REGIONS" \
    | sort -k1,1 -k2,2n > "$COUNT_TMPDIR/P1_regions.bed"
awk -v prefix="$P2_PREFIX" 'BEGIN { FS=OFS="\t" } $1 ~ "^" prefix' "$REGIONS" \
    | sort -k1,1 -k2,2n > "$COUNT_TMPDIR/P2_regions.bed"

for f in "$COUNT_TMPDIR/P1_fragments.bed" "$COUNT_TMPDIR/P2_fragments.bed" \
         "$COUNT_TMPDIR/P1_regions.bed" "$COUNT_TMPDIR/P2_regions.bed"; do
    [ -s "$f" ] || die "expected non-empty split file: $f"
done

log "Counting fragment overlaps per ATAC region"
"$BEDTOOLS" intersect "${INTERSECT_F[@]}" \
    -a "$COUNT_TMPDIR/P1_regions.bed" -b "$COUNT_TMPDIR/P1_fragments.bed" -c \
    > "$COUNT_TMPDIR/P1_region_counts.bed"
"$BEDTOOLS" intersect "${INTERSECT_F[@]}" \
    -a "$COUNT_TMPDIR/P2_regions.bed" -b "$COUNT_TMPDIR/P2_fragments.bed" -c \
    > "$COUNT_TMPDIR/P2_region_counts.bed"

log "Pairing P1/P2 region counts"
python3 - "$COUNT_TMPDIR/P1_region_counts.bed" "$COUNT_TMPDIR/P2_region_counts.bed" \
          "$REGION_PAIRS" "$OUTDIR/region_counts.tsv" <<'PYEOF'
import csv
import sys

p1_bed, p2_bed, pairs_path, output_path = sys.argv[1:5]


def read_count_bed(path):
    counts = {}
    with open(path) as source:
        for line_number, line in enumerate(source, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                raise SystemExit(f"ERROR: malformed BED count row in {path}:{line_number}")
            feature_name = fields[3]
            if feature_name in counts:
                raise SystemExit(f"ERROR: duplicate feature name in {path}: {feature_name}")
            try:
                counts[feature_name] = int(fields[-1])
            except ValueError as error:
                raise SystemExit(
                    f"ERROR: non-integer overlap count in {path}:{line_number}"
                ) from error
    return counts


def choose_column(header, candidates, path):
    for candidate in candidates:
        if candidate in header:
            return candidate
    raise SystemExit(
        f"ERROR: missing required column in {path}; expected one of: {', '.join(candidates)}"
    )


p1_counts = read_count_bed(p1_bed)
p2_counts = read_count_bed(p2_bed)

with open(pairs_path, newline="") as source:
    reader = csv.DictReader(source, delimiter="\t")
    if reader.fieldnames is None:
        raise SystemExit(f"ERROR: missing header in {pairs_path}")
    id_col = choose_column(reader.fieldnames, ("region_id", "pair_id", "feature_id"), pairs_path)
    p1_col = choose_column(reader.fieldnames, ("P1_region", "P1_peak", "P1_feature"), pairs_path)
    p2_col = choose_column(reader.fieldnames, ("P2_region", "P2_peak", "P2_feature"), pairs_path)
    rows = list(reader)

seen_ids = set()
with open(output_path, "w", newline="") as output:
    writer = csv.writer(output, delimiter="\t", lineterminator="\n")
    writer.writerow(("region_id", "P1_region", "P2_region", "P1_count", "P2_count", "total_count"))
    for row in rows:
        region_id = row[id_col]
        p1_region = row[p1_col]
        p2_region = row[p2_col]
        if not region_id or region_id in seen_ids:
            raise SystemExit(f"ERROR: blank or duplicate region ID in {pairs_path}: {region_id!r}")
        if p1_region not in p1_counts:
            raise SystemExit(f"ERROR: P1 region absent from P1 BED counts: {p1_region}")
        if p2_region not in p2_counts:
            raise SystemExit(f"ERROR: P2 region absent from P2 BED counts: {p2_region}")
        seen_ids.add(region_id)
        p1_count = p1_counts[p1_region]
        p2_count = p2_counts[p2_region]
        writer.writerow((region_id, p1_region, p2_region, p1_count, p2_count, p1_count + p2_count))
PYEOF

{
    echo -e "metric\tvalue"
    echo -e "assay\tatac"
    echo -e "sample\t${SAMPLE}"
    echo -e "rep\t${REP}"
    echo -e "library_id\t${LIBRARY_ID}"
    echo -e "bam\t${BAM}"
    echo -e "regions\t${REGIONS}"
    echo -e "region_pairs\t${REGION_PAIRS}"
    echo -e "mapq_min\t${MAPQ_MIN}"
    echo -e "min_overlap_frac\t${MIN_OVERLAP_FRAC:-default_1bp}"
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
log "ATAC region counting complete: ${OUTDIR}"
echo ""
echo "--- filter_summary.tsv ---"
sed 's/^/  /' "$OUTDIR/filter_summary.tsv"
echo ""
echo "--- region_counts.tsv ---"
show "$OUTDIR/region_counts.tsv" | sed 's/^/  /'

rm -rf "$COUNT_TMPDIR"
