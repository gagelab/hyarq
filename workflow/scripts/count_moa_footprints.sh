#!/bin/bash
# count_moa_footprints.sh
# Count MOA-seq merged-read midpoints over paired P1/P2 footprints on the
# concatenated reference.
#
# Required inputs:
#   BAM               single-end merged-read BAM from map_moa_star.sh
#   FOOTPRINTS        BED with prefixed P1/P2 footprint coordinates and names
#                     in col 4
#   FOOTPRINT_PAIRS   TSV pairing P1 and P2 footprint names
#
# Outputs (OUTDIR, default results/counts/moa/${LIBRARY_ID}):
#   footprint_counts.tsv
#   filter_summary.tsv
#
# Assignment is based only on mapped contig prefixes. A midpoint is counted only
# when it overlaps exactly one footprint on its mapped parent contig.

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
SAMPLE=${SAMPLE:-moa}
LIBRARY_ID=${LIBRARY_ID:-${SAMPLE}_rep${REP}}
BAM=${BAM:-results/mapping/moa/${LIBRARY_ID}.unique.bam}
FOOTPRINTS=${FOOTPRINTS:-moa_footprints.bed}
FOOTPRINT_PAIRS=${FOOTPRINT_PAIRS:-${REGION_PAIRS:-moa_region_pairs.tsv}}
OUTDIR=${OUTDIR:-results/counts/moa/${LIBRARY_ID}}
COUNT_TMPDIR=${COUNT_TMPDIR:-${MOA_TMPDIR:-moa_count_tmp/${LIBRARY_ID}}}
P1_PREFIX=${P1_PREFIX:-P1_}
P2_PREFIX=${P2_PREFIX:-P2_}
THREADS=${THREADS:-4}
MAPQ_MIN=${MAPQ_MIN:-255}

log "MOA footprint counting: ${LIBRARY_ID}"
echo "  BAM             = $BAM"
echo "  FOOTPRINTS      = $FOOTPRINTS"
echo "  FOOTPRINT_PAIRS = $FOOTPRINT_PAIRS"
echo "  OUTDIR          = $OUTDIR"

for f in "$BAM" "$FOOTPRINTS" "$FOOTPRINT_PAIRS"; do
    [ -f "$f" ] || die "input not found: $f"
done

mkdir -p "$OUTDIR" "$COUNT_TMPDIR"

log "Filtering BAM to unique primary alignments"
TOTAL_BAM_RECORDS=$("$SAMTOOLS" view -c "$BAM")
"$SAMTOOLS" view -@ "$THREADS" -b -q "$MAPQ_MIN" -F 0x904 "$BAM" \
    > "$COUNT_TMPDIR/filtered.bam"
FILTERED_RECORDS=$("$SAMTOOLS" view -c "$COUNT_TMPDIR/filtered.bam")

log "Converting filtered BAM to BED"
"$BEDTOOLS" bamtobed -i "$COUNT_TMPDIR/filtered.bam" \
    | sort -k1,1 -k2,2n > "$COUNT_TMPDIR/filtered_reads.bed"
[ -s "$COUNT_TMPDIR/filtered_reads.bed" ] || die "no BED records remain after MOA filtering"

log "Assigning 1 bp midpoints to paired footprints"
python3 - "$COUNT_TMPDIR/filtered_reads.bed" "$FOOTPRINTS" "$FOOTPRINT_PAIRS" \
          "$P1_PREFIX" "$P2_PREFIX" "$FILTERED_RECORDS" \
          "$COUNT_TMPDIR/midpoints.bed" "$COUNT_TMPDIR/filter_counts.tsv" \
          "$OUTDIR/footprint_counts.tsv" <<'PYEOF'
import csv
import sys

(
    reads_bed,
    footprints_bed,
    pairs_path,
    p1_prefix,
    p2_prefix,
    filtered_records,
    midpoints_path,
    summary_path,
    output_path,
) = sys.argv[1:10]
filtered_records = int(filtered_records)


def haplotype(contig):
    if contig.startswith(p1_prefix):
        return "P1"
    if contig.startswith(p2_prefix):
        return "P2"
    return None


def choose_column(header, candidates, path):
    for candidate in candidates:
        if candidate in header:
            return candidate
    raise SystemExit(
        f"ERROR: missing required column in {path}; expected one of: {', '.join(candidates)}"
    )


footprints_by_contig = {}
footprint_haplotype = {}
footprint_counts = {}
with open(footprints_bed) as source:
    for line_number, line in enumerate(source, start=1):
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 4:
            raise SystemExit(f"ERROR: malformed footprint BED row in {footprints_bed}:{line_number}")
        contig, start_text, end_text, name = fields[:4]
        try:
            start = int(start_text)
            end = int(end_text)
        except ValueError as error:
            raise SystemExit(
                f"ERROR: non-integer footprint coordinates in {footprints_bed}:{line_number}"
            ) from error
        if start < 0 or end <= start:
            raise SystemExit(f"ERROR: invalid footprint interval in {footprints_bed}:{line_number}")
        hap = haplotype(contig)
        if hap is None:
            raise SystemExit(
                f"ERROR: footprint contig lacks {p1_prefix} or {p2_prefix}: {contig}"
            )
        if name in footprint_haplotype:
            raise SystemExit(f"ERROR: duplicate footprint name in {footprints_bed}: {name}")
        footprint_haplotype[name] = hap
        footprint_counts[name] = 0
        footprints_by_contig.setdefault(contig, []).append((start, end, name))

for contig in footprints_by_contig:
    footprints_by_contig[contig].sort()

with open(pairs_path, newline="") as source:
    reader = csv.DictReader(source, delimiter="\t")
    if reader.fieldnames is None:
        raise SystemExit(f"ERROR: missing header in {pairs_path}")
    id_col = choose_column(reader.fieldnames, ("footprint_pair_id", "region_id", "pair_id", "feature_id"), pairs_path)
    p1_col = choose_column(reader.fieldnames, ("P1_footprint", "P1_feature"), pairs_path)
    p2_col = choose_column(reader.fieldnames, ("P2_footprint", "P2_feature"), pairs_path)
    pair_rows = list(reader)

pairs = []
seen_ids = set()
for row in pair_rows:
    pair_id = row[id_col]
    p1_footprint = row[p1_col]
    p2_footprint = row[p2_col]
    if not pair_id or pair_id in seen_ids:
        raise SystemExit(f"ERROR: blank or duplicate footprint pair ID in {pairs_path}: {pair_id!r}")
    if footprint_haplotype.get(p1_footprint) != "P1":
        raise SystemExit(f"ERROR: invalid P1 footprint for {pair_id}: {p1_footprint}")
    if footprint_haplotype.get(p2_footprint) != "P2":
        raise SystemExit(f"ERROR: invalid P2 footprint for {pair_id}: {p2_footprint}")
    seen_ids.add(pair_id)
    pairs.append((pair_id, p1_footprint, p2_footprint))

counts = {
    "total_midpoints": 0,
    "P1_midpoints": 0,
    "P2_midpoints": 0,
    "unrecognized_contig_midpoints": 0,
    "assigned_to_one_footprint": 0,
    "outside_all_footprints": 0,
    "multi_footprint_ambiguous": 0,
    "P1_assigned": 0,
    "P2_assigned": 0,
}

with open(reads_bed) as source, open(midpoints_path, "w") as midpoint_output:
    for line_number, line in enumerate(source, start=1):
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 4:
            raise SystemExit(f"ERROR: malformed BED row in {reads_bed}:{line_number}")
        contig, start_text, end_text, read_name = fields[:4]
        try:
            start = int(start_text)
            end = int(end_text)
        except ValueError as error:
            raise SystemExit(
                f"ERROR: non-integer BED coordinates in {reads_bed}:{line_number}"
            ) from error
        if start < 0 or end <= start:
            raise SystemExit(f"ERROR: invalid BED interval in {reads_bed}:{line_number}")

        midpoint = (start + end) // 2
        hap = haplotype(contig)
        counts["total_midpoints"] += 1
        midpoint_output.write(f"{contig}\t{midpoint}\t{midpoint + 1}\t{read_name}\t{hap or 'NA'}\n")

        if hap is None:
            counts["unrecognized_contig_midpoints"] += 1
            continue
        counts[f"{hap}_midpoints"] += 1

        candidates = [
            name for fp_start, fp_end, name in footprints_by_contig.get(contig, [])
            if fp_start <= midpoint < fp_end
        ]
        if len(candidates) == 1:
            footprint_name = candidates[0]
            footprint_counts[footprint_name] += 1
            counts["assigned_to_one_footprint"] += 1
            counts[f"{hap}_assigned"] += 1
        elif not candidates:
            counts["outside_all_footprints"] += 1
        else:
            counts["multi_footprint_ambiguous"] += 1

if counts["total_midpoints"] != filtered_records:
    raise SystemExit(
        f"ERROR: midpoint count ({counts['total_midpoints']}) does not equal "
        f"filtered BAM records ({filtered_records})"
    )

with open(output_path, "w", newline="") as output:
    writer = csv.writer(output, delimiter="\t", lineterminator="\n")
    writer.writerow(
        ("footprint_pair_id", "P1_footprint", "P2_footprint", "P1_count", "P2_count", "total_count")
    )
    for pair_id, p1_footprint, p2_footprint in pairs:
        p1_count = footprint_counts[p1_footprint]
        p2_count = footprint_counts[p2_footprint]
        writer.writerow((pair_id, p1_footprint, p2_footprint, p1_count, p2_count, p1_count + p2_count))

with open(summary_path, "w") as output:
    for metric, value in counts.items():
        output.write(f"{metric}\t{value}\n")
PYEOF

{
    echo -e "metric\tvalue"
    echo -e "assay\tmoa"
    echo -e "sample\t${SAMPLE}"
    echo -e "rep\t${REP}"
    echo -e "library_id\t${LIBRARY_ID}"
    echo -e "bam\t${BAM}"
    echo -e "footprints\t${FOOTPRINTS}"
    echo -e "footprint_pairs\t${FOOTPRINT_PAIRS}"
    echo -e "mapq_min\t${MAPQ_MIN}"
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
log "MOA footprint counting complete: ${OUTDIR}"
echo ""
echo "--- filter_summary.tsv ---"
sed 's/^/  /' "$OUTDIR/filter_summary.tsv"
echo ""
echo "--- footprint_counts.tsv ---"
show "$OUTDIR/footprint_counts.tsv" | sed 's/^/  /'

rm -rf "$COUNT_TMPDIR"
