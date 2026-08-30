#!/bin/bash
# deduplicate_atac.sh
# Remove paired-end PCR duplicates from one coordinate-sorted ATAC unique BAM.
#
# Required: BAM, OUTPUT_BAM, STATS, TMP_ROOT, and LIBRARY_ID.

set -euo pipefail

SAMTOOLS=${SAMTOOLS:-samtools}

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

require_cmd "$SAMTOOLS"

BAM=${BAM:-}
OUTPUT_BAM=${OUTPUT_BAM:-}
STATS=${STATS:-}
TMP_ROOT=${TMP_ROOT:-}
LIBRARY_ID=${LIBRARY_ID:-}
THREADS=${THREADS:-4}

[ -n "$BAM" ] || die "BAM is required"
[ -n "$OUTPUT_BAM" ] || die "OUTPUT_BAM is required"
[ -n "$STATS" ] || die "STATS is required"
[ -n "$TMP_ROOT" ] || die "TMP_ROOT is required"
[ -n "$LIBRARY_ID" ] || die "LIBRARY_ID is required"
[ -f "$BAM" ] || die "input BAM not found: $BAM"

mkdir -p "$TMP_ROOT"
mkdir -p "$(dirname "$OUTPUT_BAM")"
mkdir -p "$(dirname "$STATS")"

WORKDIR=""
cleanup() {
    if [ -n "$WORKDIR" ] && [ -d "$WORKDIR" ]; then
        rm -rf -- "$WORKDIR"
    fi
}
trap cleanup EXIT
trap 'exit 130' INT
trap 'exit 143' TERM
trap 'exit 129' HUP

WORKDIR=$(mktemp -d "${TMP_ROOT%/}/hyarq_atac_dedup.XXXXXX")
COLLATED_BAM="$WORKDIR/collated.bam"
FIXMATE_BAM="$WORKDIR/fixmate.bam"
COORDINATE_BAM="$WORKDIR/coordinate_sorted.bam"
COLLATE_TMP_PREFIX="$WORKDIR/collate_tmp"
SORT_TMP_PREFIX="$WORKDIR/sort_tmp"

log "ATAC duplicate removal: ${LIBRARY_ID}"
echo "  input  = $BAM"
echo "  output = $OUTPUT_BAM"
echo "  stats  = $STATS"
echo "  temp   = $WORKDIR"

"$SAMTOOLS" collate \
    -@ "$THREADS" \
    -T "$COLLATE_TMP_PREFIX" \
    -o "$COLLATED_BAM" \
    "$BAM"

"$SAMTOOLS" fixmate \
    -m \
    -@ "$THREADS" \
    "$COLLATED_BAM" \
    "$FIXMATE_BAM"

"$SAMTOOLS" sort \
    -@ "$THREADS" \
    -T "$SORT_TMP_PREFIX" \
    -o "$COORDINATE_BAM" \
    "$FIXMATE_BAM"

"$SAMTOOLS" markdup \
    -r \
    -s \
    -f "$STATS" \
    -@ "$THREADS" \
    "$COORDINATE_BAM" \
    "$OUTPUT_BAM"

[ -s "$OUTPUT_BAM" ] || die "samtools markdup produced no BAM for ${LIBRARY_ID}"
[ -s "$STATS" ] || die "samtools markdup produced no statistics for ${LIBRARY_ID}"

"$SAMTOOLS" index \
    -@ "$THREADS" \
    "$OUTPUT_BAM"

[ -s "${OUTPUT_BAM}.bai" ] || die "samtools produced no BAM index for ${LIBRARY_ID}"

input_records=$("$SAMTOOLS" view -c "$BAM")
output_records=$("$SAMTOOLS" view -c "$OUTPUT_BAM")

cleanup
trap - EXIT INT TERM HUP

log "ATAC duplicate removal complete: ${LIBRARY_ID}"
echo "  input alignment records: $input_records"
echo "  output alignment records: $output_records"
echo "  markdup statistics: $STATS"
