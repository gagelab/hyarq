#!/bin/bash
# map_chip_star.sh
# Map one cleaned paired-end ChIP-seq library with STAR and retain MAPQ-255
# alignments.
#
# Required: REFERENCE_ID, LIBRARY_ID, R1, R2, OUT_SAM_MULT_NMAX,
# WIN_ANCHOR_MULTIMAP_NMAX, and OUT_BAM_SORTING_BINS_N.
#
# Outputs (OUTDIR, default results/mapping/chip):
#   ${LIBRARY_ID}.Aligned.sortedByCoord.out.bam
#   ${LIBRARY_ID}.unique.bam (+ .bai)
#   ${LIBRARY_ID}.Log.final.out
#   ${LIBRARY_ID}.mapping_summary.tsv

set -euo pipefail

STAR_BIN=${STAR_BIN:-STAR}
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

for cmd in "$STAR_BIN" "$SAMTOOLS"; do
    require_cmd "$cmd"
done

REFERENCE_ID=${REFERENCE_ID:-}
LIBRARY_ID=${LIBRARY_ID:-}
R1=${R1:-}
R2=${R2:-}
OUT_SAM_MULT_NMAX=${OUT_SAM_MULT_NMAX:-}
WIN_ANCHOR_MULTIMAP_NMAX=${WIN_ANCHOR_MULTIMAP_NMAX:-}
OUT_BAM_SORTING_BINS_N=${OUT_BAM_SORTING_BINS_N:-}

[ -n "$REFERENCE_ID" ] || die "REFERENCE_ID is required"
[ -n "$LIBRARY_ID" ] || die "LIBRARY_ID is required"
[ -n "$R1" ] || die "R1 is required"
[ -n "$R2" ] || die "R2 is required"
[ -n "$OUT_SAM_MULT_NMAX" ] || die "OUT_SAM_MULT_NMAX is required"
[ -n "$WIN_ANCHOR_MULTIMAP_NMAX" ] || die "WIN_ANCHOR_MULTIMAP_NMAX is required"
[ -n "$OUT_BAM_SORTING_BINS_N" ] || die "OUT_BAM_SORTING_BINS_N is required"

INDEX_DIR=${INDEX_DIR:-resources/references/${REFERENCE_ID}/star_index}
OUTDIR=${OUTDIR:-results/mapping/chip}
THREADS=${THREADS:-4}
SAMPLE=${SAMPLE:-chip}
REP=${REP:-1}

for fastq in "$R1" "$R2"; do
    [ -f "$fastq" ] || die "input FASTQ not found: $fastq"
done
[ -d "$INDEX_DIR" ] || die "STAR index not found: $INDEX_DIR"

case "$R1:$R2" in
    *.gz:*.gz) READ_FILES_COMMAND=zcat ;;
    *) READ_FILES_COMMAND=cat ;;
esac
require_cmd "$READ_FILES_COMMAND"

mkdir -p "$OUTDIR"
PREFIX="$OUTDIR/${LIBRARY_ID}."
ALIGNED_BAM="${PREFIX}Aligned.sortedByCoord.out.bam"
UNIQUE_BAM="${PREFIX}unique.bam"
STAR_LOG="${PREFIX}Log.final.out"
STAR_TMP="${PREFIX}_STARtmp"
SUMMARY="$OUTDIR/${LIBRARY_ID}.mapping_summary.tsv"

cleanup_star_tmp() {
    rm -rf -- "$STAR_TMP"
}

# STAR refuses to reuse a stale temporary directory from an interrupted job.
cleanup_star_tmp
trap cleanup_star_tmp EXIT

log "ChIP mapping: ${LIBRARY_ID}"
echo "  R1     = $R1"
echo "  R2     = $R2"
echo "  index  = $INDEX_DIR"
echo "  prefix = $PREFIX"

"$STAR_BIN" \
    --runMode alignReads \
    --runThreadN "$THREADS" \
    --genomeDir "$INDEX_DIR" \
    --readFilesIn "$R1" "$R2" \
    --readFilesCommand "$READ_FILES_COMMAND" \
    --alignIntronMax 1 \
    --outSAMmultNmax "$OUT_SAM_MULT_NMAX" \
    --winAnchorMultimapNmax "$WIN_ANCHOR_MULTIMAP_NMAX" \
    --outBAMsortingBinsN "$OUT_BAM_SORTING_BINS_N" \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix "$PREFIX"

[ -s "$ALIGNED_BAM" ] || die "STAR produced no BAM for ${LIBRARY_ID}"

log "Filtering to unique alignments (MAPQ 255)"
"$SAMTOOLS" view \
    -b \
    -q 255 \
    -F 0x904 \
    -o "$UNIQUE_BAM" \
    "$ALIGNED_BAM"

[ -s "$UNIQUE_BAM" ] || die "samtools produced no unique BAM for ${LIBRARY_ID}"
"$SAMTOOLS" index "$UNIQUE_BAM"

alignment_records=$("$SAMTOOLS" view -c "$ALIGNED_BAM")
unique_records=$("$SAMTOOLS" view -c "$UNIQUE_BAM")
{
    echo -e "metric\tvalue"
    echo -e "reference_id\t${REFERENCE_ID}"
    echo -e "assay\tchip"
    echo -e "sample\t${SAMPLE}"
    echo -e "rep\t${REP}"
    echo -e "library_id\t${LIBRARY_ID}"
    echo -e "r1\t${R1}"
    echo -e "r2\t${R2}"
    echo -e "star_index\t${INDEX_DIR}"
    echo -e "threads\t${THREADS}"
    echo -e "read_files_command\t${READ_FILES_COMMAND}"
    echo -e "align_intron_max\t1"
    echo -e "out_sam_mult_nmax\t${OUT_SAM_MULT_NMAX}"
    echo -e "win_anchor_multimap_nmax\t${WIN_ANCHOR_MULTIMAP_NMAX}"
    echo -e "out_bam_sorting_bins_n\t${OUT_BAM_SORTING_BINS_N}"
    echo -e "mapq_min\t255"
    echo -e "alignment_records\t${alignment_records}"
    echo -e "unique_records\t${unique_records}"
    echo -e "aligned_bam\t${ALIGNED_BAM}"
    echo -e "unique_bam\t${UNIQUE_BAM}"
} > "$SUMMARY"

# Only Log.final.out is retained from STAR's auxiliary alignment files.
rm -f -- \
    "${PREFIX}Log.out" \
    "${PREFIX}Log.progress.out" \
    "${PREFIX}SJ.out.tab"
cleanup_star_tmp
trap - EXIT

echo ""
log "ChIP mapping complete: $UNIQUE_BAM"
if [ -f "$STAR_LOG" ]; then
    echo "--- STAR mapping summary (${LIBRARY_ID}) ---"
    grep -E "Number of input reads|Uniquely mapped reads %|% of reads mapped to multiple loci|% of reads unmapped" "$STAR_LOG" \
        | sed 's/^[[:space:]]*/  /'
fi
echo "  alignment records in BAM: $alignment_records"
echo "  unique (MAPQ 255) records: $unique_records"
echo "  mapping summary: $SUMMARY"
