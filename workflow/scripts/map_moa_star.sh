#!/bin/bash
# map_moa_star.sh
# Map one merged single-end MOA-seq library with STAR and retain unique primary
# alignments.
#
# Required: REFERENCE_ID, LIBRARY_ID, and MERGED.
#
# Outputs (OUTDIR, default results/mapping/moa):
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
MERGED=${MERGED:-}

[ -n "$REFERENCE_ID" ] || die "REFERENCE_ID is required"
[ -n "$LIBRARY_ID" ] || die "LIBRARY_ID is required"
[ -n "$MERGED" ] || die "MERGED is required"

INDEX_DIR=${INDEX_DIR:-resources/references/${REFERENCE_ID}/star_index}
OUTDIR=${OUTDIR:-results/mapping/moa}
THREADS=${THREADS:-4}
SAMPLE=${SAMPLE:-moa}
REP=${REP:-1}

ALIGN_INTRON_MAX=${ALIGN_INTRON_MAX:-1}
OUT_SAM_MULT_NMAX=${OUT_SAM_MULT_NMAX:-2}
WIN_ANCHOR_MULTIMAP_NMAX=${WIN_ANCHOR_MULTIMAP_NMAX:-100}
OUT_BAM_SORTING_BINS_N=${OUT_BAM_SORTING_BINS_N:-5}
MAPQ_MIN=${MAPQ_MIN:-255}

[ -f "$MERGED" ] || die "merged FASTQ not found: $MERGED"
[ -d "$INDEX_DIR" ] || die "STAR index not found: $INDEX_DIR"

case "$MERGED" in
    *.gz) READ_FILES_COMMAND=zcat ;;
    *) READ_FILES_COMMAND=cat ;;
esac
require_cmd "$READ_FILES_COMMAND"

mkdir -p "$OUTDIR"
PREFIX="$OUTDIR/${LIBRARY_ID}."
ALIGNED_BAM="${PREFIX}Aligned.sortedByCoord.out.bam"
UNIQUE_BAM="${PREFIX}unique.bam"
SUMMARY="$OUTDIR/${LIBRARY_ID}.mapping_summary.tsv"

log "MOA mapping: ${LIBRARY_ID}"
echo "  merged = $MERGED"
echo "  index  = $INDEX_DIR"
echo "  prefix = $PREFIX"

"$STAR_BIN" \
    --runMode alignReads \
    --runThreadN "$THREADS" \
    --genomeDir "$INDEX_DIR" \
    --readFilesIn "$MERGED" \
    --readFilesCommand "$READ_FILES_COMMAND" \
    --alignIntronMax "$ALIGN_INTRON_MAX" \
    --outSAMmultNmax "$OUT_SAM_MULT_NMAX" \
    --winAnchorMultimapNmax "$WIN_ANCHOR_MULTIMAP_NMAX" \
    --outBAMsortingBinsN "$OUT_BAM_SORTING_BINS_N" \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix "$PREFIX"

[ -s "$ALIGNED_BAM" ] || die "STAR produced no BAM for ${LIBRARY_ID}"

log "Filtering to unique primary alignments (MAPQ ${MAPQ_MIN})"
"$SAMTOOLS" view \
    -b \
    -q "$MAPQ_MIN" \
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
    echo -e "assay\tmoa"
    echo -e "sample\t${SAMPLE}"
    echo -e "rep\t${REP}"
    echo -e "library_id\t${LIBRARY_ID}"
    echo -e "merged_fastq\t${MERGED}"
    echo -e "star_index\t${INDEX_DIR}"
    echo -e "threads\t${THREADS}"
    echo -e "read_files_command\t${READ_FILES_COMMAND}"
    echo -e "align_intron_max\t${ALIGN_INTRON_MAX}"
    echo -e "out_sam_mult_nmax\t${OUT_SAM_MULT_NMAX}"
    echo -e "win_anchor_multimap_nmax\t${WIN_ANCHOR_MULTIMAP_NMAX}"
    echo -e "out_bam_sorting_bins_n\t${OUT_BAM_SORTING_BINS_N}"
    echo -e "mapq_min\t${MAPQ_MIN}"
    echo -e "alignment_records\t${alignment_records}"
    echo -e "unique_records\t${unique_records}"
    echo -e "aligned_bam\t${ALIGNED_BAM}"
    echo -e "unique_bam\t${UNIQUE_BAM}"
} > "$SUMMARY"

echo ""
log "MOA mapping complete: $UNIQUE_BAM"
STAR_LOG="${PREFIX}Log.final.out"
if [ -f "$STAR_LOG" ]; then
    echo "--- STAR mapping summary (${LIBRARY_ID}) ---"
    grep -E "Number of input reads|Uniquely mapped reads %|% of reads mapped to multiple loci|% of reads unmapped" "$STAR_LOG" \
        | sed 's/^[[:space:]]*/  /'
fi
echo "  alignment records in BAM: $alignment_records"
echo "  unique (MAPQ ${MAPQ_MIN}) records: $unique_records"
echo "  mapping summary: $SUMMARY"
