#!/bin/bash
# map_chip_star.sh
# Map, deduplicate, and filter one cleaned paired-end ChIP-seq library.
#
# One invocation maps one specified sample/replicate library. Each library is
# processed independently:
#   1. paired-end DNA-style alignment to the concatenated reference
#   2. name-sort and add mate tags with samtools fixmate -m
#   3. coordinate-sort and remove duplicate reads with samtools markdup -r
#   4. retain mapped, primary, non-duplicate alignments with MAPQ >= MAPQ_MIN
#
# This script only produces final BAMs, indexes, mapping logs, and simple
# summaries. Region counting and allele-specific testing are intentionally
# handled in later modules.
#
# Required input:
#   R1=/path/to/cleaned_R1.fastq.gz R2=/path/to/cleaned_R2.fastq.gz
#
# Optional convenience:
#   READS_DIR=/cleaned/fastqs LIBRARY_ID=chip_ip_rep1
# derives R1=$READS_DIR/${LIBRARY_ID}_R1.fq.gz and R2=..._R2.fq.gz.
#
# Outputs (OUTDIR):
#   ${LIBRARY_ID}.dedup.unique.bam (+ .bai)
#   ${LIBRARY_ID}.mapping_summary.tsv
#   ${LIBRARY_ID}.Log*.out
#
# Requires: STAR, samtools

set -euo pipefail

STAR_BIN=${STAR_BIN:-STAR}
SAMTOOLS=${SAMTOOLS:-samtools}

for cmd in "$STAR_BIN" "$SAMTOOLS"; do
    command -v "$cmd" >/dev/null 2>&1 || {
        echo "ERROR: required command not found in PATH: $cmd" >&2
        exit 1
    }
done

SAMPLE=${SAMPLE:-ip}
REP=${REP:-1}
LIBRARY_ID=${LIBRARY_ID:-chip_${SAMPLE}_rep${REP}}
THREADS=${THREADS:-4}
MAPQ_MIN=${MAPQ_MIN:-255}
INDEX_DIR=${INDEX_DIR:-${STAR_INDEX:-results/star_index}}
OUTDIR=${OUTDIR:-results/mapping/chip}
WORKDIR_BASE=${WORKDIR_BASE:-tmp}
R1=${R1:-}
R2=${R2:-}

[ -n "$SAMPLE" ] || {
    echo "ERROR: SAMPLE must not be empty." >&2
    exit 1
}

if [ -z "$R1" ] || [ -z "$R2" ]; then
    if [ -n "${READS_DIR:-}" ]; then
        R1=${R1:-${READS_DIR%/}/${LIBRARY_ID}_R1.fq.gz}
        R2=${R2:-${READS_DIR%/}/${LIBRARY_ID}_R2.fq.gz}
    elif [ -n "${READDIR:-}" ]; then
        R1=${R1:-${READDIR%/}/${LIBRARY_ID}_R1.fq.gz}
        R2=${R2:-${READDIR%/}/${LIBRARY_ID}_R2.fq.gz}
    else
        echo "ERROR: set R1 and R2 to cleaned FASTQ inputs, or set READS_DIR with LIBRARY_ID=${LIBRARY_ID}" >&2
        exit 1
    fi
fi

if [ -z "${READ_FILES_COMMAND:-}" ]; then
    case "$R1:$R2" in
        *.gz:*.gz) READ_FILES_COMMAND=zcat ;;
        *) READ_FILES_COMMAND=cat ;;
    esac
fi
command -v "${READ_FILES_COMMAND%% *}" >/dev/null 2>&1 || {
    echo "ERROR: required command not found in PATH: ${READ_FILES_COMMAND%% *}" >&2
    exit 1
}

[ -d "$INDEX_DIR" ] || {
    echo "ERROR: STAR index not found: $INDEX_DIR" >&2
    exit 1
}
for f in "$R1" "$R2"; do
    [ -f "$f" ] || {
        echo "ERROR: input FASTQ not found: $f" >&2
        exit 1
    }
done

mkdir -p "$OUTDIR" "$WORKDIR_BASE"
WORKDIR="${WORKDIR_BASE%/}/map_${LIBRARY_ID}"
PREFIX="$WORKDIR/${LIBRARY_ID}."
RAW_BAM="${PREFIX}Aligned.sortedByCoord.out.bam"
NAME_BAM="$WORKDIR/${LIBRARY_ID}.namesort.bam"
FIXMATE_BAM="$WORKDIR/${LIBRARY_ID}.fixmate.bam"
COORD_BAM="$WORKDIR/${LIBRARY_ID}.coordsort.bam"
DEDUP_BAM="$WORKDIR/${LIBRARY_ID}.dedup.bam"
FINAL_BAM="$OUTDIR/${LIBRARY_ID}.dedup.unique.bam"
SUMMARY="$OUTDIR/${LIBRARY_ID}.mapping_summary.tsv"

case "$(basename "$WORKDIR")" in
    map_*) ;;
    *)
        echo "ERROR: unsafe run-specific work directory: $WORKDIR" >&2
        exit 1
        ;;
esac
rm -rf "$WORKDIR"
mkdir -p "$WORKDIR"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] ChIP mapping: ${LIBRARY_ID}"
echo "  R1       = $R1"
echo "  R2       = $R2"
echo "  index    = $INDEX_DIR"
echo "  MAPQ min = $MAPQ_MIN"

"$STAR_BIN" \
    --runMode alignReads \
    --runThreadN "$THREADS" \
    --genomeDir "$INDEX_DIR" \
    --readFilesIn "$R1" "$R2" \
    --readFilesCommand "$READ_FILES_COMMAND" \
    --alignIntronMax 1 \
    --outSAMmultNmax 2 \
    --winAnchorMultimapNmax 100 \
    --outBAMsortingBinsN 5 \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix "$PREFIX"

[ -s "$RAW_BAM" ] || {
    echo "ERROR: STAR produced no BAM for ${LIBRARY_ID}" >&2
    exit 1
}

for log_name in Log.final.out Log.out Log.progress.out; do
    if [ -f "${PREFIX}${log_name}" ]; then
        cp "${PREFIX}${log_name}" "$OUTDIR/${LIBRARY_ID}.${log_name}"
    fi
done

RAW_MAPPED=$("$SAMTOOLS" view -c -F 0x4 "$RAW_BAM")

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Deduplicating ${LIBRARY_ID}"
"$SAMTOOLS" sort -n -@ "$THREADS" -o "$NAME_BAM" "$RAW_BAM"
"$SAMTOOLS" fixmate -m -@ "$THREADS" "$NAME_BAM" "$FIXMATE_BAM"
"$SAMTOOLS" sort -@ "$THREADS" -o "$COORD_BAM" "$FIXMATE_BAM"
"$SAMTOOLS" markdup -r -@ "$THREADS" "$COORD_BAM" "$DEDUP_BAM"

DEDUP_MAPPED=$("$SAMTOOLS" view -c -F 0x4 "$DEDUP_BAM")
DUPLICATES_REMOVED=$((RAW_MAPPED - DEDUP_MAPPED))

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Filtering final ${LIBRARY_ID} BAM (MAPQ ${MAPQ_MIN})"
"$SAMTOOLS" view -h -b -q "$MAPQ_MIN" -F 0xD04 "$DEDUP_BAM" > "$FINAL_BAM"
"$SAMTOOLS" index -@ "$THREADS" "$FINAL_BAM"

FINAL_UNIQUE=$("$SAMTOOLS" view -c "$FINAL_BAM")
FINAL_PAIRS=$("$SAMTOOLS" view -c -f 0x40 "$FINAL_BAM")

{
    echo -e "metric\tvalue"
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
    echo -e "out_sam_mult_nmax\t2"
    echo -e "mapq_min\t${MAPQ_MIN}"
    echo -e "raw_mapped_alignments\t${RAW_MAPPED}"
    echo -e "deduplicated_mapped_alignments\t${DEDUP_MAPPED}"
    echo -e "duplicates_removed\t${DUPLICATES_REMOVED}"
    echo -e "final_unique_alignments\t${FINAL_UNIQUE}"
    echo -e "final_unique_read_pairs\t${FINAL_PAIRS}"
    echo -e "bam\t${FINAL_BAM}"
} > "$SUMMARY"

echo ""
echo "[$(date '+%Y-%m-%d %H:%M:%S')] ChIP mapping complete: ${LIBRARY_ID}"
echo "  final BAM: $FINAL_BAM"
echo "  final unique alignments: $FINAL_UNIQUE"
echo "  final unique read pairs: $FINAL_PAIRS"
echo "  mapping summary: $SUMMARY"

echo "Cleaning temporary files for ${LIBRARY_ID}..."
rm -rf "$WORKDIR"
