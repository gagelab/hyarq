#!/bin/bash
# map_moa_star.sh
# Merge one cleaned paired-end MOA-seq library, map merged reads to the
# concatenated hybrid reference, and retain unique primary alignments.
#
# Required input:
#   R1=/path/to/cleaned_R1.fastq.gz R2=/path/to/cleaned_R2.fastq.gz
#
# Optional convenience:
#   READS_DIR=/cleaned/fastqs LIBRARY_ID=moa_rep1
# derives R1=$READS_DIR/${LIBRARY_ID}_R1.fq.gz and R2=..._R2.fq.gz.
#
# Outputs (OUTDIR, default results/mapping/moa):
#   ${LIBRARY_ID}.unique.bam (+ .bai)
#   ${LIBRARY_ID}.Aligned.sortedByCoord.out.bam
#   ${LIBRARY_ID}.Log.final.out
#   ${LIBRARY_ID}.ngmerge.log
#   ${LIBRARY_ID}.mapping_summary.tsv
#
# MOA fragments are merged before STAR and then mapped single-end with DNA-style
# STAR settings. Multimappers remain visible in STAR logs and are removed by the
# MAPQ filter before counting.

set -euo pipefail

NGMERGE=${NGMERGE:-NGmerge}
STAR_BIN=${STAR_BIN:-STAR}
SAMTOOLS=${SAMTOOLS:-samtools}

for cmd in "$NGMERGE" "$STAR_BIN" "$SAMTOOLS"; do
    command -v "$cmd" >/dev/null 2>&1 || { echo "ERROR: required command not found in PATH: $cmd" >&2; exit 1; }
done

die() {
    echo "ERROR: $*" >&2
    exit 1
}

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}

REP=${REP:-1}
SAMPLE=${SAMPLE:-moa}
LIBRARY_ID=${LIBRARY_ID:-${SAMPLE}_rep${REP}}
INDEX_DIR=${INDEX_DIR:-results/star_index}
OUTDIR=${OUTDIR:-results/mapping/moa}
MOA_TMPDIR=${MOA_TMPDIR:-moa_map_tmp/${LIBRARY_ID}}
THREADS=${THREADS:-4}
MAPQ_MIN=${MAPQ_MIN:-255}
NGMERGE_MIN_OVERLAP=${NGMERGE_MIN_OVERLAP:-15}
R1=${R1:-}
R2=${R2:-}

if [ -z "$R1" ] || [ -z "$R2" ]; then
    if [ -n "${READS_DIR:-}" ]; then
        R1=${R1:-${READS_DIR%/}/${LIBRARY_ID}_R1.fq.gz}
        R2=${R2:-${READS_DIR%/}/${LIBRARY_ID}_R2.fq.gz}
    else
        die "set R1 and R2 to cleaned FASTQ inputs, or set READS_DIR with LIBRARY_ID=${LIBRARY_ID}"
    fi
fi

[ -d "$INDEX_DIR" ] || die "STAR index not found: $INDEX_DIR"
for f in "$R1" "$R2"; do
    [ -f "$f" ] || die "input FASTQ not found: $f"
done
if [[ "$R1" == *.gz || "$R2" == *.gz ]]; then
    command -v zcat >/dev/null 2>&1 || die "required command not found in PATH: zcat"
fi

mkdir -p "$OUTDIR" "$MOA_TMPDIR"
PREFIX="$OUTDIR/${LIBRARY_ID}."
ALIGNED="${PREFIX}Aligned.sortedByCoord.out.bam"
UNIQUE_BAM="${PREFIX}unique.bam"
NGMERGE_LOG="$OUTDIR/${LIBRARY_ID}.ngmerge.log"
SUMMARY="$OUTDIR/${LIBRARY_ID}.mapping_summary.tsv"

log "MOA mapping: ${LIBRARY_ID}"
echo "  R1     = $R1"
echo "  R2     = $R2"
echo "  index  = $INDEX_DIR"
echo "  tmp    = $MOA_TMPDIR"

log "Preparing uncompressed FASTQs for NGmerge"
if [[ "$R1" == *.gz ]]; then
    zcat "$R1" > "$MOA_TMPDIR/${LIBRARY_ID}_R1.fq"
else
    cat "$R1" > "$MOA_TMPDIR/${LIBRARY_ID}_R1.fq"
fi
if [[ "$R2" == *.gz ]]; then
    zcat "$R2" > "$MOA_TMPDIR/${LIBRARY_ID}_R2.fq"
else
    cat "$R2" > "$MOA_TMPDIR/${LIBRARY_ID}_R2.fq"
fi

log "Merging overlapping paired-end reads with NGmerge"
"$NGMERGE" -1 "$MOA_TMPDIR/${LIBRARY_ID}_R1.fq" \
    -2 "$MOA_TMPDIR/${LIBRARY_ID}_R2.fq" \
    -o "$MOA_TMPDIR/merged.fastq.gz" \
    -m "$NGMERGE_MIN_OVERLAP" -n "$THREADS" -z \
    > "$MOA_TMPDIR/ngmerge.log" 2>&1 || {
        echo "ERROR: NGmerge failed (see $MOA_TMPDIR/ngmerge.log)" >&2; exit 1; }
[ -s "$MOA_TMPDIR/merged.fastq.gz" ] || {
    echo "ERROR: NGmerge produced no merged reads (see $MOA_TMPDIR/ngmerge.log)" >&2
    exit 1
}
cp "$MOA_TMPDIR/ngmerge.log" "$NGMERGE_LOG"
n_merged=$(zcat "$MOA_TMPDIR/merged.fastq.gz" | awk 'END{print NR/4}')
echo "  merged reads: $n_merged"

log "Mapping merged reads single-end with STAR"
"$STAR_BIN" \
    --runMode alignReads \
    --runThreadN "$THREADS" \
    --genomeDir "$INDEX_DIR" \
    --readFilesIn "$MOA_TMPDIR/merged.fastq.gz" \
    --readFilesCommand zcat \
    --alignIntronMax 1 \
    --outSAMmultNmax 2 \
    --winAnchorMultimapNmax 100 \
    --outBAMsortingBinsN 5 \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix "$PREFIX"

[ -s "$ALIGNED" ] || die "STAR produced no BAM for ${LIBRARY_ID}"

log "Filtering to unique primary alignments (MAPQ ${MAPQ_MIN})"
"$SAMTOOLS" view -b -q "$MAPQ_MIN" -F 0x904 "$ALIGNED" > "$UNIQUE_BAM"
"$SAMTOOLS" index "$UNIQUE_BAM"

alignment_records=$("$SAMTOOLS" view -c "$ALIGNED")
unique_records=$("$SAMTOOLS" view -c "$UNIQUE_BAM")
{
    echo -e "metric\tvalue"
    echo -e "assay\tmoa"
    echo -e "sample\t${SAMPLE}"
    echo -e "rep\t${REP}"
    echo -e "library_id\t${LIBRARY_ID}"
    echo -e "r1\t${R1}"
    echo -e "r2\t${R2}"
    echo -e "star_index\t${INDEX_DIR}"
    echo -e "threads\t${THREADS}"
    echo -e "ngmerge_min_overlap\t${NGMERGE_MIN_OVERLAP}"
    echo -e "align_intron_max\t1"
    echo -e "out_sam_mult_nmax\t2"
    echo -e "mapq_min\t${MAPQ_MIN}"
    echo -e "merged_reads\t${n_merged}"
    echo -e "alignment_records\t${alignment_records}"
    echo -e "unique_records\t${unique_records}"
    echo -e "bam\t${UNIQUE_BAM}"
    echo -e "ngmerge_log\t${NGMERGE_LOG}"
} > "$SUMMARY"

echo ""
log "MOA mapping complete: $UNIQUE_BAM"
log="${PREFIX}Log.final.out"
if [ -f "$log" ]; then
    echo "--- STAR mapping summary (${LIBRARY_ID}) ---"
    grep -E "Number of input reads|Uniquely mapped reads %|% of reads mapped to multiple loci|% of reads unmapped" "$log" \
        | sed 's/^[[:space:]]*/  /'
fi
echo "  merged reads in:        $n_merged"
echo "  unique (MAPQ ${MAPQ_MIN}) records: $unique_records"
echo "  mapping summary: $SUMMARY"

rm -rf "$MOA_TMPDIR"
