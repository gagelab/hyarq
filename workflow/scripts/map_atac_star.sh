#!/bin/bash
# map_atac_star.sh
# Map one cleaned paired-end ATAC-seq library to the concatenated hybrid
# reference with STAR DNA-style settings.
#
# Required input:
#   R1=/path/to/cleaned_R1.fastq.gz R2=/path/to/cleaned_R2.fastq.gz
#
# Optional convenience:
#   READS_DIR=/cleaned/fastqs LIBRARY_ID=atac_rep1
# derives R1=$READS_DIR/${LIBRARY_ID}_R1.fq.gz and R2=..._R2.fq.gz.
#
# Outputs (OUTDIR, default results/mapping/atac):
#   ${LIBRARY_ID}.Aligned.sortedByCoord.out.bam (+ .bai)
#   ${LIBRARY_ID}.Log.final.out
#   ${LIBRARY_ID}.mapping_summary.tsv
#
# ATAC keeps up to two STAR alignments so fragments without parent-discriminating
# sequence remain visible until downstream MAPQ and same-haplotype filtering.

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

REP=${REP:-1}
SAMPLE=${SAMPLE:-atac}
LIBRARY_ID=${LIBRARY_ID:-${SAMPLE}_rep${REP}}
INDEX_DIR=${INDEX_DIR:-results/star_index}
OUTDIR=${OUTDIR:-results/mapping/atac}
THREADS=${THREADS:-4}
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

if [ -z "${READ_FILES_COMMAND:-}" ]; then
    case "$R1:$R2" in
        *.gz:*.gz) READ_FILES_COMMAND=zcat ;;
        *) READ_FILES_COMMAND=cat ;;
    esac
fi
require_cmd "${READ_FILES_COMMAND%% *}"

[ -d "$INDEX_DIR" ] || die "STAR index not found: $INDEX_DIR"
for f in "$R1" "$R2"; do
    [ -f "$f" ] || die "input FASTQ not found: $f"
done

mkdir -p "$OUTDIR"
PREFIX="$OUTDIR/${LIBRARY_ID}."
BAM="${PREFIX}Aligned.sortedByCoord.out.bam"
SUMMARY="$OUTDIR/${LIBRARY_ID}.mapping_summary.tsv"

log "ATAC mapping: ${LIBRARY_ID}"
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
    --outSAMmultNmax 2 \
    --winAnchorMultimapNmax 100 \
    --outBAMsortingBinsN 5 \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix "$PREFIX"

[ -s "$BAM" ] || die "STAR produced no BAM for ${LIBRARY_ID}"
"$SAMTOOLS" index "$BAM"

alignment_records=$("$SAMTOOLS" view -c "$BAM")
{
    echo -e "metric\tvalue"
    echo -e "assay\tatac"
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
    echo -e "alignment_records\t${alignment_records}"
    echo -e "bam\t${BAM}"
} > "$SUMMARY"

echo ""
log "ATAC mapping complete: $BAM"
log_file="${PREFIX}Log.final.out"
if [ -f "$log_file" ]; then
    echo "--- STAR mapping summary (${LIBRARY_ID}) ---"
    grep -E "Number of input reads|Uniquely mapped reads %|% of reads mapped to multiple loci|% of reads unmapped" "$log_file" \
        | sed 's/^[[:space:]]*/  /'
fi
echo "  alignment records in BAM: $alignment_records"
echo "  mapping summary: $SUMMARY"
