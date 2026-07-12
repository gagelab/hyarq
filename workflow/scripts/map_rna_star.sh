#!/bin/bash
# map_rna_star.sh
# Map one cleaned paired-end RNA-seq library with STAR.
#
# Required: R1, R2, and REFERENCE_ID.
#
# Outputs:
#   ${LIBRARY_ID}.Aligned.sortedByCoord.out.bam (+ .bai)
#   ${LIBRARY_ID}.ReadsPerGene.out.tab
#   ${LIBRARY_ID}.Log.final.out
#   ${LIBRARY_ID}.mapping_summary.tsv
#
# Final paired parental gene counts are generated later by count_rna_genes.sh.

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
SAMPLE=${SAMPLE:-rna}
LIBRARY_ID=${LIBRARY_ID:-${SAMPLE}_rep${REP}}
if [ -z "${REFERENCE_ID:-}" ]; then
    die "REFERENCE_ID is required and must match the reference directory used by build_star_index.sh. Usual naming convention: <PARENT1_PREFIX>_<PARENT2_PREFIX> (example: B73_Mo17)"
fi
INDEX_DIR=${INDEX_DIR:-resources/references/${REFERENCE_ID}/star_index}
OUTDIR=${OUTDIR:-results/mapping/rna}
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

log "RNA mapping: ${LIBRARY_ID}"
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
    --outSAMmultNmax 1 \
    --outFilterMultimapNmax 1 \
    --winAnchorMultimapNmax 100 \
    --twopassMode Basic \
    --outFilterIntronMotifs RemoveNoncanonical \
    --outFilterType BySJout \
    --quantMode GeneCounts \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix "$PREFIX"

[ -s "$BAM" ] || die "STAR produced no BAM for ${LIBRARY_ID}"
"$SAMTOOLS" index "$BAM"

alignment_records=$("$SAMTOOLS" view -c "$BAM")
{
    echo -e "metric\tvalue"
    echo -e "reference_id\t${REFERENCE_ID}"
    echo -e "assay\trna"
    echo -e "sample\t${SAMPLE}"
    echo -e "rep\t${REP}"
    echo -e "library_id\t${LIBRARY_ID}"
    echo -e "r1\t${R1}"
    echo -e "r2\t${R2}"
    echo -e "star_index\t${INDEX_DIR}"
    echo -e "threads\t${THREADS}"
    echo -e "read_files_command\t${READ_FILES_COMMAND}"
    echo -e "out_filter_multimap_nmax\t1"
    echo -e "out_sam_mult_nmax\t1"
    echo -e "twopass_mode\tBasic"
    echo -e "alignment_records\t${alignment_records}"
    echo -e "bam\t${BAM}"
    echo -e "star_gene_counts\t${PREFIX}ReadsPerGene.out.tab"
} > "$SUMMARY"

echo ""
log "RNA mapping complete: $BAM"
log="${PREFIX}Log.final.out"
if [ -f "$log" ]; then
    echo "--- STAR mapping summary (${LIBRARY_ID}) ---"
    grep -E "Number of input reads|Uniquely mapped reads %|% of reads mapped to multiple loci|% of reads unmapped" "$log" \
        | sed 's/^[[:space:]]*/  /'
fi
echo "  alignment records in BAM: $alignment_records"
echo "  STAR gene counts (cross-check): ${PREFIX}ReadsPerGene.out.tab"
echo "  mapping summary: $SUMMARY"
