#!/bin/bash
# build_star_index.sh
# Build a STAR genome index for the concatenated reference.
#
# Required: REFERENCE_ID matching the reference directory created by
# build_concatenated_reference.sh.
# Output: resources/references/${REFERENCE_ID}/star_index/
# genomeSAindexNbases is calculated automatically because STAR scales it to
# genome length.

set -euo pipefail

STAR_BIN=${STAR_BIN:-STAR}
SAMTOOLS=${SAMTOOLS:-samtools}

for cmd in "$STAR_BIN" "$SAMTOOLS" awk; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "ERROR: required command not found in PATH: $cmd" >&2
        exit 1
    fi
done

if [ -z "${REFERENCE_ID:-}" ]; then
    echo "ERROR: REFERENCE_ID is required and must match the reference directory created by build_concatenated_reference.sh." >&2
    echo "       Usual naming convention: <PARENT1_PREFIX>_<PARENT2_PREFIX> (example: B73_Mo17)" >&2
    exit 1
fi
REF_FA=${REF_FA:-resources/references/${REFERENCE_ID}/concatenated.fa}
REF_GTF=${REF_GTF:-resources/references/${REFERENCE_ID}/concatenated.gtf}
INDEX_DIR=${INDEX_DIR:-resources/references/${REFERENCE_ID}/star_index}
THREADS=${THREADS:-4}
SJDB_OVERHANG=${SJDB_OVERHANG:-149}   # readLength - 1 (150 bp reads)
REPORT="${INDEX_DIR%/}/star_index_build_report.tsv"

for f in "$REF_FA" "$REF_GTF"; do
    [ -f "$f" ] || { echo "ERROR: input not found: $f" >&2; exit 1; }
done
[ -f "${REF_FA}.fai" ] || "$SAMTOOLS" faidx "$REF_FA"

# ---- compute genome length and genomeSAindexNbases --------------------------
GENOME_LEN=$(awk 'BEGIN{s=0} {s+=$2} END{print s}' "${REF_FA}.fai")
# floor(log2(L)/2 - 1), capped at 14, floored at 4 (STAR practical minimum)
SA_NBASES=$(awk -v L="$GENOME_LEN" 'BEGIN{
    v = int( (log(L)/log(2))/2 - 1 );
    if (v > 14) v = 14;
    if (v < 4)  v = 4;
    print v;
}')

echo "[$(date)] Genome length: ${GENOME_LEN} bp"
echo "[$(date)] Using --genomeSAindexNbases ${SA_NBASES}"

# ---- build index ------------------------------------------------------------
mkdir -p "$INDEX_DIR"
echo "[$(date)] Building STAR index..."
"$STAR_BIN" \
    --runMode genomeGenerate \
    --runThreadN "$THREADS" \
    --genomeDir "$INDEX_DIR" \
    --genomeFastaFiles "$REF_FA" \
    --sjdbGTFfile "$REF_GTF" \
    --sjdbOverhang "$SJDB_OVERHANG" \
    --genomeSAindexNbases "$SA_NBASES"

{
    echo -e "metric\tvalue"
    echo -e "reference_id\t${REFERENCE_ID}"
    echo -e "ref_fasta\t${REF_FA}"
    echo -e "ref_gtf\t${REF_GTF}"
    echo -e "index_dir\t${INDEX_DIR}"
    echo -e "genome_length\t${GENOME_LEN}"
    echo -e "genomeSAindexNbases\t${SA_NBASES}"
    echo -e "sjdbOverhang\t${SJDB_OVERHANG}"
    echo -e "threads\t${THREADS}"
    echo -e "star_bin\t${STAR_BIN}"
} > "$REPORT"

echo ""
echo "[$(date)] STAR index built in ${INDEX_DIR}/"
echo "  Report: $REPORT"
echo "Index contents:"
ls -1 "$INDEX_DIR" | sed 's/^/  /'
