#!/bin/bash
# build_star_index.sh
# Build a STAR genome index for the concatenated hybrid reference.
#
# One index is built for ALL assays. It includes the annotation (GTF) so that
# splice junctions are available for RNA-seq; assays that must not use introns
# (ATAC/MOA/ChIP) disable them at MAPPING time with --alignIntronMax 1, so the
# same index serves every assay.
#
# --genomeSAindexNbases must be scaled to genome length. STAR's rule is
#   min(14, floor(log2(genomeLength)/2 - 1))
# Small references need lower values, while large plant genomes usually use the
# default cap of 14. We compute it automatically from the .fai.
#
# Inputs:
#   results/reference/concatenated.fa   (+ .fai, created if missing)
#   results/reference/concatenated.gtf
# Output:
#   results/star_index/
#
# Requires: STAR, samtools, awk
# Run from the project directory (where results/reference/ lives).

set -euo pipefail

STAR_BIN=${STAR_BIN:-STAR}
SAMTOOLS=${SAMTOOLS:-samtools}

for cmd in "$STAR_BIN" "$SAMTOOLS" awk; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "ERROR: required command not found in PATH: $cmd" >&2
        exit 1
    fi
done

REF_FA=${REF_FA:-results/reference/concatenated.fa}
REF_GTF=${REF_GTF:-results/reference/concatenated.gtf}
INDEX_DIR=${INDEX_DIR:-results/star_index}
THREADS=${THREADS:-4}
SJDB_OVERHANG=${SJDB_OVERHANG:-149}   # readLength - 1 (150 bp reads)

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

echo ""
echo "[$(date)] STAR index built in ${INDEX_DIR}/"
echo "Index contents:"
ls -1 "$INDEX_DIR" | sed 's/^/  /'
