#!/bin/bash
set -euo pipefail

# ":?" errors on unset OR empty
PARENT1_GTF=${PARENT1_GTF:?set PARENT1_GTF}
PARENT2_GTF=${PARENT2_GTF:?set PARENT2_GTF}
PARENT1_PREFIX=${PARENT1_PREFIX:?set PARENT1_PREFIX}
PARENT2_PREFIX=${PARENT2_PREFIX:?set PARENT2_PREFIX}
REFERENCE_ID=${REFERENCE_ID:-${PARENT1_PREFIX}_${PARENT2_PREFIX}}
OUTDIR=${OUTDIR:-resources/references/${REFERENCE_ID}}

# identical prefixes would re-create the contig collision the design avoids
[ "$PARENT1_PREFIX" != "$PARENT2_PREFIX" ] || {
    echo "ERROR: PARENT1_PREFIX and PARENT2_PREFIX must differ (both '$PARENT1_PREFIX')." >&2
    exit 1
}

P1="${PARENT1_PREFIX}_"
P2="${PARENT2_PREFIX}_"
mkdir -p "$OUTDIR"
GTF="$OUTDIR/concatenated.gtf"

read_any() { case "$1" in *.gz) zcat < "$1" ;; *) cat "$1" ;; esac; }

# prefix + concatenate GTF (seqname = col 1)
echo "[$(date)] Prefixing and concatenating GTF..."
read_any "$PARENT1_GTF" | awk -v p="$P1" 'BEGIN{FS=OFS="\t"} /^#/{print;next} {$1=p$1} 1'  >  "$GTF"
read_any "$PARENT2_GTF" | awk -v p="$P2" 'BEGIN{FS=OFS="\t"} /^#/{print;next} {$1=p$1} 1'  >> "$GTF"

echo "[$(date)] Done."
echo "  GTF: $GTF"
