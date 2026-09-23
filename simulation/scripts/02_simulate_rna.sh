#!/usr/bin/env bash
set -euo pipefail

work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
mkdir -p simulated_reads

for parent in P1 P2; do
    gffread "${parent}.gtf" -g "${parent}.fa" \
        -w "$work/${parent}_transcripts.fa" > "$work/gffread_${parent}.log"
done

# Pair-table row order matches the design feature indices.
awk 'BEGIN { FS = OFS = "\t" }
    NR == FNR {
        if (FNR > 1) { p1[FNR - 1] = $1; p2[FNR - 1] = $2 }
        next
    }
    FNR == 1 {
        print "feature_index", "parent1_gene_id", "parent2_gene_id", \
              "expected_parent1", "expected_parent2", "expected_total"
        next
    }
    {
        n1 = int(1000 * $2 + 0.5)
        n2 = int(1000 * $3 + 0.5)
        print $1, p1[$1], p2[$1], n1, n2, n1 + n2
    }
' gene_pairs.tsv simulation_design.tsv > "$work/expected.tsv"

for rep in 1 2; do
    : > "$work/R1.fastq"
    : > "$work/R2.fastq"
    while IFS=$'\t' read -r feature gene1 gene2 count1 count2 total; do
        [[ "$feature" == feature_index ]] && continue
        for parent in 1 2; do
            if [[ "$parent" == 1 ]]; then
                gene=$gene1
                count=$count1
            else
                gene=$gene2
                count=$count2
            fi
            transcript="P${parent}_TX${feature}"
            awk -v id="$transcript" '
                /^>/ { keep = (substr($1, 2) == id) }
                keep { print }
            ' "$work/P${parent}_transcripts.fa" > "$work/transcript.fa"
            [[ -s "$work/transcript.fa" ]] || {
                printf 'Missing transcript: %s\n' "$transcript" >&2
                exit 1
            }
            seed=$((10000 * rep + 100 * feature + parent))
            art_illumina -p -ss HS25 -l 150 -m 300 -s 40 \
                -c "$count" -rs "$seed" -na \
                -d "P${parent}__${gene}__feature${feature}__rep${rep}__" \
                -i "$work/transcript.fa" -o "$work/art_" \
                > "$work/art.log"
            cat "$work/art_1.fq" >> "$work/R1.fastq"
            cat "$work/art_2.fq" >> "$work/R2.fastq"
        done
    done < "$work/expected.tsv"
    gzip -n -c "$work/R1.fastq" > "simulated_reads/rna_rep${rep}_R1.fastq.gz"
    gzip -n -c "$work/R2.fastq" > "simulated_reads/rna_rep${rep}_R2.fastq.gz"
    cp "$work/expected.tsv" "simulated_reads/expected_rna_rep${rep}.tsv"
done

printf '%s\n' \
    simulated_reads/rna_rep1_R1.fastq.gz \
    simulated_reads/rna_rep1_R2.fastq.gz \
    simulated_reads/rna_rep2_R1.fastq.gz \
    simulated_reads/rna_rep2_R2.fastq.gz \
    simulated_reads/expected_rna_rep1.tsv \
    simulated_reads/expected_rna_rep2.tsv
