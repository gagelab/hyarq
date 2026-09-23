#!/usr/bin/env bash
set -euo pipefail

work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
mkdir -p simulated_reads

# Keep bedtools FASTA indexes with the temporary inputs.
cp P1.fa P2.fa "$work/"

awk 'BEGIN { FS = OFS = "\t" }
    NR == 1 {
        print "feature_index", "expected_parent1", "expected_parent2", "expected_total"
        next
    }
    {
        n1 = int(1000 * $2 + 0.5)
        n2 = int(1000 * $3 + 0.5)
        print $1, n1, n2, n1 + n2
    }
' simulation_design.tsv > "$work/expected.tsv"

for rep in 1 2; do
    : > "$work/R1.fastq"
    : > "$work/R2.fastq"
    while IFS=$'\t' read -r feature count1 count2 total; do
        [[ "$feature" == feature_index ]] && continue
        for parent in 1 2; do
            if [[ "$parent" == 1 ]]; then
                count=$count1
            else
                count=$count2
            fi
            awk -v row="$((feature + 1))" -v parent="$parent" '
                BEGIN { FS = OFS = "\t" }
                NR == row {
                    column = 3 * (parent - 1) + 1
                    print $column, $(column + 1), $(column + 2)
                }
            ' chip_peak_pairs.tsv > "$work/peak.bed"
            bedtools getfasta -fi "$work/P${parent}.fa" \
                -bed "$work/peak.bed" -fo "$work/peak.fa"
            seed=$((400000 + 10000 * rep + 100 * feature + parent))
            art_illumina -p -ss HS25 -l 150 -m 250 -s 40 \
                -c "$count" -rs "$seed" -na \
                -d "P${parent}__ChIP${feature}__rep${rep}__" \
                -i "$work/peak.fa" -o "$work/art_" > "$work/art.log"
            cat "$work/art_1.fq" >> "$work/R1.fastq"
            cat "$work/art_2.fq" >> "$work/R2.fastq"
        done
    done < "$work/expected.tsv"
    gzip -n -c "$work/R1.fastq" > "simulated_reads/chip_rep${rep}_R1.fastq.gz"
    gzip -n -c "$work/R2.fastq" > "simulated_reads/chip_rep${rep}_R2.fastq.gz"
    cp "$work/expected.tsv" "simulated_reads/expected_chip_rep${rep}.tsv"
done

printf '%s\n' \
    simulated_reads/chip_rep1_R1.fastq.gz \
    simulated_reads/chip_rep1_R2.fastq.gz \
    simulated_reads/chip_rep2_R1.fastq.gz \
    simulated_reads/chip_rep2_R2.fastq.gz \
    simulated_reads/expected_chip_rep1.tsv \
    simulated_reads/expected_chip_rep2.tsv
