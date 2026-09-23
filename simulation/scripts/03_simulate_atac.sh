#!/usr/bin/env bash
set -euo pipefail

work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
mkdir -p simulated_reads

# Keep any FASTA indexes created by bedtools in the temporary directory.
cp P1.fa P2.fa "$work/"

awk 'BEGIN { FS = OFS = "\t" }
    NR == FNR {
        if (FNR > 1) peaks[FNR - 1] = $0
        next
    }
    FNR == 1 {
        print "feature_index", "parent1_chrom", "parent1_start", "parent1_end", \
              "parent2_chrom", "parent2_start", "parent2_end", \
              "expected_parent1", "expected_parent2", "expected_total"
        next
    }
    {
        n1 = int(1000 * $2 + 0.5)
        n2 = int(1000 * $3 + 0.5)
        print $1, peaks[$1], n1, n2, n1 + n2
    }
' atac_peak_pairs.tsv simulation_design.tsv > "$work/expected.tsv"

for rep in 1 2; do
    : > "$work/R1.fastq"
    : > "$work/R2.fastq"
    while IFS=$'\t' read -r feature chrom1 start1 end1 chrom2 start2 end2 count1 count2 total; do
        [[ "$feature" == feature_index ]] && continue
        for parent in 1 2; do
            if [[ "$parent" == 1 ]]; then
                printf '%s\t%s\t%s\n' "$chrom1" "$start1" "$end1" > "$work/peak.bed"
                count=$count1
            else
                printf '%s\t%s\t%s\n' "$chrom2" "$start2" "$end2" > "$work/peak.bed"
                count=$count2
            fi
            bedtools getfasta -fi "$work/P${parent}.fa" \
                -bed "$work/peak.bed" -fo "$work/peak.fa"
            seed=$((200000 + 10000 * rep + 100 * feature + parent))
            # Controlled validation fragments, not a biological ATAC size model.
            art_illumina -p -ss HS25 -l 150 -m 200 -s 30 \
                -c "$count" -rs "$seed" -na \
                -d "P${parent}__ATAC${feature}__rep${rep}__" \
                -i "$work/peak.fa" -o "$work/art_" > "$work/art.log"
            cat "$work/art_1.fq" >> "$work/R1.fastq"
            cat "$work/art_2.fq" >> "$work/R2.fastq"
        done
    done < "$work/expected.tsv"
    gzip -n -c "$work/R1.fastq" > "simulated_reads/atac_rep${rep}_R1.fastq.gz"
    gzip -n -c "$work/R2.fastq" > "simulated_reads/atac_rep${rep}_R2.fastq.gz"
    cp "$work/expected.tsv" "simulated_reads/expected_atac_rep${rep}.tsv"
done

printf '%s\n' \
    simulated_reads/atac_rep1_R1.fastq.gz \
    simulated_reads/atac_rep1_R2.fastq.gz \
    simulated_reads/atac_rep2_R1.fastq.gz \
    simulated_reads/atac_rep2_R2.fastq.gz \
    simulated_reads/expected_atac_rep1.tsv \
    simulated_reads/expected_atac_rep2.tsv
