import csv
import random


rng = random.Random(20260916)
chromosomes = ("chr1", "chr2")
genome_length = 100_000

# All internal intervals are 0-based, half-open, in the P1 genome.
starts = {
    "gene": (10_000, 35_000),
    "atac": (20_000, 50_000),
    "moa": (25_000, 60_000),
    "chip": (30_000, 70_000),
}
lengths = {"gene": 1500, "atac": 2000, "moa": 80, "chip": 400}
features = {
    assay: [(chrom, start, start + lengths[assay])
            for chrom in chromosomes for start in positions]
    for assay, positions in starts.items()
}

# Each edit replaces P1 [position, position + deleted) with inserted DNA.
# These sites lie outside every gene and peak, including gene introns.
indels = {chrom: [(5000, 0, "ACGTA"), (45_000, 7, ""),
                  (80_000, 0, "TGC")]
          for chrom in chromosomes}


def p2_coordinate(chrom, coordinate):
    return coordinate + sum(
        len(inserted) - deleted
        for position, deleted, inserted in indels[chrom]
        if position + deleted <= coordinate
    )


def write_table(filename, header, rows):
    with open(filename, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


with open("P1.fa", "w") as p1_fasta, open("P2.fa", "w") as p2_fasta:
    for chrom in chromosomes:
        p1 = rng.choices("ACGT", k=genome_length)
        # Plus-strand introns have canonical GT...AG splice sites.
        for region_chrom, start, end in features["gene"]:
            if region_chrom == chrom:
                p1[start + 500:start + 502] = "GT"
                p1[start + 998:start + 1000] = "AG"
        p1 = "".join(p1)
        p2 = list(p1)
        for assay, regions in features.items():
            for region_chrom, start, end in regions:
                if region_chrom != chrom:
                    continue
                intervals = [(start, end)]
                if assay == "gene":
                    intervals = [(start, start + 500), (start + 1000, end)]
                # Bin each exon separately; other assays retain their full region.
                # One SNP per 200 bp bin; exactly one in each short MOA peak.
                for left_bound, right_bound in intervals:
                    for left in range(left_bound, right_bound, 200):
                        position = rng.randrange(left, min(left + 200, right_bound))
                        p2[position] = rng.choice("ACGT".replace(p1[position], ""))
        # Apply edits right to left so their P1 positions remain valid.
        for position, deleted, inserted in reversed(indels[chrom]):
            p2[position:position + deleted] = list(inserted)
        for handle, sequence in ((p1_fasta, p1), (p2_fasta, "".join(p2))):
            handle.write(f">{chrom}\n")
            for offset in range(0, len(sequence), 60):
                handle.write(sequence[offset:offset + 60] + "\n")

for parent in ("P1", "P2"):
    with open(f"{parent}.gtf", "w") as handle:
        for index, (chrom, start, end) in enumerate(features["gene"], 1):
            gene_id = f"{parent}_GENE{index}"
            transcript_id = f"{parent}_TX{index}"
            records = [("gene", start, end), ("transcript", start, end),
                       ("exon", start, start + 500),
                       ("exon", start + 1000, end)]
            for kind, left, right in records:
                if parent == "P2":
                    left = p2_coordinate(chrom, left)
                    right = p2_coordinate(chrom, right)
                attributes = f'gene_id "{gene_id}";'
                if kind != "gene":
                    attributes += f' transcript_id "{transcript_id}";'
                # GTF converts [left, right) to 1-based inclusive coordinates.
                handle.write(f"{chrom}\tsynthetic\t{kind}\t{left + 1}\t{right}"
                             f"\t.\t+\t.\t{attributes}\n")

write_table("gene_pairs.tsv", ["parent1_gene_id", "parent2_gene_id"],
            [(f"P1_GENE{i}", f"P2_GENE{i}") for i in range(1, 5)])

peak_header = ["parent1_chrom", "parent1_start", "parent1_end",
               "parent2_chrom", "parent2_start", "parent2_end"]
for assay in ("atac", "moa", "chip"):
    write_table(f"{assay}_peak_pairs.tsv", peak_header,
                [(chrom, start, end, chrom, p2_coordinate(chrom, start),
                  p2_coordinate(chrom, end))
                 for chrom, start, end in features[assay]])

write_table("simulation_design.tsv",
            ["feature_index", "parent1_fraction", "parent2_fraction"],
            [(1, "0.50", "0.50"), (2, "0.80", "0.20"),
             (3, "0.20", "0.80"), (4, "0.50", "0.50")])

print("Created the nine simulation input files.")
