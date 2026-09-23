import argparse
import csv
import re
import subprocess
from collections import Counter, defaultdict
from pathlib import Path


ASSAYS = ("RNA", "MOA", "ATAC", "ChIP")
PEAK_COLUMNS = [f"parent{p}_{field}" for p in (1, 2)
                for field in ("chrom", "start", "end")]
TAGS = {
    "RNA": re.compile(r"P([12])__P\1_GENE(\d+)__feature(\d+)__rep(\d+)__"),
    **{a: re.compile(r"P([12])__" + a + r"(\d+)__rep(\d+)__")
       for a in ("MOA", "ATAC", "ChIP")},
}


def table(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def metrics(path):
    return {row["metric"]: row["value"] for row in table(path)}


def ratio(numerator, denominator):
    return numerator / denominator if denominator else None


def mean(values):
    values = [v for v in values if v is not None]
    return ratio(sum(values), len(values))


def write_table(path, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]),
                                delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: "NA" if v is None else v for k, v in row.items()})


def feature_intervals(simulation, assay, pairs):
    intervals = defaultdict(list)
    for parent in (1, 2):
        if assay == "RNA":
            ids = {row[f"parent{parent}_gene_id"]: i
                   for i, row in enumerate(pairs, 1)}
            with (simulation / f"P{parent}.gtf").open() as handle:
                for line in handle:
                    if line.startswith("#") or not line.strip():
                        continue
                    fields = line.rstrip().split("\t")
                    if fields[2] != "exon":
                        continue
                    gene = re.search(r'gene_id "([^"]+)"', fields[8]).group(1)
                    if gene in ids:
                        intervals[f"P{parent}_{fields[0]}"].append(
                            (int(fields[3]) - 1, int(fields[4]), ids[gene]))
        else:
            for i, row in enumerate(pairs, 1):
                prefix = f"parent{parent}_"
                intervals[f"P{parent}_{row[prefix + 'chrom']}"].append(
                    (int(row[prefix + "start"]), int(row[prefix + "end"]), i))
    return intervals


def assigned_feature(contig, position, cigar, intervals):
    hits = set()
    # Only aligned reference bases overlap features; skip introns and deletions.
    for length, operation in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
        length = int(length)
        if operation in "M=X":
            for start, end, feature in intervals.get(contig, []):
                if position < end and position + length > start:
                    hits.add(feature)
        if operation in "MDN=X":
            position += length
    return next(iter(hits)) if len(hits) == 1 else None


def inspect_bam(path, assay, replicate, intervals, expected):
    counts = defaultdict(Counter)
    seen = set()
    command = ["samtools", "view", "-F", "2308", "-q", "255"]
    if assay != "MOA":
        command += ["-f", "64"]
    command.append(str(path))
    with subprocess.Popen(command, stdout=subprocess.PIPE, text=True) as process:
        for line in process.stdout:
            fields = line.rstrip().split("\t")
            name = fields[0]
            match = TAGS[assay].search(name)
            if match is None:
                raise ValueError(f"Unrecognized {assay} origin tag: {name}")
            groups = list(map(int, match.groups()))
            if assay == "RNA":
                parent, gene, feature, rep = groups
                if gene != feature:
                    raise ValueError(f"Inconsistent RNA tag: {name}")
            else:
                parent, feature, rep = groups
            if rep != replicate or feature not in expected:
                raise ValueError(f"Unexpected origin tag in {path}: {name}")
            if name in seen:
                continue
            seen.add(name)
            c = counts[feature]
            c[f"retained_true_parent{parent}"] += 1
            contig = fields[2]
            assigned_parent = 1 if contig.startswith("P1_") else (
                2 if contig.startswith("P2_") else None)
            locus = assigned_feature(contig, int(fields[3]) - 1, fields[5], intervals)
            if locus != feature or assigned_parent is None:
                c["incorrect_feature_assignment"] += 1
            elif assigned_parent == parent:
                c[f"correct_parent{parent}"] += 1
            else:
                c[f"parent{parent}_to_parent{assigned_parent}"] += 1
        if process.wait():
            raise subprocess.CalledProcessError(process.returncode, command)
    return counts


def main():
    parser = argparse.ArgumentParser(description="Evaluate synthetic HyARQ libraries.")
    parser.add_argument("--hyarq-run", type=Path, required=True)
    parser.add_argument("--simulation", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    design = {int(r["feature_index"]): float(r["parent1_fraction"])
              for r in table(args.simulation / "simulation_design.tsv")}
    evaluations, libraries, assays = [], [], []
    for assay in ASSAYS:
        slug = assay.lower()
        pair_file = "gene_pairs.tsv" if assay == "RNA" else f"{slug}_peak_pairs.tsv"
        pairs = table(args.simulation / pair_file)
        keys = ["parent1_gene_id", "parent2_gene_id"] if assay == "RNA" else PEAK_COLUMNS
        intervals = feature_intervals(args.simulation, assay, pairs)
        assay_rows, assay_libraries = [], []
        for rep in (1, 2):
            library = f"{slug}_rep{rep}"
            expected = {int(r["feature_index"]): r for r in table(
                args.simulation / "simulated_reads" / f"expected_{library}.tsv")}
            mapping_dir = args.hyarq_run / "results/mapping" / slug
            count_dir = args.hyarq_run / "results/counts" / slug / library
            mapping = metrics(mapping_dir / f"{library}.mapping_summary.tsv")
            filters = metrics(count_dir / "filter_summary.tsv")
            count_file = "gene_counts.tsv" if assay == "RNA" else "peak_counts.tsv"
            observed = {tuple(r[k] for k in keys): r for r in table(count_dir / count_file)}
            suffix = ".Aligned.sortedByCoord.out.bam" if assay == "RNA" else ".unique.bam"
            counts = inspect_bam(mapping_dir / (library + suffix), assay, rep,
                                 intervals, expected)
            rows = []
            for feature in sorted(expected):
                truth, c = expected[feature], counts[feature]
                n1, n2, total = (int(truth[k]) for k in
                                 ("expected_parent1", "expected_parent2", "expected_total"))
                final = observed[tuple(pairs[feature - 1][k] for k in keys)]
                o1, o2 = int(final["parent1_count"]), int(final["parent2_count"])
                t1, t2 = c["retained_true_parent1"], c["retained_true_parent2"]
                if t1 > n1 or t2 > n2:
                    raise ValueError(f"Retained fragments exceed truth: {library}, {feature}")
                r1, r2 = ratio(t1, n1), ratio(t2, n2)
                fraction = ratio(o1, o1 + o2)
                error = fraction - design[feature] if fraction is not None else None
                row = dict(assay=assay, replicate=rep, feature_index=feature,
                           expected_parent1=n1, expected_parent2=n2, expected_total=total,
                           expected_parent1_fraction=design[feature])
                for key in ("retained_true_parent1", "retained_true_parent2", "correct_parent1",
                            "correct_parent2", "parent1_to_parent2", "parent2_to_parent1",
                            "incorrect_feature_assignment"):
                    row[key] = c[key]
                row.update(lost_parent1=n1 - t1, lost_parent2=n2 - t2,
                           retention_parent1=r1, retention_parent2=r2,
                           retention_difference=r1 - r2 if r1 is not None and r2 is not None else None,
                           misassignment_parent1_to_parent2_rate=ratio(c["parent1_to_parent2"], t1),
                           misassignment_parent2_to_parent1_rate=ratio(c["parent2_to_parent1"], t2),
                           observed_parent1=o1, observed_parent2=o2, observed_total=o1 + o2,
                           recovery=ratio(o1 + o2, total), observed_parent1_fraction=fraction,
                           fraction_error=error, absolute_fraction_error=abs(error) if error is not None else None)
                rows.append(row)
            summed = lambda key: sum(r[key] for r in rows)
            retained = summed("retained_true_parent1") + summed("retained_true_parent2")
            cross1 = summed("parent1_to_parent2")
            cross2 = summed("parent2_to_parent1")
            summary = dict(assay=assay, replicate=rep, library_id=library,
                           expected_total=summed("expected_total"), retained_unique_fragments=retained,
                           retained_true_parent1=summed("retained_true_parent1"),
                           retained_true_parent2=summed("retained_true_parent2"),
                           parent1_to_parent2=cross1, parent2_to_parent1=cross2,
                           cross_parent_misassignments=cross1 + cross2,
                           misassignment_rate=ratio(cross1 + cross2, retained),
                           final_assigned_fragments=summed("observed_total"),
                           recovery=ratio(summed("observed_total"), summed("expected_total")),
                           mapping_alignment_records=int(mapping.get("alignment_records", 0)),
                           mapping_unique_records=None, unique_fragments_before_deduplication=None,
                           duplicates_removed=None, fragments_after_deduplication=None,
                           featurecounts_assigned_fragments=None, rna_filter_retained_fragments=None)
            if assay == "RNA":
                summary["rna_filter_retained_fragments"] = int(filters["same_haplotype_pairs_kept"])
            else:
                unique = int(mapping["unique_records"])
                summary["mapping_unique_records"] = unique
                summary["featurecounts_assigned_fragments"] = int(filters["Assigned"])
                if assay == "MOA":
                    summary["retained_unique_fragments"] = unique
                else:
                    stats_path = args.hyarq_run / "results/qc/samtools" / slug / f"{library}.markdup.txt"
                    stats = {}
                    for line in stats_path.read_text().splitlines():
                        if ":" in line:
                            key, value = line.split(":", 1)
                            stats[key.strip()] = value.strip()
                    summary["unique_fragments_before_deduplication"] = unique // 2
                    summary["duplicates_removed"] = int(stats["DUPLICATE PAIR"]) // 2
                    summary["fragments_after_deduplication"] = int(stats["WRITTEN"]) // 2
            libraries.append(summary)
            assay_libraries.append(summary)
            evaluations.extend(rows)
            assay_rows.extend(rows)
        total1 = sum(r["retained_true_parent1"] for r in assay_libraries)
        total2 = sum(r["retained_true_parent2"] for r in assay_libraries)
        cross1 = sum(r["parent1_to_parent2"] for r in assay_libraries)
        cross2 = sum(r["parent2_to_parent1"] for r in assay_libraries)
        errors = [r["absolute_fraction_error"] for r in assay_rows
                  if r["absolute_fraction_error"] is not None]
        assays.append(dict(assay=assay, mean_recovery=mean([r["recovery"] for r in assay_rows]),
                           mean_absolute_fraction_error=mean(errors),
                           maximum_absolute_fraction_error=max(errors, default=None),
                           overall_parent1_to_parent2_rate=ratio(cross1, total1),
                           overall_parent2_to_parent1_rate=ratio(cross2, total2),
                           overall_cross_parent_misassignment_rate=ratio(cross1 + cross2, total1 + total2),
                           mean_absolute_parental_retention_difference=mean([
                               abs(r["retention_difference"]) for r in assay_rows
                               if r["retention_difference"] is not None])))
    args.output.mkdir(parents=True, exist_ok=True)
    for filename, rows in (("simulation_evaluation.tsv", evaluations),
                           ("library_summary.tsv", libraries), ("assay_summary.tsv", assays)):
        write_table(args.output / filename, rows)

    supplementary = []
    for row in evaluations:
        formatted = {key: row[key] for key in ("assay", "replicate", "feature_index")}
        for target, source, scale, precision in (
            ("expected_parent1_fraction", "expected_parent1_fraction", 1, 3),
            ("observed_parent1_fraction", "observed_parent1_fraction", 1, 3),
            ("absolute_fraction_error_pp", "absolute_fraction_error", 100, 2),
            ("recovery_percent", "recovery", 100, 2),
            ("retention_parent1_percent", "retention_parent1", 100, 2),
            ("retention_parent2_percent", "retention_parent2", 100, 2),
        ):
            value = row[source]
            formatted[target] = f"{value * scale:.{precision}f}" if value is not None else None
        for key in ("parent1_to_parent2", "parent2_to_parent1"):
            formatted[key] = int(row[key]) if row[key] is not None else None
        supplementary.append(formatted)
    write_table(args.output / "supplementary_simulation_validation.tsv", supplementary)


if __name__ == "__main__":
    main()
