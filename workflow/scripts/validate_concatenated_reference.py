#!/usr/bin/env python3
"""
Validate a production HyARQ concatenated reference.

The validator checks the prefixed FASTA/GTF pair produced by
build_concatenated_reference.sh and, optionally, feature tables used by
downstream counting scripts. parent1/parent2 are internal roles; FASTA and GTF
contigs are expected to use user-provided parental prefixes such as B73_ and
Mo17_.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path


def add_error(errors: list[str], message: str) -> None:
    errors.append(message)


def add_warning(warnings: list[str], message: str) -> None:
    warnings.append(message)


def metric_key(path: Path, suffix: str) -> str:
    safe_name = re.sub(r"[^A-Za-z0-9]+", "_", path.name).strip("_")
    return f"{safe_name}_{suffix}" if safe_name else suffix


def contig_prefix(parent_prefix: str | None, separator: str) -> str | None:
    if not parent_prefix:
        return None
    return parent_prefix if parent_prefix.endswith(separator) else f"{parent_prefix}{separator}"


def load_fasta_lengths(path: Path, errors: list[str]) -> dict[str, int]:
    lengths: dict[str, int] = {}
    current_name: str | None = None
    current_length = 0

    try:
        handle = path.open()
    except OSError as exc:
        add_error(errors, f"[FASTA] cannot open {path}: {exc}")
        return lengths

    with handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\n")
            if line.startswith(">"):
                if current_name is not None:
                    if current_name in lengths:
                        add_error(errors, f"[FASTA] duplicate contig name: {current_name}")
                    else:
                        lengths[current_name] = current_length
                current_name = line[1:].split()[0]
                current_length = 0
                if not current_name:
                    add_error(errors, f"[FASTA] blank contig name at {path}:{line_number}")
            elif current_name is None:
                if line.strip():
                    add_error(errors, f"[FASTA] sequence before first header at {path}:{line_number}")
            else:
                current_length += len(line.strip())

    if current_name is not None:
        if current_name in lengths:
            add_error(errors, f"[FASTA] duplicate contig name: {current_name}")
        else:
            lengths[current_name] = current_length

    if not lengths:
        add_error(errors, f"[FASTA] no contigs found in {path}")
    return lengths


def parse_gtf_gene_id(attributes: str) -> str | None:
    match = re.search(r'(?:^|;\s*)gene_id\s+"([^"]+)"', attributes)
    return match.group(1) if match else None


def validate_prefixes(
    fasta_lengths: dict[str, int],
    parent1_contig_prefix: str,
    parent2_contig_prefix: str,
    errors: list[str],
) -> dict[str, int]:
    parent1_contigs = [
        name for name in fasta_lengths if name.startswith(parent1_contig_prefix)
    ]
    parent2_contigs = [
        name for name in fasta_lengths if name.startswith(parent2_contig_prefix)
    ]

    if not parent1_contigs:
        add_error(
            errors,
            f"[FASTA] no contigs start with parent1 prefix {parent1_contig_prefix!r}",
        )
    if not parent2_contigs:
        add_error(
            errors,
            f"[FASTA] no contigs start with parent2 prefix {parent2_contig_prefix!r}",
        )
    if parent1_contig_prefix == parent2_contig_prefix:
        add_error(errors, "[prefix] parent1 and parent2 prefixes must differ")

    return {
        "parent1_contigs": len(parent1_contigs),
        "parent2_contigs": len(parent2_contigs),
        "parent1_bp": sum(fasta_lengths[name] for name in parent1_contigs),
        "parent2_bp": sum(fasta_lengths[name] for name in parent2_contigs),
    }


def validate_gtf(
    path: Path,
    fasta_lengths: dict[str, int],
    parent1_contig_prefix: str,
    parent2_contig_prefix: str,
    errors: list[str],
) -> tuple[set[str], dict[str, int], dict[str, dict[str, object]]]:
    genes: set[str] = set()
    gene_details: dict[str, dict[str, object]] = {}
    metrics = {
        "GTF_features": 0,
        "GTF_parent1_features": 0,
        "GTF_parent2_features": 0,
        "GTF_genes": 0,
    }

    try:
        handle = path.open()
    except OSError as exc:
        add_error(errors, f"[GTF] cannot open {path}: {exc}")
        return genes, metrics, gene_details

    with handle:
        for line_number, raw_line in enumerate(handle, start=1):
            if not raw_line.strip() or raw_line.startswith("#"):
                continue
            fields = raw_line.rstrip("\n").split("\t")
            if len(fields) < 9:
                add_error(errors, f"[GTF] malformed row at {path}:{line_number}")
                continue

            contig = fields[0]
            metrics["GTF_features"] += 1
            parent_role: str | None = None
            if contig.startswith(parent1_contig_prefix):
                parent_role = "parent1"
                metrics["GTF_parent1_features"] += 1
            elif contig.startswith(parent2_contig_prefix):
                parent_role = "parent2"
                metrics["GTF_parent2_features"] += 1
            else:
                add_error(
                    errors,
                    f"[GTF] contig lacks parent prefix at {path}:{line_number}: {contig}",
                )

            if contig not in fasta_lengths:
                add_error(
                    errors,
                    f"[GTF] contig not found in FASTA at {path}:{line_number}: {contig}",
                )
                continue

            try:
                start = int(fields[3])
                end = int(fields[4])
            except ValueError:
                add_error(errors, f"[GTF] non-integer coordinates at {path}:{line_number}")
                continue
            if start < 1 or end < start or end > fasta_lengths[contig]:
                add_error(
                    errors,
                    f"[GTF] interval out of bounds at {path}:{line_number}: "
                    f"{contig}:{start}-{end} length={fasta_lengths[contig]}",
                )

            gene_id = parse_gtf_gene_id(fields[8])
            if gene_id:
                genes.add(gene_id)
                detail = gene_details.setdefault(
                    gene_id,
                    {"parents": set(), "contigs": set(), "gene_feature_rows": 0},
                )
                if parent_role:
                    detail["parents"].add(parent_role)  # type: ignore[union-attr]
                detail["contigs"].add(contig)  # type: ignore[union-attr]
                if fields[2] == "gene":
                    detail["gene_feature_rows"] = int(detail["gene_feature_rows"]) + 1

    metrics["GTF_genes"] = len(genes)
    if metrics["GTF_features"] == 0:
        add_error(errors, f"[GTF] no feature rows found in {path}")
    if metrics["GTF_parent1_features"] == 0:
        add_error(errors, f"[GTF] no features on parent1 contigs")
    if metrics["GTF_parent2_features"] == 0:
        add_error(errors, f"[GTF] no features on parent2 contigs")

    return genes, metrics, gene_details


def summarize_duplicate_genes(
    gene_details: dict[str, dict[str, object]],
) -> tuple[set[str], set[str]]:
    duplicate_gene_ids: set[str] = set()
    cross_parent_duplicate_gene_ids: set[str] = set()

    for gene_id, detail in gene_details.items():
        parents = detail["parents"]
        contigs = detail["contigs"]
        gene_feature_rows = int(detail["gene_feature_rows"])

        if isinstance(parents, set) and len(parents) > 1:
            duplicate_gene_ids.add(gene_id)
            cross_parent_duplicate_gene_ids.add(gene_id)
        elif isinstance(contigs, set) and len(contigs) > 1:
            duplicate_gene_ids.add(gene_id)
        elif gene_feature_rows > 1:
            duplicate_gene_ids.add(gene_id)

    return duplicate_gene_ids, cross_parent_duplicate_gene_ids


def choose_column(header: list[str], candidates: tuple[str, ...], path: Path) -> str:
    for candidate in candidates:
        if candidate in header:
            return candidate
    raise ValueError(f"{path}: expected one of columns {', '.join(candidates)}")


def validate_gene_pairs(path: Path, gtf_genes: set[str], errors: list[str]) -> dict[str, int]:
    metrics = {"gene_pair_rows": 0}
    try:
        handle = path.open(newline="")
    except OSError as exc:
        add_error(errors, f"[gene-pairs] cannot open {path}: {exc}")
        return metrics

    with handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            add_error(errors, f"[gene-pairs] missing header in {path}")
            return metrics
        try:
            id_col = choose_column(
                reader.fieldnames,
                ("gene_pair_id", "pair_id", "gene_id"),
                path,
            )
            parent1_col = choose_column(
                reader.fieldnames,
                ("parent1_gene_id", "P1_gene_id", "P1_gene", "P1_feature"),
                path,
            )
            parent2_col = choose_column(
                reader.fieldnames,
                ("parent2_gene_id", "P2_gene_id", "P2_gene", "P2_feature"),
                path,
            )
        except ValueError as exc:
            add_error(errors, f"[gene-pairs] {exc}")
            return metrics

        seen_ids: set[str] = set()
        for row_number, row in enumerate(reader, start=2):
            pair_id = row[id_col]
            if not pair_id or pair_id in seen_ids:
                add_error(
                    errors,
                    f"[gene-pairs] blank or duplicate pair ID at {path}:{row_number}: {pair_id!r}",
                )
            seen_ids.add(pair_id)
            for label, column in (("parent1", parent1_col), ("parent2", parent2_col)):
                gene_id = row[column]
                if gene_id != "." and gene_id not in gtf_genes:
                    add_error(
                        errors,
                        f"[gene-pairs] {label} gene not found in GTF at {path}:{row_number}: {gene_id}",
                    )
            metrics["gene_pair_rows"] += 1
    return metrics


def validate_bed(
    path: Path,
    fasta_lengths: dict[str, int],
    parent1_contig_prefix: str,
    parent2_contig_prefix: str,
    errors: list[str],
) -> tuple[dict[str, int], set[str]]:
    row_key = metric_key(path, "BED_rows")
    parent1_key = metric_key(path, "BED_parent1_rows")
    parent2_key = metric_key(path, "BED_parent2_rows")
    metrics = {
        row_key: 0,
        parent1_key: 0,
        parent2_key: 0,
    }
    feature_names: set[str] = set()

    try:
        handle = path.open()
    except OSError as exc:
        add_error(errors, f"[BED] cannot open {path}: {exc}")
        return metrics, feature_names

    with handle:
        for line_number, raw_line in enumerate(handle, start=1):
            if not raw_line.strip() or raw_line.startswith("#"):
                continue
            fields = raw_line.rstrip("\n").split("\t")
            if len(fields) < 3:
                add_error(errors, f"[BED] malformed row at {path}:{line_number}")
                continue
            contig = fields[0]
            if contig.startswith(parent1_contig_prefix):
                metrics[parent1_key] += 1
            elif contig.startswith(parent2_contig_prefix):
                metrics[parent2_key] += 1
            else:
                add_error(errors, f"[BED] contig lacks parent prefix at {path}:{line_number}: {contig}")

            if contig not in fasta_lengths:
                add_error(errors, f"[BED] contig not found in FASTA at {path}:{line_number}: {contig}")
                continue
            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError:
                add_error(errors, f"[BED] non-integer coordinates at {path}:{line_number}")
                continue
            if start < 0 or end <= start or end > fasta_lengths[contig]:
                add_error(
                    errors,
                    f"[BED] interval out of bounds at {path}:{line_number}: "
                    f"{contig}:{start}-{end} length={fasta_lengths[contig]}",
                )

            if len(fields) < 4 or not fields[3]:
                add_error(errors, f"[BED] missing feature name in column 4 at {path}:{line_number}")
            else:
                name = fields[3]
                if name in feature_names:
                    add_error(errors, f"[BED] duplicate feature name in {path}: {name}")
                feature_names.add(name)
            metrics[row_key] += 1

    if metrics[row_key] == 0:
        add_error(errors, f"[BED] no feature rows found in {path}")
    if metrics[parent1_key] == 0:
        add_error(errors, f"[BED] no rows on parent1 contigs in {path}")
    if metrics[parent2_key] == 0:
        add_error(errors, f"[BED] no rows on parent2 contigs in {path}")

    return metrics, feature_names


def validate_feature_pairs(
    path: Path,
    table_kind: str,
    parent1_candidates: tuple[str, ...],
    parent2_candidates: tuple[str, ...],
    bed_feature_names: set[str],
    errors: list[str],
) -> dict[str, int]:
    metrics = {f"{table_kind}_pair_rows": 0}
    try:
        handle = path.open(newline="")
    except OSError as exc:
        add_error(errors, f"[{table_kind}-pairs] cannot open {path}: {exc}")
        return metrics

    with handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            add_error(errors, f"[{table_kind}-pairs] missing header in {path}")
            return metrics
        try:
            parent1_col = choose_column(reader.fieldnames, parent1_candidates, path)
            parent2_col = choose_column(reader.fieldnames, parent2_candidates, path)
        except ValueError as exc:
            add_error(errors, f"[{table_kind}-pairs] {exc}")
            return metrics

        for row_number, row in enumerate(reader, start=2):
            parent1_feature = row[parent1_col]
            parent2_feature = row[parent2_col]
            if not parent1_feature:
                add_error(errors, f"[{table_kind}-pairs] blank parent1 feature at {path}:{row_number}")
            if not parent2_feature:
                add_error(errors, f"[{table_kind}-pairs] blank parent2 feature at {path}:{row_number}")
            if bed_feature_names:
                if parent1_feature and parent1_feature not in bed_feature_names:
                    add_error(
                        errors,
                        f"[{table_kind}-pairs] parent1 feature not found in provided BED features "
                        f"at {path}:{row_number}: {parent1_feature}",
                    )
                if parent2_feature and parent2_feature not in bed_feature_names:
                    add_error(
                        errors,
                        f"[{table_kind}-pairs] parent2 feature not found in provided BED features "
                        f"at {path}:{row_number}: {parent2_feature}",
                    )
            metrics[f"{table_kind}_pair_rows"] += 1

    return metrics


def write_report(path: Path, metrics: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as output:
        writer = csv.writer(output, delimiter="\t", lineterminator="\n")
        writer.writerow(("metric", "value"))
        for metric, value in metrics.items():
            writer.writerow((metric, value))


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate a concatenated HyARQ reference.")
    parser.add_argument("--fasta", default="results/reference/concatenated.fa")
    parser.add_argument("--gtf", default="results/reference/concatenated.gtf")
    parser.add_argument("--parent1-prefix", "--p1-prefix", dest="parent1_prefix")
    parser.add_argument("--parent2-prefix", "--p2-prefix", dest="parent2_prefix")
    parser.add_argument("--separator", default="_")
    parser.add_argument("--gene-pairs", default=None)
    parser.add_argument("--bed", action="append", default=[], help="BED feature file to validate; may be repeated")
    parser.add_argument("--region-pairs", default=None)
    parser.add_argument("--peak-pairs", default=None)
    parser.add_argument("--footprint-pairs", default=None)
    parser.add_argument("--allow-duplicate-gene-ids", action="store_true")
    parser.add_argument("--report", default="results/reference/reference_validation_report.tsv")
    args = parser.parse_args()

    errors: list[str] = []
    warnings: list[str] = []
    metrics: dict[str, object] = {
        "validation_result": "FAIL",
        "fasta": args.fasta,
        "gtf": args.gtf,
        "parent1_prefix": args.parent1_prefix or "",
        "parent2_prefix": args.parent2_prefix or "",
    }

    if not args.parent1_prefix or not args.parent2_prefix:
        add_error(
            errors,
            "set --parent1-prefix and --parent2-prefix to real parental names "
            "(legacy aliases --p1-prefix and --p2-prefix are accepted)",
        )

    parent1_contig_prefix = contig_prefix(args.parent1_prefix, args.separator)
    parent2_contig_prefix = contig_prefix(args.parent2_prefix, args.separator)

    fasta_path = Path(args.fasta)
    gtf_path = Path(args.gtf)
    fasta_lengths = load_fasta_lengths(fasta_path, errors)
    metrics["FASTA_contigs"] = len(fasta_lengths)
    metrics["FASTA_bp"] = sum(fasta_lengths.values())

    gtf_genes: set[str] = set()
    gene_details: dict[str, dict[str, object]] = {}

    if parent1_contig_prefix and parent2_contig_prefix:
        metrics.update(
            validate_prefixes(
                fasta_lengths,
                parent1_contig_prefix,
                parent2_contig_prefix,
                errors,
            )
        )
        gtf_genes, gtf_metrics, gene_details = validate_gtf(
            gtf_path,
            fasta_lengths,
            parent1_contig_prefix,
            parent2_contig_prefix,
            errors,
        )
        metrics.update(gtf_metrics)
    else:
        metrics.update(
            {
                "parent1_contigs": 0,
                "parent2_contigs": 0,
                "parent1_bp": 0,
                "parent2_bp": 0,
                "GTF_features": 0,
                "GTF_parent1_features": 0,
                "GTF_parent2_features": 0,
                "GTF_genes": 0,
            }
        )

    duplicate_gene_ids, cross_parent_duplicate_gene_ids = summarize_duplicate_genes(gene_details)
    metrics["duplicate_gene_ids"] = len(duplicate_gene_ids)
    metrics["cross_parent_duplicate_gene_ids"] = len(cross_parent_duplicate_gene_ids)
    metrics["duplicate_gene_id_examples"] = ", ".join(sorted(duplicate_gene_ids)[:10])
    metrics["cross_parent_duplicate_gene_id_examples"] = ", ".join(
        sorted(cross_parent_duplicate_gene_ids)[:10]
    )

    if duplicate_gene_ids:
        examples = ", ".join(sorted(duplicate_gene_ids)[:10])
        message = (
            f"duplicate gene IDs detected ({len(duplicate_gene_ids)}); examples: {examples}. "
            "Duplicate gene IDs can cause featureCounts to merge parent-specific features."
        )
        if args.allow_duplicate_gene_ids:
            add_warning(warnings, message)
        else:
            add_error(errors, message)

    if cross_parent_duplicate_gene_ids:
        examples = ", ".join(sorted(cross_parent_duplicate_gene_ids)[:10])
        message = (
            f"gene IDs shared across parent1 and parent2 ({len(cross_parent_duplicate_gene_ids)}); "
            f"examples: {examples}"
        )
        if args.allow_duplicate_gene_ids:
            add_warning(warnings, message)
        else:
            add_error(errors, message)

    if args.gene_pairs:
        metrics.update(validate_gene_pairs(Path(args.gene_pairs), gtf_genes, errors))

    bed_feature_names: set[str] = set()
    if parent1_contig_prefix and parent2_contig_prefix:
        for bed_path_text in args.bed:
            bed_metrics, feature_names = validate_bed(
                Path(bed_path_text),
                fasta_lengths,
                parent1_contig_prefix,
                parent2_contig_prefix,
                errors,
            )
            metrics.update(bed_metrics)
            bed_feature_names.update(feature_names)
    elif args.bed:
        add_error(errors, "cannot validate BED files without parent1/parent2 prefixes")

    if args.region_pairs:
        metrics.update(
            validate_feature_pairs(
                Path(args.region_pairs),
                "region",
                ("parent1_region", "P1_region"),
                ("parent2_region", "P2_region"),
                bed_feature_names,
                errors,
            )
        )

    if args.peak_pairs:
        metrics.update(
            validate_feature_pairs(
                Path(args.peak_pairs),
                "peak",
                ("parent1_peak", "P1_peak"),
                ("parent2_peak", "P2_peak"),
                bed_feature_names,
                errors,
            )
        )

    if args.footprint_pairs:
        metrics.update(
            validate_feature_pairs(
                Path(args.footprint_pairs),
                "footprint",
                ("parent1_footprint", "P1_footprint"),
                ("parent2_footprint", "P2_footprint"),
                bed_feature_names,
                errors,
            )
        )

    metrics["error_count"] = len(errors)
    metrics["warning_count"] = len(warnings)
    metrics["validation_result"] = "FAIL" if errors else "PASS"

    report_path = Path(args.report)
    try:
        write_report(report_path, metrics)
    except OSError as exc:
        add_error(errors, f"[report] cannot write {report_path}: {exc}")
        metrics["error_count"] = len(errors)
        metrics["validation_result"] = "FAIL"

    print("Validation metrics:")
    for metric, value in metrics.items():
        print(f"  {metric}\t{value}")
    print(f"\nReport: {report_path}")

    if warnings:
        print(f"\nWARNINGS: {len(warnings)}")
        for warning in warnings[:50]:
            print(f"  {warning}")
        if len(warnings) > 50:
            print(f"  ... and {len(warnings) - 50} more")

    if errors:
        print(f"\nFAILURES: {len(errors)}")
        for error in errors[:50]:
            print(f"  {error}")
        if len(errors) > 50:
            print(f"  ... and {len(errors) - 50} more")
        print("\nRESULT: FAIL")
        return 1

    print("\nRESULT: PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
