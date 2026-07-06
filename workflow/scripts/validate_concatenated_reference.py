#!/usr/bin/env python3
"""
Validate a production concatenated HyARQ reference.

The validator checks structural integrity of the prefixed FASTA/GTF pair and,
optionally, feature tables used by downstream counting scripts.

Checks:
  1. FASTA contig names are unique and include both parent prefixes.
  2. GTF seqnames exist in the FASTA and feature coordinates are in bounds.
  3. Optional gene-pair tables reference gene IDs present in the GTF.
  4. Optional BED feature files reference FASTA contigs, have valid intervals,
     and use unique feature names in column 4.

Usage:
  python workflow/scripts/validate_concatenated_reference.py \
      --fasta results/reference/concatenated.fa \
      --gtf results/reference/concatenated.gtf \
      --gene-pairs gene_pairs.tsv \
      --bed atac_regions.bed \
      --bed chip_peaks.bed \
      --bed moa_footprints.bed
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path


def add_error(errors: list[str], message: str) -> None:
    errors.append(message)


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
    p1_prefix: str,
    p2_prefix: str,
    errors: list[str],
) -> dict[str, int]:
    p1_contigs = [name for name in fasta_lengths if name.startswith(p1_prefix)]
    p2_contigs = [name for name in fasta_lengths if name.startswith(p2_prefix)]
    if not p1_contigs:
        add_error(errors, f"[FASTA] no contigs start with P1 prefix {p1_prefix!r}")
    if not p2_contigs:
        add_error(errors, f"[FASTA] no contigs start with P2 prefix {p2_prefix!r}")
    if p1_prefix == p2_prefix:
        add_error(errors, "[prefix] P1 and P2 prefixes must differ")
    return {
        "P1_contigs": len(p1_contigs),
        "P2_contigs": len(p2_contigs),
        "P1_bp": sum(fasta_lengths[name] for name in p1_contigs),
        "P2_bp": sum(fasta_lengths[name] for name in p2_contigs),
    }


def validate_gtf(
    path: Path,
    fasta_lengths: dict[str, int],
    p1_prefix: str,
    p2_prefix: str,
    errors: list[str],
) -> tuple[set[str], dict[str, int]]:
    genes: set[str] = set()
    metrics = {
        "GTF_features": 0,
        "GTF_P1_features": 0,
        "GTF_P2_features": 0,
        "GTF_genes": 0,
    }

    try:
        handle = path.open()
    except OSError as exc:
        add_error(errors, f"[GTF] cannot open {path}: {exc}")
        return genes, metrics

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
            if contig.startswith(p1_prefix):
                metrics["GTF_P1_features"] += 1
            elif contig.startswith(p2_prefix):
                metrics["GTF_P2_features"] += 1
            else:
                add_error(
                    errors,
                    f"[GTF] contig lacks {p1_prefix!r} or {p2_prefix!r}: {contig}",
                )

            if contig not in fasta_lengths:
                add_error(errors, f"[GTF] contig not found in FASTA at {path}:{line_number}: {contig}")
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

    metrics["GTF_genes"] = len(genes)
    if metrics["GTF_features"] == 0:
        add_error(errors, f"[GTF] no feature rows found in {path}")
    if metrics["GTF_P1_features"] == 0:
        add_error(errors, f"[GTF] no features on {p1_prefix!r} contigs")
    if metrics["GTF_P2_features"] == 0:
        add_error(errors, f"[GTF] no features on {p2_prefix!r} contigs")
    return genes, metrics


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
            id_col = choose_column(reader.fieldnames, ("gene_pair_id", "pair_id", "gene_id"), path)
            p1_col = choose_column(reader.fieldnames, ("P1_gene_id", "P1_gene", "P1_feature"), path)
            p2_col = choose_column(reader.fieldnames, ("P2_gene_id", "P2_gene", "P2_feature"), path)
        except ValueError as exc:
            add_error(errors, f"[gene-pairs] {exc}")
            return metrics

        seen_ids: set[str] = set()
        for row_number, row in enumerate(reader, start=2):
            pair_id = row[id_col]
            if not pair_id or pair_id in seen_ids:
                add_error(errors, f"[gene-pairs] blank or duplicate pair ID at {path}:{row_number}: {pair_id!r}")
            seen_ids.add(pair_id)
            for label, column in (("P1", p1_col), ("P2", p2_col)):
                gene_id = row[column]
                if gene_id != "." and gene_id not in gtf_genes:
                    add_error(errors, f"[gene-pairs] {label} gene not found in GTF at {path}:{row_number}: {gene_id}")
            metrics["gene_pair_rows"] += 1
    return metrics


def validate_bed(
    path: Path,
    fasta_lengths: dict[str, int],
    p1_prefix: str,
    p2_prefix: str,
    errors: list[str],
) -> dict[str, int]:
    metrics = {
        f"{path.name}_rows": 0,
        f"{path.name}_P1_rows": 0,
        f"{path.name}_P2_rows": 0,
    }
    seen_names: set[str] = set()

    try:
        handle = path.open()
    except OSError as exc:
        add_error(errors, f"[BED] cannot open {path}: {exc}")
        return metrics

    with handle:
        for line_number, raw_line in enumerate(handle, start=1):
            if not raw_line.strip() or raw_line.startswith("#"):
                continue
            fields = raw_line.rstrip("\n").split("\t")
            if len(fields) < 3:
                add_error(errors, f"[BED] malformed row at {path}:{line_number}")
                continue
            contig = fields[0]
            if contig.startswith(p1_prefix):
                metrics[f"{path.name}_P1_rows"] += 1
            elif contig.startswith(p2_prefix):
                metrics[f"{path.name}_P2_rows"] += 1
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
            if len(fields) >= 4:
                name = fields[3]
                if name in seen_names:
                    add_error(errors, f"[BED] duplicate feature name in {path}: {name}")
                seen_names.add(name)
            metrics[f"{path.name}_rows"] += 1

    if metrics[f"{path.name}_rows"] == 0:
        add_error(errors, f"[BED] no feature rows found in {path}")
    return metrics


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate a concatenated HyARQ reference.")
    parser.add_argument("--fasta", default="results/reference/concatenated.fa")
    parser.add_argument("--gtf", default="results/reference/concatenated.gtf")
    parser.add_argument("--p1-prefix", default="P1_")
    parser.add_argument("--p2-prefix", default="P2_")
    parser.add_argument("--gene-pairs", default=None)
    parser.add_argument("--bed", action="append", default=[], help="BED feature file to validate; may be repeated")
    args = parser.parse_args()

    errors: list[str] = []
    metrics: dict[str, int] = {}

    fasta_path = Path(args.fasta)
    gtf_path = Path(args.gtf)
    fasta_lengths = load_fasta_lengths(fasta_path, errors)
    metrics["FASTA_contigs"] = len(fasta_lengths)
    metrics["FASTA_bp"] = sum(fasta_lengths.values())
    metrics.update(validate_prefixes(fasta_lengths, args.p1_prefix, args.p2_prefix, errors))

    gtf_genes, gtf_metrics = validate_gtf(gtf_path, fasta_lengths, args.p1_prefix, args.p2_prefix, errors)
    metrics.update(gtf_metrics)

    if args.gene_pairs:
        metrics.update(validate_gene_pairs(Path(args.gene_pairs), gtf_genes, errors))

    for bed_path in args.bed:
        metrics.update(validate_bed(Path(bed_path), fasta_lengths, args.p1_prefix, args.p2_prefix, errors))

    print("Validation metrics:")
    for metric, value in metrics.items():
        print(f"  {metric}\t{value}")

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
