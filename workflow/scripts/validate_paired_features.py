#!/usr/bin/env python3

import argparse
import csv
import sys
from pathlib import Path

GENE_HEADER = ["parent1_gene_id", "parent2_gene_id"]
COORD_HEADER = [
    "parent1_chrom", "parent1_start", "parent1_end",
    "parent2_chrom", "parent2_start", "parent2_end",
]


def read_table(path, expected_header, errors):
    # returns list of row dicts, or None if unreadable / wrong header
    try:
        fh = open(path, newline="")
    except OSError as exc:
        errors.append(f"{path}: cannot open ({exc})")
        return None
    with fh:
        reader = csv.reader(fh, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration:
            errors.append(f"{path}: empty file")
            return None
        if header != expected_header:
            errors.append(f"{path}: header must be exactly {expected_header}, got {header}")
            return None
        rows = []
        for line_number, row in enumerate(reader, start=2):
            if not row:
                continue
            if len(row) != len(expected_header):
                errors.append(
                    f"{path}:{line_number}: expected {len(expected_header)} columns, "
                    f"found {len(row)}"
                )
                continue
            rows.append(dict(zip(header, row)))
        return rows


def check_no_blanks(rows, columns, path, errors):
    for i, row in enumerate(rows, start=2):
        for col in columns:
            if not row[col].strip():
                errors.append(f"{path}:{i}: blank {col}")


def validate_gene_pairs(path, errors):
    rows = read_table(path, GENE_HEADER, errors)
    if rows is None:
        return 0
    if not rows:
        errors.append(f"{path}: no data rows")
        return 0
    check_no_blanks(rows, GENE_HEADER, path, errors)
    for col in GENE_HEADER:
        seen = set()
        duplicates = set()
        for row in rows:
            gene_id = row[col]
            if not gene_id.strip():
                continue
            if gene_id in seen:
                duplicates.add(gene_id)
            else:
                seen.add(gene_id)
        for gene_id in sorted(duplicates):
            errors.append(f"{path}: duplicate {col}: {gene_id}")
    return len(rows)


def validate_coord_pairs(path, errors):
    rows = read_table(path, COORD_HEADER, errors)
    if rows is None:
        return 0
    if not rows:
        errors.append(f"{path}: no data rows")
        return 0
    check_no_blanks(rows, COORD_HEADER, path, errors)
    for i, row in enumerate(rows, start=2):
        for side in ("parent1", "parent2"):
            try:
                start = int(row[f"{side}_start"])
                end = int(row[f"{side}_end"])
            except ValueError:
                errors.append(f"{path}:{i}: {side} coordinates not integer")
                continue
            if not (0 <= start < end):
                errors.append(f"{path}:{i}: {side} interval {start}-{end} not valid (need 0 <= start < end)")
    return len(rows)


def main():
    ap = argparse.ArgumentParser(description="Format check for paired-feature tables.")
    ap.add_argument("--gene-pairs")
    ap.add_argument("--peak-pairs")
    ap.add_argument("--footprint-pairs")
    ap.add_argument("--report", default="results/qc/paired_features_validation.tsv")
    args = ap.parse_args()

    errors = []
    metrics = {}

    if args.gene_pairs:
        metrics["gene_pair_rows"] = validate_gene_pairs(args.gene_pairs, errors)
    if args.peak_pairs:
        metrics["peak_pair_rows"] = validate_coord_pairs(args.peak_pairs, errors)
    if args.footprint_pairs:
        metrics["footprint_pair_rows"] = validate_coord_pairs(args.footprint_pairs, errors)

    if not (args.gene_pairs or args.peak_pairs or args.footprint_pairs):
        errors.append("no pair tables given (need at least one of --gene-pairs/--peak-pairs/--footprint-pairs)")

    metrics["error_count"] = len(errors)
    metrics["validation_result"] = "FAIL" if errors else "PASS"

    report = Path(args.report)
    report.parent.mkdir(parents=True, exist_ok=True)
    with report.open("w", newline="") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(("metric", "value"))
        for key, value in metrics.items():
            writer.writerow((key, value))

    for key, value in metrics.items():
        print(f"  {key}\t{value}")
    if errors:
        print(f"\nFAILURES ({len(errors)}):")
        for err in errors[:50]:
            print(f"  {err}")
        if len(errors) > 50:
            print(f"  ... and {len(errors) - 50} more")
        print("\nRESULT: FAIL")
        return 1
    print("\nRESULT: PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
