#!/usr/bin/env python3

import argparse
import sys
from pathlib import Path

import pandas as pd

from rna_gene_qc_plots import (
    plot_color,
    write_rna_gene_qc_report,
)


EXPECTED_COLUMNS = [
    "parent1_gene_id",
    "parent2_gene_id",
    "parent1_count",
    "parent2_count",
    "total_count",
]
COUNT_COLUMNS = ["parent1_count", "parent2_count", "total_count"]
SUMMARY_COLUMNS = [
    "library_id",
    "total_gene_pairs",
    "gene_pairs_with_counts",
    "parent1_only_gene_pairs",
    "parent2_only_gene_pairs",
    "total_count",
    "overall_parent1_fraction",
    "median_gene_pair_parent1_fraction",
]
PAIRED_COLUMNS = ["parent1_gene_id", "parent2_gene_id"]


class ErrorParser(argparse.ArgumentParser):
    def error(self, message):
        self.exit(2, f"ERROR: {message}\n")


def fail(errors):
    for error in errors:
        print(f"ERROR: {error}", file=sys.stderr)
    return 1


def positive_integer(value):
    try:
        parsed = int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("must be an integer greater than or equal to 1") from exc
    if parsed < 1:
        raise argparse.ArgumentTypeError("must be an integer greater than or equal to 1")
    return parsed


def read_count_table(path, errors):
    try:
        table = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    except OSError as exc:
        errors.append(f"{path}: cannot open ({exc})")
        return None
    except pd.errors.EmptyDataError:
        errors.append(f"{path}: empty file")
        return None
    except pd.errors.ParserError as exc:
        errors.append(f"{path}: cannot parse TSV ({exc})")
        return None

    if list(table.columns) != EXPECTED_COLUMNS:
        errors.append(f"{path}: header must be exactly {EXPECTED_COLUMNS}, got {list(table.columns)}")
        return None

    return table


def validate_gene_ids(table, path, errors):
    for column in PAIRED_COLUMNS:
        blank_rows = table.index[table[column].str.strip().eq("")]
        for row_index in blank_rows:
            errors.append(f"{path}:{row_index + 2}: blank {column}")


def parse_count_column(table, column, path, errors):
    values = table[column]
    invalid = ~values.str.fullmatch(r"\d+")
    for row_index in table.index[invalid]:
        errors.append(f"{path}:{row_index + 2}: {column} must be a nonnegative integer")
    parsed = pd.to_numeric(values.where(~invalid, "0"), downcast="integer")
    return parsed.astype("int64"), invalid


def validate_counts(table, path, errors):
    invalid_masks = []
    for column in COUNT_COLUMNS:
        parsed, invalid = parse_count_column(table, column, path, errors)
        table[column] = parsed
        invalid_masks.append(invalid)

    invalid_counts = invalid_masks[0] | invalid_masks[1] | invalid_masks[2]
    mismatches = (~invalid_counts) & (table["total_count"] != table["parent1_count"] + table["parent2_count"])
    for row_index in table.index[mismatches]:
        errors.append(f"{path}:{row_index + 2}: total_count does not equal parent1_count + parent2_count")


def validate_duplicate_pairs(table, path, errors):
    duplicates = table.duplicated(PAIRED_COLUMNS, keep=False)
    duplicate_pairs = table.loc[duplicates, PAIRED_COLUMNS].drop_duplicates()
    for pair in duplicate_pairs.itertuples(index=False):
        errors.append(f"{path}: duplicate gene pair: {pair.parent1_gene_id}/{pair.parent2_gene_id}")


def validate_table(path):
    errors = []
    table = read_count_table(path, errors)
    if table is None:
        return None, errors

    validate_gene_ids(table, path, errors)
    validate_counts(table, path, errors)
    validate_duplicate_pairs(table, path, errors)
    return table, errors


def pair_set(table):
    return set(map(tuple, table[PAIRED_COLUMNS].itertuples(index=False, name=None)))


def fraction_value(value):
    return "NA" if pd.isna(value) else f"{value:.6f}"


def summarize_table(library_id, table):
    counted = table["total_count"] > 0
    total_count = int(table["total_count"].sum())

    if counted.any():
        overall_parent1_fraction = table["parent1_count"].sum() / total_count
        median_parent1_fraction = (table.loc[counted, "parent1_count"] / table.loc[counted, "total_count"]).median()
    else:
        overall_parent1_fraction = pd.NA
        median_parent1_fraction = pd.NA

    return {
        "library_id": library_id,
        "total_gene_pairs": len(table),
        "gene_pairs_with_counts": int(counted.sum()),
        "parent1_only_gene_pairs": int(((table["parent1_count"] > 0) & (table["parent2_count"] == 0)).sum()),
        "parent2_only_gene_pairs": int(((table["parent1_count"] == 0) & (table["parent2_count"] > 0)).sum()),
        "total_count": total_count,
        "overall_parent1_fraction": fraction_value(overall_parent1_fraction),
        "median_gene_pair_parent1_fraction": fraction_value(median_parent1_fraction),
    }


def build_pooled_table(tables):
    pooled = pd.concat(tables, ignore_index=True)
    pooled = (
        pooled.groupby(PAIRED_COLUMNS, sort=False, as_index=False)[["parent1_count", "parent2_count"]]
        .sum()
    )
    pooled["total_count"] = pooled["parent1_count"] + pooled["parent2_count"]
    return pooled


def summarize_pooled(pooled):
    return summarize_table("pooled_all_libraries", pooled)


def write_summary(rows, path):
    summary = pd.DataFrame(rows, columns=SUMMARY_COLUMNS)
    path.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(path, sep="\t", index=False, lineterminator="\n")


def main():
    ap = ErrorParser(description="Summarize RNA gene count QC metrics.")
    ap.add_argument("--count-tables", nargs="+", required=True)
    ap.add_argument("--summary", required=True)
    ap.add_argument("--report", required=True)
    ap.add_argument(
        "--min-total-count",
        type=positive_integer,
        default=1,
        help="Minimum total count used by threshold-dependent plots; Panel A always includes pooled gene pairs with total_count > 0.",
    )
    ap.add_argument(
        "--histogram-color",
        type=plot_color,
        default="tab:blue",
        help="Color used for parental count proportion histogram bars.",
    )
    ap.add_argument(
        "--histogram-bins",
        type=positive_integer,
        default=20,
        help="Number of equal-width bins used for the Panel A parental count proportion histogram.",
    )
    ap.add_argument("--sample-histogram-columns", type=positive_integer, default=3)
    ap.add_argument("--sample-histogram-rows", type=positive_integer, default=2)
    ap.add_argument("--parent1-label", default="Parent1", help="Parent1 name displayed in the report.")
    ap.add_argument("--parent2-label", default="Parent2", help="Parent2 name displayed in the report.")
    args = ap.parse_args()

    if not args.count_tables:
        return fail(["at least one count table is required"])

    paths = [Path(path) for path in args.count_tables]
    library_ids = [path.parent.name for path in paths]
    duplicate_ids = sorted({library_id for library_id in library_ids if library_ids.count(library_id) > 1})
    if duplicate_ids:
        return fail([f"derived library IDs are not unique: {', '.join(duplicate_ids)}"])

    sorted_inputs = sorted(zip(library_ids, paths), key=lambda item: item[0])
    library_ids = [library_id for library_id, _ in sorted_inputs]
    paths = [path for _, path in sorted_inputs]

    errors = []
    tables = []
    for path in paths:
        table, table_errors = validate_table(path)
        errors.extend(table_errors)
        if table is not None:
            tables.append(table)

    if errors:
        return fail(errors)

    reference_pairs = pair_set(tables[0])
    for path, table in zip(paths[1:], tables[1:]):
        if pair_set(table) != reference_pairs:
            errors.append(f"{path}: gene pair set does not match the first count table")
    if errors:
        return fail(errors)

    pooled = build_pooled_table(tables)
    rows = [summarize_table(library_id, table) for library_id, table in zip(library_ids, tables)]
    rows.append(summarize_pooled(pooled))
    write_summary(rows, Path(args.summary))
    write_rna_gene_qc_report(
        library_ids,
        tables,
        pooled,
        Path(args.report),
        args.min_total_count,
        args.histogram_color,
        args.sample_histogram_rows,
        args.sample_histogram_columns,
        parent1_label=args.parent1_label,
        parent2_label=args.parent2_label,
        histogram_bins=args.histogram_bins,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
