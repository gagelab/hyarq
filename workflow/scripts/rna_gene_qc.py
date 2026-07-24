#!/usr/bin/env python3

import argparse
import math
import sys
from pathlib import Path

import pandas as pd

from rna_gene_qc_plots import (
    RNA_RETENTION_THRESHOLDS,
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
    "sample_id",
    "total_gene_pairs",
    "gene_pairs_with_any_count",
    "gene_pairs_with_parent1_count_only",
    "gene_pairs_with_parent2_count_only",
    "gene_pairs_with_both_parent_counts_nonzero",
    "total_parental_count",
    "overall_parent1_proportion",
    "median_gene_pair_parent1_proportion",
    *[
        f"retained_ge_{threshold}"
        for threshold in RNA_RETENTION_THRESHOLDS
    ],
]
PAIRED_COLUMNS = ["parent1_gene_id", "parent2_gene_id"]
RETENTION_COLUMNS = [
    "library_id",
    "minimum_total_count",
    "retained_gene_pairs",
    "retained_percent",
]


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


def positive_unit_float(value):
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise argparse.ArgumentTypeError(
            "must be greater than 0 and at most 1"
        ) from exc

    if math.isfinite(parsed) and 0 < parsed <= 1:
        return parsed

    raise argparse.ArgumentTypeError(
        "must be greater than 0 and at most 1"
    )


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
    gene_pairs_with_any_count = int(counted.sum())
    gene_pairs_with_parent1_count_only = int(
        ((table["parent1_count"] > 0) & (table["parent2_count"] == 0)).sum()
    )
    gene_pairs_with_parent2_count_only = int(
        ((table["parent1_count"] == 0) & (table["parent2_count"] > 0)).sum()
    )
    gene_pairs_with_both_parent_counts_nonzero = int(
        ((table["parent1_count"] > 0) & (table["parent2_count"] > 0)).sum()
    )
    total_parental_count = int(table["total_count"].sum())

    if counted.any():
        overall_parent1_proportion = (
            table["parent1_count"].sum() / total_parental_count
        )
        median_gene_pair_parent1_proportion = (
            table.loc[counted, "parent1_count"]
            / table.loc[counted, "total_count"]
        ).median()
    else:
        overall_parent1_proportion = pd.NA
        median_gene_pair_parent1_proportion = pd.NA

    retained_counts = {
        f"retained_ge_{threshold}": int(
            (table["total_count"] >= threshold).sum()
        )
        for threshold in RNA_RETENTION_THRESHOLDS
    }

    return {
        "sample_id": library_id,
        "total_gene_pairs": len(table),
        "gene_pairs_with_any_count": gene_pairs_with_any_count,
        "gene_pairs_with_parent1_count_only": (
            gene_pairs_with_parent1_count_only
        ),
        "gene_pairs_with_parent2_count_only": (
            gene_pairs_with_parent2_count_only
        ),
        "gene_pairs_with_both_parent_counts_nonzero": (
            gene_pairs_with_both_parent_counts_nonzero
        ),
        "total_parental_count": total_parental_count,
        "overall_parent1_proportion": fraction_value(
            overall_parent1_proportion
        ),
        "median_gene_pair_parent1_proportion": fraction_value(
            median_gene_pair_parent1_proportion
        ),
        **retained_counts,
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
    return summarize_table("all_samples", pooled)


def build_retention_table(library_ids, tables, thresholds):
    rows = []
    for library_id, table in zip(library_ids, tables):
        total_gene_pairs = len(table)
        for threshold in thresholds:
            retained_gene_pairs = int((table["total_count"] >= threshold).sum())
            if total_gene_pairs == 0:
                retained_percent = pd.NA
            else:
                retained_percent = 100 * retained_gene_pairs / total_gene_pairs
            rows.append(
                {
                    "library_id": library_id,
                    "minimum_total_count": threshold,
                    "retained_gene_pairs": retained_gene_pairs,
                    "retained_percent": retained_percent,
                }
            )
    return pd.DataFrame(rows, columns=RETENTION_COLUMNS)


def write_summary(rows, path):
    summary = pd.DataFrame(rows, columns=SUMMARY_COLUMNS)
    path.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(path, sep="\t", index=False, lineterminator="\n")


def write_retention(retention, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    retention.to_csv(
        path,
        sep="\t",
        index=False,
        lineterminator="\n",
        na_rep="NA",
    )


def main():
    ap = ErrorParser(description="Summarize RNA gene count QC metrics.")
    ap.add_argument("--count-tables", nargs="+", required=True)
    ap.add_argument("--summary", required=True)
    ap.add_argument("--retention", required=True)
    ap.add_argument("--report", required=True)
    ap.add_argument(
        "--min-total-count",
        type=positive_integer,
        default=1,
        help="Minimum total count used by threshold-dependent plots; Panel A always includes pooled gene pairs with total_count > 0.",
    )
    ap.add_argument(
        "--histogram-colors",
        nargs=2,
        type=plot_color,
        default=("#44C3D0", "#CD85B9"),
        metavar=("LEFT_COLOR", "RIGHT_COLOR"),
        help=(
            "Colors used for parental count proportions below 0.5 "
            "(LEFT_COLOR) and at or above 0.5 (RIGHT_COLOR)."
        ),
    )
    ap.add_argument(
        "--histogram-bins",
        type=positive_integer,
        default=20,
        help="Number of equal-width bins used for the Panel A parental count proportion histogram.",
    )
    ap.add_argument(
        "--depth-curve-color",
        type=plot_color,
        default="tab:blue",
        help="Color used for per-library cumulative gene-pair depth curves.",
    )
    ap.add_argument(
        "--depth-curve-alpha",
        type=positive_unit_float,
        default=0.15,
        help="Transparency used for per-library cumulative gene-pair depth curves; must be greater than 0 and at most 1.",
    )
    ap.add_argument(
        "--retention-colors",
        nargs=6,
        type=plot_color,
        default=(
            "tab:blue",
            "tab:orange",
            "tab:green",
            "tab:red",
            "tab:purple",
            "tab:brown",
        ),
        metavar=tuple(
            f"COLOR_{threshold}"
            for threshold in RNA_RETENTION_THRESHOLDS
        ),
        help=(
            "Colors used for retention thresholds "
            f"{', '.join(map(str, RNA_RETENTION_THRESHOLDS))}."
        ),
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
    retention = build_retention_table(
        library_ids,
        tables,
        RNA_RETENTION_THRESHOLDS,
    )
    write_retention(retention, Path(args.retention))
    write_rna_gene_qc_report(
        library_ids,
        tables,
        pooled,
        Path(args.report),
        args.min_total_count,
        args.histogram_colors,
        args.sample_histogram_rows,
        args.sample_histogram_columns,
        parent1_label=args.parent1_label,
        parent2_label=args.parent2_label,
        histogram_bins=args.histogram_bins,
        depth_curve_color=args.depth_curve_color,
        depth_curve_alpha=args.depth_curve_alpha,
        retention=retention,
        retention_colors=args.retention_colors,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
