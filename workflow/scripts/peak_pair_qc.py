#!/usr/bin/env python3

import argparse
import math
import sys
from pathlib import Path

import pandas as pd

from paired_feature_qc_plots import (
    PAIRED_FEATURE_RETENTION_THRESHOLDS,
    FeatureTerminology,
    plot_color,
    write_paired_feature_qc_report,
    write_paired_feature_qc_sample_report,
)


PEAK_PAIR_TERMINOLOGY = FeatureTerminology(
    singular="peak pair",
    plural="peak pairs",
    plural_capitalized="Peak pairs",
    hyphenated="peak-pair",
)
REGION_COLUMNS = [
    "parent1_chrom",
    "parent1_start",
    "parent1_end",
    "parent2_chrom",
    "parent2_start",
    "parent2_end",
]
COUNT_COLUMNS = [
    "parent1_count",
    "parent2_count",
    "total_count",
]
EXPECTED_COLUMNS = REGION_COLUMNS + COUNT_COLUMNS
SUMMARY_COLUMNS = [
    "sample_id",
    "total_peak_pairs",
    "peak_pairs_with_any_count",
    "peak_pairs_with_parent1_count_only",
    "peak_pairs_with_parent2_count_only",
    "peak_pairs_with_both_parent_counts_nonzero",
    "total_parental_count",
    "overall_parent1_proportion",
    "median_peak_pair_parent1_proportion",
    *[
        f"retained_ge_{threshold}"
        for threshold in PAIRED_FEATURE_RETENTION_THRESHOLDS
    ],
]
RETENTION_COLUMNS = [
    "library_id",
    "minimum_total_count",
    "retained_peak_pairs",
    "retained_percent",
]
ASSAY_REPORT_PREFIXES = {
    "moa": "MOA peak",
    "atac": "ATAC peak",
    "chip": "ChIP peak",
}


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
        table = pd.read_csv(
            path,
            sep="\t",
            dtype=str,
            keep_default_na=False,
        )
    except OSError as exc:
        errors.append(f"{path}: cannot open ({exc})")
        return None
    except UnicodeError as exc:
        errors.append(f"{path}: cannot decode ({exc})")
        return None
    except pd.errors.EmptyDataError:
        errors.append(f"{path}: empty file")
        return None
    except pd.errors.ParserError as exc:
        errors.append(f"{path}: cannot parse TSV ({exc})")
        return None

    if list(table.columns) != EXPECTED_COLUMNS:
        errors.append(
            f"{path}: header must be exactly {EXPECTED_COLUMNS}, "
            f"got {list(table.columns)}"
        )
        return None
    if not isinstance(table.index, pd.RangeIndex):
        errors.append(f"{path}: malformed TSV rows contain unexpected fields")
        return None
    if table.empty:
        errors.append(f"{path}: no data rows")
        return None

    return table


def validate_counts(table, path, errors):
    parsed_columns = {}
    valid_masks = {}

    for column in COUNT_COLUMNS:
        values = table[column]
        valid = values.str.fullmatch(r"[0-9]+")
        valid_masks[column] = valid
        for row_index in table.index[~valid]:
            errors.append(
                f"{path}:{row_index + 2}: {column} must be a "
                "nonnegative integer"
            )
        parsed_columns[column] = [
            int(value) if is_valid else None
            for value, is_valid in zip(values, valid)
        ]

    for position in range(len(table)):
        if not all(
            valid_masks[column].iloc[position]
            for column in COUNT_COLUMNS
        ):
            continue
        parent1_count = parsed_columns["parent1_count"][position]
        parent2_count = parsed_columns["parent2_count"][position]
        total_count = parsed_columns["total_count"][position]
        if total_count != parent1_count + parent2_count:
            errors.append(
                f"{path}:{position + 2}: total_count does not equal "
                "parent1_count + parent2_count"
            )

    if all(valid_masks[column].all() for column in COUNT_COLUMNS):
        for column in COUNT_COLUMNS:
            table[column] = pd.Series(
                parsed_columns[column],
                index=table.index,
                dtype="int64",
            )


def validate_reference_regions(table, path, errors):
    for row_index, row in table.iterrows():
        for parent in ("parent1", "parent2"):
            chrom_column = f"{parent}_chrom"
            start_column = f"{parent}_start"
            end_column = f"{parent}_end"

            if not row[chrom_column].strip():
                errors.append(
                    f"{path}:{row_index + 2}: blank {chrom_column}"
                )

            parsed = {}
            for column in (start_column, end_column):
                try:
                    parsed[column] = int(row[column])
                except ValueError:
                    errors.append(
                        f"{path}:{row_index + 2}: {column} must be an integer"
                    )

            if start_column not in parsed or end_column not in parsed:
                continue
            start = parsed[start_column]
            end = parsed[end_column]
            if not (0 <= start < end):
                errors.append(
                    f"{path}:{row_index + 2}: invalid {parent} interval "
                    f"{start}-{end}; expected 0 <= start < end"
                )


def compare_region_layout(
    reference_table,
    observed_table,
    library_id,
    path,
    errors,
):
    reference_rows = reference_table[REGION_COLUMNS]
    observed_rows = observed_table[REGION_COLUMNS]
    shared_row_count = min(len(reference_rows), len(observed_rows))

    for row_index in range(shared_row_count):
        reference_region = tuple(reference_rows.iloc[row_index].tolist())
        observed_region = tuple(observed_rows.iloc[row_index].tolist())
        if observed_region != reference_region:
            errors.append(
                f"library {library_id} ({path}): region layout differs at "
                f"TSV line {row_index + 2}; reference={reference_region}; "
                f"observed={observed_region}"
            )
            return

    if len(observed_rows) != len(reference_rows):
        errors.append(
            f"library {library_id} ({path}): region row count differs; "
            f"reference={len(reference_rows)}, observed={len(observed_rows)}"
        )


def fraction_value(value):
    return "NA" if pd.isna(value) else f"{value:.6f}"


def summarize_table(library_id, table):
    counted = table["total_count"] > 0
    peak_pairs_with_any_count = int(counted.sum())
    peak_pairs_with_parent1_count_only = int(
        ((table["parent1_count"] > 0) & (table["parent2_count"] == 0)).sum()
    )
    peak_pairs_with_parent2_count_only = int(
        ((table["parent1_count"] == 0) & (table["parent2_count"] > 0)).sum()
    )
    peak_pairs_with_both_parent_counts_nonzero = int(
        ((table["parent1_count"] > 0) & (table["parent2_count"] > 0)).sum()
    )
    total_parental_count = int(table["total_count"].sum())

    if counted.any():
        overall_parent1_proportion = (
            table["parent1_count"].sum() / total_parental_count
        )
        median_peak_pair_parent1_proportion = (
            table.loc[counted, "parent1_count"]
            / table.loc[counted, "total_count"]
        ).median()
    else:
        overall_parent1_proportion = pd.NA
        median_peak_pair_parent1_proportion = pd.NA

    retained_counts = {
        f"retained_ge_{threshold}": int(
            (table["total_count"] >= threshold).sum()
        )
        for threshold in PAIRED_FEATURE_RETENTION_THRESHOLDS
    }

    return {
        "sample_id": library_id,
        "total_peak_pairs": len(table),
        "peak_pairs_with_any_count": peak_pairs_with_any_count,
        "peak_pairs_with_parent1_count_only": (
            peak_pairs_with_parent1_count_only
        ),
        "peak_pairs_with_parent2_count_only": (
            peak_pairs_with_parent2_count_only
        ),
        "peak_pairs_with_both_parent_counts_nonzero": (
            peak_pairs_with_both_parent_counts_nonzero
        ),
        "total_parental_count": total_parental_count,
        "overall_parent1_proportion": fraction_value(
            overall_parent1_proportion
        ),
        "median_peak_pair_parent1_proportion": fraction_value(
            median_peak_pair_parent1_proportion
        ),
        **retained_counts,
    }


def build_pooled_table(tables):
    pooled = tables[0][REGION_COLUMNS].copy()
    pooled["parent1_count"] = sum(
        table["parent1_count"].to_numpy()
        for table in tables
    )
    pooled["parent2_count"] = sum(
        table["parent2_count"].to_numpy()
        for table in tables
    )
    pooled["total_count"] = pooled["parent1_count"] + pooled["parent2_count"]
    return pooled


def summarize_pooled(pooled):
    return summarize_table("all_samples", pooled)


def build_retention_table(library_ids, tables, thresholds):
    rows = []
    for library_id, table in zip(library_ids, tables):
        total_peak_pairs = len(table)
        for threshold in thresholds:
            retained_peak_pairs = int((table["total_count"] >= threshold).sum())
            if total_peak_pairs == 0:
                retained_percent = pd.NA
            else:
                retained_percent = 100 * retained_peak_pairs / total_peak_pairs
            rows.append(
                {
                    "library_id": library_id,
                    "minimum_total_count": threshold,
                    "retained_peak_pairs": retained_peak_pairs,
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


def resolve_sample_libraries(
    requested_library_ids,
    library_ids,
    tables,
):
    if "all" in library_ids:
        raise ValueError(
            "library ID 'all' is reserved for --sample-libraries"
        )
    if requested_library_ids == ["all"]:
        return list(library_ids), list(tables)
    if "all" in requested_library_ids:
        raise ValueError(
            "--sample-libraries cannot combine 'all' "
            "with specific library IDs"
        )

    seen = set()
    duplicate_library_ids = []
    for library_id in requested_library_ids:
        if (
            library_id in seen
            and library_id not in duplicate_library_ids
        ):
            duplicate_library_ids.append(library_id)
        seen.add(library_id)
    if duplicate_library_ids:
        raise ValueError(
            "--sample-libraries contains duplicate library IDs: "
            f"{', '.join(duplicate_library_ids)}"
        )

    unknown_library_ids = [
        library_id
        for library_id in requested_library_ids
        if library_id not in library_ids
    ]
    if unknown_library_ids:
        raise ValueError(
            "--sample-libraries contains unknown library IDs: "
            f"{', '.join(unknown_library_ids)}"
        )

    table_by_library_id = dict(zip(library_ids, tables))
    selected_tables = [
        table_by_library_id[library_id]
        for library_id in requested_library_ids
    ]
    return list(requested_library_ids), selected_tables


def main():
    ap = ErrorParser(description="Summarize paired-peak count QC metrics.")
    ap.add_argument("--count-tables", nargs="+", required=True)
    ap.add_argument(
        "--assay",
        choices=("moa", "atac", "chip"),
        required=True,
    )
    ap.add_argument("--summary", required=True)
    ap.add_argument("--retention", required=True)
    ap.add_argument("--report", required=True)
    ap.add_argument("--sample-report", default=None)
    ap.add_argument(
        "--sample-libraries",
        nargs="+",
        default=None,
    )
    ap.add_argument(
        "--sample-plot-types",
        nargs="+",
        choices=("histogram", "scatter"),
        default=None,
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
        help="Color used for per-library cumulative peak-pair depth curves.",
    )
    ap.add_argument(
        "--depth-curve-alpha",
        type=positive_unit_float,
        default=0.15,
        help="Transparency used for per-library cumulative peak-pair depth curves; must be greater than 0 and at most 1.",
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
            for threshold in PAIRED_FEATURE_RETENTION_THRESHOLDS
        ),
        help=(
            "Colors used for retention thresholds "
            f"{', '.join(map(str, PAIRED_FEATURE_RETENTION_THRESHOLDS))}."
        ),
    )
    ap.add_argument("--parent1-label", default="Parent1", help="Parent1 name displayed in the report.")
    ap.add_argument("--parent2-label", default="Parent2", help="Parent2 name displayed in the report.")
    args = ap.parse_args()
    report_prefix = ASSAY_REPORT_PREFIXES[args.assay]

    sample_option_values = (
        args.sample_report,
        args.sample_libraries,
        args.sample_plot_types,
    )
    if any(value is not None for value in sample_option_values) and not all(
        value is not None for value in sample_option_values
    ):
        ap.error(
            "--sample-report, --sample-libraries, and "
            "--sample-plot-types must be supplied together"
        )
    if args.sample_report is not None:
        if Path(args.sample_report).resolve() == Path(args.report).resolve():
            ap.error("--sample-report must differ from --report")
        duplicate_plot_types = [
            plot_type
            for index, plot_type in enumerate(
                args.sample_plot_types
            )
            if plot_type in args.sample_plot_types[:index]
        ]
        if duplicate_plot_types:
            ap.error(
                "--sample-plot-types contains duplicate values: "
                f"{', '.join(duplicate_plot_types)}"
            )
        args.sample_plot_types = [
            plot_type
            for plot_type in ("histogram", "scatter")
            if plot_type in args.sample_plot_types
        ]

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
        table = read_count_table(path, errors)
        if table is not None:
            validate_counts(table, path, errors)
        tables.append(table)

    reference_table = tables[0]
    if reference_table is not None:
        validate_reference_regions(reference_table, paths[0], errors)
    for library_id, path, table in zip(
        library_ids[1:],
        paths[1:],
        tables[1:],
    ):
        if reference_table is None or table is None:
            continue
        compare_region_layout(
            reference_table,
            table,
            library_id,
            path,
            errors,
        )
    if errors:
        return fail(errors)

    selected_sample_library_ids = None
    selected_sample_tables = None
    if args.sample_report is not None:
        try:
            (
                selected_sample_library_ids,
                selected_sample_tables,
            ) = resolve_sample_libraries(
                args.sample_libraries,
                library_ids,
                tables,
            )
        except ValueError as exc:
            ap.error(str(exc))

    pooled = build_pooled_table(tables)
    rows = [summarize_table(library_id, table) for library_id, table in zip(library_ids, tables)]
    rows.append(summarize_pooled(pooled))
    write_summary(rows, Path(args.summary))
    retention = build_retention_table(
        library_ids,
        tables,
        PAIRED_FEATURE_RETENTION_THRESHOLDS,
    )
    write_retention(retention, Path(args.retention))
    write_paired_feature_qc_report(
        tables,
        pooled,
        Path(args.report),
        args.histogram_colors,
        parent1_label=args.parent1_label,
        parent2_label=args.parent2_label,
        histogram_bins=args.histogram_bins,
        depth_curve_color=args.depth_curve_color,
        depth_curve_alpha=args.depth_curve_alpha,
        retention=retention,
        retention_colors=args.retention_colors,
        terminology=PEAK_PAIR_TERMINOLOGY,
        report_prefix=report_prefix,
    )
    if args.sample_report is not None:
        write_paired_feature_qc_sample_report(
            selected_sample_library_ids,
            selected_sample_tables,
            Path(args.sample_report),
            tuple(args.sample_plot_types),
            args.histogram_colors,
            args.parent1_label,
            args.parent2_label,
            terminology=PEAK_PAIR_TERMINOLOGY,
            report_prefix=report_prefix,
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
