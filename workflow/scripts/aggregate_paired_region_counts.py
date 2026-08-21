#!/usr/bin/env python3

import argparse
import sys
from pathlib import Path

import pandas as pd


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
MAX_REPORTED_ERRORS = 50


class ErrorParser(argparse.ArgumentParser):
    def error(self, message):
        self.exit(2, f"ERROR: {message}\n")


def fail(errors):
    for error in errors[:MAX_REPORTED_ERRORS]:
        print(f"ERROR: {error}", file=sys.stderr)
    omitted = len(errors) - MAX_REPORTED_ERRORS
    if omitted > 0:
        print(
            f"ERROR: {omitted} additional validation errors omitted",
            file=sys.stderr,
        )
    return 1


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
                dtype=object,
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


def build_matrix(library_ids, tables):
    matrix = tables[0][REGION_COLUMNS].copy()
    for library_id, table in zip(library_ids, tables):
        matrix[f"{library_id}_parent1_count"] = table[
            "parent1_count"
        ].to_numpy()
        matrix[f"{library_id}_parent2_count"] = table[
            "parent2_count"
        ].to_numpy()
    return matrix


def write_matrix(matrix, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    matrix.to_csv(path, sep="\t", index=False, lineterminator="\n")


def main():
    parser = ErrorParser(
        description="Build a study-wide paired-region count matrix."
    )
    parser.add_argument("--count-tables", nargs="+", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    paths = [Path(path) for path in args.count_tables]
    library_ids = [path.parent.name for path in paths]
    errors = []

    for library_id, path in zip(library_ids, paths):
        if not library_id:
            errors.append(
                f"{path}: cannot derive library ID from parent directory"
            )

    duplicate_ids = sorted(
        {
            library_id
            for library_id in library_ids
            if library_ids.count(library_id) > 1
        }
    )
    if duplicate_ids:
        errors.append(
            "derived library IDs are not unique: "
            f"{', '.join(duplicate_ids)}"
        )
    if errors:
        return fail(errors)

    sorted_inputs = sorted(
        zip(library_ids, paths),
        key=lambda item: item[0],
    )
    library_ids = [library_id for library_id, _ in sorted_inputs]
    paths = [path for _, path in sorted_inputs]

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

    matrix = build_matrix(library_ids, tables)
    output_path = Path(args.output)
    try:
        write_matrix(matrix, output_path)
    except OSError as exc:
        return fail([f"{output_path}: cannot write ({exc})"])

    return 0


if __name__ == "__main__":
    sys.exit(main())
