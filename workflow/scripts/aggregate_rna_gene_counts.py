#!/usr/bin/env python3

import argparse
import sys
from pathlib import Path

import pandas as pd


EXPECTED_COLUMNS = [
    "parent1_gene_id",
    "parent2_gene_id",
    "parent1_count",
    "parent2_count",
    "total_count",
]
PAIR_COLUMNS = ["parent1_gene_id", "parent2_gene_id"]
COUNT_COLUMNS = ["parent1_count", "parent2_count", "total_count"]


class ErrorParser(argparse.ArgumentParser):
    def error(self, message):
        self.exit(2, f"ERROR: {message}\n")


def fail(errors):
    for error in errors:
        print(f"ERROR: {error}", file=sys.stderr)
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

    return table


def validate_gene_ids(table, path, errors):
    for column in PAIR_COLUMNS:
        blank_rows = table.index[table[column].str.strip().eq("")]
        for row_index in blank_rows:
            errors.append(f"{path}:{row_index + 2}: blank {column}")


def validate_count_column(table, column, path, errors):
    values = table[column]
    invalid = ~values.str.fullmatch(r"[0-9]+")
    for row_index in table.index[invalid]:
        errors.append(
            f"{path}:{row_index + 2}: {column} must be a "
            "nonnegative integer"
        )

    parsed = [
        int(value) if valid else 0
        for value, valid in zip(values, ~invalid)
    ]
    table[column] = pd.Series(parsed, index=table.index, dtype=object)
    return invalid


def validate_counts(table, path, errors):
    invalid_masks = [
        validate_count_column(table, column, path, errors)
        for column in COUNT_COLUMNS
    ]
    invalid_counts = invalid_masks[0] | invalid_masks[1] | invalid_masks[2]
    mismatches = (~invalid_counts) & (
        table["total_count"]
        != table["parent1_count"] + table["parent2_count"]
    )
    for row_index in table.index[mismatches]:
        errors.append(
            f"{path}:{row_index + 2}: total_count does not equal "
            "parent1_count + parent2_count"
        )


def validate_duplicate_pairs(table, path, errors):
    duplicates = table.duplicated(PAIR_COLUMNS, keep=False)
    duplicate_pairs = table.loc[duplicates, PAIR_COLUMNS].drop_duplicates()
    for pair in duplicate_pairs.itertuples(index=False):
        errors.append(
            f"{path}: duplicate gene pair: "
            f"{pair.parent1_gene_id}/{pair.parent2_gene_id}"
        )


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
    return set(
        map(tuple, table[PAIR_COLUMNS].itertuples(index=False, name=None))
    )


def build_matrix(library_ids, tables):
    matrix = tables[0][PAIR_COLUMNS].copy()
    reference_index = pd.MultiIndex.from_frame(matrix)

    for library_id, table in zip(library_ids, tables):
        indexed = table.set_index(PAIR_COLUMNS)
        aligned = indexed.reindex(reference_index)
        matrix[f"{library_id}_parent1_count"] = aligned[
            "parent1_count"
        ].to_numpy()
        matrix[f"{library_id}_parent2_count"] = aligned[
            "parent2_count"
        ].to_numpy()

    return matrix


def write_matrix(matrix, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    matrix.to_csv(path, sep="\t", index=False, lineterminator="\n")


def main():
    parser = ErrorParser(
        description="Build a study-wide RNA paired-gene count matrix."
    )
    parser.add_argument("--count-tables", nargs="+", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    paths = [Path(path) for path in args.count_tables]
    library_ids = [path.parent.name for path in paths]
    duplicate_ids = sorted(
        {
            library_id
            for library_id in library_ids
            if library_ids.count(library_id) > 1
        }
    )
    if duplicate_ids:
        return fail(
            [
                "derived library IDs are not unique: "
                f"{', '.join(duplicate_ids)}"
            ]
        )

    sorted_inputs = sorted(
        zip(library_ids, paths),
        key=lambda item: item[0],
    )
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
            errors.append(
                f"{path}: gene pair set does not match the first count table"
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
