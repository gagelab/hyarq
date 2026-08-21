#!/usr/bin/env python3
"""Count a filtered BAM over paired parental regions with featureCounts."""

import argparse
import csv
import re
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime
from pathlib import Path


PAIRED_REGION_HEADER = [
    "parent1_chrom",
    "parent1_start",
    "parent1_end",
    "parent2_chrom",
    "parent2_start",
    "parent2_end",
]
COUNT_HEADER = PAIRED_REGION_HEADER + [
    "parent1_count",
    "parent2_count",
    "total_count",
]
SUMMARY_CATEGORIES = (
    "Assigned",
    "Unassigned_NoFeatures",
    "Unassigned_Unmapped",
    "Unassigned_MultiMapping",
)
PAIREDNESS_SAMPLE_SIZE = 4096


class CountingError(Exception):
    """An expected user input or external-tool failure."""


def log(message):
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] {message}", file=sys.stderr)


def positive_int(value):
    try:
        parsed = int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("must be a positive integer") from exc
    if parsed < 1:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return parsed


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Count a filtered single-end BAM against a six-column paired-region "
            "table using featureCounts."
        )
    )
    parser.add_argument("--bam", required=True)
    parser.add_argument("--peak-pairs", required=True)
    parser.add_argument("--parent1-prefix", required=True)
    parser.add_argument("--parent2-prefix", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--tmpdir", required=True)
    parser.add_argument("--assay", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--rep", required=True)
    parser.add_argument("--library-id", required=True)
    parser.add_argument("--threads", required=True, type=positive_int)
    parser.add_argument("--featurecounts", default="featureCounts")
    parser.add_argument("--samtools", default="samtools")
    return parser.parse_args()


def require_nonempty(value, argument):
    if not value:
        raise CountingError(f"{argument} must not be empty")


def require_command(command):
    if shutil.which(command) is None:
        raise CountingError(f"required command not found in PATH: {command}")


def validate_inputs(args):
    for value, argument in (
        (args.parent1_prefix, "--parent1-prefix"),
        (args.parent2_prefix, "--parent2-prefix"),
        (args.assay, "--assay"),
        (args.sample, "--sample"),
        (args.rep, "--rep"),
        (args.library_id, "--library-id"),
    ):
        require_nonempty(value, argument)

    bam = Path(args.bam)
    peak_pairs = Path(args.peak_pairs)
    if not bam.is_file():
        raise CountingError(f"BAM input not found: {bam}")
    if not peak_pairs.is_file():
        raise CountingError(f"paired-region table not found: {peak_pairs}")

    output = Path(args.output).resolve()
    summary = Path(args.summary).resolve()
    if output == summary:
        raise CountingError("--output and --summary must name different files")

    require_command(args.featurecounts)
    require_command(args.samtools)


def read_paired_regions(path):
    rows = []
    try:
        source = Path(path).open(newline="")
    except OSError as exc:
        raise CountingError(f"cannot open paired-region table {path}: {exc}") from exc

    with source:
        reader = csv.reader(source, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration as exc:
            raise CountingError(f"paired-region table is empty: {path}") from exc
        if header != PAIRED_REGION_HEADER:
            expected = "\t".join(PAIRED_REGION_HEADER)
            raise CountingError(
                f"paired-region header must be exactly: {expected}"
            )

        for line_number, fields in enumerate(reader, start=2):
            if not fields:
                continue
            if len(fields) != len(PAIRED_REGION_HEADER):
                raise CountingError(
                    f"{path}:{line_number}: expected six columns, found {len(fields)}"
                )

            row = dict(zip(PAIRED_REGION_HEADER, fields))
            for parent in ("parent1", "parent2"):
                chrom = row[f"{parent}_chrom"]
                if not chrom:
                    raise CountingError(
                        f"{path}:{line_number}: blank {parent}_chrom"
                    )
                try:
                    start = int(row[f"{parent}_start"])
                    end = int(row[f"{parent}_end"])
                except ValueError as exc:
                    raise CountingError(
                        f"{path}:{line_number}: {parent} coordinates must be integers"
                    ) from exc
                if start < 0 or end <= start:
                    raise CountingError(
                        f"{path}:{line_number}: invalid {parent} interval "
                        f"{start}-{end}; expected 0 <= start < end"
                    )
            rows.append(row)

    if not rows:
        raise CountingError(f"paired-region table has no data rows: {path}")
    return rows


def write_saf(rows, parent1_prefix, parent2_prefix, path):
    feature_ids = []
    used_contigs = set()
    with path.open("w", newline="") as output:
        writer = csv.writer(output, delimiter="\t", lineterminator="\n")
        writer.writerow(("GeneID", "Chr", "Start", "End", "Strand"))
        for row_number, row in enumerate(rows, start=1):
            pair_features = {}
            for parent, prefix in (
                ("parent1", parent1_prefix),
                ("parent2", parent2_prefix),
            ):
                feature_id = f"pair_{row_number:06d}_{parent}"
                contig = f"{prefix}_{row[f'{parent}_chrom']}"
                start = int(row[f"{parent}_start"]) + 1
                end = int(row[f"{parent}_end"])
                writer.writerow((feature_id, contig, start, end, "+"))
                pair_features[parent] = feature_id
                feature_ids.append(feature_id)
                used_contigs.add(contig)
            row["_feature_ids"] = pair_features
    return feature_ids, used_contigs


def bam_contigs(samtools, bam):
    try:
        result = subprocess.run(
            [samtools, "view", "-H", str(bam)],
            check=True,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
    except subprocess.CalledProcessError as exc:
        detail = exc.stderr.strip() or "no diagnostic output"
        raise CountingError(
            f"samtools view -H failed for {bam}: {detail}"
        ) from exc

    contigs = set()
    for line in result.stdout.splitlines():
        if not line.startswith("@SQ\t"):
            continue
        for field in line.split("\t")[1:]:
            if field.startswith("SN:"):
                contigs.add(field[3:])
                break
    return contigs


def check_used_contigs(used_contigs, header_contigs):
    missing = sorted(used_contigs - header_contigs)
    if missing:
        raise CountingError(
            "translated paired-region contigs absent from BAM header: "
            + ", ".join(missing)
        )


def check_single_end_bam(samtools, bam):
    process = subprocess.Popen(
        [samtools, "view", "-F", "0x4", str(bam)],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    sampled_records = 0
    paired_alignment_found = False
    terminated_intentionally = False
    stderr = ""

    try:
        while sampled_records < PAIREDNESS_SAMPLE_SIZE:
            line = process.stdout.readline()
            if not line:
                break
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                raise CountingError(
                    "samtools produced a malformed SAM record while checking "
                    "for paired alignments"
                )
            try:
                flag = int(fields[1])
            except ValueError as exc:
                raise CountingError(
                    "samtools produced a non-integer SAM flag while checking "
                    "for paired alignments"
                ) from exc
            sampled_records += 1
            if flag & 0x1:
                paired_alignment_found = True
                break

        if paired_alignment_found or sampled_records >= PAIREDNESS_SAMPLE_SIZE:
            if process.poll() is None:
                process.terminate()
                terminated_intentionally = True
        _, stderr = process.communicate()
    finally:
        if process.poll() is None:
            process.terminate()
            process.wait()
        process.stdout.close()
        process.stderr.close()

    if paired_alignment_found:
        raise CountingError(
            "count_paired_regions.py currently supports single-end BAMs only; "
            "paired-end fragment counting for assays such as ATAC/ChIP will be "
            "implemented separately"
        )
    if not terminated_intentionally and process.returncode != 0:
        detail = stderr.strip() or "no diagnostic output"
        raise CountingError(
            f"samtools mapped-alignment inspection failed for {bam}: {detail}"
        )


def get_featurecounts_version(featurecounts):
    try:
        result = subprocess.run(
            [featurecounts, "-v"],
            check=False,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
    except (OSError, subprocess.SubprocessError):
        return "NA"

    version_output = "\n".join((result.stdout, result.stderr))
    branded_match = re.search(
        r"\b(?:featureCounts|Subread)\b[^\n0-9]*"
        r"([0-9]+(?:\.[0-9]+)+(?:[-+._A-Za-z0-9]*))",
        version_output,
        flags=re.IGNORECASE,
    )
    if branded_match:
        return branded_match.group(1)
    version_match = re.search(
        r"\bv([0-9]+(?:\.[0-9]+)+(?:[-+._A-Za-z0-9]*))",
        version_output,
        flags=re.IGNORECASE,
    )
    return version_match.group(1) if version_match else "NA"


def run_featurecounts(args, saf, raw_output, workdir):
    featurecounts_log = workdir / "featureCounts.log"
    command = [
        args.featurecounts,
        "-F",
        "SAF",
        "-s",
        "0",
        "-O",
        "--minOverlap",
        "1",
        "-T",
        str(args.threads),
        "--tmpDir",
        str(workdir),
        "-a",
        str(saf),
        "-o",
        str(raw_output),
        str(args.bam),
    ]
    with featurecounts_log.open("w") as log_output:
        try:
            subprocess.run(
                command,
                check=True,
                stdout=log_output,
                stderr=subprocess.STDOUT,
            )
        except subprocess.CalledProcessError as exc:
            raise CountingError(
                f"featureCounts failed; inspect {featurecounts_log}"
            ) from exc

    if not raw_output.is_file() or raw_output.stat().st_size == 0:
        raise CountingError(
            f"featureCounts produced no count output; inspect {featurecounts_log}"
        )
    summary_path = Path(f"{raw_output}.summary")
    if not summary_path.is_file() or summary_path.stat().st_size == 0:
        raise CountingError(
            f"featureCounts produced no summary output; inspect {featurecounts_log}"
        )
    return summary_path


def parse_featurecounts(path, expected_feature_ids):
    expected = set(expected_feature_ids)
    counts = {}
    header = None

    with path.open() as source:
        for line_number, line in enumerate(source, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if header is None:
                if fields[0] != "Geneid":
                    raise CountingError(
                        f"{path}:{line_number}: expected featureCounts Geneid header"
                    )
                if len(fields) != 7:
                    raise CountingError(
                        f"{path}:{line_number}: expected one featureCounts sample column"
                    )
                header = fields
                continue

            if len(fields) != len(header):
                raise CountingError(
                    f"{path}:{line_number}: malformed featureCounts row"
                )
            feature_id = fields[0]
            if feature_id not in expected:
                raise CountingError(
                    f"unexpected feature ID in featureCounts output: {feature_id}"
                )
            if feature_id in counts:
                raise CountingError(
                    f"duplicate feature ID in featureCounts output: {feature_id}"
                )
            try:
                counts[feature_id] = int(fields[-1])
            except ValueError as exc:
                raise CountingError(
                    f"{path}:{line_number}: non-integer featureCounts value"
                ) from exc

    if header is None:
        raise CountingError(f"featureCounts Geneid header not found in {path}")
    missing = [feature_id for feature_id in expected_feature_ids if feature_id not in counts]
    if missing:
        raise CountingError(
            "feature IDs missing from featureCounts output: " + ", ".join(missing)
        )
    return counts


def parse_featurecounts_summary(path):
    selected = {}
    with path.open(newline="") as source:
        reader = csv.reader(source, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration as exc:
            raise CountingError(f"featureCounts summary is empty: {path}") from exc
        if not header or header[0] != "Status" or len(header) != 2:
            raise CountingError(
                f"unexpected featureCounts summary header in {path}"
            )
        for line_number, fields in enumerate(reader, start=2):
            if not fields:
                continue
            if len(fields) != len(header):
                raise CountingError(
                    f"{path}:{line_number}: malformed featureCounts summary row"
                )
            category = fields[0]
            if category not in SUMMARY_CATEGORIES:
                continue
            if category in selected:
                raise CountingError(
                    f"duplicate featureCounts summary category: {category}"
                )
            try:
                selected[category] = int(fields[1])
            except ValueError as exc:
                raise CountingError(
                    f"{path}:{line_number}: non-integer summary value for {category}"
                ) from exc
    return selected


def write_counts(rows, counts, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as output:
        writer = csv.writer(output, delimiter="\t", lineterminator="\n")
        writer.writerow(COUNT_HEADER)
        for row in rows:
            parent1_count = counts[row["_feature_ids"]["parent1"]]
            parent2_count = counts[row["_feature_ids"]["parent2"]]
            writer.writerow(
                [row[column] for column in PAIRED_REGION_HEADER]
                + [parent1_count, parent2_count, parent1_count + parent2_count]
            )


def write_summary(
    args,
    pair_rows,
    featurecounts_version,
    featurecounts_summary,
    path,
):
    metrics = [
        ("assay", args.assay),
        ("sample", args.sample),
        ("rep", args.rep),
        ("library_id", args.library_id),
        ("bam", args.bam),
        ("peak_pairs", args.peak_pairs),
        ("parent1_prefix", args.parent1_prefix),
        ("parent2_prefix", args.parent2_prefix),
        ("featurecounts_version", featurecounts_version),
        ("overlap_min_bp", 1),
        ("multi_feature_assignment", "all_overlaps"),
        ("peak_pair_rows", pair_rows),
        ("peak_features", pair_rows * 2),
    ]
    metrics.extend(
        (category, featurecounts_summary[category])
        for category in SUMMARY_CATEGORIES
        if category in featurecounts_summary
    )

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as output:
        writer = csv.writer(output, delimiter="\t", lineterminator="\n")
        writer.writerow(("metric", "value"))
        writer.writerows(metrics)


def run(args, workdir):
    rows = read_paired_regions(args.peak_pairs)
    saf = workdir / "paired_regions.saf"
    raw_output = workdir / "featureCounts.txt"

    feature_ids, used_contigs = write_saf(
        rows,
        args.parent1_prefix,
        args.parent2_prefix,
        saf,
    )

    log("Inspecting BAM contigs")
    check_used_contigs(
        used_contigs,
        bam_contigs(args.samtools, Path(args.bam)),
    )

    log("Checking sampled mapped alignments for paired-end flags")
    check_single_end_bam(args.samtools, Path(args.bam))

    featurecounts_version = get_featurecounts_version(args.featurecounts)
    log(f"Counting {len(feature_ids)} features across {len(rows)} paired rows")
    featurecounts_summary_path = run_featurecounts(
        args,
        saf,
        raw_output,
        workdir,
    )
    counts = parse_featurecounts(raw_output, feature_ids)
    featurecounts_summary = parse_featurecounts_summary(
        featurecounts_summary_path
    )

    write_counts(rows, counts, Path(args.output))
    write_summary(
        args,
        len(rows),
        featurecounts_version,
        featurecounts_summary,
        Path(args.summary),
    )


def main():
    args = parse_args()
    workdir = None
    try:
        validate_inputs(args)
        tmpdir = Path(args.tmpdir)
        tmpdir.mkdir(parents=True, exist_ok=True)
        if not tmpdir.is_dir():
            raise CountingError(f"temporary root is not a directory: {tmpdir}")
        workdir = Path(
            tempfile.mkdtemp(prefix="count_paired_regions.", dir=tmpdir)
        )
        log(f"Working directory: {workdir}")
        run(args, workdir)
    except (CountingError, OSError, csv.Error) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        if workdir is not None:
            print(
                f"ERROR: preserving working directory for debugging: {workdir}",
                file=sys.stderr,
            )
        return 1

    try:
        shutil.rmtree(workdir)
    except OSError as exc:
        print(
            f"ERROR: could not remove working directory {workdir}: {exc}",
            file=sys.stderr,
        )
        return 1

    log(f"Paired-region counting complete: {args.output}")
    log(f"Summary: {args.summary}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
