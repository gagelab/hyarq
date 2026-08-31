from collections.abc import Mapping

import pandas as pd


SAMPLE_COLUMNS = ["sample_id", "assay", "replicate", "r1", "r2"]

RULE_RESOURCE_DEFAULTS = {
    "build_concatenated_reference": {
        "threads": 1,
        "mem_mb": 8000,
        "runtime": 120,
    },
    "build_concatenated_gtf": {
        "threads": 1,
        "mem_mb": 8000,
        "runtime": 120,
    },
    "build_star_index": {
        "threads": 8,
        "mem_mb": 96000,
        "runtime": 720,
    },
    "validate_rna_gene_pairs": {
        "threads": 1,
        "mem_mb": 2000,
        "runtime": 30,
    },
    "validate_moa_peak_pairs": {
        "threads": 1,
        "mem_mb": 2000,
        "runtime": 30,
    },
    "validate_atac_peak_pairs": {
        "threads": 1,
        "mem_mb": 2000,
        "runtime": 30,
    },
    "validate_chip_peak_pairs": {
        "threads": 1,
        "mem_mb": 2000,
        "runtime": 30,
    },
    "raw_fastqc": {
        "threads": 2,
        "mem_mb": 4000,
        "runtime": 120,
    },
    "rna_fastp": {
        "threads": 4,
        "mem_mb": 8000,
        "runtime": 240,
    },
    "atac_fastp": {
        "threads": 4,
        "mem_mb": 8000,
        "runtime": 240,
    },
    "chip_fastp": {
        "threads": 4,
        "mem_mb": 8000,
        "runtime": 240,
    },
    "moa_seqpurge": {
        "threads": 8,
        "mem_mb": 8000,
        "runtime": 240,
    },
    "moa_flash2": {
        "threads": 4,
        "mem_mb": 8000,
        "runtime": 240,
    },
    "moa_merged_fastqc": {
        "threads": 1,
        "mem_mb": 4000,
        "runtime": 120,
    },
    "map_moa_star": {
        "threads": 8,
        "mem_mb": 64000,
        "runtime": 720,
    },
    "count_moa_peaks": {
        "threads": 1,
        "mem_mb": 8000,
        "runtime": 120,
    },
    "count_atac_peaks": {
        "threads": 1,
        "mem_mb": 8000,
        "runtime": 120,
    },
    "count_chip_peaks": {
        "threads": 1,
        "mem_mb": 8000,
        "runtime": 120,
    },
    "rna_clean_fastqc": {
        "threads": 2,
        "mem_mb": 4000,
        "runtime": 120,
    },
    "atac_clean_fastqc": {
        "threads": 2,
        "mem_mb": 4000,
        "runtime": 120,
    },
    "chip_clean_fastqc": {
        "threads": 2,
        "mem_mb": 4000,
        "runtime": 120,
    },
    "map_atac_star": {
        "threads": 8,
        "mem_mb": 64000,
        "runtime": 720,
    },
    "map_chip_star": {
        "threads": 8,
        "mem_mb": 64000,
        "runtime": 720,
    },
    "deduplicate_atac": {
        "threads": 4,
        "mem_mb": 16000,
        "runtime": 240,
    },
    "deduplicate_chip": {
        "threads": 4,
        "mem_mb": 16000,
        "runtime": 240,
    },
    "map_rna_star": {
        "threads": 8,
        "mem_mb": 64000,
        "runtime": 720,
    },
    "count_rna_genes": {
        "threads": 4,
        "mem_mb": 16000,
        "runtime": 360,
    },
    "aggregate_rna_gene_counts": {
        "threads": 1,
        "mem_mb": 8000,
        "runtime": 120,
    },
    "aggregate_moa_peak_counts": {
        "threads": 1,
        "mem_mb": 8000,
        "runtime": 120,
    },
    "aggregate_atac_peak_counts": {
        "threads": 1,
        "mem_mb": 8000,
        "runtime": 120,
    },
    "aggregate_chip_peak_counts": {
        "threads": 1,
        "mem_mb": 8000,
        "runtime": 120,
    },
    "rna_multiqc": {
        "threads": 1,
        "mem_mb": 8000,
        "runtime": 120,
    },
    "moa_multiqc": {
        "threads": 1,
        "mem_mb": 8000,
        "runtime": 120,
    },
    "rna_gene_qc": {
        "threads": 1,
        "mem_mb": 8000,
        "runtime": 120,
    },
    "moa_peak_qc": {
        "threads": 1,
        "mem_mb": 8000,
        "runtime": 120,
    },
    "atac_peak_qc": {
        "threads": 1,
        "mem_mb": 8000,
        "runtime": 120,
    },
}
RULE_RESOURCE_FIELDS = frozenset(("threads", "mem_mb", "runtime"))

configured_rule_resources = config.get("rule_resources", {})
if not isinstance(configured_rule_resources, Mapping):
    raise ValueError("rule_resources must be a mapping")

RULE_RESOURCES = {
    rule_name: resource_values.copy()
    for rule_name, resource_values in RULE_RESOURCE_DEFAULTS.items()
}

for rule_name, resource_overrides in configured_rule_resources.items():
    rule_path = f"rule_resources.{rule_name}"
    if rule_name not in RULE_RESOURCE_DEFAULTS:
        raise ValueError(f"{rule_path}: unknown rule name")
    if not isinstance(resource_overrides, Mapping):
        raise ValueError(f"{rule_path} must be a mapping")

    for resource_name, value in resource_overrides.items():
        resource_path = f"{rule_path}.{resource_name}"
        if resource_name not in RULE_RESOURCE_FIELDS:
            raise ValueError(f"{resource_path}: unknown resource field")
        if type(value) is not int or value <= 0:
            raise ValueError(f"{resource_path} must be a positive integer")
        RULE_RESOURCES[rule_name][resource_name] = value


# These helpers do not themselves defer evaluation. For optional or
# config-dependent rule inputs, call them inside an input function or lambda so
# validation occurs only when that job is instantiated.
def require_config_value(value, config_key):
    if value is None:
        raise ValueError(f"required config value is not set: {config_key}")
    validated_value = str(value).strip()
    if not validated_value:
        raise ValueError(f"required config value is not set: {config_key}")
    return validated_value


def configured_reference_path(*parts):
    reference_id = require_config_value(
        config["reference_id"],
        "reference_id",
    )
    return "/".join(
        ("resources", "references", reference_id, *parts)
    )


def _get_rule_resource(rule_name, resource_name):
    if rule_name not in RULE_RESOURCES:
        raise ValueError(
            f"unknown internal rule resource name: {rule_name}"
        )
    return RULE_RESOURCES[rule_name][resource_name]


def get_threads(rule_name):
    return _get_rule_resource(rule_name, "threads")


def get_mem_mb(rule_name):
    return _get_rule_resource(rule_name, "mem_mb")


def get_runtime(rule_name):
    return _get_rule_resource(rule_name, "runtime")


SAMPLES = pd.read_csv(
    config["samples"],
    sep="\t",
    dtype=str,
    keep_default_na=False,
)

if list(SAMPLES.columns) != SAMPLE_COLUMNS:
    raise ValueError(f"sample table header must be exactly: {SAMPLE_COLUMNS}")

if SAMPLES.empty:
    raise ValueError("sample table has no data rows")

blank_cells = SAMPLES[SAMPLE_COLUMNS].apply(lambda column: column.str.strip().eq(""))
if blank_cells.any().any():
    raise ValueError("sample table contains blank values in required columns")

ALLOWED_ASSAYS = ["rna", "atac", "chip", "moa"]
invalid_assays = sorted(set(SAMPLES.loc[~SAMPLES["assay"].isin(ALLOWED_ASSAYS), "assay"]))
if invalid_assays:
    raise ValueError(f"invalid assay values: {invalid_assays}; allowed values: {ALLOWED_ASSAYS}")

if SAMPLES["sample_id"].duplicated().any():
    duplicates = sorted(
        SAMPLES.loc[SAMPLES["sample_id"].duplicated(), "sample_id"].unique()
    )
    raise ValueError(f"duplicate sample_id values: {duplicates}")
SAMPLES["library_id"] = SAMPLES["sample_id"]

SAMPLES = SAMPLES.set_index("library_id", drop=False)

LIBRARIES = list(SAMPLES.index)
RNA_LIBRARIES = list(SAMPLES.loc[SAMPLES["assay"] == "rna"].index)
MOA_LIBRARIES = list(SAMPLES.loc[SAMPLES["assay"] == "moa"].index)
ATAC_LIBRARIES = list(SAMPLES.loc[SAMPLES["assay"] == "atac"].index)
CHIP_LIBRARIES = list(SAMPLES.loc[SAMPLES["assay"] == "chip"].index)
