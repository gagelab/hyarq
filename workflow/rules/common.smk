import pandas as pd


SAMPLE_COLUMNS = ["sample_id", "assay", "replicate", "r1", "r2"]

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

SAMPLES["library_id"] = SAMPLES["sample_id"] + "_rep" + SAMPLES["replicate"]
if SAMPLES["library_id"].duplicated().any():
    duplicates = sorted(SAMPLES.loc[SAMPLES["library_id"].duplicated(), "library_id"].unique())
    raise ValueError(f"duplicate library_id values: {duplicates}")

SAMPLES = SAMPLES.set_index("library_id", drop=False)

LIBRARIES = list(SAMPLES.index)
RNA_LIBRARIES = list(SAMPLES.loc[SAMPLES["assay"] == "rna"].index)
