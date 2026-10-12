# HyARQ Quick Start

This guide provides step-by-step instructions for preparing inputs, configuring, and running HyARQ.

All commands assume that they are run from the root directory of the HyARQ repository.

## 1. Install HyARQ

Clone the repository and create the HyARQ Conda environment:

```bash
git clone https://github.com/gagelab/hyarq.git
cd hyarq
conda env create -f environment.yaml
conda activate hyarq
```

Software required by individual workflow rules is installed automatically in separate Conda environments when the workflow is run.

## 2. Prepare inputs

HyARQ requires paired-end FASTQ files and two parental reference genomes in FASTA format. RNA-seq additionally requires a GTF annotation for each parental genome. Paired feature files defining corresponding genes or peak regions between the two parental genomes must also be provided by the user.

### Paired feature files

For RNA-seq, provide a tab-delimited file containing paired parental gene IDs:

```text
parent1_gene_id	parent2_gene_id
P1_TEST_GENE1	P2_TEST_GENE1
P1_TEST_GENE2	P2_TEST_GENE2
```

For MOA-seq, ATAC-seq, and ChIP-seq, provide a tab-delimited file containing the coordinates of paired parental peak regions:

```text
parent1_chrom	parent1_start	parent1_end	parent2_chrom	parent2_start	parent2_end
chr1	20000	20500	chr1	20005	20505
chr2	50000	50500	chr2	49998	50498
```

HyARQ validates the supplied paired feature files but does not generate the feature pairs.

### Sample table

Libraries are listed in `config/samples.tsv` using five tab-delimited columns:

```text
sample_id	assay	replicate	r1	r2
sample1	rna	1	path/to/sample1_R1.fastq.gz	path/to/sample1_R2.fastq.gz
sample2	atac	1	path/to/sample2_R1.fastq.gz	path/to/sample2_R2.fastq.gz
```

The `assay` column must be one of `rna`, `moa`, `atac`, or `chip`. The `r1` and `r2` columns specify the paired FASTQ files for each library.

## 3. Configure HyARQ

Edit `config/config.yaml` to specify the parental references, paired feature files, and sample table.

The main fields to edit in the supplied configuration are:

```yaml
reference_id: my_reference

parents:
  parent1:
    prefix: P1
    fasta: path/to/parent1.fa
    gtf: path/to/parent1.gtf
  parent2:
    prefix: P2
    fasta: path/to/parent2.fa
    gtf: path/to/parent2.gtf

rna_gene_pairs: path/to/gene_pairs.tsv
moa_peak_pairs: path/to/moa_peak_pairs.tsv
atac_peak_pairs: path/to/atac_peak_pairs.tsv
chip_peak_pairs: path/to/chip_peak_pairs.tsv

samples: config/samples.tsv
```

GTF files are required only when RNA-seq samples are included. Paired feature paths are required only for the corresponding assays being analyzed.

The remaining assay-specific parameters and computational resource settings in `config/config.yaml` can be adjusted if needed.

## 4. Choose an execution profile

HyARQ includes profiles for local execution and Slurm clusters.

### Local execution

The local profile is provided in `profiles/local/config.yaml` and can be used directly.

### Slurm execution

A Slurm profile template is provided in `profiles/slurm/config.yaml.template`. Copy the template before running HyARQ:

```bash
cp profiles/slurm/config.yaml.template profiles/slurm/config.yaml
```

Edit `profiles/slurm/config.yaml` and replace the cluster-specific placeholders, including the account, partition, and QOS where applicable. The profile controls the resources and submission settings used for individual workflow jobs.

## 5. Validate the workflow

Before starting the analysis, perform a dry run to check the configuration and planned jobs.

For local execution:

```bash
snakemake --snakefile workflow/Snakefile \
  --configfile config/config.yaml \
  --profile profiles/local \
  --dry-run all
```

For Slurm execution:

```bash
snakemake --snakefile workflow/Snakefile \
  --configfile config/config.yaml \
  --profile profiles/slurm \
  --dry-run all
```

A successful dry run should complete without configuration or missing-input errors and display the jobs that HyARQ will execute.

## 6. Run HyARQ

For local execution:

```bash
snakemake --snakefile workflow/Snakefile \
  --configfile config/config.yaml \
  --profile profiles/local \
  all
```

For Slurm execution:

```bash
snakemake --snakefile workflow/Snakefile \
  --configfile config/config.yaml \
  --profile profiles/slurm \
  all
```

On a Slurm cluster, Snakemake can also be launched from a batch job using the provided `profiles/slurm/run_hyarq.sh.template`. Copy the template, edit the cluster-specific Slurm settings, and submit it with `sbatch`.

```bash
cp profiles/slurm/run_hyarq.sh.template profiles/slurm/run_hyarq.sh
sbatch profiles/slurm/run_hyarq.sh
```

If a run is interrupted, rerunning the same Snakemake command will continue from the existing completed outputs.


## 7. Main outputs

HyARQ writes results under the `results/` directory. The main outputs are:

| Output | Location |
| --- | --- |
| Per-library parental counts | `results/counts/<assay>/<library_id>/` |
| Per-library filtering summary | `results/counts/<assay>/<library_id>/filter_summary.tsv` |
| Assay-wide paired feature count matrix | `results/counts/<assay>/` |
| Mapping summary | `results/mapping/<assay>/<library_id>.mapping_summary.tsv` |
| MultiQC report | `results/reports/<assay>/multiqc_report.html` |
| Assay-specific QC summary | `results/reports/<assay>/*_qc_summary.tsv` |
| Count-threshold retention summary | `results/reports/<assay>/*_qc_retention.tsv` |
| Assay-specific QC report | `results/reports/<assay>/*_qc_report.pdf` |

For RNA-seq, the main count files are `gene_counts.tsv` and `gene_pair_count_matrix.tsv`. For MOA-seq, ATAC-seq, and ChIP-seq, the corresponding files are `peak_counts.tsv` and `peak_pair_count_matrix.tsv`.

Additional preprocessing, alignment, and quality-control files are retained under `results/preprocessing/`, `results/mapping/`, and `results/qc/`.

## 8. Run the test dataset

HyARQ includes a small test dataset containing one RNA-seq library and one ATAC-seq library. The test inputs and configuration are provided under `test/`. The test workflow writes to the standard `results/` directory, so it is best run in a clean clone or before running HyARQ on your own data.

First perform a dry run:

```bash
snakemake --snakefile workflow/Snakefile \
  --configfile test/config.yaml \
  --profile profiles/local \
  --dry-run all
```

Then run the test workflow:

```bash
snakemake --snakefile workflow/Snakefile \
  --configfile test/config.yaml \
  --profile profiles/local \
  all
```

A successful run generates RNA-seq and ATAC-seq count files, assay-wide count matrices, MultiQC reports, and assay-specific QC reports under `results/`.

Key outputs include:

```text
results/counts/rna/test_rna/gene_counts.tsv
results/counts/rna/gene_pair_count_matrix.tsv
results/reports/rna/gene_qc_report.pdf

results/counts/atac/test_atac/peak_counts.tsv
results/counts/atac/peak_pair_count_matrix.tsv
results/reports/atac/peak_qc_report.pdf
```
