<p align="center">
  <img src="docs/images/logo.png" width="500"/>
</p>

# HyARQ

HyARQ is a Snakemake workflow for allele-specific read counting and quality control in F1 hybrids. It supports RNA-seq, MOA-seq, ATAC-seq, and ChIP-seq data and provides a common framework for preprocessing, mapping to concatenated parental genomes, counting reads over paired parental features, and generating assay-specific quality control outputs.

HyARQ requires two parental reference genomes and user-supplied paired features defining corresponding genes or regulatory regions between the parental genomes. The workflow does not generate paired features.

## Workflow overview

<p align="center">
  <img src="docs/images/workflow.png" width="900"/>
</p>

## Supported assays

| Assay | Quantified feature | Required paired feature input |
| --- | --- | --- |
| RNA-seq | Genes | Paired parental genes |
| MOA-seq | Peak regions | Paired parental peak regions |
| ATAC-seq | Peak regions | Paired parental peak regions |
| ChIP-seq | Peak regions | Paired parental peak regions |

## Installation

Clone the repository and create the HyARQ Conda environment:

```bash
git clone https://github.com/gagelab/hyarq.git
cd hyarq
conda env create -f environment.yaml
conda activate hyarq
```

Software required by individual workflow rules is installed automatically in separate Conda environments when the workflow is run.

## Inputs

HyARQ requires paired-end FASTQ files and two parental reference genomes in FASTA format. RNA-seq additionally requires a GTF annotation for each parental genome.

Paired feature files defining corresponding genes or peak regions between the parental genomes must be supplied by the user. HyARQ validates these files but does not generate the feature pairs.

For RNA-seq, the paired gene file is tab-delimited:

```text
parent1_gene_id	parent2_gene_id
P1_GENE1	P2_GENE1
P1_GENE2	P2_GENE2
```

For MOA-seq, ATAC-seq, and ChIP-seq, the paired peak file contains genomic coordinates:

```text
parent1_chrom	parent1_start	parent1_end	parent2_chrom	parent2_start	parent2_end
chr1	20000	20500	chr1	20005	20505
chr2	50000	50500	chr2	49998	50498
```

Libraries are listed in `config/samples.tsv` using five tab-delimited columns:

```text
sample_id	assay	replicate	r1	r2
sample1	rna	1	path/to/sample1_R1.fastq.gz	path/to/sample1_R2.fastq.gz
sample2	atac	1	path/to/sample2_R1.fastq.gz	path/to/sample2_R2.fastq.gz
```

The `assay` column must be one of `rna`, `moa`, `atac`, or `chip`.

## Configuration

Edit `config/config.yaml` to specify the parental references, paired feature files, and sample table. The main fields to edit are:

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

GTF files are required only when RNA-seq samples are included. Paired feature paths are required only for the corresponding assays being analyzed. Assay-specific parameters and computational resource settings can also be adjusted in `config/config.yaml`.

## Running HyARQ

HyARQ includes profiles for local execution and Slurm clusters.

For local execution, first perform a dry run:

```bash
snakemake --snakefile workflow/Snakefile \
  --configfile config/config.yaml \
  --profile profiles/local \
  --dry-run all
```

Then run HyARQ:

```bash
snakemake --snakefile workflow/Snakefile \
  --configfile config/config.yaml \
  --profile profiles/local \
  all
```

For Slurm execution, copy the supplied profile template and edit the cluster-specific account, partition, and QOS settings as needed:

```bash
cp profiles/slurm/config.yaml.template profiles/slurm/config.yaml
```

Perform a dry run:

```bash
snakemake --snakefile workflow/Snakefile \
  --configfile config/config.yaml \
  --profile profiles/slurm \
  --dry-run all
```

Then run HyARQ:

```bash
snakemake --snakefile workflow/Snakefile \
  --configfile config/config.yaml \
  --profile profiles/slurm \
  all
```

On a Slurm cluster, Snakemake can alternatively be launched from a batch job using `profiles/slurm/run_hyarq.sh.template`:

```bash
cp profiles/slurm/run_hyarq.sh.template profiles/slurm/run_hyarq.sh
sbatch profiles/slurm/run_hyarq.sh
```

If a run is interrupted, rerunning the same Snakemake command will continue from the existing completed outputs.

## Main outputs

HyARQ writes results under the `results/` directory.

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

Additional preprocessing, alignment, and quality control files are retained under `results/preprocessing/`, `results/mapping/`, and `results/qc/`.

## Test dataset

A small test dataset containing one RNA-seq library and one ATAC-seq library is provided under `test/`. It can be used to verify the workflow before analyzing user data.

Run the test with:

```bash
snakemake --snakefile workflow/Snakefile \
  --configfile test/config.yaml \
  --profile profiles/local \
  --dry-run all

snakemake --snakefile workflow/Snakefile \
  --configfile test/config.yaml \
  --profile profiles/local \
  all
```

The test workflow writes to the standard `results/` directory, so it is best run in a clean clone or before running HyARQ on user data.

## License

HyARQ is released under the MIT License. See `LICENSE` for details.
