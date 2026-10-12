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

## Quick start

For complete installation and usage instructions, including input preparation, configuration, local and Slurm execution, and the bundled test dataset, see `docs/HYARQ_quickstart.md`.

## Inputs

HyARQ is configured using `config/config.yaml` and `config/samples.tsv`. Analyses require paired-end FASTQ files, two parental genome assemblies in FASTA format, and user-supplied paired features. Parental GTF annotations are additionally required for RNA-seq.

## Main outputs

HyARQ generates per-library parental feature counts, assay-wide paired feature count matrices, mapping summaries, MultiQC reports, and assay-specific quality control summaries and reports. Outputs are written under the `results/` directory.

## Test dataset

A small test dataset containing one RNA-seq library and one ATAC-seq library is provided under `test/`. It can be used to verify the HyARQ installation and workflow before analyzing user data.

## License

HyARQ is released under the MIT License. See `LICENSE` for details.
