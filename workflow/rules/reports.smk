# QC summary and report-generation rules.

import shlex


RNA_GENE_QC_SUMMARY = "results/reports/rna/gene_qc_summary.tsv"
RNA_GENE_QC_RETENTION = "results/reports/rna/gene_qc_retention.tsv"
RNA_GENE_QC_REPORT = "results/reports/rna/gene_qc_report.pdf"
RNA_GENE_QC_SAMPLE_REPORT = (
    "results/reports/rna/gene_qc_sample_report.pdf"
)

RNA_GENE_QC_SAMPLE_REPORT_ENABLED = (
    config["rna_gene_qc"]["sample_report"]["enabled"]
)
RNA_GENE_QC_SAMPLE_REPORT_LIBRARIES = (
    config["rna_gene_qc"]["sample_report"]["libraries"]
)
RNA_GENE_QC_SAMPLE_REPORT_PLOT_TYPES = (
    config["rna_gene_qc"]["sample_report"]["plot_types"]
)

if not isinstance(RNA_GENE_QC_SAMPLE_REPORT_ENABLED, bool):
    raise ValueError(
        "rna_gene_qc.sample_report.enabled must be a Boolean"
    )
if (
    not isinstance(RNA_GENE_QC_SAMPLE_REPORT_LIBRARIES, list)
    or not RNA_GENE_QC_SAMPLE_REPORT_LIBRARIES
    or any(
        not isinstance(library_id, str) or not library_id
        for library_id in RNA_GENE_QC_SAMPLE_REPORT_LIBRARIES
    )
):
    raise ValueError(
        "rna_gene_qc.sample_report.libraries must be a nonempty list "
        "of nonempty strings"
    )
if (
    not isinstance(RNA_GENE_QC_SAMPLE_REPORT_PLOT_TYPES, list)
    or not RNA_GENE_QC_SAMPLE_REPORT_PLOT_TYPES
    or any(
        not isinstance(plot_type, str) or not plot_type
        for plot_type in RNA_GENE_QC_SAMPLE_REPORT_PLOT_TYPES
    )
):
    raise ValueError(
        "rna_gene_qc.sample_report.plot_types must be a nonempty list "
        "of nonempty strings"
    )

RNA_GENE_QC_OUTPUTS = [
    RNA_GENE_QC_SUMMARY,
    RNA_GENE_QC_RETENTION,
    RNA_GENE_QC_REPORT,
]
RNA_GENE_QC_SAMPLE_REPORT_ARGS = ""
if RNA_GENE_QC_SAMPLE_REPORT_ENABLED:
    RNA_GENE_QC_OUTPUTS.append(RNA_GENE_QC_SAMPLE_REPORT)
    RNA_GENE_QC_SAMPLE_REPORT_ARGS = shlex.join(
        [
            "--sample-report",
            RNA_GENE_QC_SAMPLE_REPORT,
            "--sample-libraries",
            *RNA_GENE_QC_SAMPLE_REPORT_LIBRARIES,
            "--sample-plot-types",
            *RNA_GENE_QC_SAMPLE_REPORT_PLOT_TYPES,
        ]
    )


rule rna_multiqc:
    input:
        raw_fastqc=expand(
            "results/qc/fastqc/raw/{library_id}",
            library_id=RNA_LIBRARIES,
        ),
        clean_fastqc=expand(
            "results/qc/fastqc/clean/rna/{library_id}",
            library_id=RNA_LIBRARIES,
        ),
        fastp_json=expand(
            "results/qc/fastp/rna/{library_id}.json",
            library_id=RNA_LIBRARIES,
        ),
        star_logs=expand(
            "results/mapping/rna/{library_id}.Log.final.out",
            library_id=RNA_LIBRARIES,
        ),
        config="config/multiqc_config.yaml",
    output:
        html="results/reports/rna/multiqc_report.html",
        data=directory("results/reports/rna/multiqc_report_data"),
    params:
        outdir="results/reports/rna",
    conda:
        "../envs/multiqc.yaml"
    log:
        "logs/reports/rna/multiqc.log"
    shell:
        r"""
        mkdir -p {params.outdir:q}
        mkdir -p "$(dirname {log:q})"
        multiqc \
            {input.raw_fastqc:q} \
            {input.clean_fastqc:q} \
            {input.fastp_json:q} \
            {input.star_logs:q} \
            --outdir {params.outdir:q} \
            --filename multiqc_report.html \
            --config {input.config:q} \
            --force \
            > {log:q} 2>&1
        """


rule rna_gene_qc:
    input:
        count_tables=expand(
            "results/counts/rna/{library_id}/gene_counts.tsv",
            library_id=RNA_LIBRARIES,
        )
    output:
        RNA_GENE_QC_OUTPUTS
    params:
        summary=RNA_GENE_QC_SUMMARY,
        retention=RNA_GENE_QC_RETENTION,
        report=RNA_GENE_QC_REPORT,
        sample_report_args=RNA_GENE_QC_SAMPLE_REPORT_ARGS,
        parent1_label=config["parents"]["parent1"]["prefix"],
        parent2_label=config["parents"]["parent2"]["prefix"],
    threads: 1
    conda:
        "../envs/rna_gene_qc.yaml"
    log:
        "logs/reports/rna/gene_qc.log"
    shell:
        r"""
        mkdir -p "$(dirname {log:q})"
        python workflow/scripts/rna_gene_qc.py \
            --count-tables {input.count_tables:q} \
            --summary {params.summary:q} \
            --retention {params.retention:q} \
            --report {params.report:q} \
            --parent1-label {params.parent1_label:q} \
            --parent2-label {params.parent2_label:q} \
            {params.sample_report_args} > {log:q} 2>&1
        """
