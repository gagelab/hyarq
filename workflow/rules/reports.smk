# QC summary and report-generation rules.

import shlex


RNA_GENE_QC_SUMMARY = "results/reports/rna/gene_qc_summary.tsv"
RNA_GENE_QC_RETENTION = "results/reports/rna/gene_qc_retention.tsv"
RNA_GENE_QC_REPORT = "results/reports/rna/gene_qc_report.pdf"
RNA_GENE_QC_SAMPLE_REPORT = (
    "results/reports/rna/gene_qc_sample_report.pdf"
)
MOA_PEAK_QC_SUMMARY = "results/reports/moa/peak_qc_summary.tsv"
MOA_PEAK_QC_RETENTION = "results/reports/moa/peak_qc_retention.tsv"
MOA_PEAK_QC_REPORT = "results/reports/moa/peak_qc_report.pdf"
MOA_PEAK_QC_SAMPLE_REPORT = (
    "results/reports/moa/peak_qc_sample_report.pdf"
)
ATAC_PEAK_QC_SUMMARY = "results/reports/atac/peak_qc_summary.tsv"
ATAC_PEAK_QC_RETENTION = "results/reports/atac/peak_qc_retention.tsv"
ATAC_PEAK_QC_REPORT = "results/reports/atac/peak_qc_report.pdf"
ATAC_PEAK_QC_SAMPLE_REPORT = (
    "results/reports/atac/peak_qc_sample_report.pdf"
)
CHIP_PEAK_QC_SUMMARY = "results/reports/chip/peak_qc_summary.tsv"
CHIP_PEAK_QC_RETENTION = "results/reports/chip/peak_qc_retention.tsv"
CHIP_PEAK_QC_REPORT = "results/reports/chip/peak_qc_report.pdf"
CHIP_PEAK_QC_SAMPLE_REPORT = (
    "results/reports/chip/peak_qc_sample_report.pdf"
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

MOA_PEAK_QC_SAMPLE_REPORT_ENABLED = (
    config["moa_peak_qc"]["sample_report"]["enabled"]
)
MOA_PEAK_QC_SAMPLE_REPORT_LIBRARIES = (
    config["moa_peak_qc"]["sample_report"]["libraries"]
)
MOA_PEAK_QC_SAMPLE_REPORT_PLOT_TYPES = (
    config["moa_peak_qc"]["sample_report"]["plot_types"]
)

if not isinstance(MOA_PEAK_QC_SAMPLE_REPORT_ENABLED, bool):
    raise ValueError(
        "moa_peak_qc.sample_report.enabled must be a Boolean"
    )
if (
    not isinstance(MOA_PEAK_QC_SAMPLE_REPORT_LIBRARIES, list)
    or not MOA_PEAK_QC_SAMPLE_REPORT_LIBRARIES
    or any(
        not isinstance(library_id, str) or not library_id
        for library_id in MOA_PEAK_QC_SAMPLE_REPORT_LIBRARIES
    )
):
    raise ValueError(
        "moa_peak_qc.sample_report.libraries must be a nonempty list "
        "of nonempty strings"
    )
if (
    not isinstance(MOA_PEAK_QC_SAMPLE_REPORT_PLOT_TYPES, list)
    or not MOA_PEAK_QC_SAMPLE_REPORT_PLOT_TYPES
    or any(
        not isinstance(plot_type, str) or not plot_type
        for plot_type in MOA_PEAK_QC_SAMPLE_REPORT_PLOT_TYPES
    )
):
    raise ValueError(
        "moa_peak_qc.sample_report.plot_types must be a nonempty list "
        "of nonempty strings"
    )

MOA_PEAK_QC_OUTPUTS = [
    MOA_PEAK_QC_SUMMARY,
    MOA_PEAK_QC_RETENTION,
    MOA_PEAK_QC_REPORT,
]
MOA_PEAK_QC_SAMPLE_REPORT_ARGS = ""
if MOA_PEAK_QC_SAMPLE_REPORT_ENABLED:
    MOA_PEAK_QC_OUTPUTS.append(MOA_PEAK_QC_SAMPLE_REPORT)
    MOA_PEAK_QC_SAMPLE_REPORT_ARGS = shlex.join(
        [
            "--sample-report",
            MOA_PEAK_QC_SAMPLE_REPORT,
            "--sample-libraries",
            *MOA_PEAK_QC_SAMPLE_REPORT_LIBRARIES,
            "--sample-plot-types",
            *MOA_PEAK_QC_SAMPLE_REPORT_PLOT_TYPES,
        ]
    )

ATAC_PEAK_QC_OUTPUTS = [
    ATAC_PEAK_QC_SUMMARY,
    ATAC_PEAK_QC_RETENTION,
    ATAC_PEAK_QC_REPORT,
]
ATAC_PEAK_QC_SAMPLE_REPORT_ARGS = ""
if ATAC_LIBRARIES:
    ATAC_PEAK_QC_SAMPLE_REPORT_ENABLED = (
        config["atac_peak_qc"]["sample_report"]["enabled"]
    )
    ATAC_PEAK_QC_SAMPLE_REPORT_LIBRARIES = (
        config["atac_peak_qc"]["sample_report"]["libraries"]
    )
    ATAC_PEAK_QC_SAMPLE_REPORT_PLOT_TYPES = (
        config["atac_peak_qc"]["sample_report"]["plot_types"]
    )

    if not isinstance(ATAC_PEAK_QC_SAMPLE_REPORT_ENABLED, bool):
        raise ValueError(
            "atac_peak_qc.sample_report.enabled must be a Boolean"
        )
    if (
        not isinstance(ATAC_PEAK_QC_SAMPLE_REPORT_LIBRARIES, list)
        or not ATAC_PEAK_QC_SAMPLE_REPORT_LIBRARIES
        or any(
            not isinstance(library_id, str) or not library_id
            for library_id in ATAC_PEAK_QC_SAMPLE_REPORT_LIBRARIES
        )
    ):
        raise ValueError(
            "atac_peak_qc.sample_report.libraries must be a nonempty "
            "list of nonempty strings"
        )
    if (
        not isinstance(ATAC_PEAK_QC_SAMPLE_REPORT_PLOT_TYPES, list)
        or not ATAC_PEAK_QC_SAMPLE_REPORT_PLOT_TYPES
        or any(
            not isinstance(plot_type, str) or not plot_type
            for plot_type in ATAC_PEAK_QC_SAMPLE_REPORT_PLOT_TYPES
        )
    ):
        raise ValueError(
            "atac_peak_qc.sample_report.plot_types must be a nonempty "
            "list of nonempty strings"
        )

    if ATAC_PEAK_QC_SAMPLE_REPORT_ENABLED:
        ATAC_PEAK_QC_OUTPUTS.append(ATAC_PEAK_QC_SAMPLE_REPORT)
        ATAC_PEAK_QC_SAMPLE_REPORT_ARGS = shlex.join(
            [
                "--sample-report",
                ATAC_PEAK_QC_SAMPLE_REPORT,
                "--sample-libraries",
                *ATAC_PEAK_QC_SAMPLE_REPORT_LIBRARIES,
                "--sample-plot-types",
                *ATAC_PEAK_QC_SAMPLE_REPORT_PLOT_TYPES,
            ]
        )


CHIP_PEAK_QC_OUTPUTS = [
    CHIP_PEAK_QC_SUMMARY,
    CHIP_PEAK_QC_RETENTION,
    CHIP_PEAK_QC_REPORT,
]
CHIP_PEAK_QC_SAMPLE_REPORT_ARGS = ""
if CHIP_LIBRARIES:
    CHIP_PEAK_QC_SAMPLE_REPORT_ENABLED = (
        config["chip_peak_qc"]["sample_report"]["enabled"]
    )
    CHIP_PEAK_QC_SAMPLE_REPORT_LIBRARIES = (
        config["chip_peak_qc"]["sample_report"]["libraries"]
    )
    CHIP_PEAK_QC_SAMPLE_REPORT_PLOT_TYPES = (
        config["chip_peak_qc"]["sample_report"]["plot_types"]
    )

    if not isinstance(CHIP_PEAK_QC_SAMPLE_REPORT_ENABLED, bool):
        raise ValueError(
            "chip_peak_qc.sample_report.enabled must be a Boolean"
        )
    if (
        not isinstance(CHIP_PEAK_QC_SAMPLE_REPORT_LIBRARIES, list)
        or not CHIP_PEAK_QC_SAMPLE_REPORT_LIBRARIES
        or any(
            not isinstance(library_id, str) or not library_id
            for library_id in CHIP_PEAK_QC_SAMPLE_REPORT_LIBRARIES
        )
    ):
        raise ValueError(
            "chip_peak_qc.sample_report.libraries must be a nonempty "
            "list of nonempty strings"
        )
    if (
        not isinstance(CHIP_PEAK_QC_SAMPLE_REPORT_PLOT_TYPES, list)
        or not CHIP_PEAK_QC_SAMPLE_REPORT_PLOT_TYPES
        or any(
            not isinstance(plot_type, str) or not plot_type
            for plot_type in CHIP_PEAK_QC_SAMPLE_REPORT_PLOT_TYPES
        )
    ):
        raise ValueError(
            "chip_peak_qc.sample_report.plot_types must be a nonempty "
            "list of nonempty strings"
        )

    if CHIP_PEAK_QC_SAMPLE_REPORT_ENABLED:
        CHIP_PEAK_QC_OUTPUTS.append(CHIP_PEAK_QC_SAMPLE_REPORT)
        CHIP_PEAK_QC_SAMPLE_REPORT_ARGS = shlex.join(
            [
                "--sample-report",
                CHIP_PEAK_QC_SAMPLE_REPORT,
                "--sample-libraries",
                *CHIP_PEAK_QC_SAMPLE_REPORT_LIBRARIES,
                "--sample-plot-types",
                *CHIP_PEAK_QC_SAMPLE_REPORT_PLOT_TYPES,
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
        config=MULTIQC_CONFIG,
    output:
        html="results/reports/rna/multiqc_report.html",
        data=directory("results/reports/rna/multiqc_report_data"),
    params:
        outdir="results/reports/rna",
    threads: get_threads("rna_multiqc")
    resources:
        mem_mb=get_mem_mb("rna_multiqc"),
        runtime=get_runtime("rna_multiqc"),
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


rule moa_multiqc:
    input:
        raw_fastqc=expand(
            "results/qc/fastqc/raw/{library_id}",
            library_id=MOA_LIBRARIES,
        ),
        merged_fastqc=expand(
            "results/qc/fastqc/merged/moa/{library_id}",
            library_id=MOA_LIBRARIES,
        ),
        star_logs=expand(
            "results/mapping/moa/{library_id}.Log.final.out",
            library_id=MOA_LIBRARIES,
        ),
        config=MULTIQC_CONFIG,
    output:
        html="results/reports/moa/multiqc_report.html",
        data=directory("results/reports/moa/multiqc_report_data"),
    params:
        outdir="results/reports/moa",
    threads: get_threads("moa_multiqc")
    resources:
        mem_mb=get_mem_mb("moa_multiqc"),
        runtime=get_runtime("moa_multiqc"),
    conda:
        "../envs/multiqc.yaml"
    log:
        "logs/reports/moa/multiqc.log"
    shell:
        r"""
        mkdir -p {params.outdir:q}
        mkdir -p "$(dirname {log:q})"
        multiqc \
            {input.raw_fastqc:q} \
            {input.merged_fastqc:q} \
            {input.star_logs:q} \
            --outdir {params.outdir:q} \
            --filename multiqc_report.html \
            --config {input.config:q} \
            --title "HyARQ MOA-seq QC" \
            --force \
            > {log:q} 2>&1
        """


rule atac_multiqc:
    input:
        raw_fastqc=expand(
            "results/qc/fastqc/raw/{library_id}",
            library_id=ATAC_LIBRARIES,
        ),
        clean_fastqc=expand(
            "results/qc/fastqc/clean/atac/{library_id}",
            library_id=ATAC_LIBRARIES,
        ),
        fastp_json=expand(
            "results/qc/fastp/atac/{library_id}.json",
            library_id=ATAC_LIBRARIES,
        ),
        star_logs=expand(
            "results/mapping/atac/{library_id}.Log.final.out",
            library_id=ATAC_LIBRARIES,
        ),
        config=MULTIQC_CONFIG,
    output:
        html="results/reports/atac/multiqc_report.html",
        data=directory("results/reports/atac/multiqc_report_data"),
    params:
        outdir="results/reports/atac",
    threads: get_threads("atac_multiqc")
    resources:
        mem_mb=get_mem_mb("atac_multiqc"),
        runtime=get_runtime("atac_multiqc"),
    conda:
        "../envs/multiqc.yaml"
    log:
        "logs/reports/atac/multiqc.log"
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
            --title "HyARQ ATAC-seq QC" \
            --force \
            > {log:q} 2>&1
        """


rule chip_multiqc:
    input:
        raw_fastqc=expand(
            "results/qc/fastqc/raw/{library_id}",
            library_id=CHIP_LIBRARIES,
        ),
        clean_fastqc=expand(
            "results/qc/fastqc/clean/chip/{library_id}",
            library_id=CHIP_LIBRARIES,
        ),
        fastp_json=expand(
            "results/qc/fastp/chip/{library_id}.json",
            library_id=CHIP_LIBRARIES,
        ),
        star_logs=expand(
            "results/mapping/chip/{library_id}.Log.final.out",
            library_id=CHIP_LIBRARIES,
        ),
        config=MULTIQC_CONFIG,
    output:
        html="results/reports/chip/multiqc_report.html",
        data=directory("results/reports/chip/multiqc_report_data"),
    params:
        outdir="results/reports/chip",
    threads: get_threads("chip_multiqc")
    resources:
        mem_mb=get_mem_mb("chip_multiqc"),
        runtime=get_runtime("chip_multiqc"),
    conda:
        "../envs/multiqc.yaml"
    log:
        "logs/reports/chip/multiqc.log"
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
            --title "HyARQ ChIP-seq QC" \
            --force \
            > {log:q} 2>&1
        """


rule rna_gene_qc:
    input:
        count_tables=expand(
            "results/counts/rna/{library_id}/gene_counts.tsv",
            library_id=RNA_LIBRARIES,
        ),
        script=RNA_GENE_QC_SCRIPT,
    output:
        RNA_GENE_QC_OUTPUTS
    params:
        summary=RNA_GENE_QC_SUMMARY,
        retention=RNA_GENE_QC_RETENTION,
        report=RNA_GENE_QC_REPORT,
        sample_report_args=RNA_GENE_QC_SAMPLE_REPORT_ARGS,
        parent1_label=lambda wildcards: require_config_value(
            config["parents"]["parent1"]["prefix"],
            "parents.parent1.prefix",
        ),
        parent2_label=lambda wildcards: require_config_value(
            config["parents"]["parent2"]["prefix"],
            "parents.parent2.prefix",
        ),
    threads: get_threads("rna_gene_qc")
    resources:
        mem_mb=get_mem_mb("rna_gene_qc"),
        runtime=get_runtime("rna_gene_qc"),
    conda:
        "../envs/rna_gene_qc.yaml"
    log:
        "logs/reports/rna/gene_qc.log"
    shell:
        r"""
        mkdir -p "$(dirname {log:q})"
        python {input.script:q} \
            --count-tables {input.count_tables:q} \
            --summary {params.summary:q} \
            --retention {params.retention:q} \
            --report {params.report:q} \
            --parent1-label {params.parent1_label:q} \
            --parent2-label {params.parent2_label:q} \
            {params.sample_report_args} > {log:q} 2>&1
        """


rule moa_peak_qc:
    input:
        count_tables=expand(
            "results/counts/moa/{library_id}/peak_counts.tsv",
            library_id=MOA_LIBRARIES,
        ),
        script=PEAK_PAIR_QC_SCRIPT,
    output:
        MOA_PEAK_QC_OUTPUTS
    params:
        summary=MOA_PEAK_QC_SUMMARY,
        retention=MOA_PEAK_QC_RETENTION,
        report=MOA_PEAK_QC_REPORT,
        sample_report_args=MOA_PEAK_QC_SAMPLE_REPORT_ARGS,
        parent1_label=lambda wildcards: require_config_value(
            config["parents"]["parent1"]["prefix"],
            "parents.parent1.prefix",
        ),
        parent2_label=lambda wildcards: require_config_value(
            config["parents"]["parent2"]["prefix"],
            "parents.parent2.prefix",
        ),
    threads: get_threads("moa_peak_qc")
    resources:
        mem_mb=get_mem_mb("moa_peak_qc"),
        runtime=get_runtime("moa_peak_qc"),
    conda:
        "../envs/rna_gene_qc.yaml"
    log:
        "logs/reports/moa/peak_qc.log"
    shell:
        r"""
        mkdir -p "$(dirname {log:q})"
        python {input.script:q} \
            --count-tables {input.count_tables:q} \
            --assay moa \
            --summary {params.summary:q} \
            --retention {params.retention:q} \
            --report {params.report:q} \
            --parent1-label {params.parent1_label:q} \
            --parent2-label {params.parent2_label:q} \
            {params.sample_report_args} > {log:q} 2>&1
        """


rule atac_peak_qc:
    input:
        count_tables=expand(
            "results/counts/atac/{library_id}/peak_counts.tsv",
            library_id=ATAC_LIBRARIES,
        ),
        script=PEAK_PAIR_QC_SCRIPT,
    output:
        ATAC_PEAK_QC_OUTPUTS
    params:
        summary=ATAC_PEAK_QC_SUMMARY,
        retention=ATAC_PEAK_QC_RETENTION,
        report=ATAC_PEAK_QC_REPORT,
        sample_report_args=ATAC_PEAK_QC_SAMPLE_REPORT_ARGS,
        parent1_label=lambda wildcards: require_config_value(
            config["parents"]["parent1"]["prefix"],
            "parents.parent1.prefix",
        ),
        parent2_label=lambda wildcards: require_config_value(
            config["parents"]["parent2"]["prefix"],
            "parents.parent2.prefix",
        ),
    threads: get_threads("atac_peak_qc")
    resources:
        mem_mb=get_mem_mb("atac_peak_qc"),
        runtime=get_runtime("atac_peak_qc"),
    conda:
        "../envs/rna_gene_qc.yaml"
    log:
        "logs/reports/atac/peak_qc.log"
    shell:
        r"""
        mkdir -p "$(dirname {log:q})"
        python {input.script:q} \
            --count-tables {input.count_tables:q} \
            --assay atac \
            --summary {params.summary:q} \
            --retention {params.retention:q} \
            --report {params.report:q} \
            --parent1-label {params.parent1_label:q} \
            --parent2-label {params.parent2_label:q} \
            {params.sample_report_args} > {log:q} 2>&1
        """


rule chip_peak_qc:
    input:
        count_tables=expand(
            "results/counts/chip/{library_id}/peak_counts.tsv",
            library_id=CHIP_LIBRARIES,
        ),
        script=PEAK_PAIR_QC_SCRIPT,
    output:
        CHIP_PEAK_QC_OUTPUTS
    params:
        summary=CHIP_PEAK_QC_SUMMARY,
        retention=CHIP_PEAK_QC_RETENTION,
        report=CHIP_PEAK_QC_REPORT,
        sample_report_args=CHIP_PEAK_QC_SAMPLE_REPORT_ARGS,
        parent1_label=lambda wildcards: require_config_value(
            config["parents"]["parent1"]["prefix"],
            "parents.parent1.prefix",
        ),
        parent2_label=lambda wildcards: require_config_value(
            config["parents"]["parent2"]["prefix"],
            "parents.parent2.prefix",
        ),
    threads: get_threads("chip_peak_qc")
    resources:
        mem_mb=get_mem_mb("chip_peak_qc"),
        runtime=get_runtime("chip_peak_qc"),
    conda:
        "../envs/rna_gene_qc.yaml"
    log:
        "logs/reports/chip/peak_qc.log"
    shell:
        r"""
        mkdir -p "$(dirname {log:q})"
        python {input.script:q} \
            --count-tables {input.count_tables:q} \
            --assay chip \
            --summary {params.summary:q} \
            --retention {params.retention:q} \
            --report {params.report:q} \
            --parent1-label {params.parent1_label:q} \
            --parent2-label {params.parent2_label:q} \
            {params.sample_report_args} > {log:q} 2>&1
        """
