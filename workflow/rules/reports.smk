# QC summary and report-generation rules.

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
        summary="results/reports/rna/gene_qc_summary.tsv",
        retention="results/reports/rna/gene_qc_retention.tsv",
        report="results/reports/rna/gene_qc_report.pdf",
    params:
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
            --summary {output.summary:q} \
            --retention {output.retention:q} \
            --report {output.report:q} \
            --parent1-label {params.parent1_label:q} \
            --parent2-label {params.parent2_label:q} \
            > {log:q} 2>&1
        """
