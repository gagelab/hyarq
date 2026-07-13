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
        data=directory("results/reports/rna/multiqc_data"),
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
