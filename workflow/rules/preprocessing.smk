rule raw_fastqc:
    input:
        r1=lambda wildcards: SAMPLES.loc[wildcards.library_id, "r1"],
        r2=lambda wildcards: SAMPLES.loc[wildcards.library_id, "r2"],
    output:
        directory("results/qc/fastqc/raw/{library_id}")
    threads: 2
    conda:
        "../envs/fastqc.yaml"
    log:
        "logs/fastqc/raw/{library_id}.log"
    shell:
        """
        mkdir -p {output:q}
        mkdir -p $(dirname {log:q})
        fastqc --threads {threads} --outdir {output:q} {input.r1:q} {input.r2:q} > {log:q} 2>&1
        """
