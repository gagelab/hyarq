rule raw_fastqc:
    input:
        r1=lambda wildcards: SAMPLES.loc[wildcards.library_id, "r1"],
        r2=lambda wildcards: SAMPLES.loc[wildcards.library_id, "r2"],
    output:
        directory("results/qc/fastqc/raw/{library_id}")
    threads: get_threads("raw_fastqc")
    resources:
        mem_mb=get_mem_mb("raw_fastqc"),
        runtime=get_runtime("raw_fastqc"),
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


rule rna_fastp:
    input:
        r1=lambda wildcards: SAMPLES.loc[wildcards.library_id, "r1"],
        r2=lambda wildcards: SAMPLES.loc[wildcards.library_id, "r2"],
    output:
        r1="results/preprocessing/rna/{library_id}/{library_id}_R1.fq.gz",
        r2="results/preprocessing/rna/{library_id}/{library_id}_R2.fq.gz",
        html="results/qc/fastp/rna/{library_id}.html",
        json="results/qc/fastp/rna/{library_id}.json",
    threads: get_threads("rna_fastp")
    resources:
        mem_mb=get_mem_mb("rna_fastp"),
        runtime=get_runtime("rna_fastp"),
    conda:
        "../envs/fastp.yaml"
    log:
        "logs/fastp/rna/{library_id}.log"
    shell:
        """
        mkdir -p $(dirname {output.r1:q})
        mkdir -p $(dirname {output.html:q})
        mkdir -p $(dirname {log:q})
        fastp \
            --in1 {input.r1:q} \
            --in2 {input.r2:q} \
            --detect_adapter_for_pe \
            --out1 {output.r1:q} \
            --out2 {output.r2:q} \
            --html {output.html:q} \
            --json {output.json:q} \
            --thread {threads} \
            > {log:q} 2>&1
        """


rule rna_clean_fastqc:
    input:
        r1="results/preprocessing/rna/{library_id}/{library_id}_R1.fq.gz",
        r2="results/preprocessing/rna/{library_id}/{library_id}_R2.fq.gz",
    output:
        directory("results/qc/fastqc/clean/rna/{library_id}")
    threads: get_threads("rna_clean_fastqc")
    resources:
        mem_mb=get_mem_mb("rna_clean_fastqc"),
        runtime=get_runtime("rna_clean_fastqc"),
    conda:
        "../envs/fastqc.yaml"
    log:
        "logs/fastqc/clean/rna/{library_id}.log"
    shell:
        """
        mkdir -p {output:q}
        mkdir -p $(dirname {log:q})
        fastqc --threads {threads} --outdir {output:q} {input.r1:q} {input.r2:q} > {log:q} 2>&1
        """
