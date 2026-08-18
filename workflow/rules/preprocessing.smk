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


rule moa_seqpurge:
    input:
        r1=lambda wildcards: SAMPLES.loc[wildcards.library_id, "r1"],
        r2=lambda wildcards: SAMPLES.loc[wildcards.library_id, "r2"],
    output:
        r1="results/preprocessing/moa/{library_id}/{library_id}_R1.fq.gz",
        r2="results/preprocessing/moa/{library_id}/{library_id}_R2.fq.gz",
        summary="results/qc/seqpurge/moa/{library_id}.summary.txt",
    threads: get_threads("moa_seqpurge")
    resources:
        mem_mb=get_mem_mb("moa_seqpurge"),
        runtime=get_runtime("moa_seqpurge"),
    conda:
        "../envs/seqpurge.yaml"
    log:
        "logs/seqpurge/moa/{library_id}.log"
    params:
        min_len=config["moa"]["seqpurge"]["min_len"],
        qcut=config["moa"]["seqpurge"]["qcut"],
    shell:
        """
        mkdir -p $(dirname {output.r1:q})
        mkdir -p $(dirname {output.summary:q})
        mkdir -p $(dirname {log:q})
        if [ {threads} -lt 4 ]; then
            echo "Error: SeqPurge requires at least 4 total threads because it may use up to three additional I/O threads." > {log:q}
            exit 1
        fi
        seqpurge_threads=$(({threads} - 3))
        SeqPurge \
          -in1 {input.r1:q} \
          -in2 {input.r2:q} \
          -out1 {output.r1:q} \
          -out2 {output.r2:q} \
          -min_len {params.min_len} \
          -qcut {params.qcut} \
          -threads "$seqpurge_threads" \
          -summary {output.summary:q} \
          > {log:q} 2>&1
        """


rule moa_ngmerge:
    input:
        r1=rules.moa_seqpurge.output.r1,
        r2=rules.moa_seqpurge.output.r2,
    output:
        merged="results/preprocessing/moa/{library_id}/{library_id}_merged.fq.gz",
    threads: get_threads("moa_ngmerge")
    resources:
        mem_mb=get_mem_mb("moa_ngmerge"),
        runtime=get_runtime("moa_ngmerge"),
    conda:
        "../envs/ngmerge.yaml"
    log:
        "logs/ngmerge/moa/{library_id}.log"
    params:
        mismatch_fraction=config["moa"]["ngmerge"]["mismatch_fraction"],
        min_overlap=config["moa"]["ngmerge"]["min_overlap"],
        dovetail_args=(
            f"-d -e {config['moa']['ngmerge']['min_dovetail_overlap']}"
            if config["moa"]["ngmerge"]["dovetail"]
            else ""
        ),
    shell:
        """
        mkdir -p $(dirname {output.merged:q})
        mkdir -p $(dirname {log:q})
        NGmerge \
          -1 {input.r1:q} \
          -2 {input.r2:q} \
          -o {output.merged:q} \
          -p {params.mismatch_fraction} \
          -m {params.min_overlap} \
          {params.dovetail_args} \
          -n {threads} \
          -z \
          -v \
          > {log:q} 2>&1
        """


rule moa_merged_fastqc:
    input:
        merged=rules.moa_ngmerge.output.merged,
    output:
        directory("results/qc/fastqc/merged/moa/{library_id}")
    threads: get_threads("moa_merged_fastqc")
    resources:
        mem_mb=get_mem_mb("moa_merged_fastqc"),
        runtime=get_runtime("moa_merged_fastqc"),
    conda:
        "../envs/fastqc.yaml"
    log:
        "logs/fastqc/merged/moa/{library_id}.log"
    shell:
        """
        mkdir -p {output:q}
        mkdir -p $(dirname {log:q})
        fastqc \
          --threads {threads} \
          --outdir {output:q} \
          {input.merged:q} \
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
