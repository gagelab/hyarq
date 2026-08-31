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


def _atac_fastp_boolean_flag(option_name, flag):
    value = config["atac"]["fastp"][option_name]
    config_key = f"atac.fastp.{option_name}"
    if not isinstance(value, bool):
        raise ValueError(f"{config_key} must be a Boolean")
    return flag if value else ""


rule atac_fastp:
    input:
        r1=lambda wildcards: SAMPLES.loc[wildcards.library_id, "r1"],
        r2=lambda wildcards: SAMPLES.loc[wildcards.library_id, "r2"],
    output:
        r1="results/preprocessing/atac/{library_id}/{library_id}_R1.fq.gz",
        r2="results/preprocessing/atac/{library_id}/{library_id}_R2.fq.gz",
        html="results/qc/fastp/atac/{library_id}.html",
        json="results/qc/fastp/atac/{library_id}.json",
    params:
        detect_adapter_for_pe=lambda wildcards: _atac_fastp_boolean_flag(
            "detect_adapter_for_pe",
            "--detect_adapter_for_pe",
        ),
        correction=lambda wildcards: _atac_fastp_boolean_flag(
            "correction",
            "--correction",
        ),
        length_required=lambda wildcards: (
            config["atac"]["fastp"]["length_required"]
        ),
    threads: get_threads("atac_fastp")
    resources:
        mem_mb=get_mem_mb("atac_fastp"),
        runtime=get_runtime("atac_fastp"),
    conda:
        "../envs/fastp.yaml"
    log:
        "logs/fastp/atac/{library_id}.log"
    shell:
        """
        mkdir -p $(dirname {output.r1:q})
        mkdir -p $(dirname {output.html:q})
        mkdir -p $(dirname {log:q})
        fastp \
            --in1 {input.r1:q} \
            --in2 {input.r2:q} \
            {params.detect_adapter_for_pe} \
            {params.correction} \
            --length_required {params.length_required} \
            --out1 {output.r1:q} \
            --out2 {output.r2:q} \
            --html {output.html:q} \
            --json {output.json:q} \
            --thread {threads} \
            > {log:q} 2>&1
        """


def _chip_fastp_boolean_flag(option_name, flag):
    value = config["chip"]["fastp"][option_name]
    config_key = f"chip.fastp.{option_name}"
    if not isinstance(value, bool):
        raise ValueError(f"{config_key} must be a Boolean")
    return flag if value else ""


rule chip_fastp:
    input:
        r1=lambda wildcards: SAMPLES.loc[wildcards.library_id, "r1"],
        r2=lambda wildcards: SAMPLES.loc[wildcards.library_id, "r2"],
    output:
        r1="results/preprocessing/chip/{library_id}/{library_id}_R1.fq.gz",
        r2="results/preprocessing/chip/{library_id}/{library_id}_R2.fq.gz",
        html="results/qc/fastp/chip/{library_id}.html",
        json="results/qc/fastp/chip/{library_id}.json",
    params:
        detect_adapter_for_pe=lambda wildcards: _chip_fastp_boolean_flag(
            "detect_adapter_for_pe",
            "--detect_adapter_for_pe",
        ),
        correction=lambda wildcards: _chip_fastp_boolean_flag(
            "correction",
            "--correction",
        ),
        length_required=lambda wildcards: (
            config["chip"]["fastp"]["length_required"]
        ),
    threads: get_threads("chip_fastp")
    resources:
        mem_mb=get_mem_mb("chip_fastp"),
        runtime=get_runtime("chip_fastp"),
    conda:
        "../envs/fastp.yaml"
    log:
        "logs/fastp/chip/{library_id}.log"
    shell:
        """
        mkdir -p $(dirname {output.r1:q})
        mkdir -p $(dirname {output.html:q})
        mkdir -p $(dirname {log:q})
        fastp \
            --in1 {input.r1:q} \
            --in2 {input.r2:q} \
            {params.detect_adapter_for_pe} \
            {params.correction} \
            --length_required {params.length_required} \
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


# FLASH2 2.2.00 settings match the validated MOA merging behavior.
# Its native minimum-overlap default is intentionally retained.
rule moa_flash2:
    input:
        r1=rules.moa_seqpurge.output.r1,
        r2=rules.moa_seqpurge.output.r2,
    output:
        merged="results/preprocessing/moa/{library_id}/{library_id}_merged.fq.gz",
    threads: get_threads("moa_flash2")
    resources:
        mem_mb=get_mem_mb("moa_flash2"),
        runtime=get_runtime("moa_flash2"),
    conda:
        "../envs/flash2.yaml"
    log:
        "logs/flash2/moa/{library_id}.log"
    shell:
        """
        mkdir -p $(dirname {output.merged:q})
        mkdir -p $(dirname {log:q})
        flash2_tmp=$(mktemp -d {resources.tmpdir:q}/hyarq_flash2.XXXXXX)
        trap 'rm -rf -- "$flash2_tmp"' EXIT
        flash2 \
          {input.r1:q} \
          {input.r2:q} \
          -t {threads} \
          --cap-mismatch-quals \
          -M 100 \
          -z \
          -d "$flash2_tmp" \
          -o flash2 \
          > {log:q} 2>&1
        mv "$flash2_tmp/flash2.extendedFrags.fastq.gz" {output.merged:q} \
          >> {log:q} 2>&1
        rm -rf -- "$flash2_tmp"
        trap - EXIT
        """


rule moa_merged_fastqc:
    input:
        merged=rules.moa_flash2.output.merged,
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


rule atac_clean_fastqc:
    input:
        r1="results/preprocessing/atac/{library_id}/{library_id}_R1.fq.gz",
        r2="results/preprocessing/atac/{library_id}/{library_id}_R2.fq.gz",
    output:
        directory("results/qc/fastqc/clean/atac/{library_id}")
    threads: get_threads("atac_clean_fastqc")
    resources:
        mem_mb=get_mem_mb("atac_clean_fastqc"),
        runtime=get_runtime("atac_clean_fastqc"),
    conda:
        "../envs/fastqc.yaml"
    log:
        "logs/fastqc/clean/atac/{library_id}.log"
    shell:
        """
        mkdir -p {output:q}
        mkdir -p $(dirname {log:q})
        fastqc --threads {threads} --outdir {output:q} {input.r1:q} {input.r2:q} > {log:q} 2>&1
        """


rule chip_clean_fastqc:
    input:
        r1="results/preprocessing/chip/{library_id}/{library_id}_R1.fq.gz",
        r2="results/preprocessing/chip/{library_id}/{library_id}_R2.fq.gz",
    output:
        directory("results/qc/fastqc/clean/chip/{library_id}")
    threads: get_threads("chip_clean_fastqc")
    resources:
        mem_mb=get_mem_mb("chip_clean_fastqc"),
        runtime=get_runtime("chip_clean_fastqc"),
    conda:
        "../envs/fastqc.yaml"
    log:
        "logs/fastqc/clean/chip/{library_id}.log"
    shell:
        """
        mkdir -p {output:q}
        mkdir -p $(dirname {log:q})
        fastqc --threads {threads} --outdir {output:q} {input.r1:q} {input.r2:q} > {log:q} 2>&1
        """
