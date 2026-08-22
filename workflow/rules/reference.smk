rule build_concatenated_reference:
    input:
        parent1_fasta=lambda wildcards: require_config_value(
            config["parents"]["parent1"]["fasta"],
            "parents.parent1.fasta",
        ),
        parent2_fasta=lambda wildcards: require_config_value(
            config["parents"]["parent2"]["fasta"],
            "parents.parent2.fasta",
        ),
        script=BUILD_CONCATENATED_REFERENCE_SCRIPT,
    output:
        fasta="resources/references/{reference_id}/concatenated.fa",
        fai="resources/references/{reference_id}/concatenated.fa.fai",
        report="resources/references/{reference_id}/reference_build_report.tsv",
    params:
        parent1_prefix=lambda wildcards: require_config_value(
            config["parents"]["parent1"]["prefix"],
            "parents.parent1.prefix",
        ),
        parent2_prefix=lambda wildcards: require_config_value(
            config["parents"]["parent2"]["prefix"],
            "parents.parent2.prefix",
        ),
        reference_id=lambda wildcards: wildcards.reference_id,
        outdir=lambda wildcards: (
            f"resources/references/{wildcards.reference_id}"
        ),
    threads: get_threads("build_concatenated_reference")
    resources:
        mem_mb=get_mem_mb("build_concatenated_reference"),
        runtime=get_runtime("build_concatenated_reference"),
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/reference/{reference_id}.log"
    shell:
        """
        mkdir -p $(dirname {log:q})
        PARENT1_FA={input.parent1_fasta:q} \
        PARENT2_FA={input.parent2_fasta:q} \
        PARENT1_PREFIX={params.parent1_prefix:q} \
        PARENT2_PREFIX={params.parent2_prefix:q} \
        REFERENCE_ID={params.reference_id:q} \
        OUTDIR={params.outdir:q} \
        bash {input.script:q} > {log:q} 2>&1
        """


rule build_concatenated_gtf:
    input:
        parent1_gtf=lambda wildcards: require_config_value(
            config["parents"]["parent1"]["gtf"],
            "parents.parent1.gtf",
        ),
        parent2_gtf=lambda wildcards: require_config_value(
            config["parents"]["parent2"]["gtf"],
            "parents.parent2.gtf",
        ),
        script=BUILD_CONCATENATED_GTF_SCRIPT,
    output:
        gtf="resources/references/{reference_id}/concatenated.gtf"
    params:
        parent1_prefix=lambda wildcards: require_config_value(
            config["parents"]["parent1"]["prefix"],
            "parents.parent1.prefix",
        ),
        parent2_prefix=lambda wildcards: require_config_value(
            config["parents"]["parent2"]["prefix"],
            "parents.parent2.prefix",
        ),
        reference_id=lambda wildcards: wildcards.reference_id,
        outdir=lambda wildcards: (
            f"resources/references/{wildcards.reference_id}"
        ),
    threads: get_threads("build_concatenated_gtf")
    resources:
        mem_mb=get_mem_mb("build_concatenated_gtf"),
        runtime=get_runtime("build_concatenated_gtf"),
    log:
        "logs/reference/{reference_id}.gtf.log"
    shell:
        """
        mkdir -p $(dirname {log:q})
        PARENT1_GTF={input.parent1_gtf:q} \
        PARENT2_GTF={input.parent2_gtf:q} \
        PARENT1_PREFIX={params.parent1_prefix:q} \
        PARENT2_PREFIX={params.parent2_prefix:q} \
        REFERENCE_ID={params.reference_id:q} \
        OUTDIR={params.outdir:q} \
        bash {input.script:q} > {log:q} 2>&1
        """


rule build_star_index:
    input:
        fasta=rules.build_concatenated_reference.output.fasta,
        fai=rules.build_concatenated_reference.output.fai,
        gtf=lambda wildcards: (
            [
                f"resources/references/{wildcards.reference_id}/concatenated.gtf"
            ]
            if RNA_LIBRARIES
            else []
        ),
    output:
        index=directory(
            "resources/references/{reference_id}/star_index"
        ),
        report="resources/references/{reference_id}/star_index_build_report.tsv",
    params:
        reference_id=lambda wildcards: wildcards.reference_id,
        sjdb_args=lambda wildcards: (
            (
                "--sjdbGTFfile "
                f"resources/references/{wildcards.reference_id}/concatenated.gtf "
                "--sjdbOverhang 149"
            )
            if RNA_LIBRARIES
            else ""
        ),
        ref_gtf=lambda wildcards: (
            f"resources/references/{wildcards.reference_id}/concatenated.gtf"
            if RNA_LIBRARIES
            else "NA"
        ),
        sjdb_overhang="149" if RNA_LIBRARIES else "NA",
    threads: get_threads("build_star_index")
    resources:
        mem_mb=get_mem_mb("build_star_index"),
        runtime=get_runtime("build_star_index"),
    conda:
        "../envs/star.yaml"
    log:
        "logs/reference/{reference_id}.star_index.log"
    shell:
        """
        mkdir -p $(dirname {log:q})
        {{
            mkdir -p {output.index:q}
            GENOME_LEN=$(awk 'BEGIN{{s=0}} {{s+=$2}} END{{print s}}' {input.fai:q})
            SA_NBASES=$(awk -v L="$GENOME_LEN" 'BEGIN{{
                v = int((log(L)/log(2))/2 - 1);
                if (v > 14) v = 14;
                if (v < 4) v = 4;
                print v;
            }}')
            STAR \
                --runMode genomeGenerate \
                --runThreadN {threads} \
                --genomeDir {output.index:q} \
                --genomeFastaFiles {input.fasta:q} \
                {params.sjdb_args} \
                --genomeSAindexNbases "$SA_NBASES"
            printf "metric\tvalue\n" > {output.report:q}
            printf "reference_id\t%s\n" "{params.reference_id}" >> {output.report:q}
            printf "ref_fasta\t%s\n" {input.fasta:q} >> {output.report:q}
            printf "ref_gtf\t%s\n" {params.ref_gtf:q} >> {output.report:q}
            printf "index_dir\t%s\n" {output.index:q} >> {output.report:q}
            printf "genome_length\t%s\n" "$GENOME_LEN" >> {output.report:q}
            printf "genomeSAindexNbases\t%s\n" "$SA_NBASES" >> {output.report:q}
            printf "sjdbOverhang\t%s\n" {params.sjdb_overhang:q} >> {output.report:q}
            printf "threads\t%s\n" "{threads}" >> {output.report:q}
            printf "star_bin\tSTAR\n" >> {output.report:q}
        }} > {log:q} 2>&1
        """


rule validate_rna_gene_pairs:
    input:
        gene_pairs=lambda wildcards: require_config_value(
            config["rna_gene_pairs"],
            "rna_gene_pairs",
        ),
        script=VALIDATE_PAIRED_FEATURES_SCRIPT,
    output:
        report="results/qc/rna/gene_pairs_validation.tsv"
    threads: get_threads("validate_rna_gene_pairs")
    resources:
        mem_mb=get_mem_mb("validate_rna_gene_pairs"),
        runtime=get_runtime("validate_rna_gene_pairs"),
    log:
        "logs/validation/rna_gene_pairs.log"
    shell:
        """
        mkdir -p $(dirname {log:q})
        python {input.script:q} \
          --gene-pairs {input.gene_pairs:q} \
          --report {output.report:q} \
          > {log:q} 2>&1
        """


rule validate_moa_peak_pairs:
    input:
        peak_pairs=lambda wildcards: require_config_value(
            config["moa_peak_pairs"],
            "moa_peak_pairs",
        ),
        script=VALIDATE_PAIRED_FEATURES_SCRIPT,
    output:
        report="results/qc/moa/peak_pairs_validation.tsv"
    threads: get_threads("validate_moa_peak_pairs")
    resources:
        mem_mb=get_mem_mb("validate_moa_peak_pairs"),
        runtime=get_runtime("validate_moa_peak_pairs"),
    log:
        "logs/validation/moa_peak_pairs.log"
    shell:
        """
        mkdir -p $(dirname {log:q})
        python {input.script:q} \
          --peak-pairs {input.peak_pairs:q} \
          --report {output.report:q} \
          > {log:q} 2>&1
        """
