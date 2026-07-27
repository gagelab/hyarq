# RNA-seq mapping and gene-level allele-specific counting rules.

rule map_rna_star:
    input:
        r1=rules.rna_fastp.output.r1,
        r2=rules.rna_fastp.output.r2,
        index=rules.build_star_index.output.index,
        script=MAP_RNA_STAR_SCRIPT,
    output:
        bam="results/mapping/rna/{library_id}.Aligned.sortedByCoord.out.bam",
        bai="results/mapping/rna/{library_id}.Aligned.sortedByCoord.out.bam.bai",
        gene_counts="results/mapping/rna/{library_id}.ReadsPerGene.out.tab",
        star_log="results/mapping/rna/{library_id}.Log.final.out",
        summary="results/mapping/rna/{library_id}.mapping_summary.tsv",
    params:
        sample=lambda wildcards: SAMPLES.loc[wildcards.library_id, "sample_id"],
        replicate=lambda wildcards: SAMPLES.loc[wildcards.library_id, "replicate"],
        reference_id=config["reference_id"],
        outdir="results/mapping/rna",
    threads: get_threads("map_rna_star")
    resources:
        mem_mb=get_mem_mb("map_rna_star"),
        runtime=get_runtime("map_rna_star"),
    conda:
        "../envs/star.yaml"
    log:
        "logs/mapping/rna/{library_id}.log"
    shell:
        r"""
        mkdir -p "$(dirname {log:q})"
        R1={input.r1:q} \
        R2={input.r2:q} \
        REFERENCE_ID={params.reference_id:q} \
        INDEX_DIR={input.index:q} \
        OUTDIR={params.outdir:q} \
        THREADS={threads} \
        SAMPLE={params.sample:q} \
        REP={params.replicate:q} \
        LIBRARY_ID={wildcards.library_id:q} \
        bash {input.script:q} > {log:q} 2>&1
        """


rule count_rna_genes:
    input:
        bam=rules.map_rna_star.output.bam,
        gtf=rules.build_concatenated_reference.output.gtf,
        gene_pairs=config["rna_gene_pairs"],
        validation=rules.validate_rna_gene_pairs.output.report,
        script=COUNT_RNA_GENES_SCRIPT,
    output:
        gene_counts="results/counts/rna/{library_id}/gene_counts.tsv",
        filter_summary="results/counts/rna/{library_id}/filter_summary.tsv",
    params:
        sample=lambda wildcards: SAMPLES.loc[wildcards.library_id, "sample_id"],
        replicate=lambda wildcards: SAMPLES.loc[wildcards.library_id, "replicate"],
        parent1_prefix=config["parents"]["parent1"]["prefix"],
        parent2_prefix=config["parents"]["parent2"]["prefix"],
        outdir="results/counts/rna/{library_id}",
        count_tmpdir="results/counts/rna/{library_id}/tmp",
    threads: get_threads("count_rna_genes")
    resources:
        mem_mb=get_mem_mb("count_rna_genes"),
        runtime=get_runtime("count_rna_genes"),
    conda:
        "../envs/rna_counting.yaml"
    log:
        "logs/counting/rna/{library_id}.log"
    shell:
        r"""
        mkdir -p "$(dirname {log:q})"
        BAM={input.bam:q} \
        GTF={input.gtf:q} \
        GENE_PAIRS={input.gene_pairs:q} \
        PARENT1_PREFIX={params.parent1_prefix:q} \
        PARENT2_PREFIX={params.parent2_prefix:q} \
        OUTDIR={params.outdir:q} \
        COUNT_TMPDIR={params.count_tmpdir:q} \
        THREADS={threads} \
        SAMPLE={params.sample:q} \
        REP={params.replicate:q} \
        LIBRARY_ID={wildcards.library_id:q} \
        bash {input.script:q} > {log:q} 2>&1
        """
