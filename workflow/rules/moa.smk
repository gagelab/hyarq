# MOA-seq mapping and peak-level allele-specific counting rules.

rule map_moa_star:
    input:
        merged=rules.moa_ngmerge.output.merged,
        index=rules.build_star_index.output.index,
        script=MAP_MOA_STAR_SCRIPT,
    output:
        aligned_bam="results/mapping/moa/{library_id}.Aligned.sortedByCoord.out.bam",
        unique_bam="results/mapping/moa/{library_id}.unique.bam",
        unique_bai="results/mapping/moa/{library_id}.unique.bam.bai",
        star_log="results/mapping/moa/{library_id}.Log.final.out",
        summary="results/mapping/moa/{library_id}.mapping_summary.tsv",
    params:
        sample=lambda wildcards: SAMPLES.loc[wildcards.library_id, "sample_id"],
        replicate=lambda wildcards: SAMPLES.loc[wildcards.library_id, "replicate"],
        reference_id=config["reference_id"],
        outdir="results/mapping/moa",
        align_intron_max=config["moa"]["star"]["align_intron_max"],
        out_sam_mult_nmax=config["moa"]["star"]["out_sam_mult_nmax"],
        win_anchor_multimap_nmax=config["moa"]["star"]["win_anchor_multimap_nmax"],
    threads: get_threads("map_moa_star")
    resources:
        mem_mb=get_mem_mb("map_moa_star"),
        runtime=get_runtime("map_moa_star"),
    conda:
        "../envs/star.yaml"
    log:
        "logs/mapping/moa/{library_id}.log"
    shell:
        r"""
        mkdir -p "$(dirname {log:q})"
        MERGED={input.merged:q} \
        REFERENCE_ID={params.reference_id:q} \
        INDEX_DIR={input.index:q} \
        OUTDIR={params.outdir:q} \
        THREADS={threads} \
        SAMPLE={params.sample:q} \
        REP={params.replicate:q} \
        LIBRARY_ID={wildcards.library_id:q} \
        ALIGN_INTRON_MAX={params.align_intron_max} \
        OUT_SAM_MULT_NMAX={params.out_sam_mult_nmax} \
        WIN_ANCHOR_MULTIMAP_NMAX={params.win_anchor_multimap_nmax} \
        bash {input.script:q} > {log:q} 2>&1
        """


rule count_moa_peaks:
    input:
        bam=rules.map_moa_star.output.unique_bam,
        peak_pairs=config["moa_peak_pairs"],
        validation=rules.validate_moa_peak_pairs.output,
        script=COUNT_PAIRED_REGIONS_SCRIPT,
    output:
        peak_counts="results/counts/moa/{library_id}/peak_counts.tsv",
        filter_summary="results/counts/moa/{library_id}/filter_summary.tsv",
    params:
        sample=lambda wildcards: SAMPLES.loc[wildcards.library_id, "sample_id"],
        replicate=lambda wildcards: SAMPLES.loc[wildcards.library_id, "replicate"],
        parent1_prefix=config["parents"]["parent1"]["prefix"],
        parent2_prefix=config["parents"]["parent2"]["prefix"],
    threads: get_threads("count_moa_peaks")
    resources:
        mem_mb=get_mem_mb("count_moa_peaks"),
        runtime=get_runtime("count_moa_peaks"),
    conda:
        "../envs/rna_counting.yaml"
    log:
        "logs/counts/moa/{library_id}.log"
    shell:
        r"""
        mkdir -p "$(dirname {log:q})"
        python {input.script:q} \
            --bam {input.bam:q} \
            --peak-pairs {input.peak_pairs:q} \
            --parent1-prefix {params.parent1_prefix:q} \
            --parent2-prefix {params.parent2_prefix:q} \
            --output {output.peak_counts:q} \
            --summary {output.filter_summary:q} \
            --tmpdir {resources.tmpdir:q} \
            --assay moa \
            --sample {params.sample:q} \
            --rep {params.replicate:q} \
            --library-id {wildcards.library_id:q} \
            --threads {threads} \
            > {log:q} 2>&1
        """


rule aggregate_moa_peak_counts:
    input:
        count_tables=expand(
            "results/counts/moa/{library_id}/peak_counts.tsv",
            library_id=MOA_LIBRARIES,
        ),
        script=AGGREGATE_PAIRED_REGION_COUNTS_SCRIPT,
    output:
        matrix="results/counts/moa/peak_pair_count_matrix.tsv",
    threads: 1
    resources:
        mem_mb=get_mem_mb("aggregate_moa_peak_counts"),
        runtime=get_runtime("aggregate_moa_peak_counts"),
    conda:
        "../envs/rna_gene_qc.yaml"
    log:
        "logs/counting/moa/aggregate_peak_counts.log"
    shell:
        r"""
        mkdir -p "$(dirname {log:q})"
        python {input.script:q} \
            --count-tables {input.count_tables:q} \
            --output {output.matrix:q} \
            > {log:q} 2>&1
        """
