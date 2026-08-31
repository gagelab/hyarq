# ATAC-seq mapping rules.

rule map_atac_star:
    input:
        r1=rules.atac_fastp.output.r1,
        r2=rules.atac_fastp.output.r2,
        index=lambda wildcards: configured_reference_path(
            "star_index"
        ),
        script=MAP_ATAC_STAR_SCRIPT,
    output:
        aligned_bam="results/mapping/atac/{library_id}.Aligned.sortedByCoord.out.bam",
        unique_bam="results/mapping/atac/{library_id}.unique.bam",
        unique_bai="results/mapping/atac/{library_id}.unique.bam.bai",
        star_log="results/mapping/atac/{library_id}.Log.final.out",
        summary="results/mapping/atac/{library_id}.mapping_summary.tsv",
    params:
        sample=lambda wildcards: SAMPLES.loc[wildcards.library_id, "sample_id"],
        replicate=lambda wildcards: SAMPLES.loc[wildcards.library_id, "replicate"],
        reference_id=lambda wildcards: require_config_value(
            config["reference_id"],
            "reference_id",
        ),
        outdir="results/mapping/atac",
        out_sam_mult_nmax=lambda wildcards: (
            config["atac"]["star"]["out_sam_mult_nmax"]
        ),
        win_anchor_multimap_nmax=lambda wildcards: (
            config["atac"]["star"]["win_anchor_multimap_nmax"]
        ),
        out_bam_sorting_bins_n=lambda wildcards: (
            config["atac"]["star"]["out_bam_sorting_bins_n"]
        ),
    threads: get_threads("map_atac_star")
    resources:
        mem_mb=get_mem_mb("map_atac_star"),
        runtime=get_runtime("map_atac_star"),
    conda:
        "../envs/star.yaml"
    log:
        "logs/mapping/atac/{library_id}.log"
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
        OUT_SAM_MULT_NMAX={params.out_sam_mult_nmax} \
        WIN_ANCHOR_MULTIMAP_NMAX={params.win_anchor_multimap_nmax} \
        OUT_BAM_SORTING_BINS_N={params.out_bam_sorting_bins_n} \
        bash {input.script:q} > {log:q} 2>&1
        """


rule deduplicate_atac:
    input:
        bam=rules.map_atac_star.output.unique_bam,
        script=DEDUPLICATE_ATAC_SCRIPT,
    output:
        bam="results/mapping/atac/{library_id}.dedup.bam",
        bai="results/mapping/atac/{library_id}.dedup.bam.bai",
        stats="results/qc/samtools/atac/{library_id}.markdup.txt",
    threads: get_threads("deduplicate_atac")
    resources:
        mem_mb=get_mem_mb("deduplicate_atac"),
        runtime=get_runtime("deduplicate_atac"),
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/deduplication/atac/{library_id}.log"
    shell:
        r"""
        mkdir -p "$(dirname {log:q})"
        BAM={input.bam:q} \
        OUTPUT_BAM={output.bam:q} \
        STATS={output.stats:q} \
        TMP_ROOT={resources.tmpdir:q} \
        THREADS={threads} \
        LIBRARY_ID={wildcards.library_id:q} \
        bash {input.script:q} > {log:q} 2>&1
        """


rule count_atac_peaks:
    input:
        bam=rules.deduplicate_atac.output.bam,
        peak_pairs=lambda wildcards: require_config_value(
            config["atac_peak_pairs"],
            "atac_peak_pairs",
        ),
        validation=rules.validate_atac_peak_pairs.output,
        script=COUNT_PAIRED_REGIONS_SCRIPT,
    output:
        peak_counts="results/counts/atac/{library_id}/peak_counts.tsv",
        filter_summary="results/counts/atac/{library_id}/filter_summary.tsv",
    params:
        sample=lambda wildcards: SAMPLES.loc[wildcards.library_id, "sample_id"],
        replicate=lambda wildcards: SAMPLES.loc[wildcards.library_id, "replicate"],
        parent1_prefix=config["parents"]["parent1"]["prefix"],
        parent2_prefix=config["parents"]["parent2"]["prefix"],
    threads: get_threads("count_atac_peaks")
    resources:
        mem_mb=get_mem_mb("count_atac_peaks"),
        runtime=get_runtime("count_atac_peaks"),
    conda:
        "../envs/rna_counting.yaml"
    log:
        "logs/counts/atac/{library_id}.log"
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
            --assay atac \
            --sample {params.sample:q} \
            --rep {params.replicate:q} \
            --library-id {wildcards.library_id:q} \
            --threads {threads} \
            --paired-end \
            > {log:q} 2>&1
        """


rule aggregate_atac_peak_counts:
    input:
        count_tables=expand(
            "results/counts/atac/{library_id}/peak_counts.tsv",
            library_id=ATAC_LIBRARIES,
        ),
        script=AGGREGATE_PAIRED_REGION_COUNTS_SCRIPT,
    output:
        matrix="results/counts/atac/peak_pair_count_matrix.tsv",
    threads: 1
    resources:
        mem_mb=get_mem_mb("aggregate_atac_peak_counts"),
        runtime=get_runtime("aggregate_atac_peak_counts"),
    conda:
        "../envs/rna_gene_qc.yaml"
    log:
        "logs/counting/atac/aggregate_peak_counts.log"
    shell:
        r"""
        mkdir -p "$(dirname {log:q})"
        python {input.script:q} \
            --count-tables {input.count_tables:q} \
            --output {output.matrix:q} \
            > {log:q} 2>&1
        """
