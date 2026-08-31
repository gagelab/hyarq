# ChIP-seq mapping rules.

rule map_chip_star:
    input:
        r1=rules.chip_fastp.output.r1,
        r2=rules.chip_fastp.output.r2,
        index=lambda wildcards: configured_reference_path(
            "star_index"
        ),
        script=MAP_CHIP_STAR_SCRIPT,
    output:
        aligned_bam="results/mapping/chip/{library_id}.Aligned.sortedByCoord.out.bam",
        unique_bam="results/mapping/chip/{library_id}.unique.bam",
        unique_bai="results/mapping/chip/{library_id}.unique.bam.bai",
        star_log="results/mapping/chip/{library_id}.Log.final.out",
        summary="results/mapping/chip/{library_id}.mapping_summary.tsv",
    params:
        sample=lambda wildcards: SAMPLES.loc[wildcards.library_id, "sample_id"],
        replicate=lambda wildcards: SAMPLES.loc[wildcards.library_id, "replicate"],
        reference_id=lambda wildcards: require_config_value(
            config["reference_id"],
            "reference_id",
        ),
        outdir="results/mapping/chip",
        out_sam_mult_nmax=lambda wildcards: (
            config["chip"]["star"]["out_sam_mult_nmax"]
        ),
        win_anchor_multimap_nmax=lambda wildcards: (
            config["chip"]["star"]["win_anchor_multimap_nmax"]
        ),
        out_bam_sorting_bins_n=lambda wildcards: (
            config["chip"]["star"]["out_bam_sorting_bins_n"]
        ),
    threads: get_threads("map_chip_star")
    resources:
        mem_mb=get_mem_mb("map_chip_star"),
        runtime=get_runtime("map_chip_star"),
    conda:
        "../envs/star.yaml"
    log:
        "logs/mapping/chip/{library_id}.log"
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


rule deduplicate_chip:
    input:
        bam=rules.map_chip_star.output.unique_bam,
        script=DEDUPLICATE_CHIP_SCRIPT,
    output:
        bam="results/mapping/chip/{library_id}.dedup.bam",
        bai="results/mapping/chip/{library_id}.dedup.bam.bai",
        stats="results/qc/samtools/chip/{library_id}.markdup.txt",
    threads: get_threads("deduplicate_chip")
    resources:
        mem_mb=get_mem_mb("deduplicate_chip"),
        runtime=get_runtime("deduplicate_chip"),
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/deduplication/chip/{library_id}.log"
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


rule count_chip_peaks:
    input:
        bam=rules.deduplicate_chip.output.bam,
        peak_pairs=lambda wildcards: require_config_value(
            config["chip_peak_pairs"],
            "chip_peak_pairs",
        ),
        validation=rules.validate_chip_peak_pairs.output,
        script=COUNT_PAIRED_REGIONS_SCRIPT,
    output:
        peak_counts="results/counts/chip/{library_id}/peak_counts.tsv",
        filter_summary="results/counts/chip/{library_id}/filter_summary.tsv",
    params:
        sample=lambda wildcards: SAMPLES.loc[wildcards.library_id, "sample_id"],
        replicate=lambda wildcards: SAMPLES.loc[wildcards.library_id, "replicate"],
        parent1_prefix=config["parents"]["parent1"]["prefix"],
        parent2_prefix=config["parents"]["parent2"]["prefix"],
    threads: get_threads("count_chip_peaks")
    resources:
        mem_mb=get_mem_mb("count_chip_peaks"),
        runtime=get_runtime("count_chip_peaks"),
    conda:
        "../envs/rna_counting.yaml"
    log:
        "logs/counts/chip/{library_id}.log"
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
            --assay chip \
            --sample {params.sample:q} \
            --rep {params.replicate:q} \
            --library-id {wildcards.library_id:q} \
            --threads {threads} \
            --paired-end \
            > {log:q} 2>&1
        """


rule aggregate_chip_peak_counts:
    input:
        count_tables=expand(
            "results/counts/chip/{library_id}/peak_counts.tsv",
            library_id=CHIP_LIBRARIES,
        ),
        script=AGGREGATE_PAIRED_REGION_COUNTS_SCRIPT,
    output:
        matrix="results/counts/chip/peak_pair_count_matrix.tsv",
    threads: 1
    resources:
        mem_mb=get_mem_mb("aggregate_chip_peak_counts"),
        runtime=get_runtime("aggregate_chip_peak_counts"),
    conda:
        "../envs/rna_gene_qc.yaml"
    log:
        "logs/counting/chip/aggregate_peak_counts.log"
    shell:
        r"""
        mkdir -p "$(dirname {log:q})"
        python {input.script:q} \
            --count-tables {input.count_tables:q} \
            --output {output.matrix:q} \
            > {log:q} 2>&1
        """
