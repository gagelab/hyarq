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
