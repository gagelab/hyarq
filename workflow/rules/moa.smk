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
