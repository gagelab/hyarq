REFERENCE_ID = config["reference_id"]
REFERENCE_OUTDIR = f"resources/references/{REFERENCE_ID}"


rule build_concatenated_reference:
    input:
        parent1_fasta=config["parents"]["parent1"]["fasta"],
        parent2_fasta=config["parents"]["parent2"]["fasta"],
        parent1_gtf=config["parents"]["parent1"]["gtf"],
        parent2_gtf=config["parents"]["parent2"]["gtf"],
    output:
        fasta=f"{REFERENCE_OUTDIR}/concatenated.fa",
        gtf=f"{REFERENCE_OUTDIR}/concatenated.gtf",
        fai=f"{REFERENCE_OUTDIR}/concatenated.fa.fai",
        report=f"{REFERENCE_OUTDIR}/reference_build_report.tsv",
    params:
        parent1_prefix=config["parents"]["parent1"]["prefix"],
        parent2_prefix=config["parents"]["parent2"]["prefix"],
        reference_id=REFERENCE_ID,
        outdir=REFERENCE_OUTDIR,
    conda:
        "../envs/samtools.yaml"
    log:
        f"logs/reference/{REFERENCE_ID}.log"
    shell:
        """
        mkdir -p $(dirname {log:q})
        PARENT1_FA={input.parent1_fasta:q} \
        PARENT2_FA={input.parent2_fasta:q} \
        PARENT1_GTF={input.parent1_gtf:q} \
        PARENT2_GTF={input.parent2_gtf:q} \
        PARENT1_PREFIX={params.parent1_prefix:q} \
        PARENT2_PREFIX={params.parent2_prefix:q} \
        REFERENCE_ID={params.reference_id:q} \
        OUTDIR={params.outdir:q} \
        bash workflow/scripts/build_concatenated_reference.sh > {log:q} 2>&1
        """


rule validate_rna_gene_pairs:
    input:
        gene_pairs=config["rna_gene_pairs"]
    output:
        report="results/qc/rna/gene_pairs_validation.tsv"
    log:
        "logs/validation/rna_gene_pairs.log"
    shell:
        """
        mkdir -p $(dirname {log:q})
        python workflow/scripts/validate_paired_features.py \
          --gene-pairs {input.gene_pairs:q} \
          --report {output.report:q} \
          > {log:q} 2>&1
        """
