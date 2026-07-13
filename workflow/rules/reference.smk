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


rule build_star_index:
    input:
        fasta=rules.build_concatenated_reference.output.fasta,
        gtf=rules.build_concatenated_reference.output.gtf,
        fai=rules.build_concatenated_reference.output.fai,
    output:
        index=directory(f"{REFERENCE_OUTDIR}/star_index"),
        report=f"{REFERENCE_OUTDIR}/star_index_build_report.tsv",
    threads: 4
    conda:
        "../envs/star.yaml"
    log:
        f"logs/reference/{REFERENCE_ID}.star_index.log"
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
                --sjdbGTFfile {input.gtf:q} \
                --sjdbOverhang 149 \
                --genomeSAindexNbases "$SA_NBASES"
            printf "metric\tvalue\n" > {output.report:q}
            printf "reference_id\t%s\n" "{REFERENCE_ID}" >> {output.report:q}
            printf "ref_fasta\t%s\n" {input.fasta:q} >> {output.report:q}
            printf "ref_gtf\t%s\n" {input.gtf:q} >> {output.report:q}
            printf "index_dir\t%s\n" {output.index:q} >> {output.report:q}
            printf "genome_length\t%s\n" "$GENOME_LEN" >> {output.report:q}
            printf "genomeSAindexNbases\t%s\n" "$SA_NBASES" >> {output.report:q}
            printf "sjdbOverhang\t149\n" >> {output.report:q}
            printf "threads\t%s\n" "{threads}" >> {output.report:q}
            printf "star_bin\tSTAR\n" >> {output.report:q}
        }} > {log:q} 2>&1
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
