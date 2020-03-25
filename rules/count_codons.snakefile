rule count_regions:
    input:
        bam = "Bowtie2_CDS/{sample}.bam",
        idx = "Bowtie2_CDS/{sample}.bam.bai",
        bed = annotation+"/selected_CDS_51b_cdsAsGenes.bed"
    output:
        counts = "counts/{sample}.CDScounts_bulk.tsv",
        bam = temp("counts/{sample}.bam.featureCounts.bam")
    params:
        filetype = "-F SAF"
    log: "logs/featurecounts_{sample}_bulk.err"
    threads: 10
    conda: CONDA_SHARED_ENV
    shell:
        "featureCounts -T {threads} -a {input.regions} {params.filetype} \
        -R BAM -s 1 --read2pos 3 -o {output.counts} {input.bam} > {log} 2>&1"

rule fcount_sort:
    input: "counts/{sample}.bam.featureCounts.bam"
    output: temp("counts/{sample}.featureCounts.bam")
    threads: 10
    conda: CONDA_SHARED_ENV
    shell:
        "samtools sort -@ {threads} -O BAM -o {output} {input}"

rule fcount_index:
    input: "counts/{sample}.featureCounts.bam"
    output: temp("counts/{sample}.featureCounts.bam.bai")
    threads: 1
    conda: CONDA_SHARED_ENV
    shell:
        "samtools index {input}"

rule count_regions_cells:
    input:
        bam = "counts/{sample}.featureCounts.bam",
        bai = "counts/{sample}.featureCounts.bam.bai"
    output: "counts/{sample}.CDScounts_per_barcode.tsv"
    log: "logs/umi_counts_{sample}.CDScounts.err"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell:
        "umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XS \
        --per-cell --method=percentile -I {input.bam} -S {output} \
        -v 4 --log2stderr --log={log}"
