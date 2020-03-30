gtf=/hpc/hub_oudenaarden/vbhardwaj/annotations/hg38_gencode31_extra/annotation/genes.gtf
idx=../annotation/STARindex/
barcodes=../whitelist.txt
r1=../fastq_reformatted/MAV-RPFv4-HEK293T_R1.fastq.gz
r2=../fastq_reformatted/MAV-RPFv4-HEK293T_R2.fastq.gz
prefix=RPFv4_
fa=

## generate genome
rule maskFasta:
    input:
        fa = genome_fasta,
        bed = mask_bed
    output: temp("annotation/genome_masked.fa")
    threads: 1
    conda: CONDA_SHARED_ENV
    shell:
        "bedtools maskfasta -fi {input.fa} -fo {output} -bed {input.bed}"

rule appendFasta:
    input:
        fa = lambda wildcards: "annotation/genome_tRNA_masked.fa" if maskFasta else genome_fasta,
        addFa = append_fasta
    output: "annotation/genome.fa"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell:
        "cat {input.fa} {input.addfa} > {output}"


### prepare annotation
rule gtfTable:
    input: genome_gtf
    output: "annotation/gtf_annotation_table.txt"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell:
        """
        cat {input} | tr -d ";\\""| awk '$3=="gene" {{ print $10, $14, $1, $12 }}' | \
        sed "1i\Geneid\tGeneSymbol\tChrom\tClass" > {output}
        """

rule gtfTable:
    input: genome_gtf
    output:
        bed = "annotation/selected_CDS.bed",
        annotation = "annotation/selected_CDS_annotation.bed"
    params:
        rscript = os.path.join(workflow.basedir, "tools", "prep_annotation.R"),
        outPrefix = "annotation/selected_CDS"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell:
        "Rscript {params.rscript} {input} {params.outPrefix} 2> {log} 2>&1"


## prepare indicies
rule STARindex:
    input:
        fa = "annotation/genome.fa",
        gtf = genome_gtf
    output:
        gidx = "annotation/STARindex/Genome",
        sidx = "annotation/STARindex/SAindex"
    params:
        outdir = "annotation/STARindex"
    threads: 10
    conda: CONDA_SHARED_ENV
    shell:
        "STAR --runThreadN {threads} \
        --runMode genomeGenerate \
        --genomeDir {params.outdir} \
        --genomeFastaFiles {input.fa} \
        --sjdbGTFfile {input.gtf} \
        --sjdbOverhang 50 2> {log} 2>&1"

rule Bowtie2index:
    input: "annotation/selected_CDS_extended.fa"
    output:
        fwd = "annotation/Bowtie2index/selected_CDS_extended.2.bt2",
        rev = "annotation/Bowtie2index/selected_CDS_extended.rev.2.bt2"
    params:
        outdir = "annotation/Bowtie2index"
    threads: 5
    conda: CONDA_SHARED_ENV
    shell:
        "bowtie2-build --seed 2020 --threads {threads} \
        {input} {params.outdir} 2> {log} 2>&1"
