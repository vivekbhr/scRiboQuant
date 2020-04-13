## generate genome
if mask_bed is not None and append_fasta is not None:
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
            fa = "annotation/genome_masked.fa",
            addfa = append_fasta
        output: "annotation/genome.fa"
        threads: 1
        conda: CONDA_SHARED_ENV
        shell:
            "cat {input.fa} {input.addfa} > {output}"
else:
    rule linkFasta:
        input:
            fa = genome_fasta
        output: "annotation/genome.fa"
        threads: 1
        conda: CONDA_SHARED_ENV
        shell:
            "ln -s {input} {output}"

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

rule filterGTF:
    input: genome_gtf
    output: temp("annotation/transcripts_filtered.gtf")
    threads: 1
    conda: CONDA_SHARED_ENV
    shell:
        """
        awk '$2 == "HAVANA" {{ print $0 }}' {input} | \
        grep 'appris_principal' | grep 'transcript_type "protein_coding"' > {output}
        """

rule prepCDSbed:
    input: "annotation/transcripts_filtered.gtf"
    output:
        bed = "annotation/selected_CDS.bed",
        annot = "annotation/selected_CDS_annotation.bed",
        exons = "annotation/selected_CDS_exons.bed"
    params:
        rscript = os.path.join(workflow.basedir, "tools", "prepareCDS.R")
    log: "logs/prepCDSbed.log"
    threads: 1
    conda: CONDA_R_ENV
    shell:
        "Rscript {params.rscript} {input} {output.bed} {output.annot} {output.exons} > {log} 2>&1"

rule prepCDSfasta:
    input:
        bed = "annotation/selected_CDS.bed",
        fasta = "annotation/genome.fa"
    output: "annotation/selected_CDS_extended.fa"
    params:
        extend = 51 #CDSextendLength
    log: "logs/prepCDSfasta.log"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell:
        "gffread -P -J -M -W -E -T --w-add {params.extend} -w {output} \
        -g {input.fasta} --in-bed {input.bed} > {log} 2>&1"

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
    log: "logs/STARindex.log"
    threads: 10
    conda: CONDA_SHARED_ENV
    shell:
        "STAR --runThreadN {threads} \
        --runMode genomeGenerate \
        --genomeDir {params.outdir} \
        --genomeFastaFiles {input.fa} \
        --sjdbGTFfile {input.gtf} \
        --sjdbOverhang 50 > {log} 2>&1"

rule Bowtie2index:
    input: "annotation/selected_CDS_extended.fa"
    output:
        fwd = "annotation/Bowtie2index/selected_CDS_extended.2.bt2",
        rev = "annotation/Bowtie2index/selected_CDS_extended.rev.2.bt2"
    params:
        outdir = "annotation/Bowtie2index/selected_CDS_extended"
    log: "logs/Bowtie2index.log"
    threads: 5
    conda: CONDA_SHARED_ENV
    shell:
        "bowtie2-build --seed 2020 --threads {threads} \
        {input} {params.outdir} > {log} 2>&1"
