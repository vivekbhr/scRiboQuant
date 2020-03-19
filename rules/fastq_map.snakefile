## optional downsampling
if downsample is None:
    rule FASTQlink:
        input:
            r1 = indir+"/{sample}"+reads[0]+ext,
            r2 = indir+"/{sample}"+reads[1]+ext
        output:
            r1 = "FASTQ/{sample}_R1.fastq.gz",
            r2 = "FASTQ/{sample}_R2.fastq.gz"
        conda: CONDA_SHARED_ENV
        shell:
          """
          ( [ -f {output.r1} ] || ln -s -r {input.r1} {output.r1} )
          ( [ -f {output.r2} ] || ln -s -r {input.r2} {output.r2} )
          """
else:
    rule FASTQdownsample:
        input:
            r1 = indir+"/{sample}"+reads[0]+ext,
            r2 = indir+"/{sample}"+reads[1]+ext
        output:
            r1 = "FASTQ/{sample}_R1.fastq.gz",
            r2 = "FASTQ/{sample}_R2.fastq.gz"
        params:
            num_reads = downsample
        threads: 10
        conda: CONDA_SHARED_ENV
        shell:
          """
          seqtk sample -s 100 {input.r1} {params.num_reads} | pigz -p {threads} -9 > {output.r1}
          seqtk sample -s 100 {input.r2} {params.num_reads} | pigz -p {threads} -9 > {output.r2}
          """

# preprocess per-sample fastqs for UMI and Cell barcode
rule preprocess:
    input:
        r1 = "FASTQ/{sample}_R1.fastq.gz",
        r2 = "FASTQ/{sample}_R2.fastq.gz"
    output:
        r1 = temp("FASTQ/trimmed/{sample}_R1.fastq.gz"),
        r2 = temp("FASTQ/trimmed/{sample}_R2.fastq.gz")
    params:
        trimFq = os.path.join(workflow.basedir, "tools/trimFastq/trimFastq")
    shell:
        "{params.trimFq} {input.r1} {input.r2} {output.r1} {output.r2}"

rule cutadapt:
    input: "FASTQ/trimmed/{sample}_R1.fastq.gz"
    output: "FASTQ/trimmed/{sample}_trimmed_R1.fastq.gz"
    log: "logs/cutadapt.{sample}.out"
    threads: 8
    conda: CONDA_SHARED_ENV
    shell:
        "cutadapt -j {threads} \
        -e 0.1 -q 30 -O 3 --trim-n -a TGGAATTCTCGG \
        -o {output} {input} > {log} 2>&1 "

rule FastQC:
    input:
        untrimmed = "FASTQ/trimmed/{sample}_R1.fastq.gz",
        trimmed = "FASTQ/trimmed/{sample}_trimmed_R1.fastq.gz"
    output:
        untrimmed = "QC/FastQC/{sample}_R1_fastqc.html",
        trimmed = "QC/FastQC/{sample}_trimmed_R1_fastqc.html"
    params:
        outdir = "QC/FastQC"
    log: "logs/FastQC.{sample}.out"
    threads: 2
    conda: CONDA_SHARED_ENV
    shell:
        "fastqc -t {threads} -o {params.outdir} {input.untrimmed} {input.trimmed} > {log} 2>&1"

rule STARsolo:
    input:
        r1 = "FASTQ/trimmed/{sample}_trimmed_R1.fastq.gz",
        r2 = "FASTQ/trimmed/{sample}_R2.fastq.gz",
        gtf = annotation+"/genes.gtf",
        bc = barcodes
    output:
        bam = "STAR/{sample}.sorted.bam",
        txbam = "STAR/{sample}.tx.sorted.bam",
        raw_counts = "STAR/{sample}/{sample}.Solo.out/Gene/raw/matrix.mtx",
        filtered_counts = "STAR/{sample}/{sample}.Solo.out/Gene/filtered/matrix.mtx",
        filtered_bc = "STAR/{sample}/{sample}.Solo.out/Gene/filtered/barcodes.tsv"
    params:
        index = annotation+"/STARindex",
        prefix = "STAR/{sample}/{sample}.",
        sample_dir = "STAR/{sample}"
    threads: 20
    shell:
        """
        ## set
        MYTEMP=$(mktemp -d ${{TMPDIR:-/tmp}}/riboseq.XXXXXXXXXX);
        ( [ -d {params.sample_dir} ] || mkdir -p {params.sample_dir} )
        ## run
        STAR --runThreadN {threads} \
          --sjdbOverhang 50 \
          --outSAMtype BAM SortedByCoordinate \
          --outSAMattributes NH HI AS nM CB UB \
          --sjdbGTFfile {input.gtf} \
          --genomeDir {params.index} \
          --readFilesIn  {input.r1} {input.r2} \
          --readFilesCommand gunzip -c \
          --outFileNamePrefix {params.prefix} \
          --soloType CB_UMI_Simple \
          --soloFeatures Gene \
          --soloUMIstart 11 \
          --soloUMIlen 10 \
          --soloCBstart 1 \
          --soloCBlen 10 \
          --soloCBwhitelist {input.bc} \
          --soloBarcodeReadLength 0 \
          --soloCBmatchWLtype 1MM \
          --soloStrand Forward \
          --soloUMIdedup Exact \
          --soloUMIfiltering MultiGeneUMI \
          --quantMode TranscriptomeSAM \
          --quanTranscriptomeBam Singleend
        ## clean
        ln -rs {params.prefix}Aligned.sortedByCoord.out.bam {output.bam}
        ln -rs {params.prefix}Aligned.toTranscriptome.out.bam  {output.txbam}
        rm -rf $MYTEMP
        """
