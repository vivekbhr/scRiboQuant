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
    input:
        R1 = "FASTQ/trimmed/{sample}_R1.fastq.gz",
        R2 = "FASTQ/trimmed/{sample}_R2.fastq.gz"
    output:
        R1 = "FASTQ/trimmed/{sample}_trimmed_R1.fastq.gz",
        R2 = "FASTQ/trimmed/{sample}_trimmed_R2.fastq.gz"
    log: "logs/cutadapt.{sample}.out"
    threads: 8
    conda: CONDA_SHARED_ENV
    shell:
        "cutadapt -j {threads} --minimum-length 18 --maximum-length 50 \
        -e 0.1 -q 20 -O 3 --trim-n -a TGGAATTCTCGG \
        -o {output.R1} -p {output.R2} {input.R1} {input.R2} > {log} 2>&1 "

rule FastQC:
    input:
        untrimmed = "FASTQ/trimmed/{sample}_R1.fastq.gz",
        trimmed = "FASTQ/trimmed/{sample}_trimmed_R1.fastq.gz"
    output:
        untrimmed = temp("QC/FastQC/{sample}_R1_fastqc.html"),
        trimmed = temp("QC/FastQC/{sample}_trimmed_R1_fastqc.html")
    params:
        outdir = "QC/FastQC",
        tmp = tempDir
    log: "logs/FastQC.{sample}.out"
    threads: 2
    conda: CONDA_SHARED_ENV
    shell:
        "fastqc -d {params.tmp} -t {threads} -o {params.outdir} \
        {input.untrimmed} {input.trimmed} > {log} 2>&1"

rule whitelist:
    input: barcodes
    output: "annotation/whitelist.txt"
    shell: "cut -f1 {input} > {output}"

rule STARsolo:
    input:
        idx = "annotation/STARindex/Genome",
        r1 = "FASTQ/trimmed/{sample}_trimmed_R1.fastq.gz",
        r2 = "FASTQ/trimmed/{sample}_trimmed_R2.fastq.gz",
        bc = "annotation/whitelist.txt"
    output:
        bam = "STAR/{sample}.sorted.bam",
        raw_counts = "STAR/{sample}/{sample}.Solo.out/Gene/raw/matrix.mtx"
    params:
        index = "annotation/STARindex",
        prefix = "STAR/{sample}/{sample}.",
        sample_dir = "STAR/{sample}"
    log: "logs/STAR.{sample}.log"
    threads: 10
    conda: CONDA_SHARED_ENV
    shell:
        """
        ## set
        MYTEMP=$(mktemp -d ${{TMPDIR:-/tmp}}/riboseq.XXXXXXXXXX);
        ( [ -d {params.sample_dir} ] || mkdir -p {params.sample_dir} )
        ## run
        STAR --runThreadN {threads} \
          --sjdbOverhang 50 \
          --seedSearchStartLmax 15 \
          --alignIntronMax 1000000 \
          --outSAMtype BAM SortedByCoordinate \
          --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
          --genomeDir {params.index} \
          --readFilesIn  {input.r1} {input.r2} \
          --readFilesCommand gunzip -c \
          --outFileNamePrefix {params.prefix} \
          --soloType CB_UMI_Simple \
          --soloFeatures Gene \
          --soloCBstart 1 \
          --soloCBlen 10 \
          --soloUMIstart 11 \
          --soloUMIlen 10 \
          --soloCBwhitelist {input.bc} \
          --soloBarcodeReadLength 0 \
          --soloCBmatchWLtype 1MM_multi \
          --soloStrand Forward \
          --soloUMIdedup 1MM_Directional \
          --soloCellFilter None > {log} 2>&1
        ## clean
        ln -rs {params.prefix}Aligned.sortedByCoord.out.bam {output.bam}
        rm -rf $MYTEMP
        """

rule idxBamSTAR:
    input: "STAR/{sample}.sorted.bam"
    output: "STAR/{sample}.sorted.bam.bai"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell: "samtools index {input}"

rule bamCoverage:
    input:
        bam = "STAR/{sample}.sorted.bam",
        bai = "STAR/{sample}.sorted.bam.bai"
    output: "bigWigs/{sample}_wholeGenome_Offset12.bw"
    params:
        ignore = "chrX chrY chrM",
        norm = '--normalizeUsing CPM',
        offset = offset
    log: "logs/bamCoverage.{sample}.log"
    threads: 10
    conda: CONDA_SHARED_ENV
    shell:
        "bamCoverage -bs 1 --Offset {params.offset} {params.norm} \
        --minMappingQuality 255 -p {threads} \
        -ignore {params.ignore}  \
        -b {input.bam} -o {output} > {log} 2>&1"
