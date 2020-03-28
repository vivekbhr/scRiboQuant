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
        r1 = "FASTQ/trimmed/{sample}_R1.fastq.gz", #temp() later
        r2 = "FASTQ/trimmed/{sample}_R2.fastq.gz"
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
        "cutadapt -j {threads} --minimum-length 18 --maximum-length 50 \
        -e 0.1 -q 20 -O 3 --trim-n -a TGGAATTCTCGG \
        -o {output} {input} > {log} 2>&1 "

rule FastQC:
    input:
        untrimmed = "FASTQ/trimmed/{sample}_R1.fastq.gz",
        trimmed = "FASTQ/trimmed/{sample}_trimmed_R1.fastq.gz"
    output:
        untrimmed = "QC/FastQC/{sample}_R1_fastqc.html",
        trimmed = "QC/FastQC/{sample}_trimmed_R1_fastqc.html"
    params:
        outdir = "QC/FastQC",
        tmp = tempDir
    log: "logs/FastQC.{sample}.out"
    threads: 2
    conda: CONDA_SHARED_ENV
    shell:
        "fastqc -d {params.tmp} -t {threads} -o {params.outdir} \
        {input.untrimmed} {input.trimmed} > {log} 2>&1"

#        txbam = "STAR/{sample}/{sample}.Aligned.toTranscriptome.out.bam",
rule STARsolo:
    input:
        r1 = "FASTQ/trimmed/{sample}_trimmed_R1.fastq.gz",
        r2 = "FASTQ/trimmed/{sample}_R2.fastq.gz",
        gtf = annotation+"/genes.gtf",
        bc = barcodes
    output:
        bam = "STAR/{sample}.sorted.bam",
        raw_counts = "STAR/{sample}/{sample}.Solo.out/Gene/raw/matrix.mtx",
        filtered_counts = "STAR/{sample}/{sample}.Solo.out/Gene/filtered/matrix.mtx",
        filtered_bc = "STAR/{sample}/{sample}.Solo.out/Gene/filtered/barcodes.tsv"
    params:
        index = annotation+"/STARindex",
        prefix = "STAR/{sample}/{sample}.",
        sample_dir = "STAR/{sample}"
    log: "logs/STAR.{sample}.log"
    threads: 20
    conda: CONDA_SHARED_ENV
    shell:
        """
        ## set
        MYTEMP=$(mktemp -d ${{TMPDIR:-/tmp}}/riboseq.XXXXXXXXXX);
        ( [ -d {params.sample_dir} ] || mkdir -p {params.sample_dir} )
        ## run
        STAR --runThreadN {threads} \
          --sjdbOverhang 50 \
          --outSAMtype BAM SortedByCoordinate \
          --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM
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
          --soloCBmatchWLtype 1MM_multi_pseudocounts \
          --soloStrand Forward \
          --soloUMIdedup Exact \
          --soloUMIfiltering MultiGeneUMI > {log} 2>&1
        ##--quantMode TranscriptomeSAM \
        ##--quantTranscriptomeBan Singleend
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


## re-map the reads mapping to tx
#rule Bam2Fq:
#    input: "STAR/{sample}/{sample}.Aligned.toTranscriptome.out.bam"
#    output: temp("STAR/{sample}.fastq.gz")
#    threads: 5
#    conda: CONDA_SHARED_ENV
#    shell: "samtools fastq -@ {threads} {input} | gzip - > {output}"

rule BamFilter:
    input: "STAR/{sample}.sorted.bam"
    output: temp("STAR/{sample}_tx.bam")
    params:
        txbed = annotation+"/selected_CDS.bed",
        tmpfile = tempDir+"/{sample}",
        mapq = 255
    threads: 5
    conda: CONDA_SHARED_ENV
    shell:
        "samtools view -q {params.mapq} -F 4 -L {params.txbed} {input} | \
         samtools sort -m 1G -n -T {params.tmpfile} -O BAM -@ {threads} -o {output} - \
         2> /dev/null"

rule CDSmap:
    input:
        bam = "STAR/{sample}_tx.bam", #"STAR/{sample}.fastq.gz",
        index = annotation+"/Bowtie2index/selected_CDS_51b.rev.2.bt2",
    output: "Bowtie2_CDS/{sample}.bam"
    params:
        idx = annotation+"/Bowtie2index/selected_CDS_51b",
        tmpfile = tempDir+"/{sample}"
    log: "logs/bowtie2_CDS.{sample}.log"
    threads: 10
    conda: CONDA_SHARED_ENV
    shell: "bowtie2 --end-to-end --preserve-tags -p {threads} -x {params.idx} -b {input.bam} 2> {log} |\
            samtools sort -m 1G -T {params.tmpfile} -@ {threads} -O BAM -o {output} 2> /dev/null"

rule idxBamBowtie:
    input: "Bowtie2_CDS/{sample}.bam"
    output: "Bowtie2_CDS/{sample}.bam.bai"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell: "samtools index {input}"

rule umi_dedup:
    input:
        bam = "Bowtie2_CDS/{sample}.bam",
        idx = "Bowtie2_CDS/{sample}.bam.bai"
    output:
        bam = "dedup/{sample}.dedup.bam",
        stats = "QC/umi_dedup/{sample}_per_umi.tsv"
    params:
        mapq = 10,
        sample = "{sample}"
    log:
        out = "logs/umi_dedup_{sample}.out",
        err = "logs/umi_dedup_{sample}.err"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell:
        """
        umi_tools dedup --mapping-quality {params.mapq} \
        --per-cell --per-gene --per-contig \
        --method unique \
        --output-stats=QC/umi_dedup/{params.sample} \
        -I {input.bam} -L {log.out} | \
        samtools view -F 4 -h | \
        awk -v sample={params.sample} \
        'OFS="\\t" {{ if($0 ~ "^@") {{print $0}} else \
        {{ split($1,a,"_"); print $0, "SM:Z:"sample"_"a[2], "CB:Z:"a[2], "UB:Z:"a[3], "MI:Z:"a[2]a[3] }} }}' | \
        samtools sort -@ {threads} -T {params.sample} -o {output.bam} > {log.out} 2> {log.err}
        """

rule idxBamDedup:
    input: "dedup/{sample}.dedup.bam"
    output: "dedup/{sample}.dedup.bam.bai"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell: "samtools index {input}"

rule bamCoverage:
    input:
        bam = "dedup/{sample}.dedup.bam",
        bai = "dedup/{sample}.dedup.bam.bai"
    output: "bigWigs/{sample}_CDS_Offset12.bw"
    params:
        ignore = "chrX chrY chrM",
        norm = '--normalizeUsing CPM',
    log: "logs/bamCoverage.{sample}.log"
    threads: 10
    conda: CONDA_SHARED_ENV
    shell:
        "bamCoverage -bs 1 --Offset -12 \
        --minMappingQuality 10 -p {threads} \
        -ignore {params.ignore}  \
        -b {input.bam} -o {output} > {log} 2>&1"
