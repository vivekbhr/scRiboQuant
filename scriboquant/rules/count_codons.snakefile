rule BamFilter:
    input:
        bam = "STAR/{sample}.sorted.bam",
        bed = "annotation/selected_CDS_exons.bed"
    output: temp("STAR/{sample}_tx.fastq")
    params:
        mapq = 255
    log: "logs/BamFilter_{sample}.log"
    threads: 5
    conda: CONDA_SHARED_ENV
    shell:
        """
        samtools view -h -q {params.mapq} -F 4 -L {input.bed} {input.bam} 2> {log} | \
        samtools fastq -@ {threads} -T "CB","UB" - | \
        awk -v RS="@" '{{ gsub("CB:Z:", "", $2);  gsub("UB:Z:", "", $3); \
        print "@"$1"_"$2"_"$3, $4, $5, $6 }}' | \
        awk 'OFS="\\n" {{ if (NF == 4) {{ print $1, $2, $3, $4}} }}' > {output} 2>> {log}
        """

rule CDSmap:
    input:
        fq = "STAR/{sample}_tx.fastq",
        index = "annotation/Bowtie2index/selected_CDS_extended.rev.2.bt2"
    output: temp("tx_bams/{sample}.bam")
    params:
        idx = "annotation/Bowtie2index/selected_CDS_extended",
        tmpfile = tempDir+"/{sample}"
    log: "logs/bowtie2_CDS.{sample}.log"
    threads: 10
    conda: CONDA_SHARED_ENV
    shell:
        """
        bowtie2 --end-to-end -p {threads} -x {params.idx} -U {input.fq} 2> {log} | \
        samtools sort -m 1G -T {params.tmpfile} -@ {threads} -O BAM -o {output} 2>> {log}
        """

rule idxBamBowtie:
    input: "tx_bams/{sample}.bam"
    output: temp("tx_bams/{sample}.bam.bai")
    threads: 1
    conda: CONDA_SHARED_ENV
    shell: "samtools index {input}"

rule umi_dedup:
    input:
        bam = "tx_bams/{sample}.bam",
        idx = "tx_bams/{sample}.bam.bai"
    output:
        bam = "deduplicated_bams/{sample}_tx.bam",
        stats = "QC/umi_dedup/{sample}_per_umi.tsv"
    params:
        mapq = 10,
        sample = "{sample}",
        splitRnameBy = "_"
    log:
        out = "logs/umi_dedup_{sample}.out",
        err = "logs/umi_dedup_{sample}.err"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell:
        """
        umi_tools dedup --mapping-quality {params.mapq} \
        --per-cell --per-gene --per-contig \
        --method=directional --umi-separator "{params.splitRnameBy}" \
        --output-stats=QC/umi_dedup/{params.sample} \
        -I {input.bam} -L {log.out} | \
        samtools view -F 4 -h | \
        awk -v sample={params.sample} \
        'OFS="\\t" {{ if($0 ~ "^@") {{print $0}} else \
        {{ split($1,a,"{params.splitRnameBy}"); print $0, "SM:Z:"sample"_"a[2], "CB:Z:"a[2], "UB:Z:"a[3], "MI:Z:"a[2]a[3] }} }}' | \
        samtools sort -@ {threads} -T {params.sample} -o {output.bam} > {log.out} 2> {log.err}
        """

rule idxBamDedup:
    input: "deduplicated_bams/{sample}_tx.bam"
    output: "deduplicated_bams/{sample}_tx.bam.bai"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell: "samtools index {input}"

rule prep_saf:
    input: "annotation/selected_CDS_annotation.bed"
    output: temp("CDS.saf")
    threads: 1
    conda: CONDA_SHARED_ENV
    shell: """
        awk 'OFS="\\t" {{ print $4, $1, $2, $3, $6 }}' {input} > {output}
        """

rule count_regions:
    input:
        bam = "tx_bams/{sample}.bam",
        idx = "tx_bams/{sample}.bam.bai",
        saf = "CDS.saf"
    output:
        counts = "counts/{sample}.CDScounts_bulk.tsv",
        bam = temp("counts/{sample}.bam.featureCounts.bam")
    params:
        filetype = "-F SAF"
    log: "logs/featurecounts_{sample}_bulk.err"
    threads: 10
    conda: CONDA_SHARED_ENV
    shell:
        "featureCounts -T {threads} -a {input.saf} {params.filetype} \
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
    params:
        mapq = 10
    log:
        out="logs/umi_counts_{sample}.CDScounts.out",
        err="logs/umi_counts_{sample}.CDScounts.err"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell:
        "umi_tools count --umi-separator '_' --mapping-quality {params.mapq} \
        --per-gene --gene-tag=XT \
        --per-cell --method=directional -I {input.bam} -S {output} \
        -v 4 --log2stderr --log={log.out} 2> {log.err}"
#        --cell-tag=CB --umi-tag=UB --extract-umi-method=tag \
if counts_codons:
    rule count_codons_cells:
        input:
            bed = "annotation/selected_CDS_annotation.bed",
            fasta = "annotation/selected_CDS_extended.fa",
            bc = barcodes,
            bam = "deduplicated_bams/{sample}_tx.bam",
            bai = "deduplicated_bams/{sample}_tx.bam.bai"
        output:
            tsv = "counts/{sample}_counts.Mtx",
            rownames = "counts/{sample}_rownames.tsv",
            colnames = "counts/{sample}_barcodes.txt"
        params:
            rscript = os.path.join(workflow.basedir, "tools", "countCodons.R"),
            offset = offset,
            prefix = "counts/{sample}"
        log:
            out="logs/codonCounts_{sample}.log"
        threads: 15
        conda: CONDA_R_ENV
        shell:
            "Rscript {params.rscript} {input.bed} {input.fasta} {input.bam} {input.bc} \
            {params.offset} {threads} {params.prefix} > {log} 2>&1"
