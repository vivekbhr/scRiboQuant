
rule split_bam:
    input:
        bam = "STAR/{sample}.sorted.bam",
        bc = barcodes
    output:
        temp(expand("split_bam/{{sample}}.CB_{barcode}.bam", barcode = bclist[0:nCellsCoverage]))
    params:
        outprefix="split_bam/{sample}",
        mapq = 255
    wildcard_constraints:
        barcode="[ATGCN]*"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell:
        "samtools view -h -b -q {params.mapq} {input.bam} | \
        bamtools split -in - -tag CB -tagPrefix '' -stub {params.outprefix}"

rule idx_split_bam:
    input: "split_bam/{sample}.CB_{barcode}.bam"
    output: temp("split_bam/{sample}.CB_{barcode}.bam.bai")
    wildcard_constraints:
        barcode="[ATGCN]*"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell:
        "samtools index {input}"


rule bamCoverage_sc:
    input:
        bam = "split_bam/{sample}.CB_{barcode}.bam",
        bai = "split_bam/{sample}.CB_{barcode}.bam.bai"
    output: "bigWigs/{sample}_{barcode}_Offset12.bw"
    wildcard_constraints:
        barcode="[ATGCN]*"
    params:
        ignore = "chrX chrY chrM",
        norm = '--normalizeUsing CPM',
        offset = offset
    log: "logs/bamCoverage.{sample}_{barcode}.log"
    threads: 4
    conda: CONDA_SHARED_ENV
    shell:
        "bamCoverage -bs 1 --Offset {params.offset} {params.norm} \
        --minMappingQuality 255 -p {threads} \
        -ignore {params.ignore}  \
        -b {input.bam} -o {output} > {log} 2>&1"
