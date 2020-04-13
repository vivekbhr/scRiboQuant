def get_multiqc_input():
    file = [
        expand("QC/FastQC/{sample}_{read}_fastqc.html", sample = samples, read=['R1', 'trimmed_R1']),
        expand("STAR/{sample}/{sample}.Solo.out/Gene/raw/matrix.mtx", sample = samples),
        expand("deduplicated_bams/{sample}.bam", sample = samples)
        ]

    return(file)

rule multiQC:
    input: get_multiqc_input()
    output: "QC/multiqc_report.html"
    params:
        outdir = "QC",
        logdir = "logs",
        star = "STAR"
    log: "logs/multiqc.out"
    threads: 1
    #conda: CONDA_SHARED_ENV
    shell:
        "multiqc -f -o {params.outdir} {params.outdir} {params.logdir} {params.star} > {log} 2>&1"
