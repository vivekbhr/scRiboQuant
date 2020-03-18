def get_multiqc_input():
    file = [
        expand("QC/FastQC/{sample}_{read}_fastqc.html", sample = samples, read=['R1', 'trimmed_R1']),
        expand("STAR/{sample}/{sample}.Solo.out/Gene/filtered/barcodes.tsv", sample = samples)
        ]

    return(file)

rule multiQC:
    input: get_multiqc_input()
    output: "QC/multiqc_report.html"
    params:
        outdir = "QC"
    log: "logs/multiqc.out"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell:
        "multiqc -f -o {params.outdir} {params.outdir} > {log} 2>&1"
