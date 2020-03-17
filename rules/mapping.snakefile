prep_annotation:


rule STAR:
    input:
    fastq_dir+"/{sample}"+reads[0]+".fastq.gz"
    output:
    bam = temp(aligner+"/{sample}.sorted.bam")
    params:
    alignerOptions = str(alignerOptions or ''),
    gtf = genes_gtf,
    index = star_index,
    prefix = aligner+"/{sample}/{sample}.",
    samsort_memory = '2G',
    sample_dir = aligner+"/{sample}",
    samtools_threads = 5
    benchmark:
    aligner+"/.benchmark/STAR.{sample}.benchmark"
    threads: 20  # 3.2G per core
    conda: CONDA_RNASEQ_ENV
    shell:
        """
        MYTEMP=$(mktemp -d ${{TMPDIR:-/tmp}}/snakepipes.XXXXXXXXXX)
        ( [ -d {params.sample_dir} ] || mkdir -p {params.sample_dir} )
        STAR --runThreadN {threads} \
            {params.alignerOptions} \
            --sjdbOverhang 100 \
            --outSAMunmapped Within \
            --outSAMtype BAM Unsorted \
            --outStd BAM_Unsorted \
            --sjdbGTFfile {params.gtf} \
            --genomeDir {params.index} \
            --readFilesIn <(gunzip -c {input}) \
            --outFileNamePrefix {params.prefix} \
            | samtools sort -m {params.samsort_memory} -T $MYTEMP/{wildcards.sample} -@ {params.samtools_threads} -O bam -o {output.bam} -
            rm -rf $MYTEMP
        """
