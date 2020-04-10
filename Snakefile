import os
import glob
import yaml

## some internal functions  ###################################################
def load_configfile(configfile):
   with open(configfile, "r") as f:
       config = yaml.load(f, Loader=yaml.FullLoader)
   return(config)

def set_condaEnv():
    return{'CONDA_SHARED_ENV': 'env.yaml'}

def get_sample_names(infiles, ext, reads):
    """
    Get sample names without file extensions
    """
    s = set()
    lext = len(ext)
    l0 = len(reads[0])
    l1 = len(reads[1])
    for x in infiles:
        x = os.path.basename(x)[:-lext]
        if x.endswith(reads[0]):
            x = x[:-l0]
        elif x.endswith(reads[1]):
            x = x[:-l1]
        else:
            continue
        s.add(x)
    return sorted(list(s))

def get_bcList(bc):
    with open(bc) as f:
        lis = [str(line.split()[0]) for line in f]
    return(lis)


# update envs
globals().update(set_condaEnv())
# load config file
globals().update(load_configfile(workflow.overwrite_configfiles[0]))

## load samples
infiles = sorted(glob.glob(os.path.join(indir, '*'+ext)))
samples = get_sample_names(infiles,ext,reads)
bclist = get_bcList(barcodes)

### include modules of other snakefiles ########################################
################################################################################
include: os.path.join(workflow.basedir, "rules", "create_annotation.snakefile")
include: os.path.join(workflow.basedir, "rules", "fastq_map.snakefile")
include: os.path.join(workflow.basedir, "rules", "count_codons.snakefile")
include: os.path.join(workflow.basedir, "rules", "QC.snakefile")
if nCellsCoverage > 0:
    include: os.path.join(workflow.basedir, "rules", "scCoverage.snakefile")

### conditional/optional rules #################################################
################################################################################

def prep_scCoverage():
    if nCellsCoverage > 0:
        out = expand("bigWigs/{sample}_{barcode}_Offset12.bw", sample = samples, barcode = bclist[0:nCellsCoverage])
    else:
        out = []
    return(out)

def count_codons():
    if counts_codons:
        out = expand("counts/{sample}_counts.Mtx", sample = samples)
    else:
        out = []
    return(out)

### main rule ##################################################################
################################################################################
def prep_annotation():
    out = ["annotation/genome.fa", "annotation/gtf_annotation_table.txt",
          "annotation/selected_CDS.bed", "annotation/selected_CDS_annotation.bed",
          "annotation/STARindex/Genome", "annotation/Bowtie2index/selected_CDS_extended.rev.2.bt2"]
    return(out)

rule all:
    input:
        prep_annotation(),
        expand("FASTQ/{sample}_{read}.fastq.gz", sample = samples, read = ['R1', 'R2']),
        expand("FASTQ/trimmed/{sample}_trimmed_R1.fastq.gz", sample = samples),
        expand("QC/FastQC/{sample}_{read}_fastqc.html", sample = samples,
                read=['R1', 'trimmed_R1']),
        expand("STAR/{sample}.sorted.bam", sample = samples),
        expand("STAR/{sample}/{sample}.Solo.out/Gene/raw/matrix.mtx", sample = samples),
        expand("bigWigs/{sample}_wholeGenome_Offset12.bw", sample = samples),
        expand("deduplicated_bams/{sample}_tx.bam", sample = samples),
        expand("counts/{sample}.CDScounts_per_barcode.tsv", sample = samples),
        prep_scCoverage(),
        count_codons(),
        "QC/multiqc_report.html"


### execute after workflow finished ############################################
################################################################################
onsuccess:
    if "verbose" in config and config["verbose"]:
        print("\n--- scRiboQuant workflow finished successfully! ------------------\n")
onerror:
    print("\n !!!! ERROR in scRiboQuant workflow! !!!!\n")
