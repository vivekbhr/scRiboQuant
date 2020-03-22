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

# update envs
globals().update(set_condaEnv())
# load config file
globals().update(load_configfile(workflow.overwrite_configfiles[0]))

## load samples
infiles = sorted(glob.glob(os.path.join(indir, '*'+ext)))
samples = get_sample_names(infiles,ext,reads)

### include modules of other snakefiles ########################################
################################################################################
include: os.path.join(workflow.basedir, "rules", "fastq_map.snakefile")
include: os.path.join(workflow.basedir, "rules", "QC.snakefile")

### conditional/optional rules #################################################
################################################################################


### main rule ##################################################################
################################################################################

rule all:
    input:
        expand("FASTQ/{sample}_{read}.fastq.gz", sample = samples, read = ['R1', 'R2']),
        expand("FASTQ/trimmed/{sample}_trimmed_R1.fastq.gz", sample = samples),
        expand("QC/FastQC/{sample}_{read}_fastqc.html", sample = samples, read=['R1', 'trimmed_R1']),
        expand("STAR/{sample}.sorted.bam", sample = samples),
        expand("STAR/{sample}/{sample}.Solo.out/Gene/filtered/barcodes.tsv", sample = samples),
        expand("bigWigs/{sample}_wholeGenome.cpm.bw", sample = samples),
        expand("Bowtie2_CDS/{sample}.bam", sample = samples),
        expand("Bowtie2_CDS/{sample}.dedup.bam", sample = samples),
        "QC/multiqc_report.html"


### execute after workflow finished ############################################
################################################################################
onsuccess:
    if "verbose" in config and config["verbose"]:
        print("\n--- scRiboQuant workflow finished successfully! ------------------\n")
onerror:
    print("\n !!!! ERROR in scRiboQuant workflow! !!!!\n")
