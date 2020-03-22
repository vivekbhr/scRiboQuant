gtf=/hpc/hub_oudenaarden/vbhardwaj/annotations/hg38_gencode31_extra/annotation/genes.gtf
idx=../annotation/STARindex/
barcodes=../whitelist.txt
r1=../fastq_reformatted/MAV-RPFv4-HEK293T_R1.fastq.gz
r2=../fastq_reformatted/MAV-RPFv4-HEK293T_R2.fastq.gz
prefix=RPFv4_
fa=

## generate genome
STAR --runThreadN 10 \
--runMode genomeGenerate \
--genomeDir STARindex \
--genomeFastaFiles $fa \
--sjdbGTFfile $gtf \
--sjdbOverhang 50

### some follow up of solo counts

dir.create("~/Desktop/riboseq_test")
setwd("riboseq_test/")
library(GenomicFeatures)
