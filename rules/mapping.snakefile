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

## preprocess fastq

~/programs/icetea/trimFastq test_R1.fastq test_R2.fastq out_R1.fastq.gz out_R2.fastq.gz

## map
STAR --runThreadN 20 \
  --sjdbOverhang 50 \
  --outSAMtype BAM SortedByCoordinate \
  --outSAMattributes NH HI AS nM CB UB \
  --sjdbGTFfile ${gtf} \
  --genomeDir ${idx} \
  --readFilesIn  ${r1} ${r2} \
  --readFilesCommand gunzip -c \
  --outFileNamePrefix ${prefix} \
  --soloType CB_UMI_Simple \
  --soloFeatures Gene \
  --soloUMIstart 11 \
  --soloUMIlen 10 \
  --soloCBstart 1 \
  --soloCBlen 10 \
  --soloCBwhitelist ${barcodes} \
  --soloBarcodeReadLength 0 \
  --soloCBmatchWLtype 1MM \
  --soloStrand Forward \
  --soloUMIdedup Exact \
  --soloUMIfiltering MultiGeneUMI \
  --quantMode TranscriptomeSAM \
  --quanTranscriptomeBam Singleend \# breaks RSEM compatibility, but gives more reads
  --clip5pNbases 10 0 # only for current reformatted fq


## umi dedup of the mapped file
bam=../STARout/Starv1-Leu6h_Aligned.toTranscriptome.sorted.bam
Aligned.sortedByCoord.out.bam
prefix=Starv1-Leu6h

# index
samtools index
# dedup

umi_tools dedup --mapping-quality 20 \
--per-cell --umi-tag=UB --cell-tag=CB --extract-umi-method=tag \
--method unique --spliced-is-unique \
--output-stats=${prefix}.stats \
-I ${bam} -L ${prefix}.log > ${prefix}.bam
