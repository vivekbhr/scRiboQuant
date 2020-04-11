#!/usr/bin/env Rscript

#setwd("~/hpc/hub_oudenaarden/vbhardwaj/2020_scRibo_seq/03_automated_workflow")
library(GenomicFeatures)
library(GenomicAlignments)
library(magrittr)
library(Matrix)
library(BiocParallel)

oArgs <- commandArgs(trailingOnly = T)
cds <- Args[1] #"annotation/selected_CDS_annotation.bed" #  #
fa <- Args[2] #"annotation/selected_CDS_extended.fa"#
bam <- Args[3] #"dedup/RPFv4-HEK293T-Starv1_S3.dedup.bam"#
barcodes <- Args[4] #"barcodeInfo.tsv" #
offset <- as.numeric(Args[5])
cores <- Args[6]
outPrefix <- Args[7]

if(offset < 0) {
  fix = "end"
} else {
  fix = "start"
}
# register cores
register(MulticoreParam(workers=cores))

## read files
cds <- rtracklayer::import.bed(cds)
fa2 <- readDNAStringSet(fa)
barcodes <- read.delim(barcodes, header = F, stringsAsFactors = F)$V1

### ------ Functions --------------- ####
ResizeReads <- function(reads, width = 1, fix = fix, trimlength = offset) {
  reads <- as(reads, "GRanges")
  stopifnot(all(GenomicRanges::strand(reads) != "*"))
  #GenomicRanges::resize(reads, width = width, fix = fix)
  GenomicRanges::resize(reads, fix = "center", width = 1)
  #reads - trimlength
}

getCounts <- function(regions, barcodes, bam) {
  BiocParallel::bplapply(barcodes, function(bc){
    o <- summarizeOverlaps(regions, bam, ignore.strand = FALSE,
                           param = ScanBamParam(tag = "CB",
                                                tagFilter = list("CB" = bc)),
                           mode = "Union",
                           preprocess.reads=ResizeReads)
    return(o)
  }) -> counts

  return(counts)
}

#### --------------------- PREPARE -------------------
## prepare 3bp annotation
fa_names <- strsplit(names(fa2), " ") %>% sapply(function(x) return(x[1]))
cds <- cds[cds$name %in% fa_names]
cds_split <- GenomicRanges::tile(cds, width=3) %>% unlist()
cds_codons <- getSeq(FaFile(fa), cds_split)
cds_split$codon <- as.character(cds_codons)
cds_split$AA <- translate(cds_codons) %>% as.character()


#### --------------------- COUNT ---------
system.time(
  regionCounts <- getCounts(cds_split, barcodes, bam)
  )
counts <- Reduce(cbind, lapply(regionCounts, assay))
counts <- Matrix::Matrix(counts, sparse = T)
#aa_counts <- Matrix.utils::aggregate.Matrix(counts, groupings=factor(cds_split$AA), fun="sum")
#codon_counts <-  Matrix.utils::aggregate.Matrix(counts, groupings=factor(cds_split$codon), fun="sum")

### --- SAVE
writeMM(counts, file = paste0(outPrefix, "_counts.Mtx"))
genemeta <- data.frame(gene = seqnames(cds_split),
                       codon = cds_split$codon,
                       aa = cds_split$AA)
rownames(genemeta) <- paste(genemeta$gene, genemeta$codon, genemeta$counts, sep = "_")
write.table(genemeta, file = paste0(outPrefix, "_rownames.tsv"), sep = "\t", quote=F)
write.table(barcodes, file = paste0(outPrefix, "_barcodes.txt"), sep = "\t", quote=F, row.names=F, col.names=F)
