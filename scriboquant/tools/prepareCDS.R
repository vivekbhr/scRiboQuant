#!/usr/bin/env Rscript

suppressMessages(library(GenomicFeatures))
library(rtracklayer)
library(magrittr)
ldply <- plyr::ldply

Args <- commandArgs(trailingOnly = T)
inGTF <- Args[1] # "raw.gtf"
CDSout <- Args[2] #  "selected_CDS.bed"
BEDout <- Args[3] # "selected_CDS_annotation.bed"
Exonsout <- Args[4] # "selected_CDS_exons.bed"

## getting the top tx from a multi-tx gene
selectmax <- function(g) {
  ## toptx = tx with longest cds
  g <- g[g$cds == max(g$cds), , drop = F]
  toptx <- g$t
  # otherwise it's with longest 5 utr
  if(length(toptx) > 1) {
    g <- g[g$utr5 == max(g$utr5), , drop = F]
    toptx <- g$t
  }
  # otherwise it's with longest 3' UTR
  if(length(toptx) > 1) {
    g <- g[g$utr3 == max(g$utr3), , drop = F]
    toptx <- g$t
  }
  # otherwise it's the first one
  if(length(toptx) > 1) {
    toptx <- g$tx_id[1]
  }
  return(toptx)
}

## prep lengths
gtf = suppressWarnings(makeTxDbFromGFF(inGTF))
#gtf <- import.gff(inGTF)
#gtf <- gtf[gtf$source == "HAVANA" &
#             gtf$transcript_type == "protein_coding" &
#             grepl("appris_principal", gtf$tag) ]
#gtf <- keepStandardChromosomes(gtf)
#gtf <- makeTxDbFromGRanges(gtf)

txlist <- list(
  uCDS = cdsBy(gtf, "tx", use.names = TRUE),
  uUTR5 = fiveUTRsByTranscript(gtf, use.names = TRUE),
  uUTR3 = threeUTRsByTranscript(gtf, use.names = TRUE)
)
## first remove CDS whose length is not multiple of 3
badCDS <- txlist$uCDS[( ((txlist$uCDS %>% width() %>% sum()) %% 3) != 0) ] %>% names()
txlist <- lapply(txlist, function(x) return(x[!(names(x) %in% badCDS)]))

message("Selecting representative Tx per gene based on longest CDS/5UTR/3UTR length")
txinfo <- lapply(txlist, function(gr) return(width(gr) %>% ldply(sum))) %>%
    Reduce(function(x,y) merge(x, y, ".id", all.x = TRUE), .)
txinfo[is.na(txinfo)] <- 0

## tx2gene map
tx <- transcriptsBy(gtf, "gene") %>% unlist()
txinfo$gene_id <- names(tx[match(txinfo$.id, tx$tx_name)])
names(txinfo) <- c("tx_id", "cds", "utr5", "utr3", "gene_id")

# remove single-tx genes
message(paste0("Total genes: ", length(unique(txinfo$gene_id)), ". ",
              "Total tx: ", length(unique(txinfo$tx_id)), ":"))

singletx <- plyr::ldply(split(txinfo, txinfo$gene_id), function(x) if(length(x$tx_id) == 1) return(x[,"tx_id"]))
txinfo <- txinfo[!(txinfo$tx_id %in% singletx$V1), ]

message(paste0(nrow(txinfo), " transcripts, belonging to ",
               length(unique(txinfo$gene_id)),
                " genes need length based filter"))

# get top tx for others
multitx <- ldply(split(txinfo, txinfo$gene_id), selectmax)
alltx <- rbind(singletx, multitx)
colnames(alltx) <- c("gene","tx")

## double check and remove duplicated entries
l <- nrow(alltx) - length(unique(alltx$tx))
if(l > 0) {
  warning(paste0(l," multi-tx genes remain after filtering! Removing these"))
  dupgene <- alltx$gene[as.data.frame(table(alltx$gene))$Freq > 1]
  alltx <- alltx[!(alltx$gene %in% dupgene),]
}

## get CDS coordinates (+- 51b seq will be extracted using gffread)
finalCDS <- txlist$uCDS[alltx$tx]
# append gene name and symbol
names(finalCDS) %<>% paste(alltx$gene, sep = "|")

# also create a bed file with CDS coord +- 51bp, for counting/plotting
message("Creating annotation bed file")
txwidths <- ldply(width(finalCDS), sum)
names(txwidths) <- c("tx", "width")
txwidths$width2 <- txwidths$width + 102 # 51b on either side

bed <- GRanges(seqnames = txwidths$tx, IRanges(52, txwidths$width2 - 51),
               strand = strand(finalCDS) %>% unique() %>% unlist(),
               name = paste0(txwidths$tx))
strand(bed) <- "+"

# select all exons from these transcripts from the gtf
selExons <- exonsBy(gtf, "tx", use.names = T)[alltx$tx]

## write outputs
export.bed(finalCDS, CDSout)
export.bed(bed, con = BEDout)
export.bed(unlist(selExons), con = Exonsout)
