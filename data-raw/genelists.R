library(readr)
library(dplyr)
library(GenomicRanges)

ref <- c(
  "hg18", "hg19", "hg38", # human
  "mm8", "mm9", "mm10", #mouse
  "rn5", "rn6" # rat
)

REFGENE_URL <- "http://hgdownload.cse.ucsc.edu/goldenPath/%s/database/refGene.txt.gz"

genes <- lapply(ref, function(refname) {
  message(refname)
  D <- read_tsv(sprintf(REFGENE_URL, refname), col_names=FALSE) %>%
    select(3,5,6,13) %>%
    magrittr::set_colnames(c("chr","start","end","name")) %>%
    mutate(chr = gsub("^chr", "", chr))
  
  # merge overlapping gene ranges in same chromosome
  ranges <- makeGRangesFromDataFrame(D, keep.extra.columns=TRUE, ignore.strand=TRUE) %>%
    split(D$name) %>%
    reduce %>%
    unlist
  
  glist <- data.frame(
    chr=as.character(seqnames(ranges)),
    start=start(ranges),
    end=end(ranges),
    name=names(ranges),
    stringsAsFactors=FALSE, row.names=NULL
  )
})
names(genes) <- ref

usethis::use_data(genes, overwrite=TRUE, internal=TRUE)
