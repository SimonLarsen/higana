#' Create mapping from genes to SNPs.
#'
#' @param snps A data frame mapping SNPs to genomic positions. Must have the columns "chromosome", "snp.name" and "position".
#' @param ref Reference genome. One of "hg18", "hg19", "hg38".
#' @param maxgap Maximum allowed distance between SNP and gene.
#' @param select How should SNPs near multiple genes by mapped? Use "all" to use all hits. Use "nearest" to use only nearest gene.
#' @param permute Permutation method. Use "none" for no permutation. "all" randomly shuffles all gene names. "binned" shuffles gene names within bins of similar gene size.
#' @param bins Number of bins to use for permute = "preserve.size".
#' @return A data frame mapping genes to SNPs.
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom GenomicRanges findOverlaps
#' @export
make_genemap <- function(snps, ref="hg19", maxgap=10e3, select="all", permute="none", bins=100) {
  if(is.null(genes[[ref]])) {
    stop("Unknown reference genome build \"", ref, "\".")
  }
  select <- match.arg(select, c("all","nearest"))

  glist <- genes[[ref]]

  permute <- match.arg(permute, c("none", "all", "binned"))
  if(permute == "all") {
    glist$name <- sample(glist$name)
  }
  else if(permute == "binned") {
    bins <- 100
    glist$size <- glist$end - glist$start + 1
    glist <- glist[order(glist$size),]

    glist$bin <- floor(seq(1, bins+1-1/nrow(glist), length.out=nrow(glist)))
    glist.bins <- split(glist, glist$bin)
    glist.bins <- lapply(glist.bins, function(b) { b$name <- sample(b$name); b })
    glist <- do.call(rbind, glist.bins)
  }

  gene_ranges <- makeGRangesFromDataFrame(glist, keep.extra.columns=TRUE, ignore.strand=TRUE)

  snp_ranges <- makeGRangesFromDataFrame(
    snps,
    start.field="position", end.field="position",
    keep.extra.columns=TRUE, ignore.strand=TRUE
  )

  hits <- findOverlaps(snp_ranges, gene_ranges, maxgap=maxgap, select="all")

  if(select == "all") {
    out <- data.frame(gene=gene_ranges$name[hits@to], snp=snp_ranges$snp.name[hits@from], stringsAsFactors=FALSE)
  }
  else if(select == "nearest") {
    snp_ranges <- snp_ranges[unique(hits@from),]
    near <- nearest(snp_ranges, gene_ranges)
    out <- data.frame(gene=gene_ranges$name[near], snp=snp_ranges$snp.name, stringsAsFactors=FALSE)
  }

  unique(out)
}
