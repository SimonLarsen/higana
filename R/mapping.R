#' Map genes to SNPs.
#' @param snps A
#' @param ref Reference genome. One of "hg18", "hg19", "hg38".
#' @param maxgap Maximum allowed distance between SNP and gene.
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom GenomicRanges findOverlaps
#' @export
get_genemap <- function(snps, ref="hg19", maxgap=10e3) {
  if(is.null(genes[[ref]])) {
    stop("Unknown reference genome build \"", ref, "\".")
  }

  snp_ranges <- makeGRangesFromDataFrame(
    snps,
    start.field="position", end.field="position",
    keep.extra.columns=TRUE, ignore.strand=TRUE
  )

  gene_ranges <- makeGRangesFromDataFrame(genes[[ref]], keep.extra.columns=TRUE, ignore.strand=TRUE)

  hits <- findOverlaps(gene_ranges, snp_ranges, maxgap=maxgap, select="all")

  out <- data.frame(gene=gene_ranges$name[hits@from], snp=snp_ranges$snp.name[hits@to])
  out <- unique(out)
  split(out$snp, out$gene)
}
