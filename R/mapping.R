#' Create mapping from genes to SNPs.
#'
#' @param snps A data frame mapping SNPs to genomic positions.
#' @param ref Reference genome. One of "hg18", "hg19", "hg38", "mm8", "mm9", "mm10", "rn5", "rn6".
#' @param maxgap Maximum allowed distance between SNP and gene.
#' @param select How should SNPs near multiple genes by mapped? Use "all" to use all hits. Use "nearest" to use only nearest gene.
#' @param chr.field A character vector of recognized names for the column in `snps` containing the chromosome name.
#' @param pos.field A character vector of recognized names for the column in `snps` containing the SNP position.
#' @param snp.field A character vector of recognized names for the column in `snps` containing the SNP name.
#' @param permute Permutation method.
#'     Use "none" for no permutation.
#'     "snps" randomly assigns SNPs to genes, while preserving the frequency of each SNP and number of SNPs assigned to each gene.
#'     "genes" randomly shuffles gene locations among all genes.
#'     "genes.binned" shuffles gene locations within bins of similar gene size.
#' @param bins Number of bins to use for permute = "genes.binned".
#' @return A data frame mapping genes to SNPs.
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom GenomicRanges findOverlaps
#' @export
make_genemap <- function(
  snps, ref="hg19", maxgap=10e3, select="all",
  chr.field=c("chr","chrom","chromosome"),
  pos.field=c("bp","position","pos"),
  snp.field=c("snp.name","snp","rsid"),
  permute="none", bins=100
) {
  if(is.null(genes[[ref]])) {
    stop("Unknown reference genome build \"", ref, "\".")
  }
  select <- match.arg(select, c("all","nearest"))
  permute <- match.arg(permute, c("none", "snps", "genes", "genes.binned"))

  glist <- genes[[ref]]

  chr.col <- .match_column(snps, chr.field)
  if(is.na(chr.col)) stop("Chromosome column not found.")
  pos.col <- .match_column(snps, pos.field)
  if(is.na(pos.col)) stop("SNP position column not found.")
  snp.col <- .match_column(snps, snp.field)
  if(is.na(snp.col)) stop("SNP name column not found.")

  snps <- data.frame(
    chr=snps[[chr.col]],
    pos=snps[[pos.col]],
    snp=snps[[snp.col]],
    stringsAsFactors=FALSE,
    row.names=NULL
  )

  if(permute == "genes") {
    glist$name <- sample(glist$name)
  }
  else if(permute == "genes.binned") {
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
    seqnames.field="chr",
    start.field="pos",
    end.field="pos",
    keep.extra.columns=TRUE, ignore.strand=TRUE
  )

  if(permute == "snps") {
    snp_ranges$snp <- sample(snp_ranges$snp)
  }

  hits <- findOverlaps(snp_ranges, gene_ranges, maxgap=maxgap, select="all")

  if(select == "all") {
    out <- data.frame(gene=gene_ranges$name[hits@to], snp=snp_ranges$snp[hits@from], stringsAsFactors=FALSE)
  }
  else if(select == "nearest") {
    snp_ranges <- snp_ranges[unique(hits@from),]
    near <- nearest(snp_ranges, gene_ranges)
    out <- data.frame(gene=gene_ranges$name[near], snp=snp_ranges$snp, stringsAsFactors=FALSE)
  }

  unique(out)
}

.match_column <- function(x, fields) {
  i <- match(tolower(fields), tolower(colnames(x)))
  na.omit(i)[1]
}
