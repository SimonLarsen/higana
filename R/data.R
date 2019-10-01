#' Open PLINK .bed file and create gene mapping.
#'
#' @param path Path to .bed file with or without extension.
#' @param ref Reference genome. One of "hg18", "hg19", "hg38", "mm8", "mm9", "mm10", "rn5", "rn6".
#' @param maxgap Maximum allowed distance between SNP and gene.
#' @param select How should SNPs near multiple genes by mapped? Use "all" to use all hits. Use "nearest" to use only nearest gene.
#' @return A \code{BEDData} object.
#' @importFrom BEDMatrix BEDMatrix
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom GenomicRanges findOverlaps
#' @export
read_bed <- function(
  path,
  ref="hg19",
  maxgap=10e3,
  select="all"
) {
  select <- match.arg(select, c("all","nearest"))
  if(is.null(genes[[ref]])) stop(sprintf("Reference genome '%s' not found.", ref))

  # Extract bed/bim/bam file prefix
  if(endsWith(tolower(path), ".bed")) {
    path <- substring(path, 1, nchar(path)-4)
  }

  # Get path to .bim file
  if(file.exists(paste0(path, ".bim"))) {
    bimpath <- paste0(path, ".bim")
  } else if (file.exists(paste0(path, ".BIM"))) {
    bimpath <- paste0(path, ".BIM")
  } else {
    stop("Could not find corresponding .bim file.")
  }

  bed <- BEDMatrix(path, simple_names=TRUE)

  message(sprintf("Creating gene map for reference genome '%s'", ref))
  map <- .make_genemap(bimpath, ref, maxgap, select)

  o <- list(snps=bed, map=map)
  class(o) <- "BEDData"
  return(o)
} 

.make_genemap <- function(path, ref, maxgap, select) {
  bim <- read.table(
    path,
    header=FALSE,
    col.names=c("chr","snp","cm","bp","a1","a2"),
    colClasses=c("character","character","integer","integer","factor","factor")
  )
  
  glist <- genes[[ref]]
  gene_ranges <- makeGRangesFromDataFrame(glist, keep.extra.columns=TRUE, ignore.strand=TRUE)

  snp_ranges <- makeGRangesFromDataFrame(
    bim,
    seqnames.field="chr",
    start.field="bp",
    end.field="bp",
    keep.extra.columns=TRUE,
    ignore.strand=TRUE
  )

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
