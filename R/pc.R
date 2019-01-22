#' Compute principal components of SNPs annotated to ontology terms.
#'
#' @param o An annotated \code{ontology} object.
#' @param geno A \code{SnpMatrix} genotype matrix.
#' @param genemap A data frame mapping genes to SNP rs numbers.
#' @param max_term_size Skip terms annotated with more than this number of genes.
#' @importFrom fastmatch "%fin%"
#' @importFrom flashpcaR flashpca
#' @export
compute_term_pcs <- function(o, geno, genemap, npcs=4, max_term_size=Inf) {
  if(class(o) != "ontology") {
    stop("o is not an ontology object.")
  }
  if(!("data.frame" %in% class(genemap))) {
    stop("genemap is not a data frame.")
  }
  if(!("SnpMatrix" %in% class(geno))) {
    stop("geno is not a SnpMatrix object.")
  }

  lapply(setNames(o$id, o$id), function(term) tryCatch({
    if(length(o$genes[[term]]) > max_term_size) return(NULL)

    snps <- unique(genemap[genemap$gene %fin% o$genes[[term]], "snp"])
    x <- geno[,snps]
    pc <- flashpca(x, npcs, "binom2")
    pc$vectors
  }, error=function(e) NULL))
}
