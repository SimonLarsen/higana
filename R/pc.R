#' Compute principal components of SNPs annotated to ontology terms.
#'
#' @param o An annotated \code{ontology} object.
#' @param geno A \code{SnpMatrix} genotype matrix.
#' @param genemap A data frame mapping genes to SNP rs numbers.
#' @param npcs Number of principal components to compute per term.
#' @param max_term_size Skip terms annotated with more than this number of genes.
#' @importFrom fastmatch "%fin%"
#' @importFrom flashpcaR flashpca
#' @importFrom pbapply pblapply
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

  terms <- o$id[lengths(o$genes) <= max_term_size]

  pcs <- pblapply(head(setNames(terms, terms), 2000), function(term) tryCatch({
    snps <- as.character(unique(genemap[genemap$gene %fin% o$genes[[term]], "snp"]))
    if(length(snps) == 0) return(NULL)

    x <- as(geno[,snps], "numeric")

    pc <- flashpca(x, npcs, "binom2")
    return(pc$vectors)
  }, error=function(e) return(NULL)))
  pcs
}
