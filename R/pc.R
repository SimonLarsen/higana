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
  message("Extracting gene SNPs.")
  term_snps <- pblapply(setNames(terms, terms), function(term) {
    as.character(unique(genemap[genemap$gene %fin% o$genes[[term]], "snp"]))
  })
  term_snps <- term_snps[lengths(term_snps) > 0]

  message("Computing PCs.")
  pc <- pblapply(term_snps, function(snps) tryCatch({
    x <- as(geno[,snps], "numeric")
    pc <- flashpca(x, npcs, "binom2", check_geno=FALSE, return_scale=FALSE, do_loadings=FALSE)
    return(pc$vectors)
  }, error=function(e) return(NULL)))

  pc[!sapply(pc, is.null)]
}
