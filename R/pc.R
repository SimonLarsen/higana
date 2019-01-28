#' Compute principal components of SNPs annotated to ontology terms.
#'
#' @param o An annotated \code{ontology} object.
#' @param geno A \code{SnpMatrix} genotype matrix.
#' @param genemap A data frame mapping genes to SNP rs numbers.
#' @param npcs Number of principal components to compute per term.
#' @param terms A character vector of terms to compute PCs for. Will use all terms if not provided.
#' @param max_term_size Skip terms annotated with more than this number of genes.
#' @param stand Which standardization method to use. One of "none", "binom" (old Eigenstrat-style), "binom2" (new Eigenstrat-style), "sd" (zero-mean unit-variance) or "center" (zero mean).
#' @param ... Further arguments to \code{\link{flashpcaR::flashpca}}.
#' @importFrom fastmatch "%fin%"
#' @importFrom flashpcaR flashpca
#' @importFrom pbapply pblapply
#' @return A list of matrices where each column corresponds to a PC.
#' @export
compute_term_pcs <- function(o, geno, genemap, npcs=4, terms=NULL, max_term_size=Inf, stand="binom2", ...) {
  if(class(o) != "ontology") {
    stop("o is not an ontology object.")
  }
  if(!("data.frame" %in% class(genemap))) {
    stop("genemap is not a data frame.")
  }
  if(!("SnpMatrix" %in% class(geno))) {
    stop("geno is not a SnpMatrix object.")
  }

  # use all terms if not provided
  if(is.null(terms)) terms <- o$id

  # restrict to terms below max_term_size
  terms <- o$id[lengths(o$genes) <= max_term_size]

  message("Extracting gene SNPs.")
  term_snps <- pblapply(setNames(terms, terms), function(term) {
    as.character(unique(genemap[genemap$gene %fin% o$genes[[term]], "snp"]))
  })
  term_snps <- term_snps[lengths(term_snps) > 0]

  message("Computing PCs.")
  pc <- pblapply(term_snps, function(snps) tryCatch({
    x <- as(geno[,snps], "numeric")
    pc <- flashpca(x, npcs, stand, check_geno=FALSE, return_scale=FALSE, do_loadings=FALSE, ...)
    return(pc$vectors)
  }, error=function(e) return(NULL)))

  pc[!sapply(pc, is.null)]
}

#' Compute principal components of SNPs annotated to genes.
#'
#' @param geno A \code{SnpMatrix} genotype matrix.
#' @param genemap A data frame mapping genes to SNP rs numbers.
#' @param npcs Number of principal components to compute per gene.
#' @param stand Which standardization method to use. One of "none", "binom" (old Eigenstrat-style), "binom2" (new Eigenstrat-style), "sd" (zero-mean unit-variance) or "center" (zero mean).
#' @param ... Further arguments to \code{\link{flashpcaR::flashpca}}.
#' @return A list of matrices where each column corresponds to a PC.
#' @export
compute_gene_pcs <- function(geno, genemap, npcs=4, stand="binom2", ...) {
  if(!("data.frame" %in% class(genemap))) {
    stop("genemap is not a data frame.")
  }
  if(!("SnpMatrix" %in% class(geno))) {
    stop("geno is not a SnpMatrix object.")
  }

  genes <- split(genemap$snp, genemap$gene)

  message("Computing PCs.")
  pc <- pblapply(genes, function(snps) tryCatch({
    x <- as(geno[,snps], "numeric")
    pc <- flashpca(x, npcs, stand, check_geno=FALSE, return_scale=FALSE, do_loadings=FALSE, ...)
    return(pc$vectors)
  }, error=function(e) return(NULL)))

  pc[!sapply(pc, is.null)]
}
