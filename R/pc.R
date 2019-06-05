#' Compute principal components of SNPs annotated to ontology terms.
#'
#' @param o An annotated \code{ontology} object.
#' @param geno A \code{SnpMatrix} genotype matrix.
#' @param genemap A data frame mapping genes to SNP rs numbers.
#' @param stand Which standardization method to use. One of "none", "binom" (old Eigenstrat-style), "binom2" (new Eigenstrat-style), "sd" (zero-mean unit-variance) or "center" (zero mean).
#' @param terms A character vector of terms to compute PCs for. Will use all terms if not provided.
#' @param max_term_size Skip terms annotated with more than this number of genes.
#' @param max_term_size Skip terms annotated with fewer than this number of genes.
#' @param explain_var Restrict number of PCs to explain at least this fraction of variance.
#' @param max_pcs Return at most this number of PCs.
#' @param rsvd_threshold Use randomized SVD when number of variables exceeds this threshold.
#' @importFrom fastmatch "%fin%"
#' @importFrom future.apply future_lapply
#' @return A list of singular value decompositions for each term. Each entry is a list with elements
#'   \describe{
#'     \item{\code{u}}{The left-singular vectors.}
#'     \item{\code{d}}{The non-zero singular values.}
#'     \item{\code{v}}{The right-singular vectors.}
#'   }
#' @export
compute_term_pcs <- function(o, geno, genemap, stand="binom2", terms=NULL, max_term_size=75, min_term_size=3, explain_var=0.95, max_pcs=50, rsvd_threshold=0) {
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
  terms <- intersect(terms, o$id[lengths(o$genes) <= max_term_size])
  # restrict to terms above min_term_size
  terms <- intersect(terms, o$id[lengths(o$genes) >= min_term_size])

  if(interactive()) message("Extracting gene SNPs.")
  term_snps <- lapply(setNames(terms, terms), function(term) {
    as.character(unique(genemap[genemap$gene %fin% o$genes[[term]], "snp"]))
  })
  term_snps <- term_snps[lengths(term_snps) > 0]

  if(interactive()) message("Computing PCs.")
  pc <- future_lapply(term_snps, function(snps) tryCatch({
    x <- as(geno[,snps], "numeric")
    cv <- apply(x, 2, var, na.rm=TRUE)
    x <- x[, cv > 1e-5, drop=FALSE]
    if(stand != "none") x <- scale2(x, stand)
    if(ncol(x) <= rsvd_threshold) {
      .get_svd(x, explain_var, max_pcs)
    } else {
      .get_rsvd(x, explain_var, max_pcs)
    }
  }, error=function(e) return(NULL)))

  pc[!sapply(pc, is.null)]
}

#' Compute principal components of SNPs annotated to genes.
#'
#' @param geno A \code{SnpMatrix} genotype matrix.
#' @param genemap A data frame mapping genes to SNP rs numbers.
#' @param stand Which standardization method to use. One of "none", "binom" (old Eigenstrat-style), "binom2" (new Eigenstrat-style), "sd" (zero-mean unit-variance) or "center" (zero mean).
#' @param explain_var Restrict number of PCs to explain at least this fraction of variance.
#' @param max_pcs Return at most this number of PCs.
#' @param rsvd_threshold Use randomized SVD when number of variables exceeds this threshold.
#' @return A list of singular value decompositions for each gene Each entry is a list with elements
#'   \describe{
#'     \item{\code{u}}{The left-singular vectors.}
#'     \item{\code{d}}{The non-zero singular values.}
#'     \item{\code{v}}{The right-singular vectors.}
#'   }
#' @importFrom future.apply future_lapply
#' @export
compute_gene_pcs <- function(geno, genemap, stand="binom2", explain_var=1, max_pcs=Inf, rsvd_threshold=0) {
  if(!("data.frame" %in% class(genemap))) {
    stop("genemap is not a data frame.")
  }
  if(!("SnpMatrix" %in% class(geno))) {
    stop("geno is not a SnpMatrix object.")
  }

  genes <- split(genemap$snp, genemap$gene)

  if(interactive()) message("Computing PCs.")
  pc <- future_lapply(genes, function(snps) tryCatch({
    x <- as(geno[,snps], "numeric")
    cv <- apply(x, 2, var, na.rm=TRUE)
    x <- x[, cv > 1e-5, drop=FALSE]
    if(stand != "none") x <- scale2(x, stand)
    if(ncol(x) <= rsvd_threshold) {
      .get_svd(x, explain_var, max_pcs)
    } else {
      .get_rsvd(x, explain_var, max_pcs)
    }
  }, error=function(e) return(NULL)))

  pc[!sapply(pc, is.null)]
}

#' @importFrom corpcor fast.svd
.get_svd <- function(x, explain_var, max_pcs) {
  sv <- fast.svd(x)

  rownames(sv$v) <- colnames(x)

  eig <- sv$d^2
  ex <- cumsum(eig) / sum(eig)
  numpc <- sum(ex < explain_var)+1
  numpc <- min(numpc, max_pcs, length(sv$d))

  sv$d <- head(sv$d, numpc)
  sv$u <- sv$u[, seq_len(numpc), drop=FALSE]
  sv$v <- sv$v[, seq_len(numpc), drop=FALSE]
  sv
}

.get_rsvd <- function(x, explain_var, max_pcs) {
  sv <- rsvd(x, max_pcs)

  rownames(sv$u) <- rownames(x)
  rownames(sv$v) <- colnames(x)

  eig <- sv$d^2
  ex <- cumsum(eig) / sum(eig)
  numpc <- sum(ex < explain_var)+1
  numpc <- min(numpc, max_pcs, length(sv$d))

  sv$d <- head(sv$d, numpc)
  sv$u <- sv$u[, seq_len(numpc), drop=FALSE]
  sv$v <- sv$v[, seq_len(numpc), drop=FALSE]
  sv
}
