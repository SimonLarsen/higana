#' Compute principal components of SNPs annotated to ontology terms.
#'
#' @param path Path to output file.
#' @param o An annotated \code{ontology} object.
#' @param data A \code{BEDData} object
#' @param stand Which standardization method to use. One of "none", "binom" (old Eigenstrat-style), "binom2" (new Eigenstrat-style), "sd" (zero-mean unit-variance) or "center" (zero mean).
#' @param terms A character vector of terms to compute PCs for. Will use all terms if not provided.
#' @param max_term_size Skip terms annotated with more than this number of genes.
#' @param max_term_size Skip terms annotated with fewer than this number of genes.
#' @param explain_var Restrict number of PCs to explain at least this fraction of variance.
#' @param max_pcs Return at most this number of PCs.
#' @param num_parts Number of parts to write output to.
#' @param exclude_snps Named list of SNP IDs to exclude from terms.
#' @param rsvd_threshold Use randomized SVD when number of variables exceeds this threshold. Set to \code{Inf} to turn off.
#' @importFrom fastmatch fmatch
#' @export
compute_term_pcs <- function(
    path,
    o,
    data,
    terms=NULL,
    stand="binom2",
    max_term_size=75,
    min_term_size=3,
    explain_var=0.95,
    max_pcs=50,
    num_parts=1,
    exclude_snps=NULL,
    rsvd_threshold=Inf
) {
  if(class(o) != "ontology") stop("'o' is not an ontology object.")
  if(class(data) != "BEDData") stop("'data' is not a BEDData object.")

  # use all terms if not provided
  if(is.null(terms)) terms <- o$id

  # restrict to terms below max_term_size
  terms <- intersect(terms, o$id[lengths(o$genes) <= max_term_size])
  # restrict to terms above min_term_size
  terms <- intersect(terms, o$id[lengths(o$genes) >= min_term_size])

  # Extract SNPs for each term
  if(interactive()) message("Extracting gene SNPs.")
  gene2snps <- split(data$map$snp, data$map$gene)
  term_snps <- lapply(setNames(terms, terms), function(term) {
    ind <- fmatch(o$genes[[term]], names(gene2snps))
    snps <- as.character(unique(unlist(gene2snps[ind])))
    setdiff(snps, exclude_snps[[term]])
  })
  term_snps <- term_snps[lengths(term_snps) > 0]

  # Separate terms into parts
  term_snps <- .split_into_parts(term_snps, num_parts)

  out <- list(num_parts=num_parts)
  out$terms <- lapply(seq_len(num_parts), function(part) {
    if(interactive()) message(sprintf("Computing PCs (%d/%d).", part, num_parts))

    pc <- .compute_snp_pcs(data, term_snps[[part]], stand, explain_var, max_pcs, rsvd_threshold)
    pc <- pc[!sapply(pc, is.null)]
    saveRDS(pc, paste0(path, ".part", part))
    names(pc)
  })

  saveRDS(out, path)
}

#' Compute principal components of ontology terms for stepdown test.
#'
#' @param path Path to output file.
#' @param o An annotated \code{ontology} object.
#' @param data A \code{BEDData} object.
#' @param result An association result computed with \code{\link{test_terms}}.
#' @param terms A character vector of terms to compute PCs for. Will use all terms if not provided.
#' @param stepdown Step-down method. One of "top.child", "top.gene".
#' @param ... Further arguments to \code{\link{compute_term_pcs}}.
#' @importFrom fastmatch "%fin%"
#' @export
compute_term_pcs_stepdown <- function(
  path,
  o,
  data,
  result,
  terms=NULL,
  stepdown="top.child",
  ...
) {
  if(class(o) != "ontology") step("'o' is not an ontology object.")
  if(!("BEDData" %in% class(data))) stop("'data' is not a BEDData object.")
  if(!("data.frame" %in% class(result$test))) stop("'result' does not contain a test statistics data frame.")
  if(!(all(c("term","p.value") %in% colnames(result$test)))) stop("malformed test statistics data frame in 'result'.")

  if(is.null(terms)) terms <- o$id
  stepdown = match.arg(stepdown, c("top.child","top.gene"))

  pvalues <- setNames(result$test$p.value, result$test$term)

  if(stepdown == "top.child") {
    exc.snps <- lapply(setNames(terms, terms), function(term) {
      top.index <- which.min(pvalues[o$children[[term]]])
      if(length(top.index) == 0) return(character(0))
      top.child <- o$children[[term]][[top.index]]
      top.genes <- o$genes[[top.child]]
      as.character(unique(data$map[data$map$gene %fin% top.genes, "snp"]))
    })
  }
  if(stepdown == "top.gene") {
    exc.snps <- lapply(setNames(terms, terms), function(term) {
      top.index <- which.min(pvalues[o$genes[[term]]])
      if(length(top.index) == 0) return(character(0))
      top.gene <- o$genes[[term]][[top.index]]
      as.character(unique(data$map[data$map$gene == top.gene, "snp"]))
    })
  }

  compute_term_pcs(path, o, data, terms=terms, exclude_snps=exc.snps, ...)
}

#' Compute principal components of SNPs annotated to genes.
#'
#' @param path Path to output file.
#' @param data A \code{BEDData} object. 
#' @param stand Which standardization method to use. One of "none", "binom" (old Eigenstrat-style), "binom2" (new Eigenstrat-style), "sd" (zero-mean unit-variance) or "center" (zero mean).
#' @param explain_var Restrict number of PCs to explain at least this fraction of variance.
#' @param max_pcs Return at most this number of PCs.
#' @param num_parts Number of parts to write output to.
#' @param rsvd_threshold Use randomized SVD when number of variables exceeds this threshold.
#' @return A list of singular value decompositions for each gene Each entry is a list with elements
#'   \describe{
#'     \item{\code{u}}{The left-singular vectors.}
#'     \item{\code{d}}{The non-zero singular values.}
#'     \item{\code{v}}{The right-singular vectors.}
#'   }
#' @export
compute_gene_pcs <- function(
  path,
  data,
  stand="binom2",
  explain_var=0.95,
  max_pcs=50,
  num_parts=1,
  rsvd_threshold=Inf
) {
  if(class(data) != "BEDData") stop("'data' is not a BEDData object.")

  gene_snps <- split(data$map$snp, data$map$gene)
  gene_snps <- .split_into_parts(gene_snps, num_parts)

  out <- list(num_parts=num_parts)
  out$terms <- lapply(seq_len(num_parts), function(part) {
    if(interactive()) message(sprintf("Computing PCs (%d/%d).", part, num_parts))

    pc <- .compute_snp_pcs(data, gene_snps[[part]], stand, explain_var, max_pcs, rsvd_threshold)
    pc <- pc[!sapply(pc, is.null)]
    saveRDS(pc, paste0(path, ".part", part))
    names(pc)
  })

  saveRDS(out, path)
}

#' @importFrom future.apply future_lapply
.compute_snp_pcs <- function(data, term_snps, stand, explain_var, max_pcs, rsvd_threshold) {
  future_lapply(term_snps, function(snps) tryCatch({
    x <- data$snps[,snps]
    cv <- apply(x, 2, var, na.rm=TRUE)
    x <- x[, cv > 1e-5, drop=FALSE]
    x <- scale2(x, stand)
    if(ncol(x) <= rsvd_threshold) {
      .get_svd(x, explain_var, max_pcs)
    } else {
      .get_rsvd(x, explain_var, max_pcs)
    }
  }, error=function(e) return(NULL)))
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
  sv$scale <- attr(x, "scaled:scale")
  sv$center <- attr(x, "scaled:center")
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
  sv$scale <- attr(x, "scaled:scale")
  sv$center <- attr(x, "scaled:center")
  sv
}

.split_into_parts <- function(sets, num_parts) {
  part <- rep(seq_len(num_parts), length.out=length(sets))
  part <- sample(part)
  split(sets, part)
}
