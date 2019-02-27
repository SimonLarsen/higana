#' Perform association test for terms.
#'
#' @param formula An object of class \code{formula}.
#' @param covars A data frame containing covariates.
#' @param pc Named list of PCs from terms. Computed with \code{\link{compute_term_pcs}}.
#' @param family Error distribution and link function to be used in the model. See \code{\link{glm.fit}} for details.
#' @param terms Character vector of terms to test. Defaults to all terms.
#' @param max_pcs Maximum number of PCs to use per term.
#' @param num_threads Number of threads.
#' @return A list with elements:
#'   \describe{
#'     \item{\code{reference}}{Reference model.}
#'     \item{\code{coef}}{Coefficients from each term model.}
#'     \item{\code{tests}}{Test objects for each term.}
#'     \item{\code{pvalue}}{p-values for each term.}
#'   }
#' @importFrom pbapply pblapply
#' @importFrom pbmcapply pbmclapply
#' @importFrom fastmatch "%fin%"
#' @export
test_terms <- function(formula, covars, pc, family=binomial("logit"), terms=NULL, max_pcs=Inf, num_threads=1) {
  if(class(formula) != "formula") stop("formula is not a formula.")
  if(!("data.frame" %in% class(covars))) stop("covars is not a data frame.")
  if(any(grepl("TERM[0-9]+", colnames(covars)))) stop("Covariates named \"TERM[n]\" where [n] is a number are not allowed.")
  if(any(sapply(pc, is.null))) stop("Principal component vector contains NULL values.")

  if(is.null(terms)) terms <- names(pc)

  reference.model <- glm(formula, covars, family=family, na.action=na.omit)

  out <- pbmclapply(pc[terms], function(term) {
    x <- term$u
    if(ncol(x) > max_pcs) x <- x[,seq_len(max_pcs), drop=FALSE]
    colnames(x) <- paste0("TERM", seq_len(ncol(x)))
    D <- cbind(covars, x, stringsAsFactors=FALSE, row.names=NULL)

    formula.term <- update(formula, paste0("~ . + ", paste0(colnames(x), collapse="+")))
    fit <- glm(formula.term, D, family=family, na.action=na.omit, model=FALSE)
    test <- anova(reference.model, fit, test="LRT")
    list(test=test, coef=coefficients(summary(fit)))
  }, mc.cores=num_threads)

  pvalues <- sapply(out, function(t) t$test$`Pr(>Chi)`[2])

  list(
    reference = reference.model,
    coef=lapply(out, `[[`, "coef"),
    tests=lapply(out, `[[`, "test"),
    pvalue=pvalues
  )
}

#' Check terms for significance against children.
#'
#' @param formula An object of class \code{formula}.
#' @param covars A data frame containing covariates.
#' @param pc Named list of PCs from terms. Computed with \code{\link{compute_term_pcs}}.
#' @param o An \code{ontology} object.
#' @param family Error distribution and link function to be used in the model. See \code{\link{glm.fit}} for details.
#' @param terms Character vector of terms to test. Defaults to all terms.
#' @param max_pcs Maximum number of PCs to use per term.
#' @param num_threads Number of threads.
#' @return A list with elements:
#'   \describe{
#'     \item{\code{tests}}{Test objects for each term.}
#'     \item{\code{pvalue}}{p-values for each term.}
#'   }
#' @importFrom fastmatch "%fin%"
#' @importFrom pbmcapply pbmclapply
#' @export
test_terms_children <- function(formula, covars, pc, o, family=binomial("logit"), terms=NULL, max_pcs=Inf, num_threads=1) {
  if(class(formula) != "formula") stop("formula is not a formula.")
  if(!("data.frame" %in% class(covars))) stop("covars is not a data frame.")
  if(any(grepl("TERM[0-9]+", colnames(covars)))) stop("Covariates named \"TERM[n]\" where [n] is a number are not allowed.")
  if(any(grepl("CHILD[0-9]+", colnames(covars)))) stop("Covariates named \"CHILD[n]\" where [n] is a number are not allowed.")
  if(any(sapply(pc, is.null))) stop("Principal component vector contains NULL values.")

  if(is.null(terms)) terms <- names(pc)
  terms <- intersect(terms, go$id[lengths(go$children) > 0])

  out <- pbmclapply(terms, function(term) {
    message(term)
    x.term <- pc[[term]]$u
    if(ncol(x.term) > max_pcs) x.term <- x.term[,seq_len(max_pcs), drop=FALSE]
    colnames(x.term) <- paste0("TERM", seq_len(ncol(x.term)))

    x.child <- lapply(pc[go$children[[term]]], `[[`, "u")
    x.child <- x.child[!sapply(x.child, is.null)]
    x.child <- lapply(x.child, function(y) {
      if(ncol(y) > max_pcs) y <- y[,seq_len(max_pcs), drop=FALSE]
      y
    })
    x.child <- do.call(cbind, x.child)

    if(!is.null(x.child)) {
      colnames(x.child) <- paste0("CHILD", seq_len(ncol(x.child)))
      D <- cbind(covars, x.child, x.term, stringsAsFactors=FALSE, row.names=NULL)
      formula.child <- update(formula, paste0("~ . + ", paste0(colnames(x.child), collapse="+")))
    } else {
      D <- cbind(covars, x.term, stringsAsFactors=FALSE, row.names=NULL)
      formula.child <- formula
    }

    fit <- glm(formula.child, D, family=family, na.action=na.omit, model=FALSE)

    formula.term <-  update(formula.child, paste0("~ . + ", paste0(colnames(x.term), collapse="+")))
    fit2 <- glm(formula.term, D, family=family, na.action=na.omit, model=FALSE)

    anova(fit, fit2, test="LRT")
  }, mc.cores=num_threads)
  names(out) <- paste0(tests$term, ".", tests$child)

  pvalues <- sapply(out, function(t) t$`Pr(>Chi)`[2])

  list(
    tests=out,
    pvalue=pvalues
  )
}

#' Perform association test for genes.
#'
#' @param formula An object of class \code{formula}.
#' @param covars A data frame containing covariates.
#' @param pc Named list of PCs from genes. Computed with \code{\link{compute_gene_pcs}}.
#' @param family Error distribution and link function to be used in the model. See \code{\link{glm.fit}} for details.
#' @param max_pcs Maximum number of PCs to use per gene.
#' @param num_threads Number of threads.
#' @return A list with elements:
#'   \describe{
#'     \item{\code{reference}}{Reference model.}
#'     \item{\code{coef}}{Coefficients from each gene model.}
#'     \item{\code{tests}}{Test objects for each gene}
#'     \item{\code{pvalue}}{p-values for each gene}
#'   }
#' @export
test_genes <- function(formula, covars, pc, family=binomial("logit"), max_pcs=Inf, num_threads=1) {
  test_terms(formula, covars, pc, family=family, max_pcs=max_pcs, num_threads=num_threads)
}
