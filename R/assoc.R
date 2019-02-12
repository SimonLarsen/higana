#' Perform association test for terms.
#'
#' @param formula An object of class \code{formula}.
#' @param covars A data frame containing covariates.
#' @param pc Named list of PCs from terms. Computed with \code{\link{compute_term_pcs}}.
#' @param o An \code{ontology} object.
#' @param family Error distribution and link function to be used in the model. See \code{\link{glm.fit}} for details.
#' @param progress Show progress bar?
#' @return A list with elements:
#'   \describe{
#'     \item{\code{test}}{Test objects for each term.}
#'     \item{\code{pvalue}}{p-values for each term.}
#'   }
#' @importFrom pbapply pblapply
#' @importFrom pbmcapply pbmclapply
#' @importFrom fastmatch "%fin%"
#' @export
test_terms <- function(formula, covars, pc, family=binomial("logit"), num_threads=1) {
  if(class(formula) != "formula") stop("formula is not a formula.")
  if(!("data.frame" %in% class(covars))) stop("covars is not a data frame.")
  if(any(grepl("TERM[0-9]+", colnames(covars)))) stop("Covariates named \"TERM[n]\" where [n] is number are not allowed.")
  if(any(grepl("PARENT[0-9]+", colnames(covars)))) stop("Covariates named \"PARENT[n]\" where [n] is number are not allowed.")
  if(any(sapply(pc, is.null))) stop("Principal component vector contains NULL values.")

  reference.model <- glm(formula, covars, family=family, na.action=na.omit)

  out <- pbmclapply(pc, function(term) {
    x <- term$u
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

#' Perform association test for genes.
#'
#' @param formula An object of class \code{formula}.
#' @param covars A data frame containing covariates.
#' @param pc Named list of PCs from genes. Computed with \code{\link{compute_gene_pcs}}.
#' @param family Error distribution and link function to be used in the model. See \code{\link{glm.fit}} for details.
#' @param progress Show progress bar?
#' @return A list with elements:
#'   \item{\code{test}}{Test objects for each gene.}
#'   \item{\code{pvalue}}{p-values for each gene.}
#' @export
test_genes <- function(formula, covars, pc, family=binomial("logit"), num_threads=1) {
  test_terms(formula, covars, pc, family=family, num_threads)
}
