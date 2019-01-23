#' Perform association test for all terms
#'
#' @param formula An object of class \code{formula}.
#' @param covars A data frame containing covariates.
#' @param pc Named list of PCs from terms. Computed with \code{\link{compute_term_pcs}}.
#' @param npcs Number of PCs to use.
#' @importFrom pbapply pblapply
#' @export
test_terms <- function(formula, covars, pc, npcs) {
  if(class(formula) != "formula") stop("formula is not a formula.")
  if(!("data.frame" %in% class(covars))) stop("covars is not a data frame.")
  if(any(grepl("TERM[0-9]+", colnames(covars)))) stop("Covariate names \"TERM[n]\" where [n] is number are not allowed.")
  if(npcs <= 0) stop("Number of principal components must be at least 1.")
  if(any(sapply(pc, is.null))) stop("Principal component vector contains NULL values.")

  fit <- glm(formula, covars, family=binomial("logit"), na.action=na.omit)

  tests <- pblapply(pc, function(term) {
    term <- term[, seq_len(min(ncol(term), npcs))]
    colnames(term) <- paste0("TERM", seq_len(ncol(term)))

    formula.term <- update(formula, paste0("~ . + ", paste0(colnames(term), collapse="+")))

    D <- data.frame(covars, term, stringsAsFactors=FALSE)
    fit2 <- glm(formula.term, D, family=binomial("logit"), na.action=na.omit)

    anova(fit, fit2, test="LRT")
  })

  pvalues <- sapply(tests, function(t) t$`Pr(>Chi)`[2])

  list(test=tests, pvalue=pvalues)
}
