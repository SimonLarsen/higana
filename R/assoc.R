#' Perform association test for all terms
#'
#' @param formula An object of class \code{formula}.
#' @param covars A data frame containing covariates.
#' @param pc Named list of PCs from terms. Computed with \code{\link{compute_term_pcs}}.
#' @param npcs Number of PCs to use.
#' @export
test_terms <- function(formula, covars, pc, npcs) {
  if(class(formula) != "formula") stop("formula is not a formula.")
  if(!("data.frame" %in% class(covars))) stop("covars is not a data frame.")
  if(any(grepl("PC[0-9]+", colnames(covars)))) stop("Covariate names \"PC[n]\" where [n] is number are not allowed.")
  if(any(grepl("TERM[0-9]+", colnames(covars)))) stop("Covariate names \"TERM[n]\" where [n] is number are not allowed.")
  if(npcs <= 0) stop("Number of principal components must be at least 1.")

  response.var <- all.vars(update(formula, .~0))
  keep <- !is.na(covars[,response.var,keep=FALSE])

  formula.pc <- update(formula, paste0("~ . + ", paste0("PC", seq_len(npcs), collapse="+")))

  # extract whole genome PCs
  pop <- pcs[["__POPULATION__"]]
  pop <- pop[, seq_len(npcs)]
  colnames(pop) <- paste0("PC", seq_len(ncol(pop)))

  D <- data.frame(covars, pop, stringsAsFactors=FALSE)

  fit <- glm(formula.pc, D[keep,], family=binomial("logit"))

  lapply(pc, function(term) {
    term <- term[, seq_len(min(ncol(term), npcs))]
    colnames(term) <- paste0("TERM", seq_len(ncol(term)))

    formula.term <- update(formula.pc, paste0("~ . + ", paste0(colnames(term), collapse="+")))
    D2 <- data.frame(D, term, stringsAsFactors=FALSE)
    fit2 <- glm(formula.term, D2[keep,], family=binomial("logit"))
    anova(fit, fit2, test="LRT")
  })
}
