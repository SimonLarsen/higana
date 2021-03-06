.pcs_from_term <- function(term) {
  if(length(term$d) > 1) term$u %*% diag(term$d)
  else term$u * term$d
}

#' Perform association test for terms.
#'
#' @param path Path to PC data.
#' @param formula An object of class \code{formula}.
#' @param covars A data frame containing covariates.
#' @param family Error distribution and link function to be used in the model. See \code{\link{glm.fit}} for details.
#' @param Test statistic for model comparison. One of "Chisq", "LRT", "Rao", "F" or "Cp". See \code{\link{stat.anova}}.
#' @param max_pcs Maximum number of PCs to use per term.
#' @return A list with elements:
#'   \describe{
#'     \item{\code{test}}{Data frame of test statistics.}
#'     \item{\code{coef}}{Coefficients from each term model.}
#'   }
#' @importFrom future.apply future_lapply
#' @importFrom fastmatch "%fin%"
#' @importFrom broom tidy
#' @importFrom dplyr bind_rows
#' @importFrom dplyr arrange
#' @export
test_terms <- function(path, formula, covars, family=binomial("logit"), test=NULL, max_pcs=Inf) {
  if(class(formula) != "formula") stop("'formula' is not a formula.")
  if(!("data.frame" %in% class(covars))) stop("'covars' is not a data frame.")
  if(any(grepl("TERM[0-9]+", colnames(covars)))) stop("Covariates named \"TERM[n]\" where [n] is a number are not allowed.")

  if(is.null(test)) {
    if(family$family %in% c("binomial","poisson")) test <- "Chisq"
    else if(family$family %in% c("gaussian","quasibinomial","quasipoisson")) test <- "F"
    else stop("Don't know which test statistic to use for family \"", family$family, "\". Please specify manually.")
  }
  if(is.null(terms)) terms <- names(pc)

  reference.model <- glm(formula, covars, family=family, na.action=na.omit)

  base <- readRDS(path)
  num_parts <- base$num_parts

  out <- lapply(seq_len(num_parts), function(part) {
    pc <- readRDS(paste0(path, ".part", part))

    future_lapply(pc, function(term) {
      x <- .pcs_from_term(term)

      if(ncol(x) > max_pcs) x <- x[,seq_len(max_pcs), drop=FALSE]

      colnames(x) <- paste0("TERM", seq_len(ncol(x)))
      D <- cbind(covars, x, stringsAsFactors=FALSE, row.names=NULL)

      formula.term <- update(formula, paste0("~ . + ", paste0(colnames(x), collapse="+")))
      fit <- glm(formula.term, data=D, family=family, na.action=na.omit, model=FALSE, x=FALSE, y=FALSE)

      list(
        test=anova(reference.model, fit, test=test),
        coef=coefficients(summary(fit))
      )
    })
  })

  out <- unlist(out, recursive=FALSE)

  tests <- lapply(out, `[[`, "test")
  tests <- suppressWarnings(lapply(tests, function(t) tidy(t)[2,]))
  tests <- cbind(term=names(out), bind_rows(tests), stringsAsFactors=FALSE)
  tests <- arrange(tests, p.value)

  list(
    coef=lapply(out, `[[`, "coef"),
    test=tests
  )
}

#' Perform association test for genes.
#'
#' @param path Path to PC data.
#' @param formula An object of class \code{formula}.
#' @param covars A data frame containing covariates.
#' @param family Error distribution and link function to be used in the model. See \code{\link{glm.fit}} for details.
#' @param test Test statistic for model comparison. One of "Chisq", "LRT", "Rao", "F" or "Cp". See \code{\link{stat.anova}}.
#' @param max_pcs Maximum number of PCs to use per gene.
#' @return A list with elements:
#'   \describe{
#'     \item{\code{test}}{Data frame of test statistics.}
#'     \item{\code{coef}}{Coefficients from each term model.}
#'   }
#' @export
test_genes <- function(path, formula, covars, family=binomial("logit"), test=NULL, max_pcs=Inf) {
  test_terms(path, formula, covars, family=family, test=test, max_pcs=max_pcs)
}
