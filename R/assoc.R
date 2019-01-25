#' Perform association test for all terms
#'
#' @param formula An object of class \code{formula}.
#' @param covars A data frame containing covariates.
#' @param pc Named list of PCs from terms. Computed with \code{\link{compute_term_pcs}}.
#' @param o An \code{ontology} object.
#' @param npcs Number of PCs to use. Set of \code{Inf} to use all available.
#' @param against Model to test against for significance. One of "none", "parent".
#' @return A list with elements:
#'   \item{\code{test}}{Test objects for each term.}
#'   \item{\code{pvalue}}{p-values for each term.}
#' @importFrom pbapply pblapply
#' @importFrom fastmatch "%fin%"
#' @export
test_terms <- function(formula, covars, pc, o=NULL, npcs=Inf, against="none") {
  if(class(formula) != "formula") stop("formula is not a formula.")
  if(!("data.frame" %in% class(covars))) stop("covars is not a data frame.")
  if(any(grepl("TERM[0-9]+", colnames(covars)))) stop("Covariates named \"TERM[n]\" where [n] is number are not allowed.")
  if(any(grepl("PARENT[0-9]+", colnames(covars)))) stop("Covariates named \"PARENT[n]\" where [n] is number are not allowed.")
  if(npcs <= 0) stop("Number of principal components must be at least 1.")
  if(any(sapply(pc, is.null))) stop("Principal component vector contains NULL values.")

  against <- match.arg(against, c("none","parent"))

  if(against == "none") {
    fit <- glm(formula, covars, family=binomial("logit"), na.action=na.omit)

    tests <- pblapply(pc, function(term) {
      term <- term[, seq_len(min(ncol(term), npcs))]
      colnames(term) <- paste0("TERM", seq_len(ncol(term)))

      formula.term <- update(formula, paste0("~ . + ", paste0(colnames(term), collapse="+")))

      D <- data.frame(covars, term, stringsAsFactors=FALSE)
      fit2 <- glm(formula.term, D, family=binomial("logit"), na.action=na.omit)

      anova(fit, fit2, test="LRT")
    })
  }
  else if(against == "parent") {
    if(is.null(o)) {
      stop("Ontology 'o' must be supplied when testing against parent.")
    }

    edges <- data.frame(parent=rep(o$id, lengths(o$children)), child=unlist(o$children), stringsAsFactors=FALSE)
    edges <- subset(edges, parent %fin% names(pc) & child %fin% names(pc))

    tests <- pblapply(seq_len(nrow(edges)), function(i) {
      pc.parent <- pc[[edges[i,"parent"]]]
      pc.parent <- pc.parent[, seq_len(min(ncol(pc.parent), npcs))]
      colnames(pc.parent) <- paste0("PARENT", seq_len(ncol(pc.parent)))

      pc.child <- pc[[edges[i,"child"]]]
      pc.child <- pc.child[, seq_len(min(ncol(pc.child), npcs))]
      colnames(pc.child) <- paste0("TERM", seq_len(ncol(pc.child)))

      formula.parent <- update(formula, paste0("~ . + ", paste0(colnames(pc.parent), collapse="+")))
      D <- data.frame(covars, pc.parent, stringsAsFactors=FALSE)
      fit <- glm(formula.parent, D, family=binomial("logit"), na.action=na.omit)

      formula.term <- update(formula.parent, paste0("~ . + ", paste0(colnames(pc.child), collapse="+")))
      D <- data.frame(D, pc.child, stringsAsFactors=FALSE)
      fit2 <- glm(formula.term, D, family=binomial("logit"), na.action=na.omit)

      anova(fit, fit2, test="LRT")
    })
    names(tests) <- apply(edges, 1, paste0, collapse=".")
  }

  pvalues <- sapply(tests, function(t) t$`Pr(>Chi)`[2])
  return(list(test=tests, pvalue=pvalues))
}
