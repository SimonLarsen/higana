#' Plot part of an ontology as a tree.
#'
#' @param o An \code{ontology} object.
#' @param terms Terms to include in plot.
#' @param include_ancestors Include the ancestors of terms in the plot.
#' @param pvalues A named numeric vector of p-values for each term. Term nodes will be colored according to -log10(p-value) if provided.
#' @param ... Further arguments to \code{\link{ontologyPlot::onto_plot}}.
#' @export
plot_subtree <- function(o, terms, include_ancestors=TRUE, pvalues=NULL, ...) UseMethod("plot_subtree")

#' @importFrom ontologyPlot onto_plot
#' @importFrom circlize colorRamp2
#' @export
plot_subtree.ontology <- function(o, terms, include_ancestors=TRUE, pvalues=NULL, ...) {
  if(class(o) != "ontology") stop("o is not an ontology object.")

  if(include_ancestors) {
    terms <- unique(do.call(c, o$ancestors[terms]))
  }

  fillcolor <- setNames(rep("lightblue", length(terms)), terms)

  if(!is.null(pvalues)) {
    if(is.null(names(pvalues))) stop("p-value vector must be named.")
    pv <- -log10(setNames(pvalues[terms], terms))
    pv[is.na(pv)] <- 0

    cramp <- colorRamp2(quantile(pv, c(0,1), na.rm=TRUE), c("lightblue", "red"))
    fillcolor <- cramp(pv)
  }

  plot.new()
  par(mar=c(0,0,0,0), bg="white")
  onto_plot(o, terms, fillcolor=fillcolor, ...)
}
