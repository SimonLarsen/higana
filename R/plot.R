#' Plot part of an ontology as a tree.
#'
#' @param path Path to save Graphviz DOT file.
#' @param o An \code{ontology} object.
#' @param terms Character vector of terms to include in the plot.
#' @param test Test result object computed with \code{\link{test_terms}}. Used to color terms according to p-value.
#' @param include_ancestors Include all ancestors of terms in `terms`.
#' @param highlight Character vector of terms to highlight.
#' @param highlight.color Border color for highlighted terms.
#' @param text_wrap Wrap term descriptions at this line width.
#' @param font Name of font to use in nodes.
#' @param color.low Gradient color for least significant terms.
#' @param color.high Gradient color for most significant terms.
#' @return Path to produced graph.
#' @export
plot_hierarchy <- function(path, o, terms, test=NULL, include_ancestors=TRUE, highlight=character(0), highlight.color="black", text_wrap=20, font="Arial", color.low="#add8e6", color.high="#ff0000") UseMethod("plot_hierarchy")

#' @importFrom stringr str_wrap
#' @importFrom stringr str_to_title
#' @importFrom circlize colorRamp2
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 scale_color_gradient
#' @importFrom ggplot2 labs
#' @export
plot_hierarchy <- function(path, o, terms, test=NULL, include_ancestors=TRUE, highlight=character(0), highlight.color="black", text_wrap=20, font="Arial", color.low="#add8e6", color.high="#ff0000") {
  if(class(o) != "ontology") stop("'o' is not an ontology object.")
  if(!is.null(pvalues) && is.null(pvalues)) stop("p-value vector must be named.")

  # include ancestors
  if(include_ancestors) {
    terms <- unique(unlist(o$ancestors[terms]))
  }

  # prepare term labels
  labels <- str_wrap(str_to_title(o$name[terms]), width = text_wrap)
  labels <- gsub("\n", "<BR/>", labels)

  # fill nodes based on p-value if supplied
  if(!is.null(test)) {
    pvalues <- setNames(test$test$p.value, test$test$term)
    pvalues2 <- setNames(pvalues[terms], terms)
    pvalues2[is.na(pvalues2)] <- 1

    color.ramp <- colorRamp2(c(0, -log10(min(pvalues2))), c(color.low, color.high))
    fillcolors <- color.ramp(-log10(pvalues2))
  } else {
    fillcolors <- setNames(rep(color.low, length(terms)), terms)
  }

  con <- file(path, "w")
  write("digraph {", con)
  write(sprintf('  node [shape=none, fixedsize=false, fontname="%s"];', font), con, append=TRUE)
  write('  edge [dir="back"];', con, append=TRUE)

  for(i in seq_along(terms)) {
    s <- sprintf('  %d [label=<
    <TABLE BORDER="0" CELLPADDING="4" CELLSPACING="0">
      <TR><TD BORDER="0" BGCOLOR="black"><FONT COLOR="white">%s</FONT></TD></TR>
      <TR><TD BORDER="%d" COLOR="%s" BGCOLOR="%s">%s</TD></TR>
    </TABLE>
  >];',
      i,
      terms[[i]],
      ifelse(terms[[i]] %in% highlight, 3, 0),
      highlight.color,
      fillcolors[[i]],
      labels[[i]]
    )
    write(s, con, append=TRUE)
  }

  for(i in seq_along(terms)) {
    pa <- intersect(o$parents[[terms[[i]]]], terms)

    pa.trans <- unique(unlist(lapply(pa, function(j) {
      setdiff(o$ancestors[[j]], j)
    })))

    pa <- setdiff(pa, pa.trans)
    pa <- na.omit(match(pa, terms))

    for(j in pa) {
      write(sprintf('  %d -> %d [dir="back"];', j, i), con, append=TRUE)
    }
  }
  write("}", con, append=TRUE)
  close(con)

  # create legend for p-value colors
  legend <- NULL
  if(!is.null(test)) {
    data <- data.frame(p=c(min(pvalues2), 1))
    legend <- ggplot(data, aes(x=1, y=1, color=-log10(p))) +
      geom_point() +
      scale_color_gradient(low=color.low, high=color.high) +
      labs(color="-log10(p-value)")
    legend <- .ggplot_get_legend(legend)

  }
  return(list(path=path, legend=legend))
}

#' @importFrom ggplot2 ggplot_gtable
#' @importFrom ggplot2 ggplot_build
.ggplot_get_legend <- function(p){ 
  tmp <- ggplot_gtable(ggplot_build(p)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  tmp$grobs[[leg]] 
} 
