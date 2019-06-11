#' Plot part of an ontology as a tree.
#'
#' @param o An \code{ontology} object.
#' @param path Path to save Graphviz DOT file.
#' @param terms Character vector of terms to include in the plot.
#' @param pvalues Named numeric vector of p-values to use for coloring nodes.
#' @param include_ancestors Include all ancestors of terms in `terms`.
#' @param highlight Character vector of terms to highlight.
#' @param text_wrap Wrap term descriptions at this line width.
#' @param font Name of font to use in nodes.
#' @return Path to produced graph.
#' @export
plot_hierarchy <- function(o, path, terms, pvalues=NULL, include_ancestors=TRUE, highlight=character(0), text_wrap=20, font="Arial") UseMethod("plot_hierarchy")

#' @importFrom stringr str_wrap
#' @importFrom stringr str_to_title
#' @importFrom circlize colorRamp2
#' @export
plot_hierarchy <- function(o, path, terms, pvalues=NULL, include_ancestors=TRUE, highlight=character(0), text_wrap=20, font="Arial") {
  if(class(o) != "ontology") stop("'o' is not an ontology object.")
  if(!is.null(pvalues) && is.null(pvalues)) stop("p-value vector must be named.")

  # include ancestors
  if(include_ancestors) {
    terms <- unique(unlist(o$ancestors[terms]))
  }

  # prepare term labels
  labels <- str_wrap(str_to_title(go$name[terms]), width = text_wrap)
  labels <- gsub("\n", "<BR/>", labels)

  # fill nodes based on p-value if supplied
  if(!is.null(pvalues)) {
    pvalues2 <- setNames(pvalues[terms], terms)
    pvalues2[is.na(pvalues2)] <- 1

    color.ramp <- colorRamp2(c(0, -log10(min(pvalues2))), c("#add8e6", "#ff0000"))
    fillcolors <- color.ramp(-log10(pvalues2))
  } else {
    fillcolors <- setNames(rep("#add8e6", length(terms)), terms)
  }

  con <- file(path, "w")
  write("digraph {", con)
  write(sprintf('  node [shape=none, fixedsize=false, fontname="%s"];', font), con, append=TRUE)
  write('  edge [dir="back"];', con, append=TRUE)

  for(i in seq_along(terms)) {
    s <- sprintf('  %d [label=<
    <TABLE BORDER="0" CELLPADDING="4" CELLSPACING="0">
      <TR><TD BORDER="0" BGCOLOR="black"><FONT COLOR="white">%s</FONT></TD></TR>
      <TR><TD BORDER="%d" BGCOLOR="%s">%s</TD></TR>
    </TABLE>
  >];',
      i,
      terms[[i]],
      ifelse(terms[[i]] %in% highlight, 2, 0),
      fillcolors[[i]],
      labels[[i]]
    )
    write(s, con, append=TRUE)
  }

  for(i in seq_along(terms)) {
    pa <- intersect(o$parents[[terms[[i]]]], terms)

    pa.trans <- unique(unlist(lapply(pa, function(j) {
      setdiff(go$ancestors[[j]], j)
    })))

    pa <- setdiff(pa, pa.trans)
    pa <- na.omit(match(pa, terms))

    for(j in pa) {
      write(sprintf('  %d -> %d [dir="back"];', j, i), con, append=TRUE)
    }
  }
  write("}", con, append=TRUE)

  close(con)
  return(path)
}
