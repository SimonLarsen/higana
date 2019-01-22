#' Read ontology from OBO file.
#'
#' @param path Path to OBO file.
#' @param relationships Character vector of relationships to include in ontology, e.g. c("is_a", "part_of").
#' @return An \code{ontology} object.
#' @importFrom ontologyIndex get_ontology
#' @export
read_obo <- function(path, relationships) {
  force(relationships)
  o <- get_ontology(path, propagate_relationships=relationships)
  class(o) <- "ontology"
  o
}

#' Read annotations from GAF file.
#'
#' @param path Path to GAF file.
#' @param aspects Character vector of aspects to include. Leave blank to include all aspects.
#' @param exclude_evidence Character vector of evidence types to exclude.
#' @importFrom readr read_delim
#' @export
read_gaf <- function(path, aspects, exclude_evidence) {
  anno <- suppressMessages(read_delim(path, comment="!", quote="", delim="\t", col_names=FALSE))
  if(!missing(aspects)) {
    anno <- anno[anno[[9]] %in% aspects,]
  }
  if(!missing(exclude_evidence)) {
    anno <- anno[!(anno[[7]] %in% exclude_evidence),]
  }
  anno <- unique(anno[,c(3,5)])
  split(anno[[1]], anno[[2]])
}

#' Filter ontology to a specific set of terms.
#'
#' @param o An \code{ontology} object.
#' @param terms Character vector of terms to keep.
#' @param keep_connected Don't remove ancestors of terms in \code{terms}.
#' @return An \code{ontology} object.
#' @importFrom fastmatch "%fin%"
#' @importFrom ontologyIndex get_ancestors
#' @export
filter <- function(o, terms, keep_connected=FALSE) UseMethod("filter")

#' @export
filter.ontology <- function(o, terms, keep_connected=TRUE) {
  if(class(o) != "ontology") stop("o is not an ontology object.")

  if(keep_connected) terms <- get_ancestors(o, terms)

  o <- lapply(o, `[`, terms)
  o$parents <- lapply(o$parents, function(x) x[x %fin% terms])
  o$children <- lapply(o$children, function(x) x[x %fin% terms])
  o$ancestors <- lapply(o$ancestors, function(x) x[x %fin% terms])
  class(o) <- "ontology"
  o
}

#' Removes obsolete terms from the ontology.
#'
#' @param o An \code{ontology} object.
#' @param keep_connected Don't remove ancestors of terms in \code{terms}.
#' @return The filtered ontology.
#' @export
filter_obsolete <- function(o, keep_connected=TRUE) UseMethod("filter_obsolete")

#' @export
filter_obsolete.ontology <- function(o, keep_connected=TRUE) {
  terms <- o$id[!o$obsolete]
  filter.ontology(o, terms)
}

#' Unites a set of roots in the ontology under a common root.
#' Primarily used to combine multiple aspects of an ontology,
#' e.g. 'biological process', 'cellular component' and 'molecular function' in Gene Ontology.
#'
#' @param o An \code{ontology} object.
#' @param roots A character vector of roots to combine.
#' @param root.name Name of the new root term.
#' @return The united \code{ontology} object.
#' @export
unite_roots <- function(o, roots, root.name="root") UseMethod("unite_roots")

#' @export
unite_roots.ontology <- function(o, roots, root.name="root") {
  if(class(o) != "ontology") stop("o is not an ontology object.")
  if(root.name %fin% o$id) {
    stop(sprintf("Root name '%s' already exists in ontology.", root.name))
  }

  for(root in roots) {
    o$parents[[root]] <- c(o$parents[[root]], root.name)
  }
  o$ancestors <- lapply(o$ancestors, c, root.name)

  o$id <- c(o$id, setNames(root.name, root.name))
  o$name <- c(o$name, setNames(root.name, root.name))
  o$parents <- c(o$parents, setNames(list(character(0)), root.name))
  o$children <- c(o$children, setNames(list(roots), root.name))
  o$ancestors <- c(o$ancestors, setNames(list(root.name), root.name))
  o$obsolete <- c(o$obsolete, setNames(FALSE, root.name))
  return(o)
}

#' Annotate ontology with gene annotations.
#'
#' @param o An \code{ontology} object.
#' @param annotations A named list of character vectors.
#' @return An \code{ontology} object with annotations.
#' @export
annotate <- function(o, annotation) UseMethod("annotate")

#' @export
annotate.ontology <- function(o, annotations) {
  if(class(o) != "ontology") stop("o is not an ontology object.")
  if(class(annotations) != "list") stop("annotations is not a list")

  o$genes <- setNames(annotations[o$id], o$id)
  o$genes[sapply(o$genes, is.null)] <- list(character(0))
  return(o)
}

#' Propagates annotations up through the hierarchy.
#' After propagation the set of annotations for a term is a superset of the
#' union of annotations for its children.
#'
#' @param o An \code{ontology} object.
#' @return A propagated \code{ontology} object.
#' @importFrom ontologyIndex get_descendants
#' @export
propagate_annotations <- function(o) UseMethod("propagate_annotations")

#' @export
propagate_annotations.ontology <- function(o) {
  if(class(o) != "ontology") stop("o is not an ontology object.")
  for(term in o$id) {
    desc <- get_descendants(o, term, exclude_roots=TRUE)
    st <- unique(unlist(o$genes[desc]))
    o$genes[[term]] <- union(o$genes[[term]], st)
  }
  return(o)
}

#' Collapse terms whose annotations are redundant with respect to their children.
#'
#' @param o An annotated \code{ontology} object.
#' @return A filtered \code{ontology} object.
#' @export
collapse_redundant_terms <- function(o) UseMethod("collapse_redundant_terms")

#' @export
collapse_redundant_terms.ontology <- function(o) {
  if(class(o) != "ontology") stop("o is not an ontology object.")
  if(is.null(o[["genes"]])) stop("Ontology is not annotated.")

  # Remove unannotated terms
  o <- filter(o, o$id[lengths(o$genes) > 0])

  # Identify redundant terms
  redundant <- sapply(o$id, function(term) {
    any(sapply(o$genes[o$children[[term]]], function(g) setequal(g, o$genes[[term]])))
  })

  # Collapse terms
  for(term in o$id[redundant]) {
    for(parent in o$parents[[term]]) {
      o$children[[parent]] <- union(o$children[[parent]], o$children[[term]])
    }
    for(child in o$children[[term]]) {
      o$parents[[child]] <- union(o$parents[[child]], o$parents[[term]])
    }
  }

  # Remove collapsed terms from ontology
  filter(o, o$id[!redundant])
}

#' @importFrom datastructures stack
#' @importFrom datastructures insert
#' @importFrom datastructures pop
#' @export
topological_order <- function(o, root) UseMethod("topological_order")

#' @export
topological_order.ontology <- function(o, root) {
  s <- stack()
  order <- stack()
  insert(s, root)
  visited <- setNames(rep(FALSE, length(o$id)), o$id)
  while(length(s) > 0) {
    v <- pop(s)
    if(!visited[[v]]) {
      visited[[v]] <- TRUE
      insert(order, v)
      insert(s, o$children[[v]])
    }
  }
}
