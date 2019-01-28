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
#' @export
read_gaf <- function(path, aspects, exclude_evidence) {
  anno <- read.table(path, header=FALSE, comment.char="!", stringsAsFactors=FALSE, sep="\t", quote="")
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
#' @param keep_connected Keep obsolete terms that are ancestors of non-obsolete terms.
#' @return An ontology with all obsolete terms removed.
#' @export
filter_obsolete <- function(o, keep_connected=TRUE) UseMethod("filter_obsolete")

#' @export
filter_obsolete.ontology <- function(o, keep_connected=TRUE) {
  if(class(o) != "ontology") stop("o is not an ontology object.")
  terms <- o$id[!o$obsolete]
  filter.ontology(o, terms, keep_connected=TRUE)
}

#' Removes unannotated terms from the ontology.
#'
#' @param o An annotated \code{ontology} object.
#' @param keep_connected Keep unannotated terms that are ancestors of annotated terms.
#' @return An ontology with all unannotated terms removed.
#' @export
filter_unannotated <- function(o, keep_connected=TRUE) UseMethod("filter_unannotated")

#' @export
filter_unannotated.ontology <- function(o, keep_connected=TRUE) {
  if(class(o) != "ontology") stop("o is not an ontology object.")
  if(is.null(o[["genes"]])) stop("Ontology is not annotated.")
  terms <- o$id[lengths(o$genes) > 0]
  filter.ontology(o, terms, keep_connected=TRUE)
}

#' Unites a set of roots in the ontology under a common root.
#'
#' Primarily used to combine multiple aspects of an ontology,
#' e.g. 'biological process', 'cellular component' and 'molecular function' in Gene Ontology.
#'
#' @param o An \code{ontology} object.
#' @param roots A character vector of roots to combine. Will unite all roots in ontology if not provided.
#' @param root.name Name of the new root term.
#' @return The united \code{ontology} object.
#' @export
unite_roots <- function(o, roots, root.name="root") UseMethod("unite_roots")

#' @export
unite_roots.ontology <- function(o, roots=NULL, root.name="root") {
  if(class(o) != "ontology") stop("o is not an ontology object.")
  if(root.name %fin% o$id) {
    stop(sprintf("Root name '%s' already exists in ontology.", root.name))
  }

  if(is.null(roots)) {
    roots <- unname(o$id[lengths(o$parents) == 0])
    message("Uniting roots ", paste0(roots, collapse=", "))
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
#' @importFrom fastmatch fmatch
#' @export
annotate <- function(o, annotation) UseMethod("annotate")

#' @export
annotate.ontology <- function(o, annotations) {
  if(class(o) != "ontology") stop("o is not an ontology object.")
  if(class(annotations) != "list") stop("annotations is not a list")

  indices <- fmatch(o$id, names(annotations))
  o$genes <- setNames(annotations[indices], o$id)
  o$genes[sapply(o$genes, is.null)] <- list(character(0))
  return(o)
}

#' Propagates annotations up through the hierarchy.
#'
#' After propagation the set of annotations for a term is a superset of the
#' union of annotations of its children.
#'
#' @param o An \code{ontology} object.
#' @return A propagated \code{ontology} object.
#' @importFrom ontologyIndex get_descendants
#' @importFrom igraph graph_from_data_frame
#' @importFrom igraph topo_sort
#' @importFrom igraph as_ids
#' @importFrom fastmatch fmatch
#' @export
propagate_annotations <- function(o) UseMethod("propagate_annotations")

#' @export
propagate_annotations.ontology <- function(o) {
  if(class(o) != "ontology") stop("o is not an ontology object.")
  if(is.null(o[["genes"]])) stop("Ontology is not annotated.")

  edges <- data.frame(from=rep(o$id, lengths(o$children)), to=unlist(o$children), stringsAsFactors=FALSE)
  g <- graph_from_data_frame(edges, directed=TRUE)
  term_order <- rev(as_ids(topo_sort(g, "out")))
  term_i <- fmatch(term_order, o$id)

  for(term in term_i) {
    ch <- o$children[[term]]
    if(length(ch) > 0) {
      st <- unique(do.call(c, o$genes[ch]))
      o$genes[[term]] <- base::union(o$genes[[term]], st)
    }
  }
  return(o)
}

#' Recomputes \code{ancestors} field for each term.
#'
#' This method can be used to repair the ontology after changing the \code{children} or \code{parents} fields.
#'
#' @param o An \code{ontology} object.
#' @importFrom igraph graph_from_data_frame
#' @importFrom igraph topo_sort
#' @importFrom igraph as_ids
#' @importFrom fastmatch fmatch
#' @export
recompute_ancestors <- function(o) UseMethod("recompute_ancestors")

#' @export
recompute_ancestors.ontology <- function(o) {
  if(class(o) != "ontology") stop("o is not an ontology object.")

  edges <- data.frame(from=rep(o$id, lengths(o$children)), to=unlist(o$children), stringsAsFactors=FALSE)
  g <- graph_from_data_frame(edges, directed=TRUE)
  term_order <- as_ids(topo_sort(g, "out"))
  term_i <- fmatch(term_order, o$id)

  o$ancestors <- setNames(rep(list(character(0)), length(o$id)), o$id)
  for(term in term_order) {
    anc <- unique(unname(do.call(c, o$ancestors[fmatch(o$parents[[term]], o$id)])))
    o$ancestors[[fmatch(term, o$id)]] <- c(o$id[fmatch(term, o$id)], anc)
  }
  o$ancestors <- o$ancestors[fmatch(o$id, names(o$ancestors))]
  return(o)
}

#' Collapse terms whose annotations are redundant with respect to their children.
#'
#' Collapses non-root and non-leaf terms in the ontology whose annotations are redundant
#' with respect to their children. Ensures that the annotations of a term is a proper superset
#' of the union of annotations of its children.
#'
#' @param o An annotated \code{ontology} object.
#' @return A filtered \code{ontology} object.
#' @importFrom fastmatch fmatch
#' @export
collapse_redundant_terms <- function(o) UseMethod("collapse_redundant_terms")

#' @export
collapse_redundant_terms.ontology <- function(o) {
  if(class(o) != "ontology") stop("o is not an ontology object.")
  if(is.null(o[["genes"]])) stop("Ontology is not annotated.")

  # Remove unannotated terms
  o <- filter.ontology(o, o$id[lengths(o$genes) > 0], keep_connected=FALSE)

  # Identify redundant terms
  setequal2 <- function(x, y) {
    length(x) == length(y) && !anyNA(fmatch(x, y)) && !anyNA(fmatch(x, y))
  }
  terms_i <- which(lengths(o$children) > 0 & lengths(o$parents) > 0)
  redundant <- sapply(terms_i, function(term) {
    ch_i <- fmatch(o$children[[term]], o$id)
    any(sapply(o$genes[ch_i], function(g) setequal2(g, o$genes[[term]])))
  })

  # Collapse terms
  for(term in which(redundant)) {
    for(parent in fmatch(o$parents[[term]], o$id)) {
      o$children[[parent]] <- base::union(o$children[[parent]], o$children[[term]])
    }
    for(child in fmatch(o$children[[term]], o$id)) {
      o$parents[[child]] <- base::union(o$parents[[child]], o$parents[[term]])
    }
  }

  # Remove collapsed terms from ontology
  filter.ontology(o, o$id[!redundant], keep_connected=FALSE)
}

#' Permute an ontology.
#'
#' @param An \code{ontology} object.
#' @param method Permutation method.
#'     "genes" assigns random genes to each term while preserving the number of annotations per term
#'     and total occurences of each gene.
#'     "sets" shuffle whole annotate sets among terms with similar number of annotations.
#' @param recompute_ancestors Should \code{ancestors} field be recomputed after perturbation?
#' @param bins Number of bins to use for method="sets".
#' @importFrom igraph graph_from_data_frame
#' @importFrom igraph rewire
#' @importFrom igraph keeping_degseq
#' @export
permute <- function(o, method="genes", recompute_ancestors=TRUE, bins=100) UseMethod("permute")

#' @export
permute.ontology <- function(o, method="genes", recompute_ancestors=TRUE, bins=100) {
  if(class(o) != "ontology") stop("o is not an ontology object.")

  method <- match.arg(method, c("genes","sets"))

  if(method == "genes") {
    m <- data.frame(term=rep(o$id, lengths(o$genes)), gene=unlist(o$genes), stringsAsFactors=FALSE, row.names=NULL)
    m$gene <- sample(m$gene)
    anno <- split(m$gene, m$term)
    return(annotate(o, anno))
  }
  else if(method == "sets") {
    sizes <- data.frame(term=o$id, genes=lengths(o$genes), stringsAsFactors=FALSE, row.names=FALSE)
    sizes <- sizes[order(sizes$genes),]
    sizes$bin <- floor(seq(1, bins+1-1/nrow(sizes), length.out=nrow(sizes)))
    term.bins <- split(sizes$term, sizes$bin)
    for(bin in term.bins) {
      if(max(lengths(o$genes[bin])) > 0) {
        o$genes[bin] <- setNames(o$genes[sample(bin)], bin)
      }
    }
    return(o)
  }
}
