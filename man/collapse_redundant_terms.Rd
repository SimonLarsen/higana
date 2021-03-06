% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ontology.R
\name{collapse_redundant_terms}
\alias{collapse_redundant_terms}
\title{Collapse terms whose annotations are redundant with respect to their children.}
\usage{
collapse_redundant_terms(o, threshold = 1)
}
\arguments{
\item{o}{An annotated \code{ontology} object.}

\item{threshold}{Filter terms with redundancy (maximum Jaccard index) above this threshold.}
}
\value{
A filtered \code{ontology} object.
}
\description{
Collapses non-root and non-leaf terms in the ontology whose annotations are redundant
with respect to their children. Ensures that the annotations of a term is a proper superset
of the union of annotations of its children.
}
