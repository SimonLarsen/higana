% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pc.R
\name{compute_term_pcs_stepdown}
\alias{compute_term_pcs_stepdown}
\title{Compute principal components of ontology terms for stepdown test.}
\usage{
compute_term_pcs_stepdown(
  path,
  o,
  data,
  result,
  terms = NULL,
  stepdown = "top.child",
  ...
)
}
\arguments{
\item{path}{Path to output file.}

\item{o}{An annotated \code{ontology} object.}

\item{data}{A \code{BEDData} object.}

\item{result}{An association result computed with \code{\link{test_terms}} or \code{\link{test_genes}}.}

\item{terms}{A character vector of terms to compute PCs for. Will use all terms if not provided.}

\item{stepdown}{Step-down method. One of "top.child", "top.gene".}

\item{...}{Further arguments to \code{\link{compute_term_pcs}}.}
}
\description{
Compute principal components of ontology terms for stepdown test.
}
