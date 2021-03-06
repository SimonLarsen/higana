% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ontology.R
\name{read_obo}
\alias{read_obo}
\title{Read ontology from OBO file.}
\usage{
read_obo(path, relationships)
}
\arguments{
\item{path}{Path to OBO file.}

\item{relationships}{Character vector of relationships to include in ontology, e.g. c("is_a", "part_of").}
}
\value{
An \code{ontology} object.
}
\description{
Read ontology from OBO file.
}
