% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ontology.R
\name{filter}
\alias{filter}
\title{Filter ontology to a specific set of terms.}
\usage{
filter(o, terms, keep_connected = TRUE)
}
\arguments{
\item{o}{An \code{ontology} object.}

\item{terms}{Character vector of terms to keep.}

\item{keep_connected}{Don't remove ancestors of terms in \code{terms}.}
}
\value{
An \code{ontology} object.
}
\description{
Filter ontology to a specific set of terms.
}
