% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/enrichment.R
\name{runEnrichment}
\alias{runEnrichment}
\title{Run enrichment}
\usage{
runEnrichment(bsdb_sets, background_sets, bp_sigs, bsdb, fdr_ths = NULL)
}
\arguments{
\item{bsdb_sets}{Vector of bsdb sets (named list).}

\item{background_sets}{Vector of background sets (named list).}

\item{bp_sigs}{bugphyzz signatures (named list).}

\item{bsdb}{BSDB data.frame.}

\item{fdr_ths}{Threshold for filtering.}
}
\value{
A data.frame.
}
\description{
\code{runEnrichment} runs the enrichment
}
