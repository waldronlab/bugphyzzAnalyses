% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{limmaVoom}
\alias{limmaVoom}
\title{Normalize assay with limma::voom}
\usage{
limmaVoom(tse)
}
\arguments{
\item{tse}{A (Tree)SummarizedExperiment.}
}
\value{
A (Tree)SummarizedExperiment.
}
\description{
\code{limmaVoom} normalizes an assay (the first one) from a
(Tree)SummarizedExperiment with lima::voom.
}
