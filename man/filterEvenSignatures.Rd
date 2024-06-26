% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filterBSDB.R
\name{filterEvenSignatures}
\alias{filterEvenSignatures}
\title{Filter even signatures from experiments}
\usage{
filterEvenSignatures(dats, rank, min_size = 10)
}
\arguments{
\item{rank}{Character string. "genus" or "species"}

\item{min_size}{Minimum number of the signatures to be included.}

\item{bsdb}{A BugSigDB dataset from importBugSigDB.}
}
\value{
A list of signatures nested by direction and condition.
}
\description{
\code{filterEvenSignatures} selects only signatures from experiments
in BSDB that have exactly two signatures: increased and decreased.
From case-control studies only.
}
