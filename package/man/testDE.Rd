% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/testDE.R
\name{testDE}
\alias{testDE}
\title{Test for DE}
\usage{
testDE(counts, sf, group1, group2, lfc = 1)
}
\arguments{
\item{counts}{Numeric matrix, count data for genes (in rows) and cells/samples (in columns).}

\item{sf}{Numeric vector, size factors for all cells.}

\item{group1}{Integer vector, identities of cells in the first group.}

\item{group2}{Integer vector, identities of cells in the second group.}

\item{lfc}{Numeric scalar, threshold on the log-fold change for \code{\link{glmTreat}}.}
}
\value{
A \code{\link{DGELRT}} object with the results of the comparison between groups.
Note that log-fold changes refer to the second group over the first group.
}
\description{
Use \pkg{edgeR} to test for differential expression between two groups.
}
\details{
This function uses the TREAT methods in \pkg{edgeR} to test for differential expression between groups.
The idea is to identify changes in the detected DEGs when the size factors change, due to more-or-less effective removal of composition biases between groups.
}
\examples{
counts <- matrix(rpois(100000, lambda=2), ncol=100)
sf <- scater::librarySizeFactors(counts)
testDE(counts, sf, 1:50, 51:100)
}
\author{
Aaron Lun
}
