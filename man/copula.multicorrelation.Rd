% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multivariance-functions.R
\name{copula.multicorrelation}
\alias{copula.multicorrelation}
\alias{CMcor}
\title{coupla versions of distance multicorrelation}
\usage{
copula.multicorrelation(x, vec = 1:ncol(x), ...)

CMcor(x, vec = 1:ncol(x), ...)
}
\arguments{
\item{x}{either a data matrix or a list of doubly centered distance matrices}

\item{vec}{if x is a matrix, then this indicates which columns are treated together as one sample; if x is a list, these are the indexes for which the multivariance is calculated. The default is all columns and all indexes, respectively.}

\item{...}{are passed to \code{\link{multicorrelation}}}
}
\description{
Formally it is nothing but distance multicorrelation applied to the Monte Carlo emprical transform of the data. Hence its values vary for repeated runs.
}
\references{
For the theoretic background see the reference [5] given on the main help page of this package: \link{multivariance-package}.
}
\seealso{
\code{\link{multicorrelation}}
}
