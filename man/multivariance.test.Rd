% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multivariance-functions.R
\name{multivariance.test}
\alias{multivariance.test}
\title{independence tests based on (total-/2-/3-) multivariance}
\usage{
multivariance.test(
  x,
  vec = 1:ncol(x),
  type = "total",
  p.value.type = "pearson_approx",
  verbose = TRUE,
  ...
)
}
\arguments{
\item{x}{matrix, each row is a sample}

\item{vec}{vector which indicates which columns are treated as one sample}

\item{type}{one of \code{"independence"}, \code{"pairwise independence"}, \code{"multi"}, \code{"total"}, \code{"m.multi.2"}, \code{"m.multi.3"}}

\item{p.value.type}{one of \code{"pearson_approx"}, \code{"distribution_free"}, \code{"resample"}, \code{"pearson_unif"}}

\item{verbose}{logical, if TRUE meaningful text output is generated.}

\item{...}{these are passed to \code{\link{cdm}}}
}
\value{
A list with class "\code{htest}" containing the following components:
\describe{
  \item{\code{statistic}}{the value of the test statistic,}
  \item{\code{p.value}}{the p-value of the test statistic,}
  \item{\code{method}}{a character string indicating the type of test performed,}
  \item{\code{data.name}}{a character string giving the name(s) of the data.}
}
}
\description{
This performs the (specified by \code{type} and \code{p.value.type}) independence test for the columns of a sample matrix.
}
\details{
For the use of \code{vec} see the examples below and the more detailed explanation of this argument for \code{\link{multivariance}}.

The types \code{"independence"} and \code{"total"} are identical: an independence test is performed.

Also the types \code{"pairwise independence"} and \code{"m.multi.2"} are identical:  a test of pairwise independence is performed.

The type \code{"m.multi.3"}, performs a test for 3-independence, assuming pairwise independence. The type  \code{"multi"} performs a test for n-independence, assuming (n-1)-independence.

There are several ways (determined by \code{p.value.type}) to estimate the p-value: The \code{"pearson_approx"} and \code{"resample"} are approximately sharp. The latter is based on a resampling approach and thus much slower. The \code{"distribution_free"} test might be very conservative, its p-value estimates are only valid for p-values lower than 0.215 - values above should be interpreted as "values larger than 0.215". Finally, \code{"pearson_unif"} uses fixed parameters in Pearson's estimate, it is only applicable for univariate uniformly distributed marginals

All tests are performed using the standard euclidean distance. Other distances can be supplied via the \code{...}, see \code{\link{cdm}} for the accepted arguments.
}
\examples{
# an independence test
multivariance.test(dep_struct_several_26_100,p.value.type = "distribution_free") # conservative
multivariance.test(dep_struct_several_26_100,p.value.type = "resample") #sharp but slow
multivariance.test(dep_struct_several_26_100,p.value.type = "pearson_approx") #

# as an example, all tests for one data set:
coins100 = coins(100)
for (ty in c("total","m.multi.2","m.multi.3","multi"))
 for (pvt in c("distribution_free","resample","pearson_approx"))
  print(multivariance.test(coins100,type=ty,p.value.type = pvt))

# using the vec argument:
x = matrix(rnorm(50*6),ncol = 10) # a 50x6 data matrix
vec = c(1,2,3,4,5,6) # each column is treated as one variable
multivariance.test(x,vec,p.value.type = "distribution_free") # is the same as the default

vec = c(1,2,2,1,3,1)
# column 1,4,6 are treated as one variable
# column 2,3 are treated as one variable
# column 5 is treated as one variable
multivariance.test(x,vec,p.value.type = "distribution_free")

}
\references{
For the theoretic background see the references given on the main help page of this package: \link{multivariance-package}.
}
