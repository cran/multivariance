% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multivariance-functions.R
\docType{data}
\name{dep_struct_iterated_13_100}
\alias{dep_struct_iterated_13_100}
\title{example dataset for \code{\link{dependence.structure}}}
\format{
\code{matrix} 13 variables (columns), 100 independent samples (rows)
}
\usage{
dep_struct_iterated_13_100
}
\description{
It was generated by \preformatted{
set.seed(532333356)
N = 100
x = matrix(sample.int(2,10*N,replace = TRUE)-1,ncol = 10)
for (i in c(2,5,9)) x = cbind(x,(rowSums(as.matrix(x[,1:(i-1)])) %% 2) == x[,i])
dep_struct_iterated_13_100 = x
save(dep_struct_iterated_13_100,file ="dep_struct_iterated_13_100.rda")
}
}
\details{
To avoid irritation, note that the seed is just a simple integer hash value of the variable name.
}
\keyword{datasets}
