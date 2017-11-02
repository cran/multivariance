#' Multivariance: detecting and measuring multivariate dependence
#'
#' The multivariance package provides basic functions to calculate distance multivariance and related quantities.
#'
# It also includes a function to perform a dependence analysis.
#'
#' Distance multivariance is a measure of dependence which can be used to detect and quantify dependence structures. The necessary functions are implemented in this packages, and examples are given. For the theoretic background we refer to the forthcoming papers [1,2].
#'
#'    The (current) code should be understood as \emph{proof of concept}, certainly there is room for improvement and development. Questions, comments and remarks are welcome: \email{bjoern.boettcher@@tu-dresden.de}
#'
#'
#' @section Multivariance:
#'
#'  \code{\link{multivariance}} computes the distance multivariance
#'
#'  \code{\link{total.multivariance}} computes the total distance multivariance
#'
#  \code{\link{k.multivariance}} computes the 'k distance multivariance'
#'
#' @section Functions to use and interpret multivariance:
#'
#'  \code{\link{rejection.level}} computes the rejection level for a given significance level. This can be used for a sensitive (conservative) interpretation of distance multivariance. The counterpart is \code{\link{multivariance.pvalue}}, which computes the pvalue for a given distance multivariance. Both use the hypothesis of independence.
#'
#'  \code{\link{independence.test}} provides the corresponding test of independence.
#'
#' \code{\link{cdm}} and \code{\link{cdms}} compute the centered distance matrix and matrices, respectively. These can be used for faster computation of distance multivariance.
#'
# @section Dependence structures:
#
#  \code{\link{dependence.structure}} computes the dependence structure as described in [3].
#'
#' @section Examples:
#'
#' \code{\link{coins}} and \code{\link{tetrahedron}} generate samples of pairwise independent random variables, with dependence of higher order.
#'
#' @references
#' [1] B. Böttcher, M. Keller-Ressel, R.L. Schilling, Detecting independence of random vectors I. Generalized distance covariance and Gaussian covariance. Preprint 2017.
#'
#' [2] B. Böttcher, M. Keller-Ressel, R.L. Schilling, Detecting independence of random vectors II. Distance multivariance and Gaussian multivariance. Preprint 2017.
#'
#'
#' @docType package
#' @name multivariance-package
NULL




################# Multivariance ###########

#' rejection level for the test statistic
#'
#' Under independence the probability for the normalized and Nscaled multivariance to be above this level is less than alpha. The same holds for the normalized, Nscaled and Escaled total multivariance.
#'
#' @param alpha level of significance
#' @details
#' The value might be very conservative. This is the counterpart to \code{\link{multivariance.pvalue}}.
#'
#' @examples
#' \dontrun{ curve(rejection.level(x),xlim = c(0.001,0.2))}
#'
#' rejection.level(0.05) #the rejection level, for comparison with the following values
#' total.multivariance(matrix(rnorm(100*3),ncol = 3)) #independent sample
#' total.multivariance(coins(100)) #dependent sample which is 2-independent
#'
#' multivariance.pvalue(total.multivariance(matrix(rnorm(100*3),ncol = 3))) #independent sample
#' multivariance.pvalue(total.multivariance(coins(100))) #dependent sample which is 2-independent
#'
#' @export
rejection.level = function(alpha) {
  return((stats::qnorm(1-alpha/2)^2))
  # identical with qchisq(1-alpha,1)
}

#' transform multivariance to pvalue
#'
#' Computes the p-value (under the hypothesis of independence) for a given multivariance/total multivariance.
#'
#' @param x value of a normalized and Nscaled \code{\link{multivariance}}
#'
#' @details
#' The p-value is conservative, i.e. it might be much smaller. This is the counterpart to \code{\link{rejection.level}}.
#'
#' @export
multivariance.pvalue = function(x) {2-2*stats::pnorm(sqrt(x))}


#' centered distance matrix
#'
#' calculates the centered distance matrix
#'
#' @param x matrix, each row of the matrix is treated as one sample
#' @param normalize logical - indicates if the matrix should be normalized
#' @param psi the continuous negative definite function to be used. If it is \code{NULL}, the euclidean distance will be used
#' @details
#' Short: If \code{normalize = TRUE} then value of multivariance is comparable and meaningful. It can be compared to the \code{\link{rejection.level}} or its p-value \code{\link{multivariance.pvalue}} can be computed.
#'
#' If \code{normalize = TRUE} the matrix is scaled such that the multivariance based on it, times the sample size, has in the limit - in the case of independence - the distribution of an L^2 norm of a Gaussian process with known expectation.
#' @examples
#' x = coins(100)
#' cdm(x) # fast euclidean distances
#' cdm(x,psi = function(x,y) sqrt(sum((x-y)^2))) # this is identical to the previous (but slower)
#' @export
cdm = function(x, normalize = TRUE, psi = NULL) {
  if (is.null(psi)) {
    dm = as.matrix(stats::dist(x,method="euclidean"))
  } else {
    x = as.matrix(x)
    n = nrow(x)
    d = ncol(x)
    dm = matrix(apply(cbind(x[rep(1:n,n),],x[rep(1:n,each = n),]), #create all combinations
                 1, # apply to each row
                 function(y) psi(y[1:d], y[(d+1):(2*d)])),nrow = n)
    # DEVELOPING NOTE: could double the speed if only the upper triangular matrix is computed.
  }
  m = mean(dm)
  # print(m)
  cdm1 = sweep(dm,1,colMeans(dm))
  cdm2 = -sweep(cdm1,2,rowMeans(dm)) - m
  if (normalize && (m != 0)) cdm2 = cdm2 / m
  return(cdm2)
}

#' computes the centered distance matrices
#' @param x matrix, each row is a sample
#' @param membership vector which indicates which columns are treated as one sample
#' @param ... these are passed to \code{\link{cdm}}
#' @return It returns an 3 dimensional array of the distance matrices. The index of the first dimension names the component for which the matrix is computed, thus it ranges from 1 to max(membership).
#' @details
#' \code{membership} has to be a vector of length \code{ncol(x)}. Furthermore it has to have values from 1 to the number of classes.
#' @export
cdms = function(x,membership = 1:ncol(x),...) {
  #n = ncol(inputmatrix)
  #rm(A,envir = pe)
  n = max(membership)
  N = nrow(x)
  array.cdm = array(,dim = c(n,N,N))
  for (i in 1:n) array.cdm[i,,] = cdm(x[,(membership == i)],...)
  return(array.cdm)
}


#' Product of distance matrices
#'
#' Computes the product of the distance matrices with the given indexes
#' @inheritParams multivariance
#'
#' @keywords internal
cdm.product = function(x,vec = if (length(dim(x)) == 2) { 1:ncol(x)} else {NA}) {
  if (length(dim(x)) == 2) { # if the input is a matrix, the distance matrices are computed
    x = cdms(x,vec)
    vec = 1:max(vec)
  }
  if (anyNA(vec)) warning("x is array, missing vec argument.")
  Aprod = x[vec[1],,]
  for (i in 2:length(vec)) Aprod = Aprod * x[vec[i],,]
  return(Aprod)
  #return(apply(array.cdm[vec,,],c(2,3),prod)) #this vector version is much much slower!!!
}

#' distance multivariance
#'
#' Computes the distance multivariance, either for given data or a given array of centered distance matrices.
#'
#' @param x either a data matrix or an array of centered distance matrices
#' @param vec if x is a matrix, then this indicates which columns are treated together as one sample; if x is an array, these are the indexes for which the multivariance is calculated. The default is all columns and all indicies, respectively.
#' @param Nscale if \code{TRUE} the multivariance is scaled up by the sample size (and thus it is exactly as required for the test of independence). This is only used if \code{correlation = FALSE}.
#' @param correlation if \code{TRUE} the multivariance is normalized by norms of their centered distance matrices.
#' @param ... these are passed to \code{\link{cdms}} (which is only invoked if \code{x} is a matrix)
#' @details
#'
#' If x is an matrix and vec is not given, then each column is treated as a separate sample. Otherwise vec has to have as many elements as x has columns and values starting from 1 up to the number of clusters. It computes the normalized multivariance, for a multivariance without normalization the argument \code{normalize = FALSE} has to be passed to \code{cdms}.
#'
#' If x is an array, then vec has to be given.
#'
#' \code{correlation = TRUE} yields values between 0 and 1. These can be interpreted similarly to classical correlations, see also \code{\link{multicorrelation}}.
#'
#' For more details see the references given in the \link[=multivariance-package]{package documentation}.
#'
#' @examples
#' multivariance(matrix(rnorm(100*3),ncol = 3)) #independent sample
#' multivariance(coins(100)) #dependent sample which is 2-independent
#'
#' x = matrix(rnorm(100*2),ncol = 2)
#' x = cbind(x,x[,2])
#' multivariance(x) #dependent sample which is not 2-independent
#' multivariance(x[,1:2]) #these are independent
#' multivariance(x[,2:3]) #these are dependent
#'
#' multivariance(x[,2:3],correlation = TRUE)
#'
#' @export
multivariance = function(x,vec = NA,Nscale = FALSE,correlation = FALSE, ...) {
  if (length(dim(x)) == 2) { # if the input is a matrix, the distance matrices are computed
    if (anyNA(vec)) vec = 1:ncol(x)
    x = cdms(x,vec,...)
    vec = 1:max(vec)
  }
  if (anyNA(vec)) vec = 1:dim(x)[1] #warning("x is array, missing vec argument.")
  if (length(vec) > dim(x)[2]) warning("More data columns than rows.")
  Aprod = x[vec[1],,]
  for (i in 2:length(vec)) Aprod = Aprod * x[vec[i],,]
  result = mean(Aprod)
  if (Nscale && !correlation) result = result *nrow(Aprod)

  if (correlation) {
    n = length(vec)
    Anorm = mean(abs(x[vec[1],,]^n))^(1/n)
    for (i in 2:length(vec)) Anorm = Anorm * mean(abs(x[vec[i],,]^n))^(1/n)
    result = result / Anorm
  }
  # DEVELOPING NOTE: The following is much slower .... It might be faster, if we store the matrices only as vectors with the upper triangular as elements.
  #ut = upper.tri(Aprod)
  #diat = diag(x[vec[1],,])
  #test = x[vec[1],,][ut]
  #for (i in 2:length(vec)) {
  #  diat = diat * diag(x[vec[i],,])
  #  test = test * x[vec[i],,][ut]
  #}
  #erg = sum(diat,2*test)/ncol(x)^2

  return(result)

  # return(mean(apply(x[vec,,],c(2,3),prod))) #this vector version is also much much slower!!!
  # DEVELOPING NOTE: mean is slow due to error correction. sum()/length() is faster, but not as precise.
}

#' distance multicorrelation
#'
#' computes the distance multicorrelation
#'
#' @inheritParams multivariance
#'
#' @details
#' This is just a wrapper for \code{\link{multivariance}(x,vec,Nscale = FALSE,correlation = TRUE,...)}.
#'
#' @export
multicorrelation = function(x,vec = NA,...) {
  multivariance(x,vec,Nscale = FALSE,correlation = TRUE,...)
}

#' total distance multivariance
#'
#' computes the total distance multivariance for the samples with the given index.
#' It uses the global array of distance matrices.
#'
#' @inheritParams multivariance
#' @param lambda a scaling parameter >0. Each k-tuple multivariance gets weight \code{lambda^(n-k)}.
#' @param Escale if \code{TRUE} then it is scaled by the number of multivariances which are theoretically summed up (in the case of independence this yields for normalized distance matrices an estimator with expectation 1)
#'
#' @details
#' For details see the references given in the \link[=multivariance-package]{package documentation}.
#'
#' @examples
#' total.multivariance(matrix(rnorm(100*3),ncol = 3)) #independent sample
#' total.multivariance(coins(100)) #dependent sample which is 2-independent
#'
#' @export
total.multivariance = function(x,vec = NA,lambda = 1, Nscale = TRUE,Escale = TRUE,...) {
  if (length(dim(x)) == 2) { # if the input is a matrix, the distance matrices are computed
    if (anyNA(vec)) vec = 1:ncol(x)
    x = cdms(x,vec,...)
    vec = 1:max(vec)
  }
  if (anyNA(vec)) warning("x is array, missing vec argument.")
  if (length(vec) > dim(x)[2]) warning("More data columns than rows.")
  Aprod = lambda + x[vec[1],,]
  for (i in 2:length(vec)) Aprod = Aprod * (lambda + x[vec[i],,])

  result = mean(Aprod)-lambda^(length(vec))
  if (Nscale) result = result *nrow(Aprod)
  if (Escale) result = result/((1+lambda)^length(vec) - length(vec)*lambda^(length(vec)-1) - lambda^length(vec))

  return(result)
}

#' test for independence
#'
#' This computes a test of independence for given columns of a sample matrix or given centered distance matrices.
#'
#' @inheritParams multivariance
#' @param alpha significance level
#'
#' @examples
#' independence.test(coins(100)) # dependent sample
#' independence.test(coins(100)[,2:3]) # independent sample
#' independence.test(coins(100)) #dependent sample which is 2-independent
#'
#' @return Returns \code{TRUE} if the hypothesis of independence is NOT rejected, otherwise \code{FALSE}.
#' @details The test might be very conservative.
#'
#' The centered distance matrices can be prepared by \code{\link{cdms}}.
#'
#' @export
independence.test = function(x,vec = NA,alpha = 0.05) {
  tm = total.multivariance(x,vec)
  R = rejection.level(alpha)
  cat("\nThe value of the test statistic is",tm,"and values above",R,"are rejected.\n")
  result = tm>R
  if (result) {cat("The hypothesis of independence is rejected.\n")
    }
  else {cat("The hypothesis of independence is NOT rejected.\n")
  }
  invisible(result)
}

################# Generating Example Data ##########

#' tetrahedron sampling
#'
#' This function creates samples of a tetrahedron-dice colored r, g, b and rgb. Each sample indicates if for the thrown dice the colors r, g and b are contained on the bottom side of the dice.
#'
#' @param N number of samples
#' @return It returns the samples of the events r, g and b as rows of a \code{N} by 3 matrix (the first column corresponds to r, the second to g,...). TRUE indicates that this color is on the bottom side of the dice. The columns are dependent but 2-independent.
#' @examples
#' tetrahedron(10)
#'
#' @export
tetrahedron = function(N = 1000) {
  # rolls the tetrahedron with sides r,g,b,rgb
  side = sample(1:4,N,replace=TRUE)
  x = (side == 1)|(side == 4)
  y = (side == 2)|(side == 4)
  z = (side == 3)|(side == 4)
  return(unname(cbind(x,y,z)))
}

#' k-independent coin sampling
#'
#' This function creates samples which are k-independent.
#' @param N number of samples
#' @param k each k-tuple will be independent
#' @return It returns the samples as rows of an \code{N} by \code{k+1} matrix. The columns are dependent but k-independent.
#' @examples
#' coins(200,4)
#'
#' @export
coins = function(N = 1000, k = 2) {
  c = matrix(sample(c(0,1),k*N,replace = TRUE),ncol = k)
  d = (rowSums(as.matrix(c[,1:(k-1)])) %% 2) == c[,k]
  return(unname(cbind(c,d)))
}

