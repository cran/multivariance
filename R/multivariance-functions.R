# Package description ####

#' multivariance: Measuring Multivariate Dependence Using Distance Multivariance
# Multivariance: detecting and measuring multivariate dependence
#'
#' The multivariance package provides basic functions to calculate distance multivariance and related quantities.
#'
# It also includes a function to perform a dependence analysis.
#'
#' Distance multivariance is a measure of dependence which can be used to detect and quantify dependence structures. The necessary functions are implemented in this packages, and examples are given. For the theoretic background we refer to the papers [1,2] and [3]. The latter includes a summary of the first two. It is the recommended starting point for users with an applied interest.
#'
#'  The (current) code should be understood as \emph{proof of concept}, certainly there is room for improvement and development. Questions, comments and remarks are welcome: \email{bjoern.boettcher@@tu-dresden.de}
#'
#' For infos on the latest changes and/or updates to the package use \code{news(package="multivariance")}.
#'
#' To cite this package use the standard citation for R packages, i.e., the output of \code{citation("multivariance")}.
#'
#' @section Multivariance:
#'
#'  \code{\link{multivariance}} computes the distance multivariance
#'
#'  \code{\link{total.multivariance}} computes the total distance multivariance
#'
#'  \code{\link{m.multivariance}} computes the m-multivariance (introduced in [3])
#'
#'  It might be convenient to compute these simultaneously using \code{\link{multivariances.all}}.
#'
#' @section Functions to use and interpret multivariance:
#'
#'  \code{\link{rejection.level}} computes the rejection level for a given significance level. This can be used for a conservative interpretation of distance multivariance. The counterpart is \code{\link{multivariance.pvalue}}, which computes a conservative p-value for a given distance multivariance. Both methods are distribution-free.
#'
#'  \code{\link{resample.rejection.level}} and \code{\link{resample.pvalue}} are the distribution dependent versions of the above. They are approximately sharp, but computational more expensive. Any resampling is done by \code{\link{resample.multivariance}}.
#'
#'  \code{\link{independence.test}} provides the corresponding tests of independence.
#'
#' \code{\link{cdm}} and \code{\link{cdms}} compute the centered distance matrix and matrices, respectively. These can be used to speed up repeated computations of distance multivariance.
#'
#' @section Dependence structures:
#'
#'  \code{\link{dependence.structure}} performs the dependence structure detection algorithm as described in [3].
#'
#'  \code{\link{find.cluster}} is the basic building block of \code{\link{dependence.structure}}. It is recommended to use \code{\link{dependence.structure}}.
#'
#' @section Examples:
#'
#' \code{\link{coins}} and \code{\link{tetrahedron}} generate samples of pairwise independent random variables, with dependence of higher order.
#'
#' \code{\link{dep_struct_iterated_13_100}}, \code{\link{dep_struct_ring_15_100}}, \code{\link{dep_struct_several_26_100}} and \code{\link{dep_struct_star_9_100}} are example data sets for the dependence structure detection. These might also serve as benchmark examples.
#'
#' @references
#' [1] B. Böttcher, M. Keller-Ressel, R.L. Schilling, Detecting independence of random vectors I. Generalized distance covariance and Gaussian covariance. Preprint 2017. \url{https://arxiv.org/abs/1711.07778}
#'
#' [2] B. Böttcher, M. Keller-Ressel, R.L. Schilling, Detecting independence of random vectors II. Distance multivariance and Gaussian multivariance. Preprint 2017. \url{https://arxiv.org/abs/1711.07775}
#'
#' [3] B. Böttcher, Dependence Structures - Estimation and Visualization Using Distance Multivariance. Preprint 2017. \url{https://arxiv.org/abs/1712.06532}
#'
#'
#' @docType package
#' @name multivariance-package
NULL




################# Multivariance ###########

#' rejection level for the test statistic
#'
#' Under independence the probability for the normalized and Nscaled multivariance to be above this level is less than \code{alpha}. The same holds for the normalized, Nscaled and Escaled total multivariance and m-multivariance.
#'
#' @param alpha level of significance
#' @details
#' This is based on a distribution-free approach. The value might be very conservative. This is the counterpart to \code{\link{multivariance.pvalue}}. For a less conservative approach see \code{\link{resample.rejection.level}}.
#'
#' The estimate is only valid for \code{alpha} smaller than 0.215.
#'
#' @examples
#' rejection.level(0.05) #the rejection level, for comparison with the following values
#' total.multivariance(matrix(rnorm(100*3),ncol = 3)) #independent sample
#' total.multivariance(coins(100)) #dependent sample which is 2-independent
#'
#' # and the p values are (to compare with alpha)
#' multivariance.pvalue(total.multivariance(matrix(rnorm(100*3),ncol = 3))) #independent sample
#' multivariance.pvalue(total.multivariance(coins(100))) #dependent sample which is 2-independent
#'
#' \dontrun{
#' # visualization of the rejection level
#' curve(rejection.level(x),xlim = c(0.001,0.215),xlab = "alpha")
#' }
#'
#' @export
rejection.level = function(alpha) {
  if (any(alpha > 0.215)) warning("alpha too large. Only valid for alpha smaller than 0.215!")
  return((stats::qnorm(1-alpha/2)^2))
  # identical with qchisq(1-alpha,1)
}

#' transform multivariance to p-value
#'
#' Computes the p-value for the hypothesis of independence for a given multivariance/total multivariance.
#'
#' @param x value of a normalized and Nscaled \code{\link{multivariance}}
#'
#' @details
#' This is based on a distribution-free approach. The p-value is conservative, i.e. it might be much smaller. This is the counterpart to \code{\link{rejection.level}}. For a less conservative approach see \code{\link{resample.pvalue}}.
#'
#' p-values larger than 0.215 might be incorrect, since the distribution-free estimate on which the computation is based only holds up to 0.215.
#'
#' @references
#' For the theoretic background see the references given on the main help page of this package: \link{multivariance-package}.
#'
#' @export
multivariance.pvalue = function(x) {
  if (any(x < 0)) print(paste("Negative multivariance = ",x[which(x<0)]))
  2-2*stats::pnorm(sqrt(x))
}


#' centered distance matrix
#'
#' computes the centered distance matrix
#'
#' @param x matrix, each row of the matrix is treated as one sample
#' @param normalize logical, indicates if the matrix should be normalized
#' @param psi a real valued function of two variables, to compute the distance of two samples based on a continuous negative definite function. If it is \code{NULL}, the euclidean distance will be used
#' @param p numeric, if it is a value between 1 and 2 then the Minkowski distance with parameter p is used.
#'
#' @details
#' The centered distance matrices are required for the computation of (total / m-) multivariance.
#'
#'If \code{normalize = TRUE} then the value of multivariance is comparable and meaningful. It can be compared to the \code{\link{rejection.level}} or its p-value \code{\link{multivariance.pvalue}} can be computed.
#'
#' More details: If \code{normalize = TRUE} the matrix is scaled such that the multivariance based on it, times the sample size, has in the limit - in the case of independence - the distribution of an L^2 norm of a Gaussian process with known expectation.
#'
#' @references
#' For the theoretic background see the references given on the main help page of this package: \link{multivariance-package}.
#'
#' @examples
#' x = coins(100)
#' cdm(x) # fast euclidean distances
#' cdm(x,psi = function(x,y) sqrt(sum((x-y)^2))) # this is identical to the previous (but slower)
#'
#' # the function cdm does the following three lines in a faster way
#' N = nrow(x)
#' C = diag(N) - matrix(1/N,nrow = N,ncol = N)
#' A = - C %*% as.matrix(stats::dist(x,method="euclidean")) %*% C #'
#' all(abs(A- cdm(x,normalize = FALSE)) < 10^(-12))
#'
#' @export
cdm = function(x, normalize = TRUE, psi = NULL, p = NULL) {
  if (is.null(psi) & is.null(p)) {
    dm = dist.to.matrix(stats::dist(x,method="euclidean"))
    #DEVELOPING NOTE: here as.matrix was slow, dist.to.matrix is faster. Instead one could just use the vector....
  } else {
    if (!is.null(p)) {
      if ((p<1) || (p>2)) warning("p is not in [1,2]")
      dm = dist.to.matrix(stats::dist(x,method="minkowski", p = p))
    } else {
      x = as.matrix(x)
      n = nrow(x)
      d = ncol(x)
      dm = matrix(apply(cbind(x[rep(1:n,n),],x[rep(1:n,each = n),]), #create all combinations
                        1, # apply to each row
                        function(y) psi(y[1:d], y[(d+1):(2*d)])),nrow = n)
      # DEVELOPING NOTE: could double the speed if only the upper triangular matrix is computed, using the idea of dist.to.matrix
    }
  }
  colm = colMeans(dm)
  m = mean(colm)  # equals mean(dm)

  if (m == 0) warning("It seems that one variable is constant. Constants are always independent.")
  if (normalize && (m != 0)) {
    return((-dm + outer(colm,colm, FUN ="+") - m)/ m)
  } else {
    return(-dm + outer(colm,colm, FUN ="+") - m)
  }
  #alternative (slower) implementations:
  #cdm1 = sweep(dm,1,colm)
  #cdm2 = -sweep(cdm1,2,rowMeans(dm)) - m
  #cdm2 = -(x - rep(colm, ncol(dm)) - rep(rowMeans(dm),each = ncol(dm))) - m # for quadratic matrix

}

#' computes the centered distance matrices
#' @param x matrix, each row is a sample
#' @param vec vector which indicates which columns are treated as one sample
#' @param membership depreciated. Now use \code{vec}.
#' @param ... these are passed to \code{\link{cdm}}
#'
#' @return It returns an 3 dimensional array of the distance matrices. The index of the first dimension names the component for which the matrix is computed, thus it ranges from 1 to max(vec).
#'
#' @export
cdms = function(x,vec = 1:ncol(x),membership = NULL,...) {
  if (!is.null(membership)) {
    vec = membership
    warning("Use 'vec' instead of 'membership' as argument to 'cdms'. 'membership' is depreciated.")
  }
  if (anyNA(vec)) vec = 1:ncol(x)
  n = max(vec)
  N = nrow(x)
  array.cdm = array(,dim = c(n,N,N))
  for (i in 1:n) array.cdm[i,,] = cdm(x[,(vec == i)],...)
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
#' @param vec if x is a matrix, then this indicates which columns are treated together as one sample; if x is an array, these are the indexes for which the multivariance is calculated. The default is all columns and all indexes, respectively.
#' @param Nscale if \code{TRUE} the multivariance is scaled up by the sample size (and thus it is exactly as required for the test of independence)
#' @param correlation if \code{TRUE} the multivariance is scaled by norms of their centered distance matrices, and \code{Nscale} will be ignored.
#' @param squared if \code{FALSE} it returns the actual multivariance, otherwise the squared multivariance (less computation)
#' @param ... these are passed to \code{\link{cdms}} (which is only invoked if \code{x} is a matrix)
#' @details
#'
#' If \code{x} is an matrix and \code{vec} is not given, then each column is treated as a separate sample. Otherwise \code{vec} has to have as many elements as \code{x} has columns and values starting from 1 up to the number of 'variables', e.g. if \code{x} is an \code{N} by 5 matrix and \code{vec = c(1,2,1,3,1)} then the multivariance of the 1-dimensional variables represented by column 2 and 4 and the 3-dimensional variable represented by the columns 1,3,5 is computed.
#'
#' As default it computes the normalized Nscaled squared multivariance, for a multivariance without normalization the argument \code{normalize = FALSE} has to be passed to \code{cdms}.
#'
#' If \code{x} is an array, then \code{vec} has to be given.
#'
#' \code{correlation = TRUE} yields values between 0 and 1. These can be interpreted similarly to classical correlations, see also \code{\link{multicorrelation}}.
#'
#' As a rough guide to interpret the value of distance multivariance note:
#' \itemize{
#' \item If the random variables are not (n-1)-independent, large values indicate dependence, but small values are meaningless. Thus in this case use \code{\link{total.multivariance}}.
#' \item If the random variables are (n-1)-independent and \code{Nscale = TRUE}, values close to 1 and smaller indicate independence, larger values indicate dependence. In fact, in the case of independence the test statistic is a Gaussian quadratic form with expectation 1 and samples of it can be generated by \code{\link{resample.multivariance}}.
#' \item If the random variables are (n-1)-independent and \code{Nscale = FALSE}, small values (close to 0) indicate independence, larger values indicate dependence.
#' }
#'
#' Finally note, that due to numerical (in)precision the value of multivariance might become negative. In these cases it is set to 0. A warning is issued, if the value is negative and further than the usual (used by \code{\link[base]{all.equal}}) tolerance away from 0.
#'
#' @references
#' For the theoretic background see the references given on the main help page of this package: \link{multivariance-package}.
#'
#' @examples
#' multivariance(matrix(rnorm(100*3),ncol = 3)) #independent sample
#' multivariance(coins(100)) #dependent sample which is 2-independent
#'
#' x = matrix(rnorm(100*2),ncol = 2)
#' x = cbind(x,x[,2])
#' multivariance(x) #dependent sample which is not 2-independent (thus small values are meaningless!)
#' multivariance(x[,1:2]) #these are independent
#' multivariance(x[,2:3]) #these are dependent
#'
#' multivariance(x[,2:3],correlation = TRUE)
#'
#' @export
multivariance = function(x,vec = NA,Nscale = TRUE,correlation = FALSE, squared = TRUE, ...) {
  if (length(dim(x)) == 2) { # if the input is a matrix, the distance matrices are computed
    if (anyNA(vec)) vec = 1:ncol(x)
    x = cdms(x,vec,...)
    vec = 1:max(vec)
  }
  if (anyNA(vec)) vec = 1:dim(x)[1] #warning("x is array, missing vec argument.")
  #if (length(vec) > dim(x)[2]) warning("More data columns than rows.")
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
  if (result < 0) {
    if (!isTRUE(all.equal(result,0))) warning(paste("Value of multivariance was negative (",result,"). This is usually due to numerical (in)precision. It was set to 0."))
    result = 0
  }

  if (squared | Nscale) { return(result)
  } else { return(sqrt(result))}


  # return(mean(apply(x[vec,,],c(2,3),prod))) #this vector version is also much much slower!!!
  # DEVELOPING NOTE: mean is slow due to error correction. sum()/length() is faster, but not as precise.
}

#' total distance multivariance
#'
#' computes the total distance multivariance
#'
#' @inheritParams multivariance
#' @param lambda a scaling parameter >0. Each k-tuple multivariance gets weight \code{lambda^(n-k)}.
#' @param Escale if \code{TRUE} then it is scaled by the number of multivariances which are theoretically summed up (in the case of independence this yields for normalized distance matrices an estimator with expectation 1)
#'
#' @details
#' Total distance multivariance is per definition the scaled sum of certain distance multivariances, and it characterize dependence.
#'
#'  As a rough guide to interpret the value of total distance multivariance note:
#' \itemize{
#' \item Large values indicate dependence.
#' \item For\code{Nscale = TRUE} values close to 1 and smaller indicate independence, larger values indicate dependence. In fact, in the case of independence the test statistic is a Gaussian quadratic form with expectation 1 and samples of it can be generated by \code{\link{resample.multivariance}}.
#' \item For \code{Nscale = FALSE} small values (close to 0) indicate independence, larger values indicate dependence.
#' }

#'
#' Finally note, that due to numerical (in)precision the value of total multivariance might become negative. In these cases it is set to 0. A warning is issued, if the value is negative and further than the usual (used by \code{\link[base]{all.equal}}) tolerance away from 0.
#'
#'@references
#' For the theoretic background see the references given on the main help page of this package: \link{multivariance-package}.
#'
#' @examples
#' x = matrix(rnorm(100*3),ncol = 3)
#' total.multivariance(x) #for an independent sample
#' # the value coincides with
#' (multivariance(x[,c(1,2)],Nscale = TRUE) + multivariance(x[,c(1,3)],Nscale = TRUE)+
#'  multivariance(x[,c(2,3)],Nscale = TRUE) + multivariance(x,Nscale = TRUE))/4
#'
#' total.multivariance(coins(100)) #value for a dependent sample which is 2-independent
#'
#' @export
total.multivariance = function(x,vec = NA,lambda = 1, Nscale = TRUE,Escale = TRUE,squared = TRUE,...) {
  if (length(dim(x)) == 2) { # if the input is a matrix, the distance matrices are computed
    if (anyNA(vec)) vec = 1:ncol(x)
    x = cdms(x,vec,...)
    vec = 1:max(vec)
  }
  if (anyNA(vec)) vec = 1:dim(x)[1] #warning("x is array, missing vec argument.")
  #if (length(vec) > dim(x)[2]) warning("More data columns than rows.")
  Aprod = lambda + x[vec[1],,]
  for (i in 2:length(vec)) Aprod = Aprod * (lambda + x[vec[i],,])

  result = mean(Aprod)-lambda^(length(vec))

  if (result < 0) {
    if (!isTRUE(all.equal(result,0))) warning(paste("Value of total multivariance was negative (",result,"). This is usually due to numerical (in)precision. It was set to 0."))
    result = 0
  }

  if (Nscale) result = result *nrow(Aprod)
  if (Escale) result = result/((1+lambda)^length(vec) - length(vec)*lambda^(length(vec)-1) - lambda^length(vec))

  if (squared | Nscale) { return(result)
  } else { return(sqrt(result))}
}

#' m distance multivariance
#'
#' Computes m distance multivariance.
#'
#' @details
#'
#' m-distance multivariance is per definition the scaled sum of certain distance multivariances, and it characterize m-dependence.
#'
#'  As a rough guide to interpret the value of total distance multivariance note:
#' \itemize{
#' \item Large values indicate dependence.
#' \item If the random variables are (m-1)-independent and \code{Nscale = TRUE}, values close to 1 and smaller indicate m-independence, larger values indicate dependence. In fact, in the case of independence the test statistic is a gaussian quadratic form with expectation 1 and samples of it can be generated by \code{\link{resample.multivariance}}.
#' \item If the random variables are (m-1)-independent and \code{Nscale = FALSE}, small values (close to 0) indicate m-independence, larger values indicate dependence.
#' }
#'
#' Since random variables are always 1-independent, the case \code{m=2} characterizes pairwise independence.
#'
#' Finally note, that due to numerical (in)precision the value of m-multivariance might become negative. In these cases it is set to 0. A warning is issued, if the value is negative and further than the usual (used by \code{\link[base]{all.equal}}) tolerance away from 0.
#'
#' @references
#' For the theoretic background see the reference [3] given on the main help page of this package: \link{multivariance-package}.
#'
#' @inheritParams multivariance
#' @param m \code{=2} or \code{3} the m-multivariance will be computed.
#' @param Escale if \code{TRUE} then it is scaled by the number of multivariances which are theoretically summed up (this yields an expectation of 1 in the case of independence)
#'
#'
#'
#' @examples
#' x = matrix(rnorm(3*30),ncol = 3)
#'
#' # the following values are identical
#' m.multivariance(x,m =2)
#' 1/choose(3,2)*(multivariance(x[,c(1,2)]) +
#'                multivariance(x[,c(1,3)]) +
#'                multivariance(x[,c(2,3)]))
#'
#' # the following values are identical
#' m.multivariance(x,m=3)
#' multivariance(x)
#'
#' # the following values are identical
#' 1/4*(3*(m.multivariance(x,m=2)) + m.multivariance(x,m=3))
#' total.multivariance(x, Nscale = TRUE)
#' 1/4*(multivariance(x[,c(1,2)], Nscale = TRUE) +
#'      multivariance(x[,c(1,3)], Nscale = TRUE) +
#'      multivariance(x[,c(2,3)], Nscale = TRUE) + multivariance(x, Nscale = TRUE))
#'
#' @export
m.multivariance = function(x, vec= NA, m = 2, Nscale = TRUE, Escale = TRUE, squared = TRUE,...) {
  if (length(dim(x)) == 2) { # if the input is a matrix, the distance matrices are computed
    if (anyNA(vec)) vec = 1:ncol(x)
    x = cdms(x,vec,...)
    vec = 1:max(vec)
  }
  if (anyNA(vec)) vec = 1:dim(x)[1] #warning("x is array, missing vec argument.")
 # if (length(vec) > dim(x)[2]) warning("More data columns than rows.")

  #k = 2
  if (m == 2) {
    Asum = x[vec[1],,]
    A2sum = Asum^2 #x[vec[1],,]^2
    for (i in 2:length(vec)) {
      tempFactor = x[vec[i],,]
      Asum = Asum + tempFactor
      A2sum = A2sum + tempFactor^2
    }

    result = mean(Asum^2 - A2sum)/2
  }

  #k more general
  #k = 3
  if (m == 3) {
    Asum = x[vec[1],,]
    A2sum = Asum^2 #x[vec[1],,]^2
    A3sum = A2sum * Asum #x[vec[1],,]^3
    for (i in 2:length(vec)) {
      tempFactor = x[vec[i],,]
      Asum = Asum + tempFactor
      summand = tempFactor^2
      A2sum = A2sum + summand #x[vec[i],,]^2
      A3sum = A3sum + summand *tempFactor #x[vec[i],,]^3
    }
    result = mean(Asum^3 - 3* Asum *A2sum + 2 * A3sum)/ 6
    #    result = mean(Asum^3 - choose(3,2)* Asum *A2sum + 2 * A3sum)/ factorial(3)
  }

  if (m > 3) { warning("m > 3, not implemented.")}

  if (result < 0) {
    if (!isTRUE(all.equal(result,0))) warning(paste("Value of m-multivariance was negative (",result,"). This is usually due to numerical (in)precision. It was set to 0."))
    result = 0
  }

  if (Nscale) result = result *nrow(Asum)
  if (Escale) result = result/(choose(length(vec),m))

  if (squared | Nscale) { return(result)
  } else { return(sqrt(result))}
}

#' simultaneous computation of total/ 2-/ 3- /... multivariance
#'
#' Computes simultaneously multivariance, total multivariance, 2-multivariance and 3-multivariance.
#'
#' @inheritParams multivariance
#'
#' @seealso \code{\link{multivariance}}, \code{\link{total.multivariance}}, \code{\link{m.multivariance}}
#'
#' @details
#' The computation is faster than the seperate computations.
#'
#' @examples
#' x = coins(100,k = 3)
#' multivariances.all(x)
#' # yields the same as:
#' multivariance(x)
#' total.multivariance(x)
#' m.multivariance(x,m=2)
#' m.multivariance(x,m=3)
#'
#'
#' @export
#'
multivariances.all = function(x, vec= NA, Nscale = TRUE, squared = TRUE,...) {
  if (length(dim(x)) == 2) { # if the input is a matrix, the distance matrices are computed
    if (anyNA(vec)) vec = 1:ncol(x)
    x = cdms(x,vec,...)
    vec = 1:max(vec)
  }
  if (anyNA(vec)) vec = 1:dim(x)[1] #warning("x is array, missing vec argument.")
  # if (length(vec) > dim(x)[2]) warning("More data columns than rows.")

  n = length(vec)

  Asum = x[vec[1],,]
  A2sum = Asum^2 # = x[vec[1],,]^2
  A3sum = A2sum * Asum # = x[vec[1],,]^3
  Aprod = Asum
  Aplusprod = 1 + Asum # = x[vec[1],,]
  for (i in 2:n) {
    tempFactor = x[vec[i],,]
    Asum = Asum + tempFactor
    Aprod = Aprod * tempFactor
    Aplusprod = Aplusprod * (1 + tempFactor)
    summand = tempFactor^2
    A2sum = A2sum + summand #x[vec[i],,]^2
    A3sum = A3sum + summand *tempFactor #x[vec[i],,]^3
  }
  m = mean(Aprod)
  mt = (mean(Aplusprod)-1)/(2^n - n - 1)
  m2 = mean(Asum^2 - A2sum)/(n *(n-1))
  m3 = mean(Asum^3 - 3* Asum *A2sum + 2 * A3sum)/ (n *(n-1)*(n-2))

  result = c(m,mt,m2,m3)

  neg.res = result<0
  if (any(neg.res)) {
    if (!isTRUE(all.equal(result[neg.res],0)))
    warning(paste("Value of ",c("","total","2","3")[neg.res],"multivariance was negative (",result[neg.res],"). This is usually due to numerical (in)precision. It was set to 0."))
    result[neg.res] = 0
  }

  if (Nscale) result = result *nrow(Asum)

  if (squared | Nscale) { return(result)
  } else { return(sqrt(result))}
}

#' distance multicorrelation
#'
#' computes the distance multicorrelation
#'
#' @inheritParams multivariance
#'
#' @details
#' This is just a wrapper for \code{\link{multivariance}(x,vec,Nscale = FALSE,correlation = TRUE,squared = squared,...)}.
#'
#'@references
#' For the theoretic background see the references given on the main help page of this package: \link{multivariance-package}.
#'
#' @export
multicorrelation = function(x,vec = NA,squared = FALSE, ...) {
  multivariance(x,vec,Nscale = FALSE,correlation = TRUE,squared = squared,...)
}

#' test for independence
#'
#' This computes a test of independence for the columns of a sample matrix (required for the resampling test) or for given centered distance matrices (only possible for the distribution-free test).
#'
#' @inheritParams multivariance
#' @param alpha significance level
#' @param type one of \code{"distribution_free","resample"}
#'
#' @return Returns \code{TRUE} if the hypothesis of independence is NOT rejected, otherwise \code{FALSE}.
#' @details The \code{"distribution_free"} test might be very conservative.
#' The centered distance matrices can be prepared by \code{\link{cdms}}. But note that for the resample test, the data matrix has to be given.
#'
#' @references
#' For the theoretic background see the references given on the main help page of this package: \link{multivariance-package}.
#'
#' @examples
#' independence.test(coins(100)) #dependent sample which is 2-independent
#' independence.test(coins(100),type = "resample") #dependent sample which is 2-independent
#'
#' independence.test(coins(100)[,2:3]) # independent sample
#' independence.test(coins(100)[,2:3],type = "resample") # independent sample
#'
#' independence.test(coins(10),type = "resample") #dependent sample which is 2-independent
#' independence.test(coins(10)[,2:3],type = "resample") #dependent sample which is 2-independent
#'
#' @export
independence.test = function(x,vec = 1:ncol(x),alpha = 0.05,type = "distribution_free") {
  tm = total.multivariance(x,vec)
  if (type == "distribution_free") {
    R = rejection.level(alpha)
    cat("\nThe value of the test statistic is",tm,"and values above",R,"are rejected.\n")
    result = tm>R
  } else {
    p.value = resample.pvalue(tm,x=x,vec=vec,times = 300,type="total")
  cat("\nThe value of the test statistic is",tm,"and (by resampling) its p-value is",p.value,"\n")
  result = p.value<alpha
  }

  if (result) {cat("The hypothesis of independence is rejected.\n")
  }
  else {cat("The hypothesis of independence is NOT rejected.\n")
  }
  invisible(result)
}

################# Example data ##########

#' dependence example: tetrahedron sampling
#'
#' This function creates samples of a tetrahedron-dice colored r, g, b and rgb. Each sample indicates if for the thrown dice the colors r, g and b are contained on the bottom side of the dice.
#'
#' @param N number of samples
#' @return It returns the samples of the events r, g and b as rows of a \code{N} by 3 matrix (the first column corresponds to r, the second to g,...). TRUE indicates that this color is on the bottom side of the dice. The columns are dependent but 2-independent.
#' @examples
#' tetrahedron(10)
#'
#' @references
#' For the theoretic background see the reference [3] given on the main help page of this package: \link{multivariance-package}.
#'
#' @export
tetrahedron = function(N = 1000) {
  # rolls the tetrahedron with sides r,g,b,rgb
  side = sample.int(4,N,replace=TRUE)
  x = (side == 1)|(side == 4)
  y = (side == 2)|(side == 4)
  z = (side == 3)|(side == 4)
  return(unname(cbind(x,y,z)))
}

#' dependence example: k-independent coin sampling
#'
#' This function creates samples which are dependent but k-independent.
#' @param N number of samples
#' @param k each k-tuple will be independent
#' @return It returns the samples as rows of an \code{N} by \code{k+1} matrix. The columns are dependent but k-independent.
#'
#' @references
#' For the theoretic background see the reference [3] given on the main help page of this package: \link{multivariance-package}.
#'
#' @examples
#' coins(200,4)
#'
#' @export
coins = function(N = 1000, k = 2) {
  c = matrix(sample.int(2,k*N,replace = TRUE)-1,ncol = k)
  d = (rowSums(as.matrix(c[,1:(k-1)])) %% 2) == c[,k]
  return(unname(cbind(c,d)))
}


# resample ######


#' resample the columns of a matrix
#' @param x matrix
#' @param vec vector, indicates which columns belong together
#' @param replace boolean, sampling with or without replacement
#'
#' @return Returns a matrix with the same dimensions as \code{x}. The columns are resampled from the original columns. The resampling is done with replacement (\code{replace = TRUE}) or without (\code{replace = FALSE}). Columns which belong together (indicated by vec) are resampled identically, i.e., all values in rows of these are kept together.
#'
#' @examples
#' sample.cols(matrix(1:15,nrow = 5),vec = c(1,1,2))
#'
#' @export
sample.cols = function(x,vec = 1:ncol(x),replace = TRUE) {
  if (anyNA(vec)) vec = 1:ncol(x)
  N = nrow(x)
  xnew = x
  for (i in 1:max(vec)) {
    xnew[,vec == i] = x[sample.int(N,replace = replace),vec == i]
  }
  return(xnew)
}


#' resampling (total /m-) multivariance
#'
#' The distribution of the test statistic under the hypothesis of independence is required for the independence tests. This function generates approximate samples of this distribution either by sampling without replacement (permutations) or by sampling with replacement (bootstrap).
#'
#' @details
#' The resampling is done by sampling from the original data either without replacement (\code{"permutation"}) or with replacement (\code{"bootstrap"}).
#'
#' For convenience also the actual (total /m-) multivariance is computed and its p-value.
#'
#' @param x matrix, the rows should be iid samples
#' @param vec vector, which indicates which columns of \code{x} are treated together as one sample
#' @param times integer, number of samples to generate
#' @param type one of \code{"multi","total","m.multi.2","m.multi.3"}
#' @param resample.type one of \code{"permutation", "bootstrap"}. The samples are generated without replacement (permutations) or with replacement (bootstrap).
#' @param ... is passed to \code{\link{multivariance}, \link{total.multivariance}, \link{m.multivariance}}, respectively.
#'
#' @return A list with elements
#' \describe{
#'   \item{\code{resampled}}{the (total/m-)multivariances of the resampled data,}
#'   \item{\code{original}}{the (total/m-)multivariance of the original data,}
#'   \item{\code{p.value}}{the p-value of the original data, computed using the resampled data}
#' }
#'
#' @references
#' For the theoretic background see the reference [3] given on the main help page of this package: \link{multivariance-package}.
#'
#' @examples
#' re.m = resample.multivariance(matrix(rnorm(30*2),nrow = 30),
#'                         Nscale = TRUE,type= "multi",times = 300)$resampled
#' curve(ecdf(re.m)(x), xlim = c(0,4),main = "empirical distribution of the test statistic under H_0")
#' @export
resample.multivariance = function(x,vec = 1:ncol(x),times = 300,type = "multi",resample.type = "permutation",...) {
  # formerly required arguments: Nscale = TRUE,
  N = length(x[,1])
  n = length(x[1,])
  doAll = FALSE
  switch(resample.type,
         #distinct.permutation = {resample = function() matrix(x[derangements.without.fixpoint(N,n, distinctcols,vec)],ncol = n)},
         bootstrap = {resample = function() sample.cols(x,vec)},
         permutation = {resample = function() sample.cols(x,vec,replace =FALSE)}
  )

  switch(type,
         multi = {fun = function (x) multivariance(x,vec = vec,...)},
         total = {fun = function (x) total.multivariance(x,vec = vec,...)},
         m.multi.2 = {fun = function (x) m.multivariance(x,vec = vec,...)},
         m.multi.3 = {fun = function (x) m.multivariance(x,vec = vec,m = 3,...)},
         all = {fun = function (x) multivariances.all(x,vec = vec,...)}#doAll = TRUE}
  )

  results = matrix(,nrow = times, ncol = 1 + (type == "all")*3)
  for (i in 1:times) {# we use a for loop instead of replicate to prevent trouble with '...'
    results[i,] = fun(resample())
  }

  multi = fun(x)
  invisible(list(resampled = results,
                 original = multi,
                 p.value = rowSums(t(results) >= multi)/times ))

}


#' rejection level via resampling
#'
#' Uses the resample method to sample from the test statistic under the hypothesis of independence. The alpha quantile of these samples is returned.
#'
#' @param alpha numeric, the significance value
#' @param ... passed to \code{\link{resample.multivariance}}. Required is the data matrix \code{x}.
#'
#'@references
#' For the theoretic background see the reference [3] given on the main help page of this package: \link{multivariance-package}.
#'
#' @examples
#' resample.rejection.level(0.05,matrix(rnorm(30*2),nrow = 30))
#' resample.rejection.level(0.05,matrix(rnorm(30*3),nrow = 30),vec = c(1,1,2))
#'
#' @export
resample.rejection.level = function(alpha = 0.05,...) {
  samples = resample.multivariance(...)$resampled
  stats::quantile(samples,probs= 1-alpha)
}

#' p-value via resampling
#'
#' Use a resampling method to generate samples of the test statistic under the hypothesis of independence. Based on these the p.value of a given value of a test statistic is computed.
#'
#' @return It returns 1 minus the value of the empirical distribution function of the resampling samples evaluated at the given value is returned.

#' @param value numeric, the value of (total-/m-)multivariance for which the p-value shall be computed
#' @param ... passed to \code{\link{resample.multivariance}}. Required is the data matrix \code{x}.
#'
#' @references
#' For the theoretic background see the reference [3] given on the main help page of this package: \link{multivariance-package}.
#'
#' @examples
#' x = coins(100)
#' resample.pvalue(multivariance(x),x=x,times = 300)
#'
#' @export
resample.pvalue = function(value,...){
  samples = resample.multivariance(...)$resampled
  sum(samples >= value) / length(samples)
#slower:  1-stats::ecdf(samples)(value)
}



##### Dependence structure ####



# * Data ####

#' example dataset for \code{\link{dependence.structure}}
#'
#' It was generated by \preformatted{
#' set.seed(1348879148)
#' N = 100
#' dep_struct_several_26_100 = cbind(coins(N,2),tetrahedron(N),coins(N,4),
#'     tetrahedron(N),tetrahedron(N),coins(N,3),coins(N,3),rnorm(N))
#'save(dep_struct_several_26_100,file ="dep_struct_several_26_100.rda")
#'}
#'
#' To avoid irritation, note that the seed is just a simple integer hash value of the variable name.
#'
#' @format \code{matrix} 26 variables (columns), 100 independent samples (rows)
#'
"dep_struct_several_26_100"

#' example dataset for \code{\link{dependence.structure}}
#'
#' It was generated by \preformatted{
#' set.seed(222454572)
#' N = 100
#' y = coins(N,2)
#' dep_struct_star_9_100 = cbind(y,y,y)
#' save(dep_struct_star_9_100,file ="dep_struct_star_9_100.rda")
#'}
#'
#' To avoid irritation, note that the seed is just a simple integer hash value of the variable name.
#'
#' @format \code{matrix} 9 variables (columns), 100 independent samples (rows)
#'
"dep_struct_star_9_100"

#' example dataset for \code{\link{dependence.structure}}
#'
#' It was generated by \preformatted{
#' set.seed(532333356)
#' N = 100
#' x = matrix(sample.int(2,10*N,replace = TRUE)-1,ncol = 10)
#' for (i in c(2,5,9)) x = cbind(x,(rowSums(as.matrix(x[,1:(i-1)])) %% 2) == x[,i])
#' dep_struct_iterated_13_100 = x
#' save(dep_struct_iterated_13_100,file ="dep_struct_iterated_13_100.rda")
#'}
#'
#' To avoid irritation, note that the seed is just a simple integer hash value of the variable name.
#'
#' @format \code{matrix} 13 variables (columns), 100 independent samples (rows)
#'
"dep_struct_iterated_13_100"

#' example dataset for \code{\link{dependence.structure}}
#'
#' It was generated by \preformatted{
#' set.seed(436646700)
#' N = 100
#' n= 15
#' x=matrix(sample.int(2,N*n,replace = TRUE)-1,nrow =N)
#' x[,4] = rowSums(x[,1:3]) %% 2
#' x[,7] = rowSums(x[,4:6]) %% 2
#' x[,10] = rowSums(x[,7:9]) %% 2
#' x[,13] = rowSums(x[,10:12]) %% 2
#' x[,15] = rowSums(x[,c(13,14,1)]) %% 2
#' dep_struct_ring_15_100 = x
#' save(dep_struct_ring_15_100,file ="dep_struct_ring_15_100.rda")
#'}
#'
#' To avoid irritation, note that the seed is just a simple integer hash value of the variable name.
#'
#' @format \code{matrix} 13 variables (columns), 100 independent samples (rows)
#'
"dep_struct_ring_15_100"


# * detection ####



#' determines the dependence structure
#'
#' Determines the dependence structure as described in [3].
#'
#' @param x matrix, each row of the matrix is treated as one sample
#' @param vec vector, it indicates which columns are initially treated together as one sample
#' @param verbose boolean, if \code{TRUE} details are printed during the detection and whenever a cluster is newly detected the (so far) detected dependence structure is plotted.
#' @param detection.aim \code{=NULL} or a list of vectors which indicate the expected detection, see below for more details
#' @param ... these are passed to \code{\link{find.cluster}}
#'
#' @details
#' Performs the detection of the dependence structure as described in [3].
#'
#' If \code{fixed.rejection.level} is not provided, the significance level \code{alpha} is used to determine which multivariances are significant using the distribution-free rejection level. As default the Holm method is used for p-value correction corresponding to multiple testing.
#'
#' The resulting graph can be simplified (pairwise dependence can be represented by edges instead of vertices) using \code{\link{clean.graph}}.
#'
#' Advanced:
#' The argument \code{detection.aim} can be used to check, if an expected dependence structure was detected. This might be useful for simulation studies to determine the empirical power of the detection algorithm. Hereto  \code{detection.aim} is set to a list of vectors which indicate the expected detected dependence structures (one for each run of \code{\link{find.cluster}}). The vector has as first element the \code{k} for which k-tuples are detected (for this aim the detection stops without success if no k-tuple is found), and the other elements, indicate to which clusters all present vertices belong after the detection, e.g. \code{c(3,2,2,1,2,1,1,2,1)} expects that 3-tuples are detected and in the graph are 8 vertices (including those representing the detected 3 dependences), the order of the 2's and 1's indicate which vertices belong to which cluster. If \code{detection.aim} is provided, the vector representing the actual detection is printed, thus one can use the output with copy-paste to fix successively the expected detection aims.
#'
#' Note that a failed detection might invoce the warning:
#' \preformatted{
#' run$mem == detection.aim[[k]][-1] :
#' longer object length is not a multiple of shorter object length
#' }
#'
#'
#'
#' @return returns a list with elements:
#' \describe{
#'   \item{\code{multivariances}}{calculated multivariances,}
#'   \item{\code{cdms}}{calculated centered distance matrices,}
#'   \item{\code{graph}}{graph representing the dependence structure.}
#'   \item{\code{detected}}{boolean, this is only included if a \code{detection.aim} is given.}
#' }
#'
#' @references
#' For the theoretic background see the reference [3] given on the main help page of this package: \link{multivariance-package}.
#'
#' @example inst/examples/dependence-structures.R
#' @export
#'
dependence.structure = function(x, vec = 1:ncol(x), verbose = TRUE, detection.aim = NULL,  ...) {

  array.cdm = cdms(x,vec = vec) # creates the distance matrices

  all.multivariances = numeric(0) # vector which will contain all distance multivariances which are calculated

  mem = as.numeric(1:max(vec))
  #its length is the number of vertices, its content is the number of the corresponding cluster for the current iteration!!!
  # has to be numeric, since otherwise 'identical' fails to end the loop (in the case of
  # no detected clusters in the first run)

  cluster.to.vertex = 1:max(mem) # cluster to vertex relation - gets renewed each iteration (since the names of the clusters change)

  vertex.to.cdm = 1:max(mem) # vertex to A (the centered distance matrices) relation - gets appended each iteration


  previous.n.o.cdms = rep(0,max(mem)) # number of As in the previous iteration


  n = max(mem) # number of clusters
  #n = length(array.cdm[,1,1]) # number of clusters

  g = igraph::graph.empty(,directed=FALSE)
  g = igraph::add.vertices(g,n,label = sapply(1:max(mem),function(r) paste(colnames(x,do.NULL = FALSE,prefix = "")[vec == r],collapse = ",")),shape = "circle")

  # Loop through the tuples
  detected = TRUE
  k = 1
  while (detected) {
    if (!is.null(detection.aim)) {
      run = find.cluster(x,vec,array.cdm,mem,cluster.to.vertex,vertex.to.cdm,previous.n.o.cdms,all.multivariances,g,kvec = 2:detection.aim[[k]][1], verbose = verbose, ...)
      if (verbose) {
        cat("last detected structure (in detection.aim format): ")
        dput(c(run$k,run$mem))
      }
      success = all(run$mem == detection.aim[[k]][-1])
      k = k+1
    } else {
      run = find.cluster(x,vec,array.cdm,mem,cluster.to.vertex,vertex.to.cdm,previous.n.o.cdms,all.multivariances,g,...)
    }

    detected = run$detected
    array.cdm = run$array.cdm
    mem = run$mem
    cluster.to.vertex = run$cluster.to.vertex
    vertex.to.cdm = run$vertex.to.cdm
    previous.n.o.cdms = run$previous.n.o.cdms
    all.multivariances = run$all.multivariances
    g = run$g

    if (!is.null(detection.aim)) if (!success) break
  }

  if (!is.null(detection.aim)) {
    return(invisible(list(cdms = run$array.cdm, multivariances = run$all.multivariances, graph = run$g,detected = success)))
  } else {
    return(invisible(list(cdms = run$array.cdm, multivariances = run$all.multivariances, graph = run$g)))
  }
}

#' cluster detection
#'
#' Performs the detection of dependence structures algorithm until a cluster is found. This function is the basic building block \code{\link{dependence.structure}}. Advanced users, might use it directly.
#'
#' @param x matrix with the samples
#'
#' @param vec vector, it indicates which columns are initially treated together as one sample
#' @param array.cdm array of centered distance matrices
#' @param mem numeric vector, its length is the number of vertices, its content is the number of the corresponding cluster for the current iteration, i.e., vertex \code{i} belongs to cluster \code{mem[i]}
#' @param cluster.to.vertex vector, contains the cluster to vertex relations, i.e., \code{cluster.to.vertex[i]} is the index of the vertex which represents cluster \code{i}
#' @param vertex.to.cdm vector, contains the vertex to centered distance matrix relations, i.e., \code{vertex.to.cdm[i]} is the index centered distance matrix in \code{array.cdm} which corresponds to vertex \code{i}
#' @param previous.n.o.cdms vector, number of centered distance matrices in the previous iteration (it is used to ensure that previously check tuples are not checked again)
#' @param all.multivariances vector, which contains all distance multivariances which have been calculated so far. Only used to finally return all distance multivariances which have been calculated.
#' @param g dependence structure graph
#'                         fixed.rejection.level = NA, alpha=0.05,method = "holm",explore = FALSE, verbose = TRUE, kvec = 2:max(mem)
#' @param alpha numeric, significance level used for the (distribution-free) tests
#' @param fixed.rejection.level vector, if not \code{NA} the \code{fixed.rejection.level[k]} is used for the k-tuples, instead of a level derived from the significance level \code{alpha}
#' @param p.adjust.method name of the method used to adjust the p-values for multiple testing, see \code{\link[stats]{p.adjust}} for all possible options.
#' @param verbose boolean, if \code{TRUE} details during the detection are printed and whenever a cluster is newly detected the (so far) detected dependence structure is plotted.
#' @param kvec vector, k-tuples are only checked for each k in \code{kvec}, i.e., for \code{kvec = 2:4} only 2,3 and 4-tuples would be check and then the algorithm stops.
#'
#' @details
#' For further details see \code{\link{dependence.structure}}.
#'
find.cluster = function(x,
                        vec = 1:ncol(x), # which columns should be treated as one sample
                        array.cdm = cdms(x,vec = vec), # creates the distance matrices
                        mem = as.numeric(1:max(vec)),
                        #its length is the number of vertices, its content is the number of the corresponding cluster for the current iteration!!!
                        # has to be numeric, since otherwise 'identical' fails to end the loop (in the case of
                        # no detected clusters in the first run)
                        cluster.to.vertex = 1:max(mem), # cluster to vertex relation - gets renewed each iteration (since the names of the clusters change)
                        vertex.to.cdm = 1:max(mem), # vertex to A (the centered distance matrices) relation - gets appended each iteration
                        previous.n.o.cdms = rep(0,max(mem)), # number of As in the iteration before. it is used to speed up the detection.
                        all.multivariances = numeric(0), # vector which will contain all distance multivariances which are calculated
                        g = igraph::add.vertices(igraph::graph.empty(,directed=FALSE),max(mem),label = sapply(1:max(mem),function(r) paste(colnames(x,do.NULL = FALSE,prefix = "")[vec == r],collapse = ",")),shape = "circle"), #the graph
                        fixed.rejection.level = NA, alpha=0.05,p.adjust.method = "holm", verbose = TRUE, kvec = 2:max(mem)) {
  explore = FALSE # undocumented option, which would provide some more graphs during the detection
  if (verbose) graphics::plot(g)

  n = max(mem) # number of clusters
  nV = length(igraph::V(g)) #length(mem) # number of vertices at the start of the iteration
  n.o.cdm = length(array.cdm[,1,1]) # number of As at the start of the iteration

  cluster.to.cdm = vertex.to.cdm[cluster.to.vertex] #cluster to A relation - gets renewed each iteration
  # Each cluster is represented by its 'largest' vertex

  cdm.to.vertex = NA # A to vertex relation
  for (i in 1:n.o.cdm) cdm.to.vertex[i] = which(vertex.to.cdm == i)

  for (k in 2:min(max(kvec),max(mem))) { # look at the k-tuples of the n variables.

    tuples = utils::combn(n,k) # all k-tuples of 1,..,n

    tuples = matrix(cluster.to.cdm[tuples],ncol = k,byrow = TRUE)
    #transform tuples into the A indexes

    tuples = tuples[apply(tuples,1,max) > previous.n.o.cdms[k],]
    # to speed up, we only consider those with new A

    if (length(tuples) == k) dim(tuples) = c(1,k)
    # make sure that it is a matrix

    multivariances = apply(tuples,1,function(x) multivariance(array.cdm,x,Nscale = TRUE)) #calculates all distance multivariances

    all.multivariances = c(all.multivariances,multivariances)
    # print(multivariances)
    if (explore) {
      graphics::plot(multivariances,main = paste(k,"tuple multivariances"))
      graphics::abline(h = c( rejection.level(alpha/choose(n,k)), rejection.level(alpha)), col = c("red","green"))
      readline("continue with [enter]")
    }

    for (i in which(((anyNA(fixed.rejection.level) & (stats::p.adjust(multivariance.pvalue(multivariances),method = p.adjust.method) < alpha))) | (multivariances > fixed.rejection.level[k]) )) {
      # for each tuple, with adjusted p value less than the significance level (or multivariance less than a prescribed fixed.rejection.level, if given) we add a vertex and the edges

      new = length(igraph::V(g))+1
      g = igraph::add.vertices(g,1,label = signif(multivariances[i],4), shape = "none", level=k)

      g = igraph::add_edges(g, as.vector(t(cbind(new,cdm.to.vertex[tuples[i,]]))), weight= NA, color = k)
      #  }
    }
    if (verbose) cat(paste(k,"-tuples: max. multivariance: ",max(multivariances),"; min. p-value: ",multivariance.pvalue(max(multivariances)),"\n",sep =""))
    #readline(paste("level",k,", press [Enter] for next level"))

    previous.n.o.cdms[k] = n.o.cdm
    if ((anyNA(fixed.rejection.level) && (stats::p.adjust(multivariance.pvalue(max(multivariances)),method = p.adjust.method,n = length(multivariances)) < alpha)) || (!anyNA(fixed.rejection.level) && (max(multivariances) > fixed.rejection.level[k]))) {
      #if a cluster was found exit the loop
      break
    }
  } # end loop k


  previous.n.o.cdms[(k+1):length(previous.n.o.cdms)] = 0 # ensures that for all higher tuples all combinations are checked. (Otherwise those tuples which were checked previously are skipped)

  if (verbose) graphics::plot(g)

  #print(paste("mem: ",paste(mem,collapse = " "),"clust",paste(igraph::clusters(g)$membership,collapse = " ")))
  if (identical(mem, igraph::clusters(g)$membership) || (igraph::clusters(g)$no == 1) ) {
    detected = FALSE
    if (verbose) {
      if (igraph::clusters(g)$no == 1) {cat("All vertices are in one cluster.\n")} else { cat("No new cluster detected.\n")}
    }    # will end recursion if the membership did not change or there is only one cluster
  } else {
    if (verbose) cat("New cluster detected and not all vertices are in this cluster.\n")
    detected = TRUE
  }
  mem = igraph::clusters(g)$membership


  cluster.to.vertex = NA
  for (i in 1:max(mem)) cluster.to.vertex[i] = max(which(mem == i))
  # cluster to vertex relation. Each cluster is represented by its 'largest' vertex

  previous.n.o.cdms[2] = n.o.cdm

  for (i in cluster.to.vertex[cluster.to.vertex > nV]) {
    # for every vertex which represents a new cluster the distance matrix representing the cluster is calculated.
    #  the vertex i gets the cdm with index vertex.to.cdm[i]
    n.o.cdm = n.o.cdm + 1
    vertex.to.cdm[i] = n.o.cdm
    # print(which(mem[1:ncol(x)] == which(cluster.to.vertex== i)))
    #      array.cdm = abind::abind(array.cdm, cdm(x[,which(mem[1:ncol(x)] == which(cluster.to.vertex== i))]),along=1)
    array.cdm = abind::abind(array.cdm, cdm(x[,which(vec %in% which(mem[1:ncol(x)] == which(cluster.to.vertex== i)))]),along=1)
    # which(cluster.to.vertex== i) is the number of the cluster represented by vertex i
    # which(mem ...) gives the "original" vertices which belong to the same cluster as vertex i
    # which(vec ...) finally gives the columns of x which belong to the same cluster as vertex i

  } # at the end of this for loop, n.o.cdm contains the new number of As


  invisible(list(detected = detected,array.cdm = array.cdm,mem = mem,cluster.to.vertex = cluster.to.vertex,vertex.to.cdm = vertex.to.cdm,previous.n.o.cdms = previous.n.o.cdms,all.multivariances = all.multivariances,g = g, k = k))
}

#' cleanup dependence structure graph
#'
#' Given a dependence structure graph: vertices representing the multivariances of only two verticies become an edge labeled with the label of the vertex.
#'
#' @param g graph, created by \code{\link{dependence.structure}}
#' @return graph
#'
#' @examples
#'
#' N = 200
#' y = coins(N,2)
#' x = cbind(y,y,y)
#' ds = dependence.structure(x)
#' plot(clean.graph(ds$graph))
#' @export
clean.graph = function(g) {
  #  g = ds$graph
  vert = which(igraph::V(g)$level == 2) # might still have more neighbors!!!
  only.two.neighbors = logical(igraph::vcount(g))
  for (i in vert) {
    only.two.neighbors[i] = length(igraph::neighbors(g,igraph::V(g)[i]))==2
    if (only.two.neighbors[i]) g = igraph::add_edges(g, igraph::neighbors(g,i), weight= NA, label = igraph::V(g)[i]$label, color = 2)
  }

  return(igraph::delete.vertices(g,only.two.neighbors))
}

# Utility functions ####


#' Transforms a distance matrix to a matrix
#'
#' Does for a distance matrix generated via \code{dist} the same as \code{as.matrix} only slightly faster.
#'
#' @param ds a distance matrix object, e.g. generated by \code{\link{dist}}
#'
#' @keywords internal
dist.to.matrix = function(ds) {
  N=attr(ds,"Size")
  m = matrix(0, nrow = N, ncol = N)
  m[outer(1:N,1:N,">")] = ds
  m+t(m)
}
