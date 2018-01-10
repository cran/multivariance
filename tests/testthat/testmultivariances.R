set.seed(123412)
N = 50
n = 3
x = matrix(rnorm(N*n),nrow = N)

y = matrix(rnorm(N*4),nrow = N)
vec = c(1,1,2,3)

context("centered distance matrices")
test_that("cdms",{
  expect_warning(cdm(matrix(1,nrow= N,ncol = n)),"constant")
  expect_equal(cdm(x),cdm(x,psi = function(x,y) sqrt(sum((x-y)^2))))
})

context("definition of multivariances")
test_that("multivariance, total.multivariance, m.multivariance", {
  expect_warning(multivariance(matrix(1,nrow= N,ncol = n)),"constant")

  expect_equal(multivariance(x), multivariance(x[,c(2,3,1)]))
  expect_equal(multivariance(x), m.multivariance(x,m=3))
  expect_equal(multivariance(x[,c(1,2)]), m.multivariance(x[,c(1,2)],m=2))
  expect_equal((multivariance(x[,c(1,2)]) + multivariance(x[,c(1,3)]) + multivariance(x[,c(3,2)]))/3, m.multivariance(x,m=2))
  expect_equal(multivariance(x[,c(1,2)]), total.multivariance(x[,c(1,2)]))
  expect_equal(total.multivariance(x), (m.multivariance(x,m=3)+m.multivariance(x,m=2)*3)/4)

  expect_equal(multivariances.all(x),c(multivariance(x),total.multivariance(x),m.multivariance(x,m=2),m.multivariance(x,m=3)))
  expect_equal(multivariances.all(y,vec),c(multivariance(y,vec),total.multivariance(y,vec),m.multivariance(y,vec,m=2),m.multivariance(y,vec,m=3)))
})

context("resampling")
set.seed(1)
quant = quantile(resample.multivariance(x)$resampled,0.95)
pval = sum(resample.multivariance(x)$resampled>=multivariance(x))/300

set.seed(1)
test_that("resampling p-values and quantiles",{
  expect_equal(resample.rejection.level(0.05,x), quant)
  expect_equal(resample.pvalue(multivariance(x),x), pval)
})
