### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true"); Sys.setenv(RUN_VIS_TESTS="true")

source("settings.r")

context("Checking misc: deltamethod() function")

test_that("deltamethod() works correctly for a one-dimensional input with one transformation", {

   mu     <- 0.7
   sigma2 <- 0.2

   res1 <- deltamethod(x=mu, vcov=sigma2, fun=function(x) x^3, order=1)
   res2 <- deltamethod(x=mu, vcov=sigma2, fun=function(x) x^3, order=2)
   res3 <- deltamethod(x=mu, vcov=sigma2, fun=function(x) x^3, order=3)

   expected1 <- matrix(9 * mu^4 * sigma2)
   expected2 <- matrix(9 * mu^4 * sigma2 + 18 * mu^2 * sigma2^2)
   expected3 <- matrix(9 * mu^4 * sigma2 + 36 * mu^2 * sigma2^2 + 15 * sigma2^3)

   expect_equivalent(vcov(res1), expected1)
   expect_equivalent(vcov(res2), expected2)
   expect_equivalent(vcov(res3), expected3)

   if (FALSE) {
      xsim <- rnorm(10^8, mean=mu, sd=sqrt(sigma2))
      vsim <- var(xsim^3)
      round(vcov(res3) - vsim, digits=3)
      round(expected3 - vsim, digits=3)
   }

})

test_that("deltamethod() works correctly for a one-dimensional input with two transformations", {

   mu     <- 0.7
   sigma2 <- 0.2

   res1 <- deltamethod(x=mu, vcov=sigma2, fun=function(x) c(x^2, x^3), order=1)
   res2 <- deltamethod(x=mu, vcov=sigma2, fun=function(x) c(x^2, x^3), order=2)
   res3 <- deltamethod(x=mu, vcov=sigma2, fun=function(x) c(x^2, x^3), order=3)

   var11 <- 4 * mu^2 * sigma2
   var12 <- 6 * mu^3 * sigma2
   var22 <- 9 * mu^4 * sigma2
   expected1 <- matrix(c(var11, var12, var12, var22), nrow=2)
   var11 <- 4 * mu^2 * sigma2 + 2 * sigma2^2
   var12 <- 6 * mu^3 * sigma2 + 6 * mu * sigma2^2
   var22 <- 9 * mu^4 * sigma2 + 18 * mu^2 * sigma2^2
   expected2 <- matrix(c(var11, var12, var12, var22), nrow=2)
   var11 <- 4 * mu^2 * sigma2 + 2 * sigma2^2
   var12 <- 6 * mu^3 * sigma2 + 12 * mu * sigma2^2
   var22 <- 9 * mu^4 * sigma2 + 36 * mu^2 * sigma2^2 + 15 * sigma2^3
   expected3 <- matrix(c(var11, var12, var12, var22), nrow=2)

   expect_equivalent(vcov(res1), expected1)
   expect_equivalent(vcov(res2), expected2)
   expect_equivalent(vcov(res3), expected3)

   if (FALSE) {
      xsim <- rnorm(10^8, mean=mu, sd=sqrt(sigma2))
      vsim <- var(cbind(xsim^2, xsim^3))
      round(vcov(res3) - vsim, digits=3)
      round(expected3 - vsim, digits=3)
      round(cov2cor(vcov(res3)) - cov2cor(vsim), digits=3)
      round(cov2cor(expected3) - cov2cor(vsim), digits=3)
   }

})

test_that("deltamethod() works correctly for a two-dimensional input with one transformation", {

   mu1 <- 0.7
   mu2 <- -0.4
   var1 <- 0.20
   var2 <- 0.35
   covar12 <- 0.08
   mu <- c(mu1, mu2)
   V <- matrix(c(var1, covar12, covar12, var2), nrow=2)

   res1 <- deltamethod(x=mu, vcov=V, fun=function(x) x[1]^2 * x[2], order=1)
   res2 <- deltamethod(x=mu, vcov=V, fun=function(x) x[1]^2 * x[2], order=2)
   res3 <- deltamethod(x=mu, vcov=V, fun=function(x) x[1]^2 * x[2], order=3)

   expected1 <- matrix(mu1^4 * var2 + 4 * mu1^3 * mu2 * covar12 + 4 * mu1^2 * mu2^2 * var1)
   expected2 <- matrix(mu1^4 * var2 + 4 * mu1^3 * mu2 * covar12 + 4 * mu1^2 * mu2^2 * var1 + 4 * mu1^2 * covar12^2 + 8 * mu1 * mu2 * covar12 * var1 + 4 * mu1^2 * var1 * var2 + 2 * mu2^2 * var1^2)
   expected3 <- matrix(mu1^4 * var2 + 4 * mu1^3 * mu2 * covar12 + 4 * mu1^2 * mu2^2 * var1 + 8 * mu1^2 * covar12^2 + 6 * mu1^2 * var1 * var2 + 20 * mu1 * mu2 * covar12 * var1 + 2 * mu2^2 * var1^2 + 12 * covar12^2 * var1 + 3 * var1^2 * var2)

   expect_equivalent(vcov(res1), expected1)
   expect_equivalent(vcov(res2), expected2)
   expect_equivalent(vcov(res3), expected3)

   res1 <- deltamethod(x=mu, vcov=V, fun=function(x1, x2) x1^2 * x2, order=1)
   res2 <- deltamethod(x=mu, vcov=V, fun=function(x1, x2) x1^2 * x2, order=2)
   res3 <- deltamethod(x=mu, vcov=V, fun=function(x1, x2) x1^2 * x2, order=3)

   expect_equivalent(vcov(res1), expected1)
   expect_equivalent(vcov(res2), expected2)
   expect_equivalent(vcov(res3), expected3)

   if (FALSE) {
      xsim <- metafor:::.mvrnorm(10^8, mu, V)
      vsim <- var(xsim[,1]^2 * xsim[,2])
      round(vcov(res3) - vsim, digits=3)
      round(expected3 - vsim, digits=3)
   }

})

test_that("deltamethod() works correctly for a two-dimensional input with two transformations", {

   mu1 <- 0.7
   mu2 <- -0.4
   var1 <- 0.20
   var2 <- 0.35
   covar12 <- 0.08
   mu <- c(mu1, mu2)
   V <- matrix(c(var1, covar12, covar12, var2), nrow=2)

   res1 <- deltamethod(x=mu, vcov=V, fun=function(x) c(x[1]^2 * x[2], x[1] * x[2]^2), order=1)
   res2 <- deltamethod(x=mu, vcov=V, fun=function(x) c(x[1]^2 * x[2], x[1] * x[2]^2), order=2)
   res3 <- deltamethod(x=mu, vcov=V, fun=function(x) c(x[1]^2 * x[2], x[1] * x[2]^2), order=3)

   var11 <- mu1^4 * var2 + 4 * mu1^3 * mu2 * covar12+ 4 * mu1^2 * mu2^2 * var1
   var22 <- mu2^4 * var1 + 4 * mu1 * mu2^3 * covar12+ 4 * mu1^2 * mu2^2 * var2
   var12 <- 2 * mu1^3 * mu2 * var2 + 2 * mu1 * mu2^3 * var1 + 5 * mu1^2 * mu2^2 * covar12
   expected1 <- matrix(c(var11, var12, var12, var22), nrow=2)
   var11 <- mu1^4 * var2 + 4 * mu1^3 * mu2 * covar12+ 4 * mu1^2 * mu2^2 * var1 + 4 * mu1^2 * covar12^2 + 8 * mu1 * mu2 * var1 * covar12+ 4 * mu1^2 * var1 * var2 + 2 * mu2^2 * var1^2
   var22 <- mu2^4 * var1 + 4 * mu1 * mu2^3 * covar12+ 4 * mu1^2 * mu2^2 * var2 + 4 * mu2^2 * covar12^2 + 8 * mu1 * mu2 * var2 * covar12+ 4 * mu2^2 * var1 * var2 + 2 * mu1^2 * var2^2
   var12 <- 2 * mu1^3 * mu2 * var2 + 2 * mu1 * mu2^3 * var1 + 5 * mu1^2 * mu2^2 * covar12+ 6 * mu1 * mu2 * covar12^2 + 4 * mu1^2 * var2 * covar12+ 4 * mu2^2 * var1 * covar12+ 4 * mu1 * mu2 * var1 * var2
   expected2 <- matrix(c(var11, var12, var12, var22), nrow=2)
   var11 <- mu1^4 * var2 + 4 * mu1^3 * mu2 * covar12+ 4 * mu1^2 * mu2^2 * var1 + 8 * mu1^2 * covar12^2 + 6 * mu1^2 * var1 * var2 + 20 * mu1 * mu2 * var1 * covar12+ 2 * mu2^2 * var1^2 + 12 * var1 * covar12^2 + 3 * var1^2 * var2
   var22 <- mu2^4 * var1 + 4 * mu1 * mu2^3 * covar12+ 4 * mu1^2 * mu2^2 * var2 + 8 * mu2^2 * covar12^2 + 6 * mu2^2 * var1 * var2 + 20 * mu1 * mu2 * var2 * covar12+ 2 * mu1^2 * var2^2 + 12 * var2 * covar12^2 + 3 * var1 * var2^2
   var12 <- 2 * mu1^3 * mu2 * var2 + 2 * mu1 * mu2^3 * var1 + 5 * mu1^2 * mu2^2 * covar12+ 7 * mu1^2 * var2 * covar12+ 7 * mu2^2 * var1 * covar12+ 8 * mu1 * mu2 * var1 * var2 + 14 * mu1 * mu2 * covar12^2 + 9 * var1 * var2 * covar12+ 6 * covar12^3
   expected3 <- matrix(c(var11, var12, var12, var22), nrow=2)

   expect_equivalent(vcov(res1), expected1)
   expect_equivalent(vcov(res2), expected2)
   expect_equivalent(vcov(res3), expected3)

   res1 <- deltamethod(x=mu, vcov=V, fun=function(x1, x2) c(x1^2 * x2, x1 * x2^2), order=1)
   res2 <- deltamethod(x=mu, vcov=V, fun=function(x1, x2) c(x1^2 * x2, x1 * x2^2), order=2)
   res3 <- deltamethod(x=mu, vcov=V, fun=function(x1, x2) c(x1^2 * x2, x1 * x2^2), order=3)

   expect_equivalent(vcov(res1), expected1)
   expect_equivalent(vcov(res2), expected2)
   expect_equivalent(vcov(res3), expected3)

   if (FALSE) {
      xsim <- metafor:::.mvrnorm(10^8, mu, V)
      vsim <- var(cbind(xsim[,1]^2 * xsim[,2], xsim[,1] * xsim[,2]^2))
      round(vcov(res3) - vsim, digits=3)
      round(expected3 - vsim, digits=3)
      round(cov2cor(vcov(res3)) - cov2cor(vsim), digits=3)
      round(cov2cor(expected3) - cov2cor(vsim), digits=3)
   }

})

test_that("deltamethod() works correctly for a two-dimensional input with three transformations", {

   mu1 <- 0.7
   mu2 <- -0.4
   var1 <- 0.20
   var2 <- 0.35
   covar12 <- 0.08
   mu <- c(mu1, mu2)
   V <- matrix(c(var1, covar12, covar12, var2), nrow=2)

   res1 <- deltamethod(x=mu, vcov=V, fun=function(x) c(x[1]^2 * x[2], x[1] * x[2]^2, x[1]^3), order=1)
   res2 <- deltamethod(x=mu, vcov=V, fun=function(x) c(x[1]^2 * x[2], x[1] * x[2]^2, x[1]^3), order=2)
   res3 <- deltamethod(x=mu, vcov=V, fun=function(x) c(x[1]^2 * x[2], x[1] * x[2]^2, x[1]^3), order=3)

   var11 <- mu1^4 * var2 + 4 * mu1^3 * mu2 * covar12+ 4 * mu1^2 * mu2^2 * var1
   var22 <- mu2^4 * var1 + 4 * mu1 * mu2^3 * covar12+ 4 * mu1^2 * mu2^2 * var2
   V33 <- 9 * mu1^4 * var1
   var12 <- 2 * mu1^3 * mu2 * var2 + 2 * mu1 * mu2^3 * var1 + 5 * mu1^2 * mu2^2 * covar12
   var13 <- 3 * mu1^4 * covar12+ 6 * mu1^3 * mu2 * var1
   var23 <- 6 * mu1^3 * mu2 * covar12+ 3 * mu1^2 * mu2^2 * var1
   expected1 <- matrix(c(var11, var12, var13, var12, var22, var23, var13, var23, V33), nrow=3)
   var11 <- mu1^4 * var2 + 4 * mu1^3 * mu2 * covar12+ 4 * mu1^2 * mu2^2 * var1 + 4 * mu1^2 * covar12^2 + 8 * mu1 * mu2 * var1 * covar12+ 4 * mu1^2 * var1 * var2 + 2 * mu2^2 * var1^2
   var22 <- mu2^4 * var1 + 4 * mu1 * mu2^3 * covar12+ 4 * mu1^2 * mu2^2 * var2 + 4 * mu2^2 * covar12^2 + 8 * mu1 * mu2 * var2 * covar12+ 4 * mu2^2 * var1 * var2 + 2 * mu1^2 * var2^2
   V33 <- 9 * mu1^4 * var1 + 18 * mu1^2 * var1^2
   var12 <- 2 * mu1^3 * mu2 * var2 + 2 * mu1 * mu2^3 * var1 + 5 * mu1^2 * mu2^2 * covar12+ 6 * mu1 * mu2 * covar12^2 + 4 * mu1^2 * var2 * covar12+ 4 * mu2^2 * var1 * covar12+ 4 * mu1 * mu2 * var1 * var2
   var13 <- 3 * mu1^4 * covar12+ 6 * mu1^3 * mu2 * var1 + 12 * mu1^2 * var1 * covar12+ 6 * mu1 * mu2 * var1^2
   var23 <- 6 * mu1^3 * mu2 * covar12+ 3 * mu1^2 * mu2^2 * var1 + 6 * mu1^2 * covar12^2 + 12 * mu1 * mu2 * var1 * covar12
   expected2 <- matrix(c(var11, var12, var13, var12, var22, var23, var13, var23, V33), nrow=3)
   var11 <- mu1^4 * var2 + 4 * mu1^3 * mu2 * covar12+ 4 * mu1^2 * mu2^2 * var1 + 8 * mu1^2 * covar12^2 + 6 * mu1^2 * var1 * var2 + 20 * mu1 * mu2 * var1 * covar12+ 2 * mu2^2 * var1^2 + 12 * var1 * covar12^2 + 3 * var1^2 * var2
   var22 <- mu2^4 * var1 + 4 * mu1 * mu2^3 * covar12+ 4 * mu1^2 * mu2^2 * var2 + 8 * mu2^2 * covar12^2 + 6 * mu2^2 * var1 * var2 + 20 * mu1 * mu2 * var2 * covar12+ 2 * mu1^2 * var2^2 + 12 * var2 * covar12^2 + 3 * var1 * var2^2
   V33 <- 9 * mu1^4 * var1 + 36 * mu1^2 * var1^2 + 15 * var1^3
   var12 <- 2 * mu1^3 * mu2 * var2 + 2 * mu1 * mu2^3 * var1 + 5 * mu1^2 * mu2^2 * covar12+ 7 * mu1^2 * var2 * covar12+ 7 * mu2^2 * var1 * covar12+ 8 * mu1 * mu2 * var1 * var2 + 14 * mu1 * mu2 * covar12^2 + 9 * var1 * var2 * covar12+ 6 * covar12^3
   var13 <- 3 * mu1^4 * covar12+ 6 * mu1^3 * mu2 * var1 + 24 * mu1^2 * var1 * covar12+ 12 * mu1 * mu2 * var1^2 + 15 * var1^2 * covar12
   var23 <- 6 * mu1^3 * mu2 * covar12+ 3 * mu1^2 * mu2^2 * var1 + 3 * mu1^2 * var1 * var2 + 12 * mu1^2 * covar12^2 + 18 * mu1 * mu2 * var1 * covar12+ 3 * mu2^2 * var1^2 + 12 * var1 * covar12^2 + 3 * var1^2 * var2
   expected3 <- matrix(c(var11, var12, var13, var12, var22, var23, var13, var23, V33), nrow=3)

   expect_equal(vcov(res1), expected1)
   expect_equal(vcov(res2), expected2)
   expect_equal(vcov(res3), expected3)

   res1 <- deltamethod(x=mu, vcov=V, fun=function(x1, x2) c(x1^2 * x2, x1 * x2^2, x1^3), order=1)
   res2 <- deltamethod(x=mu, vcov=V, fun=function(x1, x2) c(x1^2 * x2, x1 * x2^2, x1^3), order=2)
   res3 <- deltamethod(x=mu, vcov=V, fun=function(x1, x2) c(x1^2 * x2, x1 * x2^2, x1^3), order=3)

   expect_equal(vcov(res1), expected1)
   expect_equal(vcov(res2), expected2)
   expect_equal(vcov(res3), expected3)

   if (FALSE) {
      xsim <- metafor:::.mvrnorm(10^8, mu, V)
      vsim <- var(cbind(xsim[,1]^2 * xsim[,2], xsim[,1] * xsim[,2]^2, xsim[,1]^3))
      round(vcov(res3) - vsim, digits=3)
      round(expected3 - vsim, digits=3)
      round(cov2cor(vcov(res3)) - cov2cor(vsim), digits=3)
      round(cov2cor(expected3) - cov2cor(vsim), digits=3)
   }

})

rm(list=ls())
