### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

context("Comparing rma.uni() against direct computations")

test_that("results match (FE model).", {

   data(dat.bcg, package="metafor")
   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

   res <- rma(yi, vi, mods = ~ ablat + year, data=dat, method="FE")

   X <- cbind(1, dat$ablat, dat$year)
   W <- diag(1/dat$vi)
   y <- cbind(dat$yi)

   b  <- solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% y
   vb <- solve(t(X) %*% W %*% X)

   expect_that(c(res$b), equals(c(b)))
   expect_that(unname(res$vb), equals(vb))

   yhat <- c(X %*% b)

   expect_that(as.vector(fitted(res)), equals(yhat))

   H <- X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W

   expect_that(unname(hatvalues(res, type="matrix")), equals(H))

   ei <- c((diag(res$k) - H) %*% y)

   expect_that(as.vector(resid(res)), equals(ei))

})
