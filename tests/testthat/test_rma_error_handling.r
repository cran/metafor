### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

context("Checking that rma() handles various errors properly")

test_that("rma() handles NAs correctly.", {

   data(dat.bcg, package="metafor")
   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

   dat$yi[1] <- NA
   dat$yi[2] <- NA

   res <- rma(yi, vi, data=dat, digits=3)

   expect_that(res$k, equals(11))
   expect_that(res$k.f, equals(13))
   expect_that(length(res$yi), equals(11))
   expect_that(length(res$yi.f), equals(13))
   expect_that(res$not.na, equals(rep(c(FALSE,TRUE),times=c(2,11))))

   dat$ablat[3] <- NA

})
