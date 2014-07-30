context("Comparing rma.uni() against metan with 'dat.bcg'")

### library(metafor); library(testthat)

test_that("results match (FE model, measure='RR').", {

   data(dat.bcg, package="metafor")
   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

   ### compare results with: metan tpos tneg cpos cneg, fixedi nograph rr log

   res <- rma(yi, vi, data=dat, method="FE")

   expect_that(round(c(res$b),  digits=3), equals(-0.430))
   expect_that(round(res$ci.lb, digits=3), equals(-0.510))
   expect_that(round(res$ci.ub, digits=3), equals(-0.351))
   expect_that(round(res$zval,  digits=2), equals(-10.62)) ### 10.62 in Stata
   expect_that(round(res$QE,    digits=2), equals(152.23))

   ### compare results with: metan tpos tneg cpos cneg, fixedi nograph rr

   sav <- predict(res, transf=exp)

   expect_that(round(c(sav$pred),  digits=3), equals(0.650))
   expect_that(round(c(as.vector(sav$ci.lb)), digits=3), equals(0.601))
   expect_that(round(c(as.vector(sav$ci.ub)), digits=3), equals(0.704))

})

test_that("results match (RE model w/ DL estimator, measure='RR').", {

   data(dat.bcg, package="metafor")
   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

   ### compare results with: metan tpos tneg cpos cneg, randomi nograph rr log

   res <- rma(yi, vi, data=dat, method="DL")

   expect_that(round(c(res$b),  digits=3), equals(-0.714))
   expect_that(round(res$ci.lb, digits=3), equals(-1.064))
   expect_that(round(res$ci.ub, digits=3), equals(-0.364))
   expect_that(round(res$zval,  digits=2), equals(-4.00)) ### 4.00 in Stata
   expect_that(round(res$tau2,  digits=4), equals(.3088))
   expect_that(round(res$I2,    digits=1), equals(92.1))

   ### compare results with: metan tpos tneg cpos cneg, randomi nograph rr

   sav <- predict(res, transf=exp)

   expect_that(round(c(sav$pred),  digits=3), equals(0.490))
   expect_that(round(c(as.vector(sav$ci.lb)), digits=3), equals(0.345))
   expect_that(round(c(as.vector(sav$ci.ub)), digits=3), equals(0.695))

})

test_that("results match (FE model, measure='OR').", {

   data(dat.bcg, package="metafor")
   dat <- escalc(measure="OR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

   ### compare results with: metan tpos tneg cpos cneg, fixedi nograph or log

   res <- rma(yi, vi, data=dat, method="FE")

   expect_that(round(c(res$b),  digits=3), equals(-0.436))
   expect_that(round(res$ci.lb, digits=3), equals(-0.519))
   expect_that(round(res$ci.ub, digits=3), equals(-0.353))
   expect_that(round(res$zval,  digits=2), equals(-10.32)) ### 10.32 in Stata
   expect_that(round(res$QE,    digits=2), equals(163.16))

   ### compare results with: metan tpos tneg cpos cneg, fixedi nograph or

   sav <- predict(res, transf=exp)

   expect_that(round(c(sav$pred),  digits=3), equals(0.647))
   expect_that(round(c(as.vector(sav$ci.lb)), digits=3), equals(0.595))
   expect_that(round(c(as.vector(sav$ci.ub)), digits=3), equals(0.702))

})

test_that("results match (RE model w/ DL estimator, measure='OR').", {

   data(dat.bcg, package="metafor")
   dat <- escalc(measure="OR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

   ### compare results with: metan tpos tneg cpos cneg, randomi nograph or log

   res <- rma(yi, vi, data=dat, method="DL")

   expect_that(round(c(res$b),  digits=3), equals(-0.747))
   expect_that(round(res$ci.lb, digits=3), equals(-1.124))
   expect_that(round(res$ci.ub, digits=3), equals(-0.371))
   expect_that(round(res$zval,  digits=2), equals(-3.89)) ### 3.89 in Stata
   expect_that(round(res$tau2,  digits=4), equals(.3663))
   expect_that(round(res$I2,    digits=1), equals(92.6))

   ### compare results with: metan tpos tneg cpos cneg, randomi nograph or

   sav <- predict(res, transf=exp)

   expect_that(round(c(sav$pred),  digits=3), equals(0.474))
   expect_that(round(c(as.vector(sav$ci.lb)), digits=3), equals(0.325))
   expect_that(round(c(as.vector(sav$ci.ub)), digits=3), equals(0.690))

})

test_that("results match (FE model, measure='RD').", {

   data(dat.bcg, package="metafor")
   dat <- escalc(measure="RD", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

   ### compare results with: metan tpos tneg cpos cneg, fixedi nograph rd

   res <- rma(yi, vi, data=dat, method="FE")

   expect_that(round(c(res$b),  digits=3), equals(-0.001))
   expect_that(round(res$ci.lb, digits=3), equals(-0.001))
   expect_that(round(res$ci.ub, digits=3), equals(-0.000))
   expect_that(round(res$zval,  digits=2), equals(-4.04)) ### 4.04 in Stata
   expect_that(round(res$QE,    digits=2), equals(276.47))

})

test_that("results match (RE model w/ DL estimator, measure='RD').", {

   data(dat.bcg, package="metafor")
   dat <- escalc(measure="RD", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

   ### compare results with: metan tpos tneg cpos cneg, randomi nograph rd

   res <- rma(yi, vi, data=dat, method="DL")

   expect_that(round(c(res$b),  digits=3), equals(-0.007))
   expect_that(round(res$ci.lb, digits=3), equals(-0.010))
   expect_that(round(res$ci.ub, digits=3), equals(-0.004))
   expect_that(round(res$zval,  digits=2), equals(-4.51)) ### 4.51 in Stata
   expect_that(round(res$tau2,  digits=4), equals(.0000))
   expect_that(round(res$I2,    digits=1), equals(95.7))

})

#expect_that(rma(yi ~ ablat, vi, data=dat, subset=1:2), throws_error("Number of parameters to be estimated is larger than the number of observations."))
