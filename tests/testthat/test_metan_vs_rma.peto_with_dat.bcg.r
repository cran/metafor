context("Comparing rma.peto() against metan with 'dat.bcg'")

### library(metafor); library(testthat)

test_that("results match (FE model, measure='OR').", {

   data(dat.bcg, package="metafor")

   ### compare results with: metan tpos tneg cpos cneg, peto nograph or log

   res <- rma.peto(ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

   expect_that(round(as.vector(res$b),  digits=3), equals(-0.474))
   expect_that(round(res$ci.lb, digits=3), equals(-0.554))
   expect_that(round(res$ci.ub, digits=3), equals(-0.395))
   expect_that(round(res$zval,  digits=2), equals(-11.67)) ### 11.67 in Stata
   expect_that(round(res$QE,    digits=2), equals(167.73))

   ### compare results with: metan tpos tneg cpos cneg, peto nograph or

   sav <- predict(res, transf=exp)

   expect_that(round(c(sav$pred),  digits=3), equals(0.622))
   expect_that(round(c(as.vector(sav$ci.lb)), digits=3), equals(0.575))
   expect_that(round(c(as.vector(sav$ci.ub)), digits=3), equals(0.674))

})
