context("Checking that results for rma() match up with those from lm()")

### this is essentially checking the equivalence of the results as explained here:
### http://www.metafor-project.org/doku.php/tips:regression_with_rma

### library(metafor); library(testthat)

test_that("results for rma() and lm() match for method='FE'.", {

   stackloss$vi <- 0

   res.lm  <- lm(stack.loss ~ Air.Flow + Water.Temp + Acid.Conc., data=stackloss)
   res.rma <- rma(stack.loss, vi, mods =  ~ Air.Flow + Water.Temp + Acid.Conc., data=stackloss, knha=TRUE, control=list(REMLf=FALSE))

   ### log likelihood (REML) should be the same
   expect_that(logLik(res.lm, REML=TRUE), equals(logLik(res.rma)))

   ### coefficients should be the same (need to strip names)
   expect_that(as.vector(coef(res.lm)), equals(as.vector(coef(res.rma))))

   ### var-cov matrix should be the same (need to strip names)
   expect_that(matrix(vcov(res.lm), nrow=4, ncol=4), equals(matrix(vcov(res.rma), nrow=4, ncol=4)))

   ### fitted values should be the same
   expect_that(fitted(res.lm), equals(fitted(res.rma)))

   ### standardized residuals should be the same
   expect_that(as.vector(rstandard(res.lm)), equals(rstandard(res.rma)$z))

   ### studentized residuals should be the same
   expect_that(as.vector(rstudent(res.lm)), equals(rstudent(res.rma)$z))

   ### hat values should be the same
   expect_that(as.vector(hatvalues(res.lm)), equals(as.vector(hatvalues(res.rma))))

   ### dffits should be the same
   expect_that(as.vector(dffits(res.lm)), equals(influence(res.rma)$inf$dffits))

   ### covratios should be the same
   expect_that(as.vector(covratio(res.lm)), equals(influence(res.rma)$inf$cov.r))

   ### dfbetas should be the same
   expect_that(unname(as.matrix(dfbetas(res.lm))), equals(unname(as.matrix(dfbetas(res.rma)))))

   ### Cook's distancs should differ by a factor of p
   expect_that(as.vector(cooks.distance(res.lm)), equals(as.vector(cooks.distance(res.rma)/res.rma$p)))

})
