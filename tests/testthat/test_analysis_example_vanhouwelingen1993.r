### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

### see also: http://www.metafor-project.org/doku.php/analyses:vanhouwelingen1993

context("Checking analysis example vanhouwelingen1993")

### load data
dat <- get(data(dat.collins1985a, package="metafor"))

test_that("the log likelihood plot can be created.", {

   skip_on_cran()

   llplot(measure="OR", ai=b.xci, n1i=nci, ci=b.xti, n2i=nti, data=dat,
          xlim=c(-4,4), lwd=1, col="black", refline=NA, drop00=FALSE)
   dev.off()

})

test_that("results of the fixed-effects conditional logistic model are correct.", {

   skip_on_cran()

   res <- rma.glmm(measure="OR", ai=b.xci, n1i=nci, ci=b.xti, n2i=nti, data=dat, model="CM.EL", method="FE")

   ### compare with results on page 2275 (in text)
   expect_that(round(coef(res),3), is_equivalent_to(0.122))
   expect_that(round(res$se,3), is_equivalent_to(0.100))
   expect_that(round(res$ci.lb,3), is_equivalent_to(-0.074))
   expect_that(round(res$ci.ub,3), is_equivalent_to(0.317)) ### 0.31 in paper
   expect_that(round(c(logLik(res)), 3), is_equivalent_to(-53.679))

})

test_that("results of the random-effects conditional logistic model are correct.", {

   skip_on_cran()

   res <- rma.glmm(measure="OR", ai=b.xci, n1i=nci, ci=b.xti, n2i=nti, data=dat, model="CM.EL", method="ML")

   ### compare with results on page 2277 (in text)
   expect_that(round(coef(res),3), is_equivalent_to(0.175))
   expect_that(round(res$se,3), is_equivalent_to(0.134))
   expect_that(round(res$ci.lb,3), is_equivalent_to(-0.088))
   expect_that(round(res$ci.ub,3), is_equivalent_to(0.437))
   expect_that(round(c(logLik(res)), 3), is_equivalent_to(-52.989))
   expect_that(round(res$tau2,3), is_equivalent_to(0.119))

})
