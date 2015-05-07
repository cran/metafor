### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

### see also: http://www.metafor-project.org/doku.php/analyses:rothman2008

context("Checking analysis example rothman2008")

############################################################################

### create dataset (Table 15-1)
dat <- data.frame(
age = c("Age <55", "Age 55+"),
ai = c(8,22),
bi = c(98,76),
ci = c(5,16),
di = c(115,69))

test_that("the to.table() function works.", {

   tmp <- to.table(ai=ai, bi=bi, ci=ci, di=di, data=dat, measure="OR", slab=age, rows=c("Tolbutamide", "Placebo"), cols=c("Dead", "Surviving"))

   expected <- structure(c(8, 5, 98, 115, 22, 16, 76, 69), .Dim = c(2L, 2L, 2L), .Dimnames = list(c("Tolbutamide", "Placebo"), c("Dead", "Surviving"), c("Age <55", "Age 55+")))

   ### compare with data in Table 15-1
   expect_that(tmp, is_identical_to(expected))

})

test_that("the stratum-specific and crude risk differences are computed correctly.", {

   ### stratum-specific risk differences
   tmp <- summary(escalc(ai=ai, bi=bi, ci=ci, di=di, data=dat, measure="RD", digits=3, append=FALSE))
   tmp <- round(as.matrix(tmp[1:4]), 4)

   expected <- structure(c(0.0338, 0.0363, 0.001, 0.0036, 0.0315, 0.0598, 1.0738, 0.6064), .Dim = c(2L, 4L), .Dimnames = list(NULL, c("yi", "vi", "sei", "zi")))

   ### compare with data in Table 15-1
   expect_that(tmp, is_identical_to(expected))

   ### crude risk difference
   tmp <- summary(escalc(ai=sum(ai), bi=sum(bi), ci=sum(ci), di=sum(di), data=dat, measure="RD", digits=3, append=FALSE))
   tmp <- round(as.matrix(tmp[1:4]), 4)

   expected <- structure(c(0.0446, 0.0011, 0.0326, 1.3683), .Dim = c(1L, 4L), .Dimnames = list(NULL, c("yi", "vi", "sei", "zi")))

   ### compare with data in Table 15-1
   expect_that(tmp, is_identical_to(expected))

})

test_that("the stratum-specific and crude risk ratios are computed correctly.", {

   ### stratum-specific risk ratios
   tmp <- summary(escalc(ai=ai, bi=bi, ci=ci, di=di, data=dat, measure="RR", digits=2), transf=exp, append=FALSE)
   tmp <- round(as.matrix(tmp[1:4]), 2)

   expected <- structure(c(1.81, 1.19, 0.31, 0.09, 0.55, 0.29, 1.07, 0.6), .Dim = c(2L, 4L), .Dimnames = list(NULL, c("yi", "vi", "sei", "zi")))

   ### compare with data in Table 15-1
   expect_that(tmp, is_identical_to(expected))

   ### crude risk ratio
   tmp <- summary(escalc(ai=sum(ai), bi=sum(bi), ci=sum(ci), di=sum(di), data=dat, measure="RR", digits=2, append=FALSE), transf=exp)
   tmp <- round(as.matrix(tmp[1:4]), 2)

   expected <- structure(c(1.44, 0.07, 0.27, 1.36), .Dim = c(1L, 4L), .Dimnames = list(NULL, c("yi", "vi", "sei", "zi")))

   ### compare with data in Table 15-1
   expect_that(tmp, is_identical_to(expected))

})

test_that("results are correct for Mantel-Haenszel method.", {

   ### Mantel-Haenszel method with risk differences
   res <- rma.mh(ai=ai, bi=bi, ci=ci, di=di, data=dat, measure="RD", digits=3, level=90)

   expect_that(round(coef(res),3), is_equivalent_to( 0.035))
   expect_that(round(res$ci.lb,3), is_equivalent_to(-0.018))
   expect_that(round(res$ci.ub,3), is_equivalent_to( 0.087))
   expect_that(round(res$QE,3), is_equivalent_to(0.002)) ### 0.001 in chapter
   expect_that(round(res$QEp,3), is_equivalent_to(0.967))

   ### Mantel-Haenszel method with risk ratios
   res <- rma.mh(ai=ai, bi=bi, ci=ci, di=di, data=dat, measure="RR", digits=2, level=90)

   expect_that(round(coef(res),2), is_equivalent_to( 0.28))
   expect_that(round(res$ci.lb,2), is_equivalent_to(-0.14))
   expect_that(round(res$ci.ub,2), is_equivalent_to( 0.71))
   expect_that(round(res$QE,3), is_equivalent_to(0.447))
   expect_that(round(res$QEp,3), is_equivalent_to(0.504))

   tmp <- c(round(confint(res, transf=exp)$fixed, 2))
   expect_that(tmp, is_equivalent_to(c(1.33, 0.87, 2.03)))

   ### Mantel-Haenszel method with odds ratios
   res <- rma.mh(ai=ai, bi=bi, ci=ci, di=di, data=dat, measure="OR", correct=FALSE, digits=2, level=90)

   expect_that(round(coef(res),2), is_equivalent_to( 0.34))
   expect_that(round(res$ci.lb,2), is_equivalent_to(-0.17))
   expect_that(round(res$ci.ub,2), is_equivalent_to( 0.85))
   expect_that(round(res$QE,3), is_equivalent_to(0.347))
   expect_that(round(res$QEp,3), is_equivalent_to(0.556))
   expect_that(round(res$CO,3), is_equivalent_to(1.198))
   expect_that(round(res$COp,3), is_equivalent_to(0.274))
   expect_that(round(res$MH,3), is_equivalent_to(1.191))
   expect_that(round(res$MHp,3), is_equivalent_to(0.275))
   expect_that(round(res$TA,3), is_equivalent_to(0.349))
   expect_that(round(res$TAp,3), is_equivalent_to(0.555))

   tmp <- c(round(confint(res, transf=exp)$fixed, 2))
   expect_that(tmp, is_equivalent_to(c(1.40, 0.84, 2.34)))

   ### conditional MLE of the odds ratio
   res <- rma.glmm(ai=ai, bi=bi, ci=ci, di=di, data=dat, measure="OR", model="CM.EL", method="FE")

   expect_that(round(coef(res),3), is_equivalent_to( 0.338))
   expect_that(round(res$ci.lb,3), is_equivalent_to(-0.271))
   expect_that(round(res$ci.ub,3), is_equivalent_to( 0.947))
   expect_that(round(res$QE.Wld,3), is_equivalent_to(0.346))
   expect_that(round(res$QEp.Wld,3), is_equivalent_to(0.556))
   expect_that(round(res$QE.LRT,3), is_equivalent_to(0.350))
   expect_that(round(res$QEp.LRT,3), is_equivalent_to(0.554))

   tmp <- predict(res, transf=exp)
   expect_that(round(tmp$pred,3), is_equivalent_to(1.402))
   expect_that(round(tmp$ci.lb,3), is_equivalent_to(0.763))
   expect_that(round(tmp$ci.ub,3), is_equivalent_to(2.578))

})

############################################################################

### create dataset (Table 15-2)
dat <- data.frame(
age = c("35-44", "45-54", "55-64", "65-74", "75-84"),
x1i = c(32, 104, 206, 186, 102),
t1i = c(52407, 43248, 28612, 12663, 5317) / 10000,
x2i = c(2, 12, 28, 28, 31),
t2i = c(18790, 10673, 5710, 2585, 1462) / 10000)

test_that("the to.table() function works.", {

   tmp <- to.table(x1i=x1i, x2i=x2i, t1i=t1i, t2i=t2i, data=dat, measure="IRR", slab=age, rows=c("Smokers", "Nonsmokers"), cols=c("Deaths", "Years"))

   expected <- structure(c(32, 2, 5.2407, 1.879, 104, 12, 4.3248, 1.0673, 206, 28, 2.8612, 0.571, 186, 28, 1.2663, 0.2585, 102, 31, 0.5317, 0.1462), .Dim = c(2L, 2L, 5L), .Dimnames = list(c("Smokers", "Nonsmokers"), c("Deaths", "Years"), c("35-44", "45-54", "55-64", "65-74", "75-84")))

   ### compare with data in Table 15-2
   expect_that(tmp, is_identical_to(expected))

})

test_that("the stratum-specific and crude rate differences are computed correctly.", {

   ### stratum-specific rate differences
   tmp <- summary(escalc(x1i=x1i, x2i=x2i, t1i=t1i, t2i=t2i, data=dat, measure="IRD", digits=1, append=FALSE))
   tmp <- round(as.matrix(tmp[1:4]), 4)

   expected <- structure(c(5.0417, 12.804, 22.961, 38.5674, -20.2008, 1.7316, 16.0947, 111.0423, 535.0172, 1811.1307, 1.3159, 4.0118, 10.5377, 23.1304, 42.5574, 3.8313, 3.1916, 2.1789, 1.6674, -0.4747), .Dim = c(5L, 4L), .Dimnames = list(NULL, c("yi", "vi", "sei", "zi")))

   ### compare with data in Table 15-2
   expect_that(tmp, is_identical_to(expected))

   ### crude rate difference
   tmp <- summary(escalc(x1i=sum(x1i), x2i=sum(x2i), t1i=sum(t1i), t2i=sum(t2i), data=dat, measure="IRD", digits=1, append=FALSE))
   tmp <- round(as.matrix(tmp[1:4]), 4)

   expected <- structure(c(18.537, 9.6796, 3.1112, 5.9581), .Dim = c(1L, 4L), .Dimnames = list(NULL, c("yi", "vi", "sei", "zi")))

   ### compare with data in Table 15-2
   expect_that(tmp, is_identical_to(expected))

})

test_that("the stratum-specific and crude rate ratios are computed correctly.", {

   ### stratum-specific rate ratios
   tmp <- summary(escalc(x1i=x1i, x2i=x2i, t1i=t1i, t2i=t2i, data=dat, measure="IRR", digits=1, append=FALSE), transf=exp)
   tmp <- round(as.matrix(tmp[1:4]), 2)

   expected <- structure(c(5.74, 2.14, 1.47, 1.36, 0.9, 0.53, 0.09, 0.04, 0.04, 0.04, 0.73, 0.3, 0.2, 0.2, 0.21, 2.4, 2.49, 1.91, 1.5, -0.49), .Dim = c(5L, 4L), .Dimnames = list(NULL, c("yi", "vi", "sei", "zi")))

   ### compare with data in Table 15-2
   expect_that(tmp, is_identical_to(expected))

   ### crude rate ratio
   tmp <- summary(escalc(x1i=sum(x1i), x2i=sum(x2i), t1i=sum(t1i), t2i=sum(t2i), data=dat, measure="IRR", digits=1, append=FALSE), transf=exp)
   tmp <- round(as.matrix(tmp[1:4]), 2)

   expected <- structure(c(1.72, 0.01, 0.11, 5.06), .Dim = c(1L, 4L), .Dimnames = list(NULL, c("yi", "vi", "sei", "zi")))

   ### compare with data in Table 15-2
   expect_that(tmp, is_identical_to(expected))

})

test_that("results are correct for Mantel-Haenszel method.", {

   ### Mantel-Haenszel method with rate differences
   res <- rma.mh(x1i=x1i, x2i=x2i, t1i=t1i, t2i=t2i, data=dat, measure="IRD", digits=2, level=90)

   expect_that(round(coef(res),2), is_equivalent_to(11.44))
   expect_that(round(res$ci.lb,2), is_equivalent_to( 6.35))
   expect_that(round(res$ci.ub,2), is_equivalent_to(16.53))
   expect_that(round(res$QE,3), is_equivalent_to(26.876))
   expect_that(round(res$QEp,3), is_equivalent_to(0))

   ### Mantel-Haenszel method with rate ratios
   res <- rma.mh(x1i=x1i, x2i=x2i, t1i=t1i, t2i=t2i, data=dat, measure="IRR", digits=2, level=90)

   expect_that(round(coef(res),2), is_equivalent_to(0.35))
   expect_that(round(res$ci.lb,2), is_equivalent_to(0.18))
   expect_that(round(res$ci.ub,2), is_equivalent_to(0.53))
   expect_that(round(res$QE,3), is_equivalent_to(10.412))
   expect_that(round(res$QEp,3), is_equivalent_to(0.034))
   expect_that(round(res$MH,3), is_equivalent_to(10.702))
   expect_that(round(res$MHp,3), is_equivalent_to(0.001))

   tmp <- c(round(confint(res, transf=exp)$fixed, 2))
   expect_that(tmp, is_equivalent_to(c(1.42, 1.19, 1.70)))

   ### Mantel-Haenszel test without continuity correction
   res <- rma.mh(x1i=x1i, x2i=x2i, t1i=t1i, t2i=t2i, data=dat, measure="IRR", level=90, correct=FALSE)

   expect_that(round(res$MH,3), is_equivalent_to(11.016))
   expect_that(round(res$MHp,3), is_equivalent_to(0.001))

   ### unconditional MLE of the rate ratio
   res <- rma.glmm(x1i=x1i, x2i=x2i, t1i=t1i, t2i=t2i, data=dat, measure="IRR", digits=2, level=90, model="UM.FS", method="FE")

   expect_that(round(coef(res),3), is_equivalent_to(0.355))
   expect_that(round(res$ci.lb,3), is_equivalent_to(0.178))
   expect_that(round(res$ci.ub,3), is_equivalent_to(0.531))
   expect_that(round(res$QE.Wld,3), is_equivalent_to(10.199))
   expect_that(round(res$QEp.Wld,3), is_equivalent_to(0.037))
   expect_that(round(res$QE.LRT,3), is_equivalent_to(12.132))
   expect_that(round(res$QEp.LRT,3), is_equivalent_to(0.016))

   tmp <- predict(res, transf=exp)
   expect_that(round(tmp$pred,3), is_equivalent_to(1.426))
   expect_that(round(tmp$ci.lb,3), is_equivalent_to(1.195))
   expect_that(round(tmp$ci.ub,3), is_equivalent_to(1.701))

   ### conditional MLE of the rate ratio
   res <- rma.glmm(x1i=x1i, x2i=x2i, t1i=t1i, t2i=t2i, data=dat, measure="IRR", digits=2, level=90, model="CM.EL", method="FE")

   expect_that(round(coef(res),3), is_equivalent_to(0.355))
   expect_that(round(res$ci.lb,3), is_equivalent_to(0.178))
   expect_that(round(res$ci.ub,3), is_equivalent_to(0.531))
   expect_that(round(res$QE.Wld,3), is_equivalent_to(10.199))
   expect_that(round(res$QEp.Wld,3), is_equivalent_to(0.037))
   expect_that(round(res$QE.LRT,3), is_equivalent_to(12.132))
   expect_that(round(res$QEp.LRT,3), is_equivalent_to(0.016))

   tmp <- predict(res, transf=exp)
   expect_that(round(tmp$pred,3), is_equivalent_to(1.426))
   expect_that(round(tmp$ci.lb,3), is_equivalent_to(1.195))
   expect_that(round(tmp$ci.ub,3), is_equivalent_to(1.701))

})

############################################################################

### create dataset (Table 15-5)
dat <- data.frame(
age = c("<35", "35+"),
ai = c(3,1),
bi = c(9,3),
ci = c(104,5),
di = c(1059,86))

test_that("the to.table() function works.", {

   tmp <- to.table(ai=ai, bi=bi, ci=ci, di=di, data=dat, measure="OR", slab=age, rows=c("Down Syndrome", "Control"), cols=c("Spermicide Use", "No Spermicide"))

   expected <- structure(c(3, 104, 9, 1059, 1, 5, 3, 86), .Dim = c(2L, 2L, 2L), .Dimnames = list(c("Down Syndrome", "Control"), c("Spermicide Use", "No Spermicide"), c("<35", "35+")))

   ### compare with data in Table 15-5
   expect_that(tmp, is_identical_to(expected))

})

test_that("results are correct for Mantel-Haenszel method.", {

   ### Mantel-Haenszel method with odds ratios
   res <- rma.mh(ai=ai, bi=bi, ci=ci, di=di, data=dat, measure="OR", digits=2, level=90, correct=FALSE)

   expect_that(round(coef(res),3), is_equivalent_to(1.330))
   expect_that(round(res$ci.lb,3), is_equivalent_to(0.358))
   expect_that(round(res$ci.ub,3), is_equivalent_to(2.302))
   expect_that(round(res$QE,3), is_equivalent_to(0.138))
   expect_that(round(res$QEp,3), is_equivalent_to(0.711))
   expect_that(round(res$CO,3), is_equivalent_to(5.825))
   expect_that(round(res$COp,3), is_equivalent_to(0.016))
   expect_that(round(res$MH,3), is_equivalent_to(5.809))
   expect_that(round(res$MHp,3), is_equivalent_to(0.016))
   expect_that(round(res$TA,3), is_equivalent_to(0.139))
   expect_that(round(res$TAp,3), is_equivalent_to(0.709))

   tmp <- c(round(confint(res, transf=exp)$fixed, 2))
   expect_that(tmp, is_equivalent_to(c(3.78, 1.43, 10.00)))

   ### unconditional MLE of the odds ratio
   res <- rma.glmm(ai=ai, bi=bi, ci=ci, di=di, data=dat, measure="OR", digits=2, level=90, model="UM.FS", method="FE")

   expect_that(round(coef(res),3), is_equivalent_to(1.332))
   expect_that(round(res$ci.lb,3), is_equivalent_to(0.358))
   expect_that(round(res$ci.ub,3), is_equivalent_to(2.305))
   expect_that(round(res$QE.Wld,3), is_equivalent_to(0.137))
   expect_that(round(res$QEp.Wld,3), is_equivalent_to(0.711))
   expect_that(round(res$QE.LRT,3), is_equivalent_to(0.132))
   expect_that(round(res$QEp.LRT,3), is_equivalent_to(0.716))

   tmp <- predict(res, transf=exp)
   expect_that(round(tmp$pred,3), is_equivalent_to(3.788))
   expect_that(round(tmp$ci.lb,3), is_equivalent_to(1.431))
   expect_that(round(tmp$ci.ub,3), is_equivalent_to(10.028))

   ### conditional MLE of the odds ratio
   res <- rma.glmm(ai=ai, bi=bi, ci=ci, di=di, data=dat, measure="OR", digits=2, level=90, model="CM.EL", method="FE")

   expect_that(round(coef(res),3), is_equivalent_to(1.326))
   expect_that(round(res$ci.lb,3), is_equivalent_to(0.356))
   expect_that(round(res$ci.ub,3), is_equivalent_to(2.296))
   expect_that(round(res$QE.Wld,3), is_equivalent_to(0.123))
   expect_that(round(res$QEp.Wld,3), is_equivalent_to(0.726))
   expect_that(round(res$QE.LRT,3), is_equivalent_to(0.119))
   expect_that(round(res$QEp.LRT,3), is_equivalent_to(0.730))

   tmp <- predict(res, transf=exp)
   expect_that(round(tmp$pred,3), is_equivalent_to(3.765))
   expect_that(round(tmp$ci.lb,3), is_equivalent_to(1.427))
   expect_that(round(tmp$ci.ub,3), is_equivalent_to(9.932))

})

############################################################################
