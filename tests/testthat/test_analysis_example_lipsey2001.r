### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

### see also: https://www.metafor-project.org/doku.php/analyses:lipsey2001

context("Checking analysis example: lipsey2001")

source("settings.r")

### create dataset
dat <- data.frame(
id = c(100, 308, 1596, 2479, 9021, 9028, 161, 172, 537, 7049),
yi = c(-0.33, 0.32, 0.39, 0.31, 0.17, 0.64, -0.33, 0.15, -0.02, 0.00),
vi = c(0.084, 0.035, 0.017, 0.034, 0.072, 0.117, 0.102, 0.093, 0.012, 0.067),
random = c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1),
intensity = c(7, 3, 7, 5, 7, 7, 4, 4, 5, 6))

test_that("results are correct for the equal-effects model.", {

   res <- rma(yi, vi, data=dat, method="EE")

   ### compare with results on page 133 (Exhibit 7.3)
   expect_equivalent(c(as.matrix(coef(summary(res)))), c(0.1549, 0.0609, 2.5450, 0.0109, 0.0356, 0.2742), tolerance=.tol[["misc"]])
   expect_equivalent(res$QE, 14.7640, tolerance=.tol[["test"]])
   expect_equivalent(res$QEp, 0.0976, tolerance=.tol[["pval"]])

})

test_that("results are correct for the random-effects model.", {

   res <- rma(yi, vi, data=dat, method="DL")

   ### compare with results on page 133 (Exhibit 7.3)
   expect_equivalent(c(as.matrix(coef(summary(res)))), c(0.1534, 0.0858, 1.7893, 0.0736, -0.0146, 0.3215), tolerance=.tol[["misc"]])
   expect_equivalent(res$tau2, 0.025955, tolerance=.tol[["var"]])

})

test_that("results are correct for the ANOVA-type analysis.", {

   res <- rma(yi, vi, mods = ~ random, data=dat, method="FE")

   res0 <- rma(yi, vi, data=dat, method="EE", subset=random==0)
   res1 <- rma(yi, vi, data=dat, method="EE", subset=random==1)

   tmp <- predict(res, newmods=c(0,1))
   tmp <- do.call(cbind, unclass(tmp)[1:4])

   ### compare with results on page 138 (Exhibit 7.4)
   expect_equivalent(tmp[1,], c( 0.2984, 0.0813,  0.1390, 0.4578), tolerance=.tol[["pred"]])
   expect_equivalent(tmp[2,], c(-0.0277, 0.0917, -0.2075, 0.1521), tolerance=.tol[["se"]])
   expect_equivalent(res$QM,   7.0739, tolerance=.tol[["test"]]) ### 7.0738 in chapter
   expect_equivalent(res$QMp,  0.0078, tolerance=.tol[["pval"]])
   expect_equivalent(res$QE,   7.6901, tolerance=.tol[["test"]]) ### 7.6902 in chapter
   expect_equivalent(res$QEp,  0.4643, tolerance=.tol[["pval"]])
   expect_equivalent(res0$QE,  6.4382, tolerance=.tol[["test"]]) ### 6.4383 in chapter
   expect_equivalent(res0$QEp, 0.2659, tolerance=.tol[["pval"]])
   expect_equivalent(res1$QE,  1.2519, tolerance=.tol[["test"]])
   expect_equivalent(res1$QEp, 0.7406, tolerance=.tol[["pval"]])

})

test_that("results are correct for the meta-regression analysis (fixed-effects with moderators model).", {

   res <- rma(yi, vi, mods = ~ random + intensity, data=dat, method="FE")

   expected <- structure(list(estimate = c(0.32233263, -0.32978043, -0.00408559), se = c(0.29977632, 0.13041815, 0.04928185), zval = c(1.0752438, -2.52863907, -0.08290246), pval = c(0.28226559, 0.01145057, 0.9339291), ci.lb = c(-0.26521816, -0.58539531, -0.10067623), ci.ub = c(0.90988342, -0.07416555, 0.09250506)), row.names = c("intrcpt", "random", "intensity"), class = "data.frame")

   ### compare with results on page 141 (Exhibit 7.6)

   expect_equivalent(coef(summary(res)), expected, tolerance=.tol[["misc"]])
   expect_equivalent(res$QM,  7.0807, tolerance=.tol[["test"]])
   expect_equivalent(res$QMp, 0.0290, tolerance=.tol[["pval"]])
   expect_equivalent(res$QE,  7.6832, tolerance=.tol[["test"]]) ### 7.6833 in chapter
   expect_equivalent(res$QEp, 0.3614, tolerance=.tol[["pval"]]) ### 0.3613 in chapter

})

test_that("results are correct for the meta-regression analysis (mixed-effects model).", {

   res <- rma(yi, vi, mods = ~ random + intensity, data=dat, method="DL")

   expected <- structure(list(estimate = c(0.33106915, -0.32691858, -0.00682302), se = c(0.31983925, 0.1439395, 0.0528008), zval = c(1.03511109, -2.2712222, -0.12922184), pval = c(0.30061703, 0.02313353, 0.89718211), ci.lb = c(-0.29580425, -0.60903481, -0.11031068), ci.ub = c(0.95794255, -0.04480235, 0.09666464)), row.names = c("intrcpt", "random", "intensity"), class = "data.frame")

   ### compare with results on page 141 (Exhibit 7.7)

   expect_equivalent(coef(summary(res)), expected, tolerance=.tol[["misc"]])
   expect_equivalent(res$QM,   5.5711, tolerance=.tol[["test"]]) ### 5.5709 in chapter
   expect_equivalent(res$QMp,  0.0617, tolerance=.tol[["pval"]])
   expect_equivalent(res$tau2, 0.00488, tolerance=.tol[["var"]])

})

test_that("results are correct for the comutation of R^2 via the anova() function.", {

   res.ME <- rma(yi, vi, mods = ~ random + intensity, data=dat, method="DL")
   res.RE <- rma(yi, vi, data=dat, method="DL")
   expect_warning(tmp <- anova(res.RE, res.ME))

   expect_equivalent(tmp$R2, 81.2023, tolerance=.tol[["r2"]])

})

rm(list=ls())
