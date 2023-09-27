### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

### see also: https://www.metafor-project.org/doku.php/analyses:viechtbauer2007b

context("Checking analysis example: viechtbauer2007b")

source("settings.r")

### create dataset for example
dat <- escalc(measure="RR", ai=ai, ci=ci, n1i=n1i, n2i=n2i, data=dat.linde2005)
dat <- dat[c(7:10,13:25), c(13:16,18:19,11,6,7,9)]
dat$dosage <- (dat$dosage * 7) / 1000

test_that("results are correct for the CIs.", {

   sav <- summary(dat, transf=exp)[c(13,17),]

   ### compare with results on page 106
   tmp <- sav$ci.lb
   expect_equivalent(tmp, c(.7397, 1.0039), tolerance=.tol[["ci"]]) ### 1.01 in article
   tmp <- sav$ci.ub
   expect_equivalent(tmp, c(1.2793, 1.5434), tolerance=.tol[["ci"]])

})

test_that("results are correct for the equal-effects model.", {

   res <- rma(yi, vi, data=dat, method="EE")
   sav <- predict(res, transf=exp)
   tmp <- c(sav$pred, sav$ci.lb, sav$ci.ub)

   ### compare with results on page 107
   expect_equivalent(tmp, c(1.3840, 1.2599, 1.5204), tolerance=.tol[["pred"]]) ### 1.39 in article
   expect_equivalent(res$QE, 51.5454, tolerance=.tol[["test"]]) ### 55.54 in article

})

test_that("results are correct for the random-effects model.", {

   res <- rma(yi, vi, data=dat, method="DL")
   sav <- predict(res, transf=exp)

   ### compare with results on page 109
   tmp <- c(sav$pred, sav$ci.lb, sav$ci.ub)
   expect_equivalent(tmp, c(1.5722, 1.3103, 1.8864), tolerance=.tol[["pred"]]) ### 1.90 in article
   tmp <- c(sav$pi.lb, sav$pi.ub)
   expect_equivalent(tmp, c(.8488, 2.9120), tolerance=.tol[["ci"]]) ### .87, 2.83 in article (but this was calculated without taking Var[hat(mu)] into consideration)
   expect_equivalent(res$tau2, .0903, tolerance=.tol[["var"]]) ### .091 in article

})

test_that("results are correct for the mixed-effects model.", {

   dat$dosage <- dat$dosage * dat$duration
   res <- rma(yi, vi, mods = ~ I(dosage-34) * I(baseline-20), data=dat, method="DL")

   ### compare with results on page 112
   expect_equivalent(res$tau2, .0475, tolerance=.tol[["var"]])
   expect_equivalent(res$R2, 47.3778, tolerance=.tol[["r2"]]) ### 48% in article

   sav <- structure(list(estimate = c(0.47625885, -0.0058448, -0.06722782, -0.00156996), se = c(0.08764097, 0.00999872, 0.03522283, 0.00344659), zval = c(5.43420301, -0.58455444, -1.9086436, -0.45551255), pval = c(6e-08, 0.55884735, 0.05630808, 0.64874054)), row.names = c("intrcpt", "I(dosage - 34)", "I(baseline - 20)", "I(dosage - 34):I(baseline - 20)"), class = "data.frame")

   ### compare with results in Table II on page 113
   expect_equivalent(coef(summary(res))[,1:4], sav, tolerance=.tol[["misc"]])

   ### compare with results on page 113
   sav <- predict(res, newmods=c(34-34, 12.5-20, (34-34)*(12.5-20)), transf=exp)
   tmp <- c(sav$pred, sav$ci.lb, sav$ci.ub)
   expect_equivalent(tmp, c(2.6657, 1.4560, 4.8806), tolerance=.tol[["pred"]]) ### 2.66, 1.46, 4.90 in article
   sav <- predict(res, newmods=c(34-34, 23.6-20, (34-34)*(23.6-20)), transf=exp)
   tmp <- c(sav$pred, sav$ci.lb, sav$ci.ub)
   expect_equivalent(tmp, c(1.2639, 0.9923, 1.6099), tolerance=.tol[["pred"]]) ### 1.61 in article

   skip_on_cran()

   png(filename="images/test_analysis_example_viechtbauer2007b_test.png", res=200, width=1800, height=1600, type="cairo")

   par(mar=c(4,4,1,1))

   xvals   <- seq(12, 24, by=0.1) - 20
   modvals <- cbind(0, cbind(xvals, 0))
   preds   <- predict(res, modvals)

   regplot(res, mod=3, pred=preds, xvals=xvals,
           shade=FALSE, bty="l", las=1, digits=1, transf=exp,
           xlim=c(12,24)-20, ylim=c(0.5,4), xaxt="n",
           xlab="Baseline HRSD Score", ylab="Relative Rate")
   axis(side=1, at=seq(12, 24, by=2) - 20, labels=seq(12, 24, by=2))

   dev.off()

   expect_true(.vistest("images/test_analysis_example_viechtbauer2007b_test.png", "images/test_analysis_example_viechtbauer2007b.png"))

   ### check results for all tau^2 estimators

   res.HS   <- rma(yi, vi, mods = ~ I(dosage-34) * I(baseline-20), data=dat, method="HS")
   res.HE   <- rma(yi, vi, mods = ~ I(dosage-34) * I(baseline-20), data=dat, method="HE")
   res.DL   <- rma(yi, vi, mods = ~ I(dosage-34) * I(baseline-20), data=dat, method="DL")
   res.GENQ <- rma(yi, vi, mods = ~ I(dosage-34) * I(baseline-20), data=dat, method="GENQ", weights = n1i + n2i)
   res.SJ   <- rma(yi, vi, mods = ~ I(dosage-34) * I(baseline-20), data=dat, method="SJ")
   res.DLIT <- rma(yi, vi, mods = ~ I(dosage-34) * I(baseline-20), data=dat, method="DLIT", control=list(maxiter=500))
   res.SJIT <- rma(yi, vi, mods = ~ I(dosage-34) * I(baseline-20), data=dat, method="SJIT")
   res.PM   <- rma(yi, vi, mods = ~ I(dosage-34) * I(baseline-20), data=dat, method="PM")
   res.ML   <- rma(yi, vi, mods = ~ I(dosage-34) * I(baseline-20), data=dat, method="ML")
   res.REML <- rma(yi, vi, mods = ~ I(dosage-34) * I(baseline-20), data=dat, method="REML")
   res.EB   <- rma(yi, vi, mods = ~ I(dosage-34) * I(baseline-20), data=dat, method="EB")

   res <- list(res.HS, res.HE, res.DL, res.GENQ, res.SJ, res.DLIT, res.SJIT, res.PM, res.ML, res.REML, res.EB)

   res <- data.frame(method=sapply(res, function(x) x$method),
                     tau2=sapply(res, function(x) x$tau2),
                     se.tau2=sapply(res, function(x) x$se.tau2))

   expect_equivalent(res$tau2,    c(0.0253, 0.0388, 0.0475, 0.06, 0.0912, 0.0299, 0.0633, 0.0633, 0.024, 0.0558, 0.0633), tolerance=.tol[["var"]])
   expect_equivalent(res$se.tau2, c(0.0197, 0.0764, 0.0376, 0.0528, 0.0436, 0.0437, 0.046, 0.046, 0.0222, 0.0409, 0.046), tolerance=.tol[["sevar"]])

})

rm(list=ls())
