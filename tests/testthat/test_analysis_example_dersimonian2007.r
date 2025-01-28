### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true"); Sys.setenv(RUN_VIS_TESTS="true")

### see: https://www.metafor-project.org/doku.php/analyses:dersimonian2007

source("settings.r")

context("Checking analysis example: dersimonian2007")

### data for the CLASP example
n1i <- c(156, 303, 565, 1570, 103, 4659)
n2i <- c( 74, 303, 477, 1565, 105, 4650)
ai  <- c(  5,   5,  12,   69,   9,  313)
ci  <- c(  8,  17,   9,   94,  11,  352)

test_that("results are correct for the CLASP example.", {

   skip_on_cran()

   ### calculate log(OR)s and corresponding sampling variances
   dat <- escalc(measure="OR", ai=ai, n1i=n1i, ci=ci, n2i=n2i)

   ### fit RE model with various tau^2 estimators
   res.PM    <- rma(yi, vi, method="PM", data=dat)
   res.CA    <- rma(yi, vi, method="HE", data=dat)
   res.DL    <- rma(yi, vi, method="DL", data=dat)

   res.CA2   <- rma(yi, vi, method="GENQ", weights=1/(vi + res.CA$tau2), data=dat)
   res.DL2   <- rma(yi, vi, method="GENQ", weights=1/(vi + res.DL$tau2), data=dat)
   res.CA2   <- rma(yi, vi, tau2=res.CA2$tau2, data=dat)
   res.DL2   <- rma(yi, vi, tau2=res.DL2$tau2, data=dat)

   res.EB    <- rma(yi, vi, method="EB",   data=dat)
   res.ML    <- rma(yi, vi, method="ML",   data=dat)
   res.REML  <- rma(yi, vi, method="REML", data=dat)
   res.HS    <- rma(yi, vi, method="HS",   data=dat)
   res.SJ    <- rma(yi, vi, method="SJ",   data=dat)
   res.SJ2   <- rma(yi, vi, method="SJ",   data=dat, control=list(tau2.init=res.CA$tau2))

   ### some extra ones
   res.HSk   <- rma(yi, vi, method="HSk",   data=dat)
   res.GENQM <- rma(yi, vi, method="GENQM", weights=1/vi, data=dat)
   res.PMM   <- rma(yi, vi, method="PMM",   data=dat)

   ### combine results into one long list of fitted models
   res.all <- list(res.PM, res.CA, res.DL, res.CA2, res.DL2, res.EB, res.ML, res.REML, res.HS, res.SJ, res.SJ2, res.HSk, res.GENQM, res.PMM)

   ### create table with estimate of tau, mu, and standard error
   results <- rbind(
   tau = sapply(res.all, function(x) sqrt(x$tau2)),
   mu  = sapply(res.all, coef),
   se  = sapply(res.all, se))
   colnames(results) <- c("PM", "CA", "DL", "CA2", "DL2", "EB", "ML", "REML", "HS", "SJ", "SJ2", "HSk", "GENQM", "PMM")
   tmp <- t(results)

   ### compare with results on page 111-112 (Tables 3 and 4)
   expected <- structure(c( 0.3681,  0.4410,  0.2323,  0.3831,  0.3254,  0.3681,  0.0023,  0.1843,  0.1330,  0.4572,  0.4084,  0.1644,  0.2929,  0.4341,
                           -0.3811, -0.4035, -0.3240, -0.3861, -0.3655, -0.3811, -0.1974, -0.2980, -0.2666, -0.4079, -0.3941, -0.2863, -0.1973, -0.4016,
                            0.2060,  0.2327,  0.1540,  0.2115,  0.1901,  0.2060,  0.0694,  0.1343,  0.1125,  0.2386,  0.2208,  0.1259,  0.2342,  0.2302),
                          .Dim = c(14L, 3L), .Dimnames = list(c("PM", "CA", "DL", "CA2", "DL2", "EB", "ML", "REML", "HS", "SJ", "SJ2", "HSk", "GENQM", "PMM"), c("tau", "mu", "se")))
   expect_equivalent(tmp[,1], expected[,1], tolerance=.tol[["var"]])
   expect_equivalent(tmp[,2], expected[,2], tolerance=.tol[["coef"]])
   expect_equivalent(tmp[,3], expected[,3], tolerance=.tol[["se"]])

})

rm(list=ls())
