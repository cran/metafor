### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true"); Sys.setenv(RUN_VIS_TESTS="true")

context("Checking misc: handling of edge cases due to zeros")

source("settings.r")

test_that("rma.peto(), rma.mh(), and rma.glmm() handle outcome1 never occurring properly.", {

   ai <- c(0,0,0,0)
   bi <- c(10,15,20,25)
   ci <- c(0,0,0,0)
   di <- c(10,10,30,20)

   expect_that(suppressWarnings(rma.peto(ai=ai, bi=bi, ci=ci, di=di)), throws_error())

   expect_warning(res <- rma.mh(measure="OR", ai=ai, bi=bi, ci=ci, di=di))
   expect_true(is.na(res$beta))
   expect_warning(res <- rma.mh(measure="RR", ai=ai, bi=bi, ci=ci, di=di))
   expect_true(is.na(res$beta))
   expect_warning(res <- rma.mh(measure="RD", ai=ai, bi=bi, ci=ci, di=di))
   expect_equivalent(res$beta, 0)

   skip_on_cran()

   expect_error(suppressWarnings(rma.glmm(measure="OR", ai=ai, bi=bi, ci=ci, di=di)))

})

test_that("rma.peto(), rma.mh(), and rma.glmm() handle outcome2 never occurring properly.", {

   ai <- c(10,15,20,25)
   bi <- c(0,0,0,0)
   ci <- c(10,10,30,20)
   di <- c(0,0,0,0)

   expect_error(suppressWarnings(rma.peto(ai=ai, bi=bi, ci=ci, di=di)))

   expect_warning(res <- rma.mh(measure="OR", ai=ai, bi=bi, ci=ci, di=di))
   expect_true(is.na(res$beta))
   expect_warning(res <- rma.mh(measure="RR", ai=ai, bi=bi, ci=ci, di=di))
   expect_equivalent(res$beta, 0)
   expect_warning(res <- rma.mh(measure="RD", ai=ai, bi=bi, ci=ci, di=di))
   expect_equivalent(res$beta, 0)

   skip_on_cran()

   expect_error(suppressWarnings(rma.glmm(measure="OR", ai=ai, bi=bi, ci=ci, di=di)))

})

rm(list=ls())
