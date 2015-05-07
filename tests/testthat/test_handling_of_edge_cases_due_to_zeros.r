### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

context("Checking that edge cases due to zeros are handled properly")

test_that("rma.peto(), rma.mh(), and rma.glmm() handle outcome1 never occurring properly.", {

   ai <- c(0,0,0,0)
   bi <- c(10,15,20,25)
   ci <- c(0,0,0,0)
   di <- c(10,10,30,20)

   expect_that(rma.peto(ai=ai, bi=bi, ci=ci, di=di), throws_error())
   
   res <- rma.mh(measure="OR", ai=ai, bi=bi, ci=ci, di=di)
   expect_that(is.na(res$b), is_true())
   res <- rma.mh(measure="RR", ai=ai, bi=bi, ci=ci, di=di)
   expect_that(is.na(res$b), is_true())
   res <- rma.mh(measure="RD", ai=ai, bi=bi, ci=ci, di=di)
   expect_that(unname(res$b), equals(0))

   expect_that(rma.glmm(measure="OR", ai=ai, bi=bi, ci=ci, di=di), throws_error())

})

test_that("rma.peto(), rma.mh(), and rma.glmm() handle outcome2 never occurring properly.", {

   ai <- c(10,15,20,25)
   bi <- c(0,0,0,0)
   ci <- c(10,10,30,20)
   di <- c(0,0,0,0)

   expect_that(rma.peto(ai=ai, bi=bi, ci=ci, di=di), throws_error())
   
   res <- rma.mh(measure="OR", ai=ai, bi=bi, ci=ci, di=di)
   expect_that(is.na(res$b), is_true())
   res <- rma.mh(measure="RR", ai=ai, bi=bi, ci=ci, di=di)
   expect_that(unname(res$b), equals(0))
   res <- rma.mh(measure="RD", ai=ai, bi=bi, ci=ci, di=di)
   expect_that(unname(res$b), equals(0))

   expect_that(rma.glmm(measure="OR", ai=ai, bi=bi, ci=ci, di=di), throws_error())

})
