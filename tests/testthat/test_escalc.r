context("Testing escalc()")

test_that("escalc() works correctly for measure='RR'", {

   data(dat.bcg, package="metafor")

   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

   expect_that(round(dat$yi[1],4), equals(-0.8893))
   expect_that(round(dat$vi[1],4), equals(0.3256))

})
