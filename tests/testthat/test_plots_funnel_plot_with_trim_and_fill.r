### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

### see also: https://www.metafor-project.org/doku.php/plots:funnel_plot_with_trim_and_fill

source("settings.r")

context("Checking plots example: funnel plot with trim and fill")

test_that("plot can be drawn.", {

   expect_equivalent(TRUE, TRUE) # avoid 'Empty test' message

   skip_on_cran()

   png("images/test_plots_funnel_plot_with_trim_and_fill_test.png", res=200, width=1800, height=1500, type="cairo")

   ### decrease margins so the full space is used
   par(mar=c(5,4,1,2))

   ### fit random-effects model
   res <- rma(yi, vi, data=dat.hackshaw1998, measure="OR")

   ### carry out trim-and-fill analysis
   taf <- trimfill(res)

   ### draw funnel plot with missing studies filled in
   funnel(taf, legend=list(show="cis"))

   dev.off()

   expect_true(.vistest("images/test_plots_funnel_plot_with_trim_and_fill_test.png", "images/test_plots_funnel_plot_with_trim_and_fill.png"))

   out <- capture.output(print(taf))

})

rm(list=ls())
