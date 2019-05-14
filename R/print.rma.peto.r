print.rma.peto <- function(x, digits, showfit=FALSE, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   if (!inherits(x, "rma.peto"))
      stop(mstyle$stop("Argument 'x' must be an object of class \"rma.peto\"."))

   if (missing(digits)) {
      digits <- .get.digits(xdigits=x$digits, dmiss=TRUE)
   } else {
      digits <- .get.digits(digits=digits, xdigits=x$digits, dmiss=FALSE)
   }

   if (!exists(".rmspace"))
      cat("\n")

   cat(mstyle$section("Fixed-Effects Model"))
   cat(mstyle$section(paste0(" (k = ", x$k, ")")))

   if (showfit) {
      cat("\n")
      if (anyNA(x$fit.stats$ML)) {
         fs <- x$fit.stats$ML
      } else {
         fs <- .fcf(x$fit.stats$ML, digits[["fit"]])
      }
      names(fs) <- c("logLik", "deviance", "AIC", "BIC", "AICc")
      cat("\n")
      tmp <- capture.output(print(fs, quote=FALSE, print.gap=2))
      .print.table(tmp, mstyle)
      cat("\n")
   } else {
      cat("\n\n")
   }

   if (!is.na(x$QE)) {
      cat(mstyle$section("Test for Heterogeneity:"), "\n")
      cat(mstyle$result(paste0("Q(df = ", x$k.pos-1, ") = ", .fcf(x$QE, digits[["test"]]), ", p-val ", .pval(x$QEp, digits=digits[["pval"]], showeq=TRUE, sep=" "))))
   }

   res.table <- c(estimate=.fcf(unname(x$beta), digits[["est"]]), se=.fcf(x$se, digits[["se"]]), zval=.fcf(x$zval, digits[["test"]]), pval=.pval(x$pval, digits[["pval"]]), ci.lb=.fcf(x$ci.lb, digits[["ci"]]), ci.ub=.fcf(x$ci.ub, digits[["ci"]]))
   res.table.exp <- c(estimate=.fcf(exp(unname(x$beta)), digits[["est"]]), ci.lb=.fcf(exp(x$ci.lb), digits[["ci"]]), ci.ub=.fcf(exp(x$ci.ub), digits[["ci"]]))

   cat("\n\n")
   cat(mstyle$section("Model Results (log scale):"))
   cat("\n\n")
   tmp <- capture.output(.print.vector(res.table))
   .print.table(tmp, mstyle)

   cat("\n")
   cat(mstyle$section("Model Results (OR scale):"))
   cat("\n\n")
   tmp <- capture.output(.print.vector(res.table.exp))
   .print.table(tmp, mstyle)

   if (!exists(".rmspace"))
      cat("\n")

   invisible()

}
