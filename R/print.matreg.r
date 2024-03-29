print.matreg <- function(x, digits=x$digits, signif.stars=getOption("show.signif.stars"), signif.legend=signif.stars, ...) {

   mstyle <- .get.mstyle()

   .chkclass(class(x), must="matreg")

   if (missing(digits)) {
      digits <- .get.digits(xdigits=x$digits, dmiss=TRUE)
   } else {
      digits <- .get.digits(digits=digits, xdigits=x$digits, dmiss=FALSE)
   }

   .space()

   if (x$test == "t") {
      res.table <- data.frame(estimate=fmtx(c(x$tab$beta), digits[["est"]]), se=fmtx(x$tab$se, digits[["se"]]), tval=fmtx(x$tab$tval, digits[["test"]]), df=round(x$tab$df,2), pval=fmtp(x$tab$pval, digits[["pval"]]), ci.lb=fmtx(x$tab$ci.lb, digits[["ci"]]), ci.ub=fmtx(x$tab$ci.ub, digits[["ci"]]), stringsAsFactors=FALSE)
   } else {
      res.table <- data.frame(estimate=fmtx(c(x$tab$beta), digits[["est"]]), se=fmtx(x$tab$se, digits[["se"]]), zval=fmtx(x$tab$zval, digits[["test"]]), pval=fmtp(x$tab$pval, digits[["pval"]]), ci.lb=fmtx(x$tab$ci.lb, digits[["ci"]]), ci.ub=fmtx(x$tab$ci.ub, digits[["ci"]]), stringsAsFactors=FALSE)
   }
   rownames(res.table) <- rownames(x$tab)
   signif <- symnum(x$tab$pval, corr=FALSE, na=FALSE, cutpoints=c(0, 0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", ".", " "))
   if (signif.stars) {
      res.table <- cbind(res.table, signif)
      colnames(res.table)[ncol(res.table)] <- ""
   }

   tmp <- capture.output(print(res.table, quote=FALSE, right=TRUE, print.gap=2))
   #tmp[1] <- paste0(tmp[1], "\u200b")
   .print.table(tmp, mstyle)

   if (signif.legend) {
      cat("\n")
      cat(mstyle$legend("---"))
      cat("\n")
      cat(mstyle$legend("Signif. codes: "), mstyle$legend(attr(signif, "legend")))
      cat("\n")
   }

   .space()

   invisible()

}
