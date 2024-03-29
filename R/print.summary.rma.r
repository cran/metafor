print.summary.rma <- function(x, digits=x$digits, showfit=TRUE, signif.stars=getOption("show.signif.stars"), signif.legend=signif.stars, ...) {

   mstyle <- .get.mstyle()

   .chkclass(class(x), must="summary.rma")

   digits <- .get.digits(digits=digits, xdigits=x$digits, dmiss=FALSE)

   ### strip summary.rma class from object (otherwise get recursion)

   class(x) <- class(x)[-1]

   ### print with showfit=TRUE

   print(x, digits=digits, showfit=showfit, signif.stars=signif.stars, signif.legend=signif.legend, ...)

   invisible()

}
