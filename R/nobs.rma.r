nobs.rma <- function(object, ...) {

   mstyle <- .get.mstyle()

   .chkclass(class(object), must="rma")

   n.obs <- object$k.eff - ifelse(object$method == "REML", 1, 0) * object$p.eff

   return(n.obs)

}
