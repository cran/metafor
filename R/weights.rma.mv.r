weights.rma.mv <- function(object, type="diagonal", ...) {

   mstyle <- .get.mstyle()

   .chkclass(class(object), must="rma.mv")

   if (is.null(object$not.na))
      stop(mstyle$stop("Information needed to compute the weights is not available in the model object."))

   na.act <- getOption("na.action")

   if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
      stop(mstyle$stop("Unknown 'na.action' specified under options()."))

   type <- match.arg(type, c("diagonal", "matrix", "rowsum"))

   x <- object

   #########################################################################

   if (is.null(x$W)) {
      W <- chol2inv(chol(x$M))
   } else {
      W <- x$W
   }

   #########################################################################

   if (type == "diagonal") {

      wi <- as.vector(diag(W))
      weight <- rep(NA_real_, x$k.f)
      weight[x$not.na] <- wi / sum(wi) * 100
      names(weight) <- x$slab

      if (na.act == "na.omit")
         weight <- weight[x$not.na]

      if (na.act == "na.fail" && any(!x$not.na))
         stop(mstyle$stop("Missing values in weights."))

      return(weight)

   }

   if (type == "matrix") {

      Wfull <- matrix(NA_real_, nrow=x$k.f, ncol=x$k.f)
      Wfull[x$not.na, x$not.na] <- as.matrix(W) # as.matrix() needed when sparse=TRUE

      rownames(Wfull) <- x$slab
      colnames(Wfull) <- x$slab

      if (na.act == "na.omit")
         Wfull <- Wfull[x$not.na, x$not.na, drop=FALSE]

      if (na.act == "na.fail" && any(!x$not.na))
         stop(mstyle$stop("Missing values in results."))

      return(Wfull)

   }

   if (type == "rowsum") {

      if (!x$int.only)
         stop("Row-sum weights are only meaningful for intercept-only models.")

      wi <- rowSums(W)
      weight <- rep(NA_real_, x$k.f)
      weight[x$not.na] <- wi / sum(wi) * 100
      names(weight) <- x$slab

      if (na.act == "na.omit")
         weight <- weight[x$not.na]

      if (na.act == "na.fail" && any(!x$not.na))
         stop(mstyle$stop("Missing values in weights."))

      return(weight)

   }

}
