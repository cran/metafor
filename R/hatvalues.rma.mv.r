hatvalues.rma.mv <- function(model, type="diagonal", ...) {

   mstyle <- .get.mstyle()

   .chkclass(class(model), must="rma.mv")

   na.act <- getOption("na.action")

   if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
      stop(mstyle$stop("Unknown 'na.action' specified under options()."))

   type <- match.arg(type, c("diagonal", "matrix"))

   #########################################################################

   x <- model

   if (is.null(x$W)) {
      W     <- chol2inv(chol(x$M))
      stXWX <- chol2inv(chol(as.matrix(t(x$X) %*% W %*% x$X)))
      H     <- as.matrix(x$X %*% stXWX %*% crossprod(x$X,W))
      #H <- as.matrix(x$X %*% x$vb %*% crossprod(x$X,W)) # x$vb may have been changed through robust()
   } else {
      A     <- x$W
      stXAX <- chol2inv(chol(as.matrix(t(x$X) %*% A %*% x$X)))
      H     <- as.matrix(x$X %*% stXAX %*% crossprod(x$X,A))
   }


   #########################################################################

   if (type == "diagonal") {

      hii <- rep(NA_real_, x$k.f)
      hii[x$not.na] <- as.vector(diag(H))
      hii[hii > 1 - 10 * .Machine$double.eps] <- 1 # as in lm.influence()
      names(hii) <- x$slab

      if (na.act == "na.omit")
         hii <- hii[x$not.na]

      if (na.act == "na.fail" && any(!x$not.na))
         stop(mstyle$stop("Missing values in results."))

      return(hii)

   }

   if (type == "matrix") {

      Hfull <- matrix(NA_real_, nrow=x$k.f, ncol=x$k.f)
      Hfull[x$not.na, x$not.na] <- H

      rownames(Hfull) <- x$slab
      colnames(Hfull) <- x$slab

      if (na.act == "na.omit")
         Hfull <- Hfull[x$not.na, x$not.na, drop=FALSE]

      if (na.act == "na.fail" && any(!x$not.na))
         stop(mstyle$stop("Missing values in results."))

      return(Hfull)

   }

}
