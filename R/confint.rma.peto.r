confint.rma.peto <- function(object, parm, level, digits, transf, targs, ...) {

   mstyle <- .get.mstyle()

   .chkclass(class(object), must="rma.peto")

   if (!missing(parm))
      warning(mstyle$warning("Argument 'parm' (currently) ignored."), call.=FALSE)

   x <- object

   if (missing(level))
      level <- x$level

   if (missing(digits)) {
      digits <- .get.digits(xdigits=x$digits, dmiss=TRUE)
   } else {
      digits <- .get.digits(digits=digits, xdigits=x$digits, dmiss=FALSE)
   }

   if (missing(transf))
      transf <- FALSE

   if (missing(targs))
      targs <- NULL

   ddd <- list(...)

   .chkdots(ddd, c("time"))

   if (.isTRUE(ddd$time))
      time.start <- proc.time()

   #########################################################################

   level <- .level(level)
   crit  <- qnorm(level/2, lower.tail=FALSE)

   beta  <- x$beta
   ci.lb <- beta - crit * x$se
   ci.ub <- beta + crit * x$se

   ### if requested, apply transformation function

   if (.isTRUE(transf)) # if transf=TRUE, apply exp transformation to ORs
      transf <- exp

   if (is.function(transf)) {
      if (is.null(targs)) {
         beta  <- sapply(beta, transf)
         ci.lb <- sapply(ci.lb, transf)
         ci.ub <- sapply(ci.ub, transf)
      } else {
         beta  <- sapply(beta, transf, targs)
         ci.lb <- sapply(ci.lb, transf, targs)
         ci.ub <- sapply(ci.ub, transf, targs)
      }
   }

   ### make sure order of intervals is always increasing

   tmp <- .psort(ci.lb, ci.ub)
   ci.lb <- tmp[,1]
   ci.ub <- tmp[,2]

   #########################################################################

   res <- cbind(estimate=beta, ci.lb, ci.ub)
   res <- list(fixed=res)
   rownames(res$fixed) <- ""

   res$digits <- digits

   if (.isTRUE(ddd$time)) {
      time.end <- proc.time()
      .print.time(unname(time.end - time.start)[3])
   }

   class(res) <- "confint.rma"
   return(res)

}
