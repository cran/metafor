plot.profile.rma <- function(x, xlim, ylim, pch=19, xlab, ylab, main, refline=TRUE, cline=FALSE, ...) {

   #########################################################################

   mstyle <- .get.mstyle("crayon" %in% .packages())

   .chkclass(class(x), must="profile.rma")

   .start.plot()

   if (dev.cur() == 1) {
      par(mfrow=n2mfrow(x$comps))
      #on.exit(par(mfrow=c(1,1)), add=TRUE)
   }

   missing.xlim <- missing(xlim)
   missing.ylim <- missing(ylim)
   missing.xlab <- missing(xlab)
   missing.ylab <- missing(ylab)
   missing.main <- missing(main)

   ### filter out some arguments for the plot() function

   lplot   <- function(..., time, LB, startmethod, sub1) plot(...)
   lpoints <- function(..., time, LB, startmethod, sub1, log) points(...) # need 'log' here so profile(res, log="x") doesn't throw a warning

   #########################################################################

   if (x$comps == 1) {

      if (missing.xlim)
         xlim <- x$xlim

      if (missing.ylim)
         ylim <- x$ylim

      if (missing.xlab)
         xlab <- x$xlab

      if (missing.ylab)
         ylab <- paste(ifelse(x$method=="REML", "Restricted ", ""), "Log-Likelihood", sep="")

      if (missing.main)
         main <- x$title

      if (min(x[[1]]) <= x$vc && max(x[[1]]) >= x$vc) {
         pos <- which(x[[1]] >= x$vc)[1]
         x[[1]] <- c(x[[1]][seq_len(pos-1)], x$vc,    x[[1]][pos:length(x[[1]])])
         x[[2]] <- c(x[[2]][seq_len(pos-1)], x$maxll, x[[2]][pos:length(x[[2]])])
      }

      lplot(x[[1]], x[[2]], type="n", xlab=xlab, ylab=ylab, main=main, bty="l", xlim=xlim, ylim=ylim, ...)

      if (refline) {
         abline(v=x$vc, lty="dotted")
         abline(h=x$maxll, lty="dotted")
      }

      if (cline)
         abline(h=x$maxll - qchisq(0.95, df=1)/2, lty="dotted")

      lpoints(x[[1]], x[[2]], type="o", pch=pch, ...)

   } else {

      for (j in seq_len(x$comps)) {

         if (missing.xlim)
            xlim <- x[[j]]$xlim

         if (missing.ylim)
            ylim <- x[[j]]$ylim

         if (missing.xlab) {
            xlab <- x[[j]]$xlab
         } else {
            if (length(xlab) == 1L) {
               xlab <- rep(xlab, x$comps)
            }
         }

         if (missing.ylab) {
            ylab <- paste(ifelse(x[[j]]$method=="REML", "Restricted ", ""), "Log-Likelihood", sep="")
         } else {
            if (length(ylab) == 1L) {
               ylab <- rep(ylab, x$comps)
            }
         }

         if (missing.main) {
            main <- x[[j]]$title
         } else {
            if (length(main) == 1L) {
               main <- rep(main, x$comps)
            }
         }

         lplot(x[[j]], xlim=xlim, ylim=ylim, pch=pch,
               xlab=if (missing.xlab) xlab else xlab[j],
               ylab=if (missing.ylab) ylab else ylab[j],
               main=if (missing.main) main else main[j],
               cline=cline, ...)

      }

   }

}
