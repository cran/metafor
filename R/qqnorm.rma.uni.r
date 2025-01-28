qqnorm.rma.uni <- function(y, type="rstandard", pch=21, col, bg, grid=FALSE,
envelope=TRUE, level=y$level, bonferroni=FALSE, reps=1000, smooth=TRUE, bass=0,
label=FALSE, offset=0.3, pos=13, lty, ...) {

   mstyle <- .get.mstyle()

   .chkclass(class(y), must="rma.uni", notav=c("rma.gen", "rma.uni.selmodel"))

   na.act <- getOption("na.action")
   on.exit(options(na.action=na.act), add=TRUE)

   x <- y

   type <- match.arg(type, c("rstandard", "rstudent"))

   if (x$k == 1L)
      stop(mstyle$stop("Stopped because k = 1."))

   if (length(label) != 1L)
      stop(mstyle$stop("Argument 'label' should be of length 1."))

   .start.plot()

   envelopecol <- .coladj(par("bg","fg"), dark=0.15, light=-0.15)

   if (label == "out" && is.logical(envelope))
      envelope <- TRUE

   if (is.logical(envelope))
      draw.envelope <- envelope

   if (is.character(envelope)) {
      envelopecol <- envelope
      draw.envelope <- TRUE
   }

   if (missing(col))
      col <- par("fg")

   if (missing(bg))
      bg <- .coladj(par("bg","fg"), dark=0.35, light=-0.35)

   if (is.logical(grid))
      gridcol <- .coladj(par("bg","fg"), dark=c(0.2,-0.6), light=c(-0.2,0.6))

   if (is.character(grid)) {
      gridcol <- grid
      grid <- TRUE
   }

   if (missing(lty)) {
      lty <- c("solid", "dotted") # 1st value = diagonal line, 2nd value = pseudo confidence envelope
   } else {
      if (length(lty) == 1L)
         lty <- c(lty, lty)
   }

   ddd <- list(...)

   lqqnorm  <- function(..., seed) qqnorm(...)
   lpoints  <- function(..., seed) points(...)
   labline  <- function(..., seed) abline(...)
   lpolygon <- function(..., seed) polygon(...)
   llines   <- function(..., seed) lines(...)
   lbox     <- function(..., seed) box(...)
   ltext    <- function(..., seed) text(...)

   #########################################################################

   if (type == "rstandard") {
      res    <- rstandard(x)
      not.na <- !is.na(res$z)
      zi     <- res$z[not.na]
      slab   <- res$slab[not.na]
      ord    <- order(zi)
      slab   <- slab[ord]
   } else {
      res    <- rstudent(x)
      not.na <- !is.na(res$z)
      zi     <- res$z[not.na]
      slab   <- res$slab[not.na]
      ord    <- order(zi)
      slab   <- slab[ord]
   }

   sav <- lqqnorm(zi, pch=pch, col=col, bg=bg, bty="l", ...)

   #########################################################################

   ### construct simulation based pseudo confidence envelope

   if (draw.envelope) {

      level <- .level(level)

      if (!is.null(ddd$seed))
         set.seed(ddd$seed)

      dat <- matrix(rnorm(x$k*reps), nrow=x$k, ncol=reps)

      options(na.action="na.omit")
      H <- hatvalues(x, type="matrix")
      options(na.action = na.act)

      ImH <- diag(x$k) - H
      ei  <- ImH %*% dat
      ei  <- apply(ei, 2, sort)
      if (bonferroni) {
         lb <- apply(ei, 1, quantile,   (level/2)/x$k) # consider using rowQuantiles() from matrixStats package
         ub <- apply(ei, 1, quantile, 1-(level/2)/x$k) # consider using rowQuantiles() from matrixStats package
      } else {
         lb <- apply(ei, 1, quantile,   (level/2)) # consider using rowQuantiles() from matrixStats package
         ub <- apply(ei, 1, quantile, 1-(level/2)) # consider using rowQuantiles() from matrixStats package
      }

      temp.lb <- qqnorm(lb, plot.it=FALSE)
      temp.ub <- qqnorm(ub, plot.it=FALSE)

      if (smooth) {
         temp.lb <- supsmu(temp.lb$x, temp.lb$y, bass=bass)
         temp.ub <- supsmu(temp.ub$x, temp.ub$y, bass=bass)
      }

      if (draw.envelope) {

         lpolygon(c(temp.lb$x,rev(temp.ub$x)), c(temp.lb$y,rev(temp.ub$y)), col=envelopecol, border=NA, ...)
         llines(temp.lb$x, temp.lb$y, lty=lty[2], ...)
         llines(temp.ub$x, temp.ub$y, lty=lty[2], ...)

      }

   }

   ### add grid (and redraw box)

   if (.isTRUE(grid)) {
      grid(col=gridcol)
      lbox(..., bty="l")
   }

   ### draw the diagonal line

   labline(a=0, b=1, lty=lty[1], ...)
   #qqline(zi, ...)
   #abline(h=0, lty="dotted", ...)
   #abline(v=0, lty="dotted", ...)

   ### add the points

   lpoints(sav$x, sav$y, pch=pch, col=col, bg=bg, ...)

   #########################################################################

   ### labeling of points

   if ((is.character(label) && label=="none") || .isFALSE(label))
      return(invisible(sav))

   if ((is.character(label) && label=="all") || .isTRUE(label))
      label <- x$k

   if (is.numeric(label)) {

      label <- round(label)

      if (label < 1 | label > x$k)
         stop(mstyle$stop("Out of range value for 'label' argument."))

      pos.x <- sav$x[ord]
      pos.y <- sav$y[ord]

      dev <- abs(pos.x - pos.y)

      for (i in seq_len(x$k)) {

         if (sum(dev > dev[i]) < label) {
            if (pos <= 4)
               ltext(pos.x[i], pos.y[i], slab[i], pos=pos, offset=offset, ...)
            if (pos == 13)
               ltext(pos.x[i], pos.y[i], slab[i], pos=ifelse(pos.x[i]-pos.y[i] >= 0, 1, 3), offset=offset, ...)
            if (pos == 24)
               ltext(pos.x[i], pos.y[i], slab[i], pos=ifelse(pos.x[i]-pos.y[i] <= 0, 2, 4), offset=offset, ...)
               #ltext(pos.x[i], pos.y[i], slab[i], pos=ifelse(pos.x[i] >= 0, 2, 4), offset=offset, ...)
         }

      }

   } else {

      pos.x <- sav$x[ord]
      pos.y <- sav$y[ord]

      for (i in seq_len(x$k)) {

         if (pos.y[i] < temp.lb$y[i] || pos.y[i] > temp.ub$y[i]) {
            if (pos <= 4)
               ltext(pos.x[i], pos.y[i], slab[i], pos=pos, offset=offset, ...)
            if (pos == 13)
               ltext(pos.x[i], pos.y[i], slab[i], pos=ifelse(pos.x[i]-pos.y[i] >= 0, 1, 3), offset=offset, ...)
            if (pos == 24)
               ltext(pos.x[i], pos.y[i], slab[i], pos=ifelse(pos.x[i]-pos.y[i] <= 0, 2, 4), offset=offset, ...)
         }

      }

   }

   #########################################################################

   #if (envelope) {
   #   invisible(list(pts=sav, ci.lb=temp.lb, ci.ub=temp.ub))
   #} else {
   #   invisible(sav)
   #}

   invisible(sav)

}
