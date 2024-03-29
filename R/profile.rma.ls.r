profile.rma.ls <- function(fitted, alpha,
   xlim, ylim, steps=20, lltol=1e-03, progbar=TRUE, parallel="no", ncpus=1, cl, plot=TRUE, ...) {

   mstyle <- .get.mstyle()

   .chkclass(class(fitted), must="rma.ls")

   x <- fitted

   if (x$optbeta)
      stop(mstyle$stop("Profiling not yet implemented for models fitted with 'optbeta=TRUE'."))

   if (anyNA(steps))
      stop(mstyle$stop("No missing values allowed in 'steps' argument."))

   if (length(steps) >= 2L) {
      if (missing(xlim))
         xlim <- range(steps)
      stepseq <- TRUE
   } else {
      if (steps < 2)
         stop(mstyle$stop("Argument 'steps' must be >= 2."))
      stepseq <- FALSE
   }

   parallel <- match.arg(parallel, c("no", "snow", "multicore"))

   if (parallel == "no" && ncpus > 1)
      parallel <- "snow"

   if (missing(cl))
      cl <- NULL

   if (!is.null(cl) && inherits(cl, "SOCKcluster")) {
      parallel <- "snow"
      ncpus <- length(cl)
   }

   if (parallel == "snow" && ncpus < 2)
      parallel <- "no"

   if (parallel == "snow" || parallel == "multicore") {

      if (!requireNamespace("parallel", quietly=TRUE))
         stop(mstyle$stop("Please install the 'parallel' package for parallel processing."))

      ncpus <- as.integer(ncpus)

      if (ncpus < 1L)
         stop(mstyle$stop("Argument 'ncpus' must be >= 1."))

   }

   if (!progbar) {
      pbo <- pbapply::pboptions(type="none")
      on.exit(pbapply::pboptions(pbo), add=TRUE)
   }

   ddd <- list(...)

   if (.isTRUE(ddd$time))
      time.start <- proc.time()

   #########################################################################

   ### check if user has not specified alpha argument

   if (missing(alpha)) {

      mc <- match.call()

      ### total number of non-fixed components

      comps <- sum(!x$alpha.fix)

      if (comps == 0)
         stop(mstyle$stop("No components in the model for which a profile likelihood can be constructed."))

      if (plot) {
         if (dev.cur() == 1L) { # if only the 'null device' is currently open, set mfrow
            par(mfrow=n2mfrow(comps))
            #on.exit(par(mfrow=c(1,1)), add=TRUE)
         }
      }

      sav <- list()
      j <- 0

      if (any(!x$alpha.fix)) {
         for (pos in seq_len(x$alphas)[!x$alpha.fix]) {
            j <- j + 1
            mc.vc <- mc
            mc.vc$alpha <- pos
            mc.vc$time <- FALSE
            #mc.vc$fitted <- quote(x)
            if (progbar)
               cat(mstyle$verbose(paste("Profiling alpha =", pos, "\n")))
            sav[[j]] <- eval(mc.vc, envir=parent.frame())
         }
      }

      ### if there is just one component, turn the list of lists into a simple list

      if (comps == 1)
         sav <- sav[[1]]

      sav$comps <- comps

      if (.isTRUE(ddd$time)) {
         time.end <- proc.time()
         .print.time(unname(time.end - time.start)[3])
      }

      class(sav) <- "profile.rma"
      return(invisible(sav))

   }

   ### round and take unique values

   if (!missing(alpha) && is.numeric(alpha))
      alpha <- unique(round(alpha))

   ### check if model actually contains (at least one) such a component and that it was actually estimated

   if (!missing(alpha) && all(x$alpha.fix))
      stop(mstyle$stop("Model does not contain any estimated 'alpha' components."))

   ### check if user specified more than one alpha component

   if (!missing(alpha) && (length(alpha) > 1L))
      stop(mstyle$stop("Can only specify one 'alpha' component."))

   ### check if user specified a logical

   if (!missing(alpha) && is.logical(alpha))
      stop(mstyle$stop("Must specify a number for the 'alpha' component."))

   ### check if user specified a component that does not exist

   if (!missing(alpha) && (alpha > x$alphas || alpha <= 0))
      stop(mstyle$stop("No such 'alpha' component in the model."))

   ### check if user specified a component that was fixed

   if (!missing(alpha) && x$alpha.fix[alpha])
      stop(mstyle$stop("Specified 'alpha' component was fixed."))

   ### if everything is good so far, get value of the component and set 'comp'

   alpha.pos <- NA_integer_

   if (!missing(alpha)) {
      vc <- x$alpha[alpha]
      comp <- "alpha"
      alpha.pos <- alpha
   }

   #return(list(comp=comp, vc=vc))

   #########################################################################

   if (missing(xlim) || is.null(xlim)) {

      ### if the user has not specified xlim, set it automatically

      if (comp == "alpha") {
         if (is.na(x$se.alpha[alpha])) {
            vc.lb <- vc - 4 * abs(vc)
            vc.ub <- vc + 4 * abs(vc)
         } else {
            vc.lb <- vc - qnorm(0.995) * x$se.alpha[alpha]
            vc.ub <- vc + qnorm(0.995) * x$se.alpha[alpha]
         }
      }

      ### if that fails, throw an error

      if (is.na(vc.lb) || is.na(vc.ub))
         stop(mstyle$stop("Cannot set 'xlim' automatically. Please set this argument manually."))

      ### apply alpha.min/alpha.max limits (if they exist) on vc.lb/vc.ub as well

      if (!is.null(x$control$alpha.min)) {
         if (length(x$control$alpha.min) == 1L)
            x$control$alpha.min <- rep(x$control$alpha.min, x$q)
         vc.lb <- max(vc.lb, x$con$alpha.min[alpha])
      }
      if (!is.null(x$control$alpha.max)) {
         if (length(x$control$alpha.max) == 1L)
            x$control$alpha.max <- rep(x$control$alpha.max, x$q)
         vc.ub <- min(vc.ub, x$con$alpha.max[alpha])
      }

      xlim <- sort(c(vc.lb, vc.ub))

   } else {

      if (length(xlim) != 2L)
         stop(mstyle$stop("Argument 'xlim' should be a vector of length 2."))

      xlim <- sort(xlim)

   }

   if (stepseq) {
      vcs <- steps
   } else {
      vcs <- seq(xlim[1], xlim[2], length.out=steps)
   }

   #return(vcs)

   if (length(vcs) <= 1L)
      stop(mstyle$stop("Cannot set 'xlim' automatically. Please set this argument manually."))

   if (parallel == "no")
      res <- pbapply::pblapply(vcs, .profile.rma.ls, obj=x, comp=comp, alpha.pos=alpha.pos, parallel=parallel, profile=TRUE)

   if (parallel == "multicore")
      res <- pbapply::pblapply(vcs, .profile.rma.ls, obj=x, comp=comp, alpha.pos=alpha.pos, parallel=parallel, profile=TRUE, cl=ncpus)
      #res <- parallel::mclapply(vcs, .profile.rma.ls, obj=x, comp=comp, alpha.pos=alpha.pos, parallel=parallel, profile=TRUE, mc.cores=ncpus)

   if (parallel == "snow") {
      if (is.null(cl)) {
         cl <- parallel::makePSOCKcluster(ncpus)
         on.exit(parallel::stopCluster(cl), add=TRUE)
      }
      if (.isTRUE(ddd$LB)) {
         res <- parallel::parLapplyLB(cl, vcs, .profile.rma.ls, obj=x, comp=comp, alpha.pos=alpha.pos, parallel=parallel, profile=TRUE)
         #res <- parallel::clusterApplyLB(cl, vcs, .profile.rma.ls, obj=x, comp=comp, alpha.pos=alpha.pos, parallel=parallel, profile=TRUE)
         #res <- parallel::clusterMap(cl, .profile.rma.ls, vcs, MoreArgs=list(obj=x, comp=comp, alpha.pos=alpha.pos, parallel=parallel, profile=TRUE), .scheduling = "dynamic")
      } else {
         res <- pbapply::pblapply(vcs, .profile.rma.ls, obj=x, comp=comp, alpha.pos=alpha.pos, parallel=parallel, profile=TRUE, cl=cl)
         #res <- parallel::parLapply(cl, vcs, .profile.rma.ls, obj=x, comp=comp, alpha.pos=alpha.pos, parallel=parallel, profile=TRUE)
         #res <- parallel::clusterApply(cl, vcs, .profile.rma.ls, obj=x, comp=comp, alpha.pos=alpha.pos, parallel=parallel, profile=TRUE)
         #res <- parallel::clusterMap(cl, .profile.rma.ls, vcs, MoreArgs=list(obj=x, comp=comp, alpha.pos=alpha.pos, parallel=parallel, profile=TRUE))
      }
   }

   lls   <- sapply(res, function(x) x$ll)
   beta  <- do.call(rbind, lapply(res, function(x) t(x$beta)))
   ci.lb <- do.call(rbind, lapply(res, function(x) t(x$ci.lb)))
   ci.ub <- do.call(rbind, lapply(res, function(x) t(x$ci.ub)))

   beta  <- data.frame(beta)
   ci.lb <- data.frame(ci.lb)
   ci.ub <- data.frame(ci.ub)
   names(beta)  <- rownames(x$beta)
   names(ci.lb) <- rownames(x$beta)
   names(ci.ub) <- rownames(x$beta)

   #########################################################################

   maxll <- c(logLik(x))

   if (any(lls >= maxll + lltol, na.rm=TRUE))
      warning(mstyle$warning("At least one profiled log-likelihood value is larger than the log-likelihood of the fitted model."), call.=FALSE)

   if (all(is.na(lls)))
      warning(mstyle$warning("All model fits failed. Cannot draw profile likelihood plot."), call.=FALSE)

   if (.isTRUE(ddd$exp)) {
      lls <- exp(lls)
      maxll <- exp(maxll)
   }

   if (missing(ylim)) {

      if (any(is.finite(lls))) {
         if (xlim[1] <= vc && xlim[2] >= vc) {
            ylim <- range(c(maxll,lls[is.finite(lls)]), na.rm=TRUE)
         } else {
            ylim <- range(lls[is.finite(lls)], na.rm=TRUE)
         }
      } else {
         ylim <- rep(maxll, 2L)
      }

      if (!.isTRUE(ddd$exp))
         ylim <- ylim + c(-0.1, 0.1)

   } else {

      if (length(ylim) != 2L)
         stop(mstyle$stop("Argument 'ylim' should be a vector of length 2."))

      ylim <- sort(ylim)

   }

   if (comp == "alpha") {
      if (x$alphas == 1L) {
         xlab  <- expression(paste(alpha, " Value"))
         title <- expression(paste("Profile Plot for ", alpha))
      } else {
         if (.isTRUE(ddd$sub1))
            alpha <- alpha - 1
         xlab  <- bquote(alpha[.(alpha)] ~ "Value")
         title <- bquote("Profile Plot for" ~ alpha[.(alpha)])
      }
   }

   sav <- list(alpha=vcs, ll=lls, beta=beta, ci.lb=ci.lb, ci.ub=ci.ub, comps=1, ylim=ylim, method=x$method, vc=vc, maxll=maxll, xlab=xlab, title=title, exp=ddd$exp)
   class(sav) <- "profile.rma"

   #########################################################################

   if (plot)
      plot(sav, ...)

   #########################################################################

   if (.isTRUE(ddd$time)) {
      time.end <- proc.time()
      .print.time(unname(time.end - time.start)[3])
   }

   invisible(sav)

}
