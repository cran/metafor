rma <- rma.uni <- function(yi, vi, sei, weights, ai, bi, ci, di, n1i, n2i, x1i, x2i, t1i, t2i, m1i, m2i, sd1i, sd2i, xi, mi, ri, ti, fi, pi, sdi, r2i, ni, mods, scale,
measure="GEN", intercept=TRUE,
data, slab, subset,
add=1/2, to="only0", drop00=FALSE, vtype="LS",
method="REML", weighted=TRUE,
test="z", level=95, btt, att, tau2, verbose=FALSE, digits, control, ...) {

   #########################################################################

   ###### setup

   mstyle <- .get.mstyle()

   ### check argument specifications
   ### (arguments "to" and "vtype" are checked inside escalc function)

   if (!is.element(measure, c("RR","OR","PETO","RD","AS","PHI","ZPHI","YUQ","YUY","RTET","ZTET", # 2x2 table measures
                              "PBIT","OR2D","OR2DN","OR2DL",                                     # 2x2 table transformations to SMDs
                              "MPRD","MPRR","MPOR","MPORC","MPPETO","MPORM",                     # 2x2 table measures for matched pairs / pre-post data
                              "IRR","IRD","IRSD",                                                # two-group person-time data (incidence) measures
                              "MD","SMD","SMDH","SMD1","SMD1H","ROM",                            # two-group mean/SD measures
                              "CVR","VR",                                                        # coefficient of variation ratio, variability ratio
                              "RPB","ZPB","RBIS","ZBIS","D2OR","D2ORN","D2ORL",                  # two-group mean/SD transformations to r_pb, r_bis, and log(OR)
                              "COR","UCOR","ZCOR",                                               # correlations (raw and r-to-z transformed)
                              "PCOR","ZPCOR","SPCOR","ZSPCOR",                                   # partial and semi-partial correlations
                              "R2","ZR2",                                                        # coefficient of determination / R^2 (raw and r-to-z transformed)
                              "PR","PLN","PLO","PAS","PFT",                                      # single proportions (and transformations thereof)
                              "IR","IRLN","IRS","IRFT",                                          # single-group person-time (incidence) data (and transformations thereof)
                              "MN","SMN","MNLN","CVLN","SDLN",                                   # mean, single-group standardized mean, log(mean), log(CV), log(SD),
                              "MC","SMCC","SMCR","SMCRH","SMCRP","SMCRPH","ROMC","CVRC","VRC",   # raw/standardized mean change, log(ROM), CVR, and VR for dependent samples
                              "ARAW","AHW","ABT",                                                # alpha (and transformations thereof)
                              "REH",                                                             # relative excess heterozygosity
                              "HR","HD",                                                         # hazard (rate) ratios and differences
                              "GEN")))
      stop(mstyle$stop("Unknown 'measure' specified."))

   if (!is.element(method[1], c("FE","EE","CE","HS","HSk","HE","DL","DLIT","GENQ","GENQM","SJ","SJIT","PM","MP","PMM","ML","REML","EB")))
      stop(mstyle$stop("Unknown 'method' specified."))

   ### in case user specifies more than one add/to value (as one can do with rma.mh() and rma.peto())
   ### (any kind of continuity correction is directly applied to the outcomes, which are then analyzed as such)

   if (length(add) > 1L)
      add <- add[1]

   if (length(to) > 1L)
      to <- to[1]

   na.act <- getOption("na.action")

   if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
      stop(mstyle$stop("Unknown 'na.action' specified under options()."))

   if (missing(tau2))
      tau2 <- NULL

   if (missing(control))
      control <- list()

   time.start <- proc.time()

   ### get ... argument and check for extra/superfluous arguments

   ddd <- list(...)

   .chkdots(ddd, c("knha", "onlyo1", "addyi", "addvi", "i2def", "r2def", "skipr2", "abbrev", "dfs", "time", "outlist", "link", "optbeta", "alpha", "beta", "skiphes", "retopt", "pleasedonotreportI2thankyouverymuch"))

   ### handle 'knha' argument from ... (note: overrides test argument)

   if (.isFALSE(ddd$knha))
      test <- "z"
   if (.isTRUE(ddd$knha))
      test <- "knha"

    test <- tolower(test)

   if (!is.element(test, c("z", "t", "knha", "hksj", "adhoc")))
      stop(mstyle$stop("Invalid option selected for 'test' argument."))

   if (test == "hksj")
      test <- "knha"

   if (missing(scale)) {
      model <- "rma.uni"
   } else {
      model <- "rma.ls"
   }

   ### set defaults or get onlyo1, addyi, and addvi arguments

   onlyo1 <- .chkddd(ddd$onlyo1, FALSE)
   addyi  <- .chkddd(ddd$addyi,  TRUE)
   addvi  <- .chkddd(ddd$addvi,  TRUE)

   ### set defaults for i2def and r2def

   i2def <- .chkddd(ddd$i2def, "1")
   r2def <- .chkddd(ddd$r2def, "1")

   ### handle arguments for location-scale models

   link <- .chkddd(ddd$link, "log", match.arg(ddd$link, c("log", "identity")))
   optbeta <- .chkddd(ddd$optbeta, FALSE, .isTRUE(ddd$optbeta))

   if (optbeta && !weighted)
      stop(mstyle$stop("Must use 'weighted=TRUE' when 'optbeta=TRUE'."))

   alpha <- .chkddd(ddd$alpha, NA_real_)
   beta  <- .chkddd(ddd$beta,  NA_real_)

   if (model == "rma.uni" && !missing(att))
      warning(mstyle$warning("Argument 'att' only relevant for location-scale models and hence ignored."), call.=FALSE)

   ### set defaults for digits

   if (missing(digits)) {
      digits <- .set.digits(dmiss=TRUE)
   } else {
      digits <- .set.digits(digits, dmiss=FALSE)
   }

   ### set defaults for formulas

   formula.yi    <- NULL
   formula.mods  <- NULL
   formula.scale <- NULL

   ### set options(warn=1) if verbose > 2

   if (verbose > 2) {
      opwarn <- options(warn=1)
      on.exit(options(warn=opwarn$warn), add=TRUE)
   }

   #########################################################################

   if (verbose) .space()

   if (verbose > 1)
      message(mstyle$message("Extracting/computing yi/vi values ..."))

   ### check if data argument has been specified

   if (missing(data))
      data <- NULL

   if (is.null(data)) {
      data <- sys.frame(sys.parent())
   } else {
      if (!is.data.frame(data))
         data <- data.frame(data)
   }

   mf <- match.call()

   ### for certain measures, set add=0 by default unless user explicitly sets the add argument

   addval <- mf[[match("add", names(mf))]]

   if (is.element(measure, c("AS","PHI","ZPHI","RTET","ZTET","IRSD","PAS","PFT","IRS","IRFT")) && is.null(addval))
      add <- 0

   ### extract yi (either NULL if not specified, a vector, a formula, or an escalc object)

   yi <- .getx("yi", mf=mf, data=data)

   ### if yi is not NULL and it is an escalc object, then use that object in place of the data argument

   if (!is.null(yi) && inherits(yi, "escalc"))
      data <- yi

   ### extract weights, slab, subset, mods, and scale values, possibly from the data frame specified via data or yi (arguments not specified are NULL)

   weights <- .getx("weights", mf=mf, data=data, checknumeric=TRUE)
   slab    <- .getx("slab",    mf=mf, data=data)
   subset  <- .getx("subset",  mf=mf, data=data)
   mods    <- .getx("mods",    mf=mf, data=data)
   scale   <- .getx("scale",   mf=mf, data=data)

   ai <- bi <- ci <- di <- x1i <- x2i <- t1i <- t2i <- NA_real_

   if (!is.null(weights) && optbeta)
      stop(mstyle$stop("Cannot use custom weights when 'optbeta=TRUE'."))

   if (!is.null(yi)) {

      ### if yi is not NULL, then yi now either contains the yi values, a formula, or an escalc object

      ### if yi is a formula, extract yi and X (this overrides anything specified via the mods argument further below)

      if (inherits(yi, "formula")) {
         formula.yi <- yi
         formula.mods <- formula.yi[-2]
         options(na.action = "na.pass")                   # set na.action to na.pass, so that NAs are not filtered out (we'll do that later)
         mods <- model.matrix(yi, data=data)              # extract model matrix (now mods is no longer a formula, so [a] further below is skipped)
         attr(mods, "assign") <- NULL                     # strip assign attribute (not needed at the moment)
         attr(mods, "contrasts") <- NULL                  # strip contrasts attribute (not needed at the moment)
         yi <- model.response(model.frame(yi, data=data)) # extract yi values from model frame
         options(na.action = na.act)                      # set na.action back to na.act
         names(yi) <- NULL                                # strip names (1:k) from yi (so res$yi is the same whether yi is a formula or not)
         intercept <- FALSE                               # set to FALSE since formula now controls whether the intercept is included or not
      }                                                   # note: code further below ([b]) actually checks whether intercept is included or not

      ### if yi is an escalc object, try to extract yi and vi (note that moderators must then be specified via the mods argument)

      if (inherits(yi, "escalc")) {

         if (!is.null(attr(yi, "yi.names"))) { # if yi.names attributes is available
            yi.name <- attr(yi, "yi.names")[1] # take the first entry to be the yi variable
         } else {                              # if not, see if 'yi' is in the object and assume that is the yi variable
            if (!is.element("yi", names(yi)))
               stop(mstyle$stop("Cannot determine name of the 'yi' variable."))
            yi.name <- "yi"
         }

         if (!is.null(attr(yi, "vi.names"))) { # if vi.names attributes is available
            vi.name <- attr(yi, "vi.names")[1] # take the first entry to be the vi variable
         } else {                              # if not, see if 'vi' is in the object and assume that is the vi variable
            if (!is.element("vi", names(yi)))
               stop(mstyle$stop("Cannot determine name of the 'vi' variable."))
            vi.name <- "vi"
         }

         ### get vi and yi variables from the escalc object (vi first, then yi, since yi is overwritten)

         vi <- yi[[vi.name]]
         yi <- yi[[yi.name]]

         ### could still be NULL if attributes do not match up with actual contents of the escalc object

         if (is.null(yi))
            stop(mstyle$stop(paste0("Cannot find variable '", yi.name, "' in the object.")))
         if (is.null(vi))
            stop(mstyle$stop(paste0("Cannot find variable '", vi.name, "' in the object.")))

         yi.escalc <- TRUE

      } else {

         yi.escalc <- FALSE

      }

      ### in case user passed a data frame to yi, convert it to a vector (if possible)

      if (is.data.frame(yi)) {
         if (ncol(yi) == 1L) {
            yi <- yi[[1]]
         } else {
            stop(mstyle$stop("The object/variable specified for the 'yi' argument is a data frame with multiple columns."))
         }
      }

      ### in case user passed a matrix to yi, convert it to a vector (if possible)

      if (.is.matrix(yi)) {
         if (nrow(yi) == 1L || ncol(yi) == 1L) {
            yi <- as.vector(yi)
         } else {
            stop(mstyle$stop("The object/variable specified for the 'yi' argument is a matrix with multiple rows/columns."))
         }
      }

      ### check if yi is numeric

      if (!is.numeric(yi))
         stop(mstyle$stop("The object/variable specified for the 'yi' argument is not numeric."))

      ### number of outcomes before subsetting

      k <- length(yi)
      k.all <- k

      ### if the user has specified 'measure' to be something other than "GEN", then use that for the measure argument
      ### otherwise, if yi has a 'measure' attribute, use that to set the 'measure' argument

      if (measure == "GEN" && !is.null(attr(yi, "measure")))
         measure <- attr(yi, "measure")

      ### add measure attribute (back) to the yi vector

      attr(yi, "measure") <- measure

      ### extract vi and sei values (but only if yi wasn't an escalc object)

      if (!yi.escalc) {
         vi  <- .getx("vi",  mf=mf, data=data, checknumeric=TRUE)
         sei <- .getx("sei", mf=mf, data=data, checknumeric=TRUE)
      }

      ### extract ni argument

      ni <- .getx("ni", mf=mf, data=data, checknumeric=TRUE)

      ### if neither vi nor sei is specified, then throw an error
      ### if only sei is specified, then square those values to get vi
      ### if vi is specified, use those values

      if (is.null(vi)) {
         if (is.null(sei)) {
            stop(mstyle$stop("Must specify 'vi' or 'sei' argument."))
         } else {
            vi <- sei^2
         }
      }

      ### check 'vi' argument for potential misuse

      .chkviarg(mf$vi)

      ### in case user passes a matrix to vi, convert it to a vector
      ### note: only a row or column matrix with the right dimensions will have the right length

      if (.is.matrix(vi)) {
         if (nrow(vi) == 1L || ncol(vi) == 1L) {
            vi <- as.vector(vi)
         } else {
            if (.is.square(vi) && isSymmetric(unname(vi))) {
               vi <- as.matrix(vi) # in case vi is sparse
               if (any(vi[!diag(nrow(vi))] != 0))
                  warning(mstyle$warning("Using only the diagonal elements from 'vi' argument as the sampling variances."), call.=FALSE)
               vi <- diag(vi)
            } else {
               stop(mstyle$stop("The object/variable specified for the 'vi' argument is a matrix with multiple rows/columns."))
            }
         }
      }

      ### check if user constrained vi to 0

      if ((length(vi) == 1L && vi == 0) || (length(vi) == k && !anyNA(vi) && all(vi == 0))) {
         vi0 <- TRUE
      } else {
         vi0 <- FALSE
      }

      ### allow easy setting of vi to a single value

      if (length(vi) == 1L)
         vi <- rep(vi, k) # note: k is number of outcomes before subsetting

      ### check length of yi and vi

      if (length(vi) != k)
         stop(mstyle$stop("Length of 'yi' and 'vi' (or 'sei') is not the same."))

      ### if ni has not been specified, try to get it from the attributes of yi

      if (is.null(ni))
         ni <- attr(yi, "ni")

      ### check length of yi and ni (only if ni is not NULL)
      ### if there is a mismatch, then ni cannot be trusted, so set it to NULL

      if (!is.null(ni) && length(ni) != k)
         ni <- NULL

      ### if ni is now available, add it (back) as an attribute to yi

      if (!is.null(ni))
         attr(yi, "ni") <- ni

      ### note: one or more yi/vi pairs may be NA/NA (also a corresponding ni value may be NA)

      ### if slab has not been specified but is an attribute of yi, get it

      if (is.null(slab)) {

         slab <- attr(yi, "slab") # will be NULL if there is no slab attribute

         ### check length of yi and slab (only if slab is now not NULL)
         ### if there is a mismatch, then slab cannot be trusted, so set it to NULL

         if (!is.null(slab) && length(slab) != k)
            slab <- NULL

      }

      ### subsetting of yi/vi/ni values (note: mods and slab are subsetted further below)

      if (!is.null(subset)) {

         subset <- .chksubset(subset, k)
         yi <- .getsubset(yi, subset)
         vi <- .getsubset(vi, subset)
         ni <- .getsubset(ni, subset)

         attr(yi, "measure") <- measure # add measure attribute back
         attr(yi, "ni")      <- ni      # add ni attribute back

      }

   } else {

      ### if yi is NULL, try to compute yi/vi based on specified measure and supplied data

      if (is.element(measure, c("RR","OR","PETO","RD","AS","PHI","ZPHI","YUQ","YUY","RTET","ZTET","PBIT","OR2D","OR2DN","OR2DL","MPRD","MPRR","MPOR","MPORC","MPPETO","MPORM"))) {

         ai  <- .getx("ai",  mf=mf, data=data, checknumeric=TRUE)
         bi  <- .getx("bi",  mf=mf, data=data, checknumeric=TRUE)
         ci  <- .getx("ci",  mf=mf, data=data, checknumeric=TRUE)
         di  <- .getx("di",  mf=mf, data=data, checknumeric=TRUE)
         n1i <- .getx("n1i", mf=mf, data=data, checknumeric=TRUE)
         n2i <- .getx("n2i", mf=mf, data=data, checknumeric=TRUE)
         ri  <- .getx("ri",  mf=mf, data=data, checknumeric=TRUE)
         pi  <- .getx("pi",  mf=mf, data=data, checknumeric=TRUE)

         if (is.null(bi)) bi <- n1i - ai
         if (is.null(di)) di <- n2i - ci

         k <- length(ai) # number of outcomes before subsetting
         k.all <- k

         if (!is.null(subset)) {
            subset <- .chksubset(subset, k)
            ai <- .getsubset(ai, subset)
            bi <- .getsubset(bi, subset)
            ci <- .getsubset(ci, subset)
            di <- .getsubset(di, subset)
            ri <- .getsubset(ri, subset)
            pi <- .getsubset(pi, subset)
         }

         args <- list(measure=measure, ai=ai, bi=bi, ci=ci, di=di, ri=ri, pi=pi, add=add, to=to, drop00=drop00, vtype=vtype, onlyo1=onlyo1, addyi=addyi, addvi=addvi)

      }

      if (is.element(measure, c("IRR","IRD","IRSD"))) {

         x1i <- .getx("x1i", mf=mf, data=data, checknumeric=TRUE)
         x2i <- .getx("x2i", mf=mf, data=data, checknumeric=TRUE)
         t1i <- .getx("t1i", mf=mf, data=data, checknumeric=TRUE)
         t2i <- .getx("t2i", mf=mf, data=data, checknumeric=TRUE)

         k <- length(x1i) # number of outcomes before subsetting
         k.all <- k

         if (!is.null(subset)) {
            subset <- .chksubset(subset, k)
            x1i <- .getsubset(x1i, subset)
            x2i <- .getsubset(x2i, subset)
            t1i <- .getsubset(t1i, subset)
            t2i <- .getsubset(t2i, subset)
         }

         args <- list(measure=measure, x1i=x1i, x2i=x2i, t1i=t1i, t2i=t2i, add=add, to=to, drop00=drop00, vtype=vtype, addyi=addyi, addvi=addvi)

      }

      if (is.element(measure, c("MD","SMD","SMDH","SMD1","SMD1H","ROM","RPB","ZPB","RBIS","ZBIS","D2OR","D2ORN","D2ORL","CVR","VR"))) {

         m1i  <- .getx("m1i",  mf=mf, data=data, checknumeric=TRUE)
         m2i  <- .getx("m2i",  mf=mf, data=data, checknumeric=TRUE)
         sd1i <- .getx("sd1i", mf=mf, data=data, checknumeric=TRUE)
         sd2i <- .getx("sd2i", mf=mf, data=data, checknumeric=TRUE)
         n1i  <- .getx("n1i",  mf=mf, data=data, checknumeric=TRUE)
         n2i  <- .getx("n2i",  mf=mf, data=data, checknumeric=TRUE)
         di   <- .getx("di",   mf=mf, data=data, checknumeric=TRUE)
         ti   <- .getx("ti",   mf=mf, data=data, checknumeric=TRUE)
         pi   <- .getx("pi",   mf=mf, data=data, checknumeric=TRUE)

         if (is.element(measure, c("SMD","RPB","ZPB","RBIS","ZBIS","D2OR","D2ORN","D2ORL"))) {

            if (!.equal.length(m1i, m2i, sd1i, sd2i, n1i, n2i, di, ti, pi))
               stop(mstyle$stop("Supplied data vectors are not all of the same length."))

            ti <- replmiss(ti, .convp2t(pi, df=n1i+n2i-2))
            di <- replmiss(di, ti * sqrt(1/n1i + 1/n2i))

            m1i[!is.na(di)]  <- di[!is.na(di)]
            m2i[!is.na(di)]  <- 0
            sd1i[!is.na(di)] <- 1
            sd2i[!is.na(di)] <- 1

         }

         k <- length(n1i) # number of outcomes before subsetting
         k.all <- k

         if (!is.null(subset)) {
            subset <- .chksubset(subset, k)
            m1i  <- .getsubset(m1i,  subset)
            m2i  <- .getsubset(m2i,  subset)
            sd1i <- .getsubset(sd1i, subset)
            sd2i <- .getsubset(sd2i, subset)
            n1i  <- .getsubset(n1i,  subset)
            n2i  <- .getsubset(n2i,  subset)
         }

         args <- list(measure=measure, m1i=m1i, m2i=m2i, sd1i=sd1i, sd2i=sd2i, n1i=n1i, n2i=n2i, vtype=vtype)

      }

      if (is.element(measure, c("COR","UCOR","ZCOR"))) {

         ri <- .getx("ri", mf=mf, data=data, checknumeric=TRUE)
         ni <- .getx("ni", mf=mf, data=data, checknumeric=TRUE)
         ti <- .getx("ti", mf=mf, data=data, checknumeric=TRUE)
         pi <- .getx("pi", mf=mf, data=data, checknumeric=TRUE)

         if (!.equal.length(ri, ni, ti, pi))
            stop(mstyle$stop("Supplied data vectors are not all of the same length."))

         ti <- replmiss(ti, .convp2t(pi, df=ni-2))
         ri <- replmiss(ri, ti / sqrt(ti^2 + ni-2))

         k <- length(ri) # number of outcomes before subsetting
         k.all <- k

         if (!is.null(subset)) {
            subset <- .chksubset(subset, k)
            ri <- .getsubset(ri, subset)
            ni <- .getsubset(ni, subset)
         }

         args <- list(measure=measure, ri=ri, ni=ni, vtype=vtype)

      }

      if (is.element(measure, c("PCOR","ZPCOR","SPCOR","ZSPCOR"))) {

         ri  <- .getx("ri",  mf=mf, data=data, checknumeric=TRUE)
         ti  <- .getx("ti",  mf=mf, data=data, checknumeric=TRUE)
         mi  <- .getx("mi",  mf=mf, data=data, checknumeric=TRUE)
         ni  <- .getx("ni",  mf=mf, data=data, checknumeric=TRUE)
         pi  <- .getx("pi",  mf=mf, data=data, checknumeric=TRUE)
         r2i <- .getx("r2i", mf=mf, data=data, checknumeric=TRUE)

         if (!.equal.length(ri, ti, mi, ni, pi, r2i))
            stop(mstyle$stop("Supplied data vectors are not all of the same length."))

         ti <- replmiss(ti, .convp2t(pi, df=ni-mi-1))

         if (is.element(measure, c("PCOR","ZPCOR")))
            ri <- replmiss(ri, ti / sqrt(ti^2 + ni-mi-1))

         if (is.element(measure, c("SPCOR","ZSPCOR")))
            ri <- replmiss(ri, ti * sqrt(1-r2i) / sqrt(ni-mi-1))

         k <- length(ri) # number of outcomes before subsetting
         k.all <- k

         if (!is.null(subset)) {
            subset <- .chksubset(subset, k)
            ri  <- .getsubset(ri,  subset)
            mi  <- .getsubset(mi,  subset)
            ni  <- .getsubset(ni,  subset)
            r2i <- .getsubset(r2i, subset)
         }

         args <- list(measure=measure, ri=ri, mi=mi, ni=ni, r2i=r2i, vtype=vtype)

      }

      if (is.element(measure, c("R2","ZR2"))) {

         r2i <- .getx("r2i", mf=mf, data=data, checknumeric=TRUE)
         mi  <- .getx("mi",  mf=mf, data=data, checknumeric=TRUE)
         ni  <- .getx("ni",  mf=mf, data=data, checknumeric=TRUE)
         fi  <- .getx("fi",  mf=mf, data=data, checknumeric=TRUE)
         pi  <- .getx("pi",  mf=mf, data=data, checknumeric=TRUE)

         if (!.equal.length(r2i, mi, ni))
            stop(mstyle$stop("Supplied data vectors are not all of the same length."))

         fi  <- replmiss(fi, .convp2f(pi, df1=mi, df2=ni-mi-1))
         r2i <- replmiss(r2i, mi*fi / (mi*fi + (ni-mi-1)))

         k <- length(r2i) # number of outcomes before subsetting
         k.all <- k

         if (!is.null(subset)) {
            subset <- .chksubset(subset, k)
            r2i <- .getsubset(r2i,  subset)
            mi  <- .getsubset(mi,  subset)
            ni  <- .getsubset(ni,  subset)
         }

         args <- list(measure=measure, r2i=r2i, mi=mi, ni=ni, vtype=vtype)

      }

      if (is.element(measure, c("PR","PLN","PLO","PAS","PFT"))) {

         xi <- .getx("xi", mf=mf, data=data, checknumeric=TRUE)
         mi <- .getx("mi", mf=mf, data=data, checknumeric=TRUE)
         ni <- .getx("ni", mf=mf, data=data, checknumeric=TRUE)

         if (is.null(mi)) mi <- ni - xi

         k <- length(xi) # number of outcomes before subsetting
         k.all <- k

         if (!is.null(subset)) {
            subset <- .chksubset(subset, k)
            xi <- .getsubset(xi, subset)
            mi <- .getsubset(mi, subset)
         }

         args <- list(measure=measure, xi=xi, mi=mi, add=add, to=to, vtype=vtype, addyi=addyi, addvi=addvi)

      }

      if (is.element(measure, c("IR","IRLN","IRS","IRFT"))) {

         xi <- .getx("xi", mf=mf, data=data, checknumeric=TRUE)
         ti <- .getx("ti", mf=mf, data=data, checknumeric=TRUE)

         k <- length(xi) # number of outcomes before subsetting
         k.all <- k

         if (!is.null(subset)) {
            subset <- .chksubset(subset, k)
            xi <- .getsubset(xi, subset)
            ti <- .getsubset(ti, subset)
         }

         args <- list(measure=measure, xi=xi, ti=ti, add=add, to=to, vtype=vtype, addyi=addyi, addvi=addvi)

      }

      if (is.element(measure, c("MN","SMN","MNLN","CVLN","SDLN"))) {

         mi  <- .getx("mi",  mf=mf, data=data, checknumeric=TRUE)
         sdi <- .getx("sdi", mf=mf, data=data, checknumeric=TRUE)
         ni  <- .getx("ni",  mf=mf, data=data, checknumeric=TRUE)

         k <- length(ni) # number of outcomes before subsetting
         k.all <- k

         if (!is.null(subset)) {
            subset <- .chksubset(subset, k)
            mi  <- .getsubset(mi,  subset)
            sdi <- .getsubset(sdi, subset)
            ni  <- .getsubset(ni,  subset)
         }

         args <- list(measure=measure, mi=mi, sdi=sdi, ni=ni, vtype=vtype)

      }

      if (is.element(measure, c("MC","SMCC","SMCR","SMCRH","SMCRP","SMCRPH","ROMC","CVRC","VRC"))) {

         m1i  <- .getx("m1i",  mf=mf, data=data, checknumeric=TRUE)
         m2i  <- .getx("m2i",  mf=mf, data=data, checknumeric=TRUE)
         sd1i <- .getx("sd1i", mf=mf, data=data, checknumeric=TRUE)
         sd2i <- .getx("sd2i", mf=mf, data=data, checknumeric=TRUE)
         ri   <- .getx("ri",   mf=mf, data=data, checknumeric=TRUE)
         ni   <- .getx("ni",   mf=mf, data=data, checknumeric=TRUE)
         di   <- .getx("di",   mf=mf, data=data, checknumeric=TRUE)
         ti   <- .getx("ti",   mf=mf, data=data, checknumeric=TRUE)
         pi   <- .getx("pi",   mf=mf, data=data, checknumeric=TRUE)

         if (measure == "SMCC") {

            if (!.equal.length(m1i, m2i, sd1i, sd2i, ri, ni, di, ti, pi))
               stop(mstyle$stop("Supplied data vectors are not all of the same length."))

            ti <- replmiss(ti, .convp2t(pi, df=ni-1))
            di <- replmiss(di, ti * sqrt(1/ni))

            m1i[!is.na(di)]  <- di[!is.na(di)]
            m2i[!is.na(di)]  <- 0
            sd1i[!is.na(di)] <- 1
            sd2i[!is.na(di)] <- 1
            ri[!is.na(di)]   <- 0.5

         }

         k <- length(m1i) # number of outcomes before subsetting
         k.all <- k

         if (!is.null(subset)) {
            subset <- .chksubset(subset, k)
            m1i  <- .getsubset(m1i,  subset)
            m2i  <- .getsubset(m2i,  subset)
            sd1i <- .getsubset(sd1i, subset)
            sd2i <- .getsubset(sd2i, subset)
            ni   <- .getsubset(ni,   subset)
            ri   <- .getsubset(ri,   subset)
         }

         args <- list(measure=measure, m1i=m1i, m2i=m2i, sd1i=sd1i, sd2i=sd2i, ri=ri, ni=ni, vtype=vtype)

      }

      if (is.element(measure, c("ARAW","AHW","ABT"))) {

         ai <- .getx("ai", mf=mf, data=data, checknumeric=TRUE)
         mi <- .getx("mi", mf=mf, data=data, checknumeric=TRUE)
         ni <- .getx("ni", mf=mf, data=data, checknumeric=TRUE)

         k <- length(ai) # number of outcomes before subsetting
         k.all <- k

         if (!is.null(subset)) {
            subset <- .chksubset(subset, k)
            ai <- .getsubset(ai, subset)
            mi <- .getsubset(mi, subset)
            ni <- .getsubset(ni, subset)
         }

         args <- list(measure=measure, ai=ai, mi=mi, ni=ni, vtype=vtype)

      }

      if (measure == "REH") {

         ai <- .getx("ai", mf=mf, data=data, checknumeric=TRUE)
         bi <- .getx("bi", mf=mf, data=data, checknumeric=TRUE)
         ci <- .getx("ci", mf=mf, data=data, checknumeric=TRUE)

         k <- length(ai) # number of outcomes before subsetting
         k.all <- k

         if (!is.null(subset)) {
            subset <- .chksubset(subset, k)
            ai <- .getsubset(ai, subset)
            bi <- .getsubset(bi, subset)
            ci <- .getsubset(ci, subset)
         }

         args <- list(measure=measure, ai=ai, bi=bi, ci=ci, vtype=vtype)

      }

      dat <- .do.call(escalc, args)

      if (is.element(measure, "GEN"))
         stop(mstyle$stop("Specify the desired outcome measure via the 'measure' argument."))

      ### note: these values are already subsetted

      yi <- dat$yi         # one or more yi/vi pairs may be NA/NA
      vi <- dat$vi         # one or more yi/vi pairs may be NA/NA
      ni <- attr(yi, "ni") # unadjusted total sample sizes (ni.u in escalc)

   }

   #########################################################################

   ### allow easy setting of weights to a single value

   if (length(weights) == 1L)
      weights <- rep(weights, k) # note: k is number of outcomes before subsetting

   ### check length of yi and weights (only if weights is not NULL)

   if (!is.null(weights) && (length(weights) != k))
      stop(mstyle$stop("Length of 'yi' and 'weights' is not the same."))

   ### subsetting of weights

   if (!is.null(subset))
      weights <- .getsubset(weights, subset)

   #########################################################################

   if (verbose > 1)
      message(mstyle$message("Creating model matrix ..."))

   ### convert mods formula to X matrix and set intercept equal to FALSE
   ### skipped if formula has already been specified via yi argument, since mods is then no longer a formula (see [a])

   if (inherits(mods, "formula")) {
      formula.mods <- mods
      if (isTRUE(all.equal(formula.mods, ~ 1))) { # needed so 'mods = ~ 1' without 'data' specified works
         mods <- matrix(1, nrow=k, ncol=1)
         intercept <- FALSE
      } else {
         options(na.action = "na.pass")        # set na.action to na.pass, so that NAs are not filtered out (we'll do that later)
         mods <- model.matrix(mods, data=data) # extract model matrix
         attr(mods, "assign") <- NULL          # strip assign attribute (not needed at the moment)
         attr(mods, "contrasts") <- NULL       # strip contrasts attribute (not needed at the moment)
         options(na.action = na.act)           # set na.action back to na.act
         intercept <- FALSE                    # set to FALSE since formula now controls whether the intercept is included or not
      }                                        # note: code further below ([b]) actually checks whether intercept is included or not
   }

   ### turn a vector for mods into a column vector

   if (.is.vector(mods))
      mods <- cbind(mods)

   ### turn a mods data frame into a matrix

   if (is.data.frame(mods))
      mods <- as.matrix(mods)

   ### check if model matrix contains character variables

   if (is.character(mods))
      stop(mstyle$stop("Model matrix contains character variables."))

   ### check if mods matrix has the right number of rows

   if (!is.null(mods) && nrow(mods) != k)
      stop(mstyle$stop(paste0("Number of rows in the model matrix (", nrow(mods), ") does not match length of the outcome vector (", k, ").")))

   ### for rma.ls models, get model matrix for the scale part

   if (model == "rma.ls") {
      if (inherits(scale, "formula")) {
         formula.scale <- scale
         if (isTRUE(all.equal(formula.scale, ~ 1))) { # needed so 'scale = ~ 1' without 'data' specified works
            Z <- matrix(1, nrow=k, ncol=1)
            colnames(Z) <- "intrcpt"
         } else {
            options(na.action = "na.pass")
            Z <- model.matrix(scale, data=data)
            colnames(Z)[grep("(Intercept)", colnames(Z), fixed=TRUE)] <- "intrcpt"
            attr(Z, "assign") <- NULL
            attr(Z, "contrasts") <- NULL
            options(na.action = na.act)
         }
      } else {
         Z <- scale
         if (.is.vector(Z))
            Z <- cbind(Z)
         if (is.data.frame(Z))
            Z <- as.matrix(Z)
         if (is.character(Z))
            stop(mstyle$stop("Scale model matrix contains character variables."))
      }
      if (nrow(Z) != k)
         stop(mstyle$stop(paste0("Number of rows in the model matrix specified via the 'scale' argument (", nrow(Z), ") does not match length of the outcome vector (", k, ").")))
   } else {
      Z <- NULL
   }

   ### generate study labels if none are specified (or none have been found in yi)

   if (verbose > 1)
      message(mstyle$message("Generating/extracting study labels ..."))

   ### study ids (1:k sequence before subsetting)

   ids <- seq_len(k)

   if (is.null(slab)) {

      slab.null <- TRUE
      slab      <- ids

   } else {

      if (anyNA(slab))
         stop(mstyle$stop("NAs in study labels."))

      if (length(slab) != k)
         stop(mstyle$stop("Study labels not of same length as data."))

      slab.null <- FALSE

   }

   ### if a subset of studies is specified

   if (!is.null(subset)) {

      if (verbose > 1)
         message(mstyle$message("Subsetting ..."))

      mods <- .getsubset(mods, subset)
      slab <- .getsubset(slab, subset)
      ids  <- .getsubset(ids,  subset)
      Z    <- .getsubset(Z,    subset)

   }

   ### check if study labels are unique; if not, make them unique

   if (anyDuplicated(slab))
      slab <- .make.unique(slab)

   ### add slab attribute back

   attr(yi, "slab") <- slab

   ### number of outcomes after subsetting

   k <- length(yi)

   ### check for negative/infinite weights

   if (any(weights < 0, na.rm=TRUE))
      stop(mstyle$stop("Negative weights not allowed."))

   if (any(is.infinite(weights)))
      stop(mstyle$stop("Infinite weights not allowed."))

   ### save full data (including potential NAs in yi/vi/weights/ni/mods/Z.f)

   outdat.f <- list(ai=ai, bi=bi, ci=ci, di=di, x1i=x1i, x2i=x2i, t1i=t1i, t2i=t2i)

   yi.f      <- yi
   vi.f      <- vi
   weights.f <- weights
   ni.f      <- ni
   mods.f    <- mods
   Z.f       <- Z

   k.f <- k # total number of observed outcomes including all NAs

   ### check for NAs and act accordingly

   has.na <- is.na(yi) | is.na(vi) | (if (is.null(mods)) FALSE else apply(is.na(mods), 1, any)) | (if (is.null(Z)) FALSE else apply(is.na(Z), 1, any)) | (if (is.null(weights)) FALSE else is.na(weights))
   not.na <- !has.na

   if (any(has.na)) {

      if (verbose > 1)
         message(mstyle$message("Handling NAs ..."))

      if (na.act == "na.omit" || na.act == "na.exclude" || na.act == "na.pass") {

         yi      <- yi[not.na]
         vi      <- vi[not.na]
         weights <- weights[not.na]
         ni      <- ni[not.na]
         mods    <- mods[not.na,,drop=FALSE]
         Z       <- Z[not.na,,drop=FALSE]
         k       <- length(yi)
         warning(mstyle$warning(paste(sum(has.na), ifelse(sum(has.na) > 1, "studies", "study"), "with NAs omitted from model fitting.")), call.=FALSE)

         attr(yi, "measure") <- measure # add measure attribute back
         attr(yi, "ni")      <- ni      # add ni attribute back

         ### note: slab is always of the same length as the full yi vector (after subsetting), so missings are not removed and slab is not added back to yi

      }

      if (na.act == "na.fail")
         stop(mstyle$stop("Missing values in data."))

   }

   ### at least one study left?

   if (k < 1L)
      stop(mstyle$stop("Processing terminated since k = 0."))

   ### check for non-positive sampling variances (and set negative values to 0)
   ### note: done after removing NAs since only the included studies are relevant

   if (any(vi <= 0)) {
      allvipos <- FALSE
      if (!vi0)
         warning(mstyle$warning("There are outcomes with non-positive sampling variances."), call.=FALSE)
      vi.neg <- vi < 0
      if (any(vi.neg)) {
         vi[vi.neg] <- 0
         warning(mstyle$warning("Negative sampling variances constrained to zero."), call.=FALSE)
      }
   } else {
      allvipos <- TRUE
   }

   ### but even in vi.f, constrain negative sampling variances to 0 (not needed)
   #vi.f[vi.f < 0] <- 0

   ### if k=1 and test != "z", set test="z" (other methods cannot be used)

   if (k == 1L && test != "z") {
      warning(mstyle$warning("Setting argument test=\"z\" since k=1."), call.=FALSE)
      test <- "z"
   }

   ### make sure that there is at least one column in X ([b])

   if (is.null(mods) && !intercept) {
      warning(mstyle$warning("Must either include an intercept and/or moderators in model.\nCoerced intercept into the model."), call.=FALSE)
      intercept <- TRUE
   }

   if (!is.null(mods) && ncol(mods) == 0L) {
      warning(mstyle$warning("Cannot fit model with an empty model matrix. Coerced intercept into the model."), call.=FALSE)
      intercept <- TRUE
   }

   ### add vector of 1s to the X matrix for the intercept (if intercept=TRUE)

   if (intercept) {
      X   <- cbind(intrcpt=rep(1,k), mods)
      X.f <- cbind(intrcpt=rep(1,k.f), mods.f)
   } else {
      X   <- mods
      X.f <- mods.f
   }

   ### drop redundant predictors
   ### note: need to save coef.na for functions that modify the data/model and then refit the model (regtest() and the
   ### various function that leave out an observation); so we can check if there are redundant/dropped predictors then

   tmp <- try(lm(yi ~ X - 1), silent=TRUE)
   if (inherits(tmp, "lm")) {
      coef.na <- is.na(coef(tmp))
   } else {
      coef.na <- rep(FALSE, NCOL(X))
   }
   if (any(coef.na)) {
      warning(mstyle$warning("Redundant predictors dropped from the model."), call.=FALSE)
      X   <- X[,!coef.na,drop=FALSE]
      X.f <- X.f[,!coef.na,drop=FALSE]
   }

   ### check whether intercept is included and if yes, move it to the first column (NAs already removed, so na.rm=TRUE for any() not necessary)

   is.int <- apply(X, 2, .is.intercept)
   if (any(is.int)) {
      int.incl <- TRUE
      int.indx <- which(is.int, arr.ind=TRUE)
      X        <- cbind(intrcpt=1,   X[,-int.indx, drop=FALSE]) # this removes any duplicate intercepts
      X.f      <- cbind(intrcpt=1, X.f[,-int.indx, drop=FALSE]) # this removes any duplicate intercepts
      intercept <- TRUE # set intercept appropriately so that the predict() function works
   } else {
      int.incl <- FALSE
   }

   p <- NCOL(X) # number of columns in X (including the intercept if it is included)

   ### make sure variable names in X and Z are unique

   colnames(X) <- colnames(X.f) <- .make.unique(colnames(X))
   colnames(Z) <- colnames(Z.f) <- .make.unique(colnames(Z))

   ### check whether this is an intercept-only model

   if ((p == 1L) && .is.intercept(X)) {
      int.only <- TRUE
   } else {
      int.only <- FALSE
   }

   ### check if there are too many parameters for given k (TODO: what about rma.ls models?)

   if (!(int.only && k == 1L)) {
      if (is.element(method[1], c("FE","EE","CE"))) { # have to estimate p parms
         if (p > k)
            stop(mstyle$stop("Number of parameters to be estimated is larger than the number of observations."))
      } else {
         if (!is.null(tau2) && !is.na(tau2)) {        # have to estimate p parms (tau2 is fixed at value specified)
            if (p > k)
               stop(mstyle$stop("Number of parameters to be estimated is larger than the number of observations."))
         } else {
            if ((p+1) > k)                            # have to estimate p+1 parms
               stop(mstyle$stop("Number of parameters to be estimated is larger than the number of observations."))
         }
      }
   }

   ### set/check 'btt' argument

   btt <- .set.btt(btt, p, int.incl, colnames(X))
   m <- length(btt) # number of betas to test (m = p if all betas are tested)

   #########################################################################

   ### set default control parameters

   con <- list(verbose = FALSE,
               evtol = 1e-07,                  # lower bound for eigenvalues to determine if model matrix is positive definite (also for checking if vimaxmin >= 1/con$evtol)
               REMLf = TRUE)                   # should |X'X| term be included in the REML log-likelihood?

   if (model == "rma.uni") {

      con <- c(con,
          list(tau2.init = NULL,               # initial value for iterative estimators (ML, REML, EB, SJ, SJIT, DLIT)
               tau2.min = 0,                   # lower bound for tau^2 value (passed down to confint.rma.uni())
               tau2.max = 100,                 # upper bound for tau^2 value (for PM/PMM/GENQM estimators) but see [c]
               threshold = 10^-5,              # convergence threshold (for ML, REML, EB, SJIT, DLIT)
               tol = .Machine$double.eps^0.25, # convergence tolerance for uniroot() as used for PM, PMM, GENQM (also used in 'll0 - ll > con$tol' check for ML/REML)
               ll0check = TRUE,                # should the 'll0 - ll > con$tol' check be conducted for ML/REML?
               maxiter = 100,                  # maximum number of iterations (for ML, REML, EB, SJIT, DLIT)
               stepadj = 1))                   # step size adjustment for Fisher scoring algorithm (for ML, REML, EB)

      ### [c] for some applications, tau2.max = 100 may not be enough; use an adaptive max instead

      con$tau2.max <- max(con$tau2.max, 10*mad(yi)^2)

   }

   if (model == "rma.ls") {

      con <- c(con,
          list(beta.init = NULL,               # initial values for location parameters (only relevant when optbeta=TRUE)
               hesspack = "numDeriv",          # package for computing the Hessian (numDeriv or pracma)
               optimizer = "nlminb",           # optimizer to use ("optim","nlminb","uobyqa","newuoa","bobyqa","nloptr","nlm","hjk","nmk","mads","ucminf","lbfgsb3c","subplex","BBoptim","optimParallel","constrOptim","solnp","alabama"/"constrOptim.nl","Rcgmin","Rvmmin")
               optmethod = "BFGS",             # argument 'method' for optim() ("Nelder-Mead" and "BFGS" are sensible options)
               parallel = list(),              # parallel argument for optimParallel() (note: 'cl' argument in parallel is not passed; this is directly specified via 'cl')
               cl = NULL,                      # arguments for optimParallel()
               ncpus = 1L,                     # arguments for optimParallel()
               tau2.min = 0,                   # lower bound for tau^2 values (can be used to constrain tau^2 values but see [d])
               tau2.max = Inf,                 # upper bound for tau^2 values (can be used to constrain tau^2 values but see [d])
               alpha.init = NULL,              # initial values for scale parameters
               alpha.min = -Inf,               # min possible value(s) for scale parameter(s)
               alpha.max = Inf,                # max possible value(s) for scale parameter(s)
               hessianCtrl=list(r=8),          # arguments passed on to 'method.args' of hessian()
               scaleZ = TRUE))                 # rescale Z matrix (only if Z.int.incl, is.na(alpha[1]), all(is.infinite(con$alpha.min)), all(is.infinite(con$alpha.max)), !optbeta)

      ### [d] can constrain the tau^2 values in location-scale models, but this is done in a very crude way
      ### in the optimization (by returning Inf when any tau^2 value falls outside the bounds) and this is
      ### not recommended/documented (instead, one can constrain the alpha values via alpha.min/alpha.max);
      ### note: the tau^2 bounds are only in effect when tau2.min or tau2.max are actually used in 'control'
      ### (if not, tau2.min and tau2.max are set to 0 and Inf, respectively)

   }

   ### replace defaults with any user-defined values

   con.pos <- pmatch(names(control), names(con))
   con[c(na.omit(con.pos))] <- control[!is.na(con.pos)]

   if (verbose)
      con$verbose <- verbose

   verbose <- con$verbose

   if (model == "rma.ls") {

      con$hesspack <- match.arg(con$hesspack, c("numDeriv","pracma"))

      if (!isTRUE(ddd$skiphes) && !requireNamespace(con$hesspack, quietly=TRUE))
         stop(mstyle$stop(paste0("Please install the '", con$hesspack, "' package to compute the Hessian.")))

   }

   if (model == "rma.uni") {

      ### constrain a negative tau2.min value to -min(vi) (to ensure that the marginal variance is always >= 0)

      if (con$tau2.min < 0 && (con$tau2.min < -min(vi))) {
         con$tau2.min <- -min(vi) # + .Machine$double.eps^0.25 # to force tau2.min just above -min(vi)
         warning(mstyle$warning(paste0("Value of 'tau2.min' constrained to -min(vi) = ", fmtx(-min(vi), digits[["est"]]), ".")), call.=FALSE)
      }

   } else {

      ### constrain a negative tau2.min value to 0 for ls models

      if (is.element("tau2.min", names(control)))
         con$tau2.min[con$tau2.min < 0] <- 0

   }

   ### check whether model matrix is of full rank

   if (!.chkpd(crossprod(X), tol=con$evtol))
      stop(mstyle$stop("Model matrix not of full rank. Cannot fit model."))

   ### check ratio of largest to smallest sampling variance
   ### note: need to exclude some special cases (0/0 = NaN, max(vi)/0 = Inf)
   ### TODO: use the condition number of diag(vi) here instead?

   vimaxmin <- max(vi) / min(vi)

   if (is.finite(vimaxmin) && vimaxmin >= 1/con$evtol)
      warning(mstyle$warning("Ratio of largest to smallest sampling variance extremely large. May not be able to obtain stable results."), call.=FALSE)

   ### set some defaults

   se.tau2 <- I2 <- H2 <- QE <- QEp <- NA_real_
   s2w <- 1
   level <- .level(level)

   Y <- as.matrix(yi)

   ### mean center yi for some calculations to increase the stability of the computations

   ymci <- scale(yi, center=TRUE, scale=FALSE)
   Ymc  <- as.matrix(ymci)

   #########################################################################

   ###### heterogeneity estimation for the standard normal-normal model (rma.uni)

   tau2.inf <- FALSE

   if (model == "rma.uni") {

      if (!is.null(tau2) && !is.na(tau2) && !is.element(method[1], c("FE","EE","CE"))) { # if user has fixed the tau2 value
         tau2.fix <- TRUE
         tau2.arg <- tau2
         tau2.inf <- identical(tau2, Inf)
      } else {
         tau2.fix <- FALSE
         tau2.arg <- NA_real_
      }

      if (verbose > 1 && !tau2.fix && !is.element(method[1], c("FE","EE","CE")))
         message(mstyle$message("Estimating tau^2 value ...\n"))

      if (k == 1L) {
         method.sav <- method[1]
         method <- "k1" # set method to k1 so all of the stuff below is skipped
         if (!tau2.fix)
            tau2 <- 0
      }

      conv <- FALSE

      while (!conv && !tau2.inf) {

         ### convergence indicator and change variable

         conv   <- TRUE # assume TRUE for now unless things go wrong below
         change <- con$threshold + 1

         ### iterations counter for iterative estimators (i.e., DLIT, SJIT, ML, REML, EB)
         ### (note: PM, PMM, and GENQM are also iterative, but uniroot() handles that)

         iter <- 0

         ### Hunter & Schmidt (HS) estimator (or k-corrected HS estimator (HSk))

         if (is.element(method[1], c("HS","HSk"))) {

            if (!allvipos)
               stop(mstyle$stop(paste0(method[1], " estimator cannot be used when there are non-positive sampling variances in the data.")))

            wi    <- 1/vi
            W     <- diag(wi, nrow=k, ncol=k)
            stXWX <- .invcalc(X=X, W=W, k=k)
            P     <- W - W %*% X %*% stXWX %*% crossprod(X,W)
            RSS   <- crossprod(Ymc,P) %*% Ymc
            if (method[1] == "HS") {
               tau2 <- ifelse(tau2.fix, tau2.arg, (RSS - k) / sum(wi))
            } else {
               tau2 <- ifelse(tau2.fix, tau2.arg, (k/(k-p)*RSS - k) / sum(wi))
            }

         }

         ### Hedges (HE) estimator (or initial value for ML, REML, EB)

         if (is.element(method[1], c("HE","ML","REML","EB"))) {

            stXX  <- .invcalc(X=X, W=diag(k), k=k)
            P     <- diag(k) - X %*% tcrossprod(stXX,X)
            RSS   <- crossprod(Ymc,P) %*% Ymc
            V     <- diag(vi, nrow=k, ncol=k)
            PV    <- P %*% V # note: this is not symmetric
            trPV  <- .tr(PV) # since PV needs to be computed anyway, can use .tr()
            tau2  <- ifelse(tau2.fix, tau2.arg, (RSS - trPV) / (k-p))

         }

         ### DerSimonian-Laird (DL) estimator

         if (method[1] == "DL") {

            if (!allvipos)
               stop(mstyle$stop("DL estimator cannot be used when there are non-positive sampling variances in the data."))

            wi    <- 1/vi
            W     <- diag(wi, nrow=k, ncol=k)
            stXWX <- .invcalc(X=X, W=W, k=k)
            P     <- W - W %*% X %*% stXWX %*% crossprod(X,W)
            RSS   <- crossprod(Ymc,P) %*% Ymc
            trP   <- .tr(P)
            tau2  <- ifelse(tau2.fix, tau2.arg, (RSS - (k-p)) / trP)

         }

         ### DerSimonian-Laird (DL) estimator with iteration (when this converges, same as PM)

         if (method[1] == "DLIT") {

            if (is.null(con$tau2.init)) {
               tau2 <- 0
            } else {
               tau2 <- con$tau2.init
            }

            while (change > con$threshold) {

               if (verbose)
                  cat(mstyle$verbose(paste("Iteration", formatC(iter, width=5, flag="-", format="f", digits=0), "tau^2 =", fmtx(tau2, digits[["var"]]), "\n")))

               iter <- iter + 1
               old2 <- tau2
               wi   <- 1/(vi + tau2)
               if (any(tau2 + vi < 0))
                  stop(mstyle$stop("Some marginal variances are negative."))
               if (any(is.infinite(wi)))
                  stop(mstyle$stop("Division by zero when computing the inverse variance weights."))
               W     <- diag(wi, nrow=k, ncol=k)
               stXWX <- .invcalc(X=X, W=W, k=k)
               P     <- W - W %*% X %*% stXWX %*% crossprod(X,W)
               V     <- diag(vi, nrow=k, ncol=k)
               trP   <- .tr(P)
               trPV  <- .tr(P %*% V)
               RSS   <- crossprod(Ymc,P) %*% Ymc
               tau2  <- ifelse(tau2.fix, tau2.arg, (RSS - trPV) / trP)
               tau2[tau2 < con$tau2.min] <- con$tau2.min
               change <- abs(old2 - tau2)

               if (iter > con$maxiter) {
                  conv <- FALSE
                  break
               }

            }

            if (!conv) {
               if (length(method) == 1L) {
                  stop(mstyle$stop("Iterative DL estimator did not converge."))
               } else {
                  if (verbose)
                     warning(mstyle$warning("Iterative DL estimator did not converge."), call.=FALSE)
               }
            }

         }

         ### generalized Q-statistic estimator

         if (method[1] == "GENQ") {

            #if (!allvipos)
            #   stop(mstyle$stop("GENQ estimator cannot be used when there are non-positive sampling variances in the data."))

            if (is.null(weights))
               stop(mstyle$stop("Must specify 'weights' when method='GENQ'."))

            A     <- diag(weights, nrow=k, ncol=k)
            stXAX <- .invcalc(X=X, W=A, k=k)
            P     <- A - A %*% X %*% stXAX %*% crossprod(X,A)
            V     <- diag(vi, nrow=k, ncol=k)
            PV    <- P %*% V # note: this is not symmetric
            trP   <- .tr(P)
            trPV  <- .tr(PV)
            RSS   <- crossprod(Ymc,P) %*% Ymc
            tau2  <- ifelse(tau2.fix, tau2.arg, (RSS - trPV) / trP)

         }

         ### generalized Q-statistic estimator (median unbiased version)

         if (method[1] == "GENQM") {

            if (is.null(weights))
               stop(mstyle$stop("Must specify 'weights' when method='GENQM'."))

            A     <- diag(weights, nrow=k, ncol=k)
            stXAX <- .invcalc(X=X, W=A, k=k)
            P     <- A - A %*% X %*% stXAX %*% crossprod(X,A)
            V     <- diag(vi, nrow=k, ncol=k)
            PV    <- P %*% V # note: this is not symmetric
            trP   <- .tr(P)

            if (!tau2.fix) {

               RSS <- crossprod(Ymc,P) %*% Ymc

               if (.GENQ.func(con$tau2.min, P=P, vi=vi, Q=RSS, level=0, k=k, p=p, getlower=TRUE) > 0.5) {

                  ### if GENQ.tau2.min is > 0.5, then estimate < tau2.min

                  tau2 <- con$tau2.min

               } else {

                  if (.GENQ.func(con$tau2.max, P=P, vi=vi, Q=RSS, level=0, k=k, p=p, getlower=TRUE) < 0.5) {

                     ### if GENQ.tau2.max is < 0.5, then estimate > tau2.max

                     conv <- FALSE

                     if (length(method) == 1L) {
                        stop(mstyle$stop("Value of 'tau2.max' too low. Try increasing 'tau2.max' or switch to another 'method'."))
                     } else {
                        if (verbose)
                           warning(mstyle$warning("Value of 'tau2.max' too low. Try increasing 'tau2.max' or switch to another 'method'."), call.=FALSE)
                     }

                  } else {

                     tau2 <- try(uniroot(.GENQ.func, interval=c(con$tau2.min, con$tau2.max), tol=con$tol, maxiter=con$maxiter, P=P, vi=vi, Q=RSS, level=0.5, k=k, p=p, getlower=FALSE, verbose=verbose, digits=digits, extendInt="no")$root, silent=TRUE)

                     if (inherits(tau2, "try-error")) {
                        conv <- FALSE
                        if (length(method) == 1L) {
                           stop(mstyle$stop("Error in iterative search for tau^2 using uniroot()."))
                        } else {
                           if (verbose)
                              warning(mstyle$warning("Error in iterative search for tau^2 using uniroot()."), call.=FALSE)
                        }
                     }

                  }

               }

            } else {

               tau2 <- tau2.arg

            }

         }

         ### Sidik-Jonkman (SJ) estimator

         if (method[1] == "SJ") {

            if (is.null(con$tau2.init)) {
               tau2.0 <- c(var(ymci) * (k-1)/k)
            } else {
               tau2.0 <- con$tau2.init
            }

            wi    <- 1/(vi + tau2.0)
            W     <- diag(wi, nrow=k, ncol=k)
            stXWX <- .invcalc(X=X, W=W, k=k)
            P     <- W - W %*% X %*% stXWX %*% crossprod(X,W)
            RSS   <- crossprod(Ymc,P) %*% Ymc
            V     <- diag(vi, nrow=k, ncol=k)
            PV    <- P %*% V # note: this is not symmetric
            tau2  <- ifelse(tau2.fix, tau2.arg, tau2.0 * RSS / (k-p))

         }

         ### Sidik-Jonkman (SJ) estimator with iteration

         if (method[1] == "SJIT") {

            if (is.null(con$tau2.init)) {
               tau2 <- c(var(ymci) * (k-1)/k)
            } else {
               tau2 <- con$tau2.init
            }

            tau2.0 <- tau2

            while (change > con$threshold) {

               if (verbose)
                  cat(mstyle$verbose(paste("Iteration", formatC(iter, width=5, flag="-", format="f", digits=0), "tau^2 =", fmtx(tau2, digits[["var"]]), "\n")))

               iter <- iter + 1
               old2 <- tau2

               wi     <- 1/(vi + tau2)
               W      <- diag(wi, nrow=k, ncol=k)
               stXWX  <- .invcalc(X=X, W=W, k=k)
               P      <- W - W %*% X %*% stXWX %*% crossprod(X,W)
               RSS    <- crossprod(Ymc,P) %*% Ymc
               V      <- diag(vi, nrow=k, ncol=k)
               PV     <- P %*% V # note: this is not symmetric
               tau2   <- ifelse(tau2.fix, tau2.arg, tau2 * RSS / (k-p))
               change <- abs(old2 - tau2)

               if (iter > con$maxiter) {
                  conv <- FALSE
                  break
               }

            }

            if (!conv) {
               if (length(method) == 1L) {
                  stop(mstyle$stop("Iterative SJ estimator did not converge."))
               } else {
                  if (verbose)
                     warning(mstyle$warning("Iterative SJ estimator did not converge."), call.=FALSE)
               }
            }

         }

         ### Paule-Mandel (PM) estimator (regular and median unbiased version)

         if (is.element(method[1], c("PM","MP","PMM"))) {

            if (!allvipos)
               stop(mstyle$stop(method[1], " estimator cannot be used when there are non-positive sampling variances in the data."))

            if (method[1] == "PMM") {
               target <- qchisq(0.5, df=k-p)
            } else {
               target <- k-p
            }

            if (!tau2.fix) {

               if (.QE.func(con$tau2.min, Y=Ymc, vi=vi, X=X, k=k, objective=0) < target) {

                  tau2 <- con$tau2.min

               } else {

                  if (.QE.func(con$tau2.max, Y=Ymc, vi=vi, X=X, k=k, objective=0) > target) {

                     conv <- FALSE

                     if (length(method) == 1L) {
                        stop(mstyle$stop("Value of 'tau2.max' too low. Try increasing 'tau2.max' or switch to another 'method'."))
                     } else {
                        if (verbose)
                           warning(mstyle$warning("Value of 'tau2.max' too low. Try increasing 'tau2.max' or switch to another 'method'."), call.=FALSE)
                     }

                  } else {

                     tau2 <- try(uniroot(.QE.func, interval=c(con$tau2.min, con$tau2.max), tol=con$tol, maxiter=con$maxiter, Y=Ymc, vi=vi, X=X, k=k, objective=target, verbose=verbose, digits=digits, extendInt="no")$root, silent=TRUE)

                     if (inherits(tau2, "try-error")) {
                        conv <- FALSE
                        if (length(method) == 1L) {
                           stop(mstyle$stop("Error in iterative search for tau^2 using uniroot()."))
                        } else {
                           if (verbose)
                              warning(mstyle$warning("Error in iterative search for tau^2 using uniroot()."), call.=FALSE)
                        }
                     }

                  }

               }

               #W <- diag(wi, nrow=k, ncol=k)
               #stXWX <- .invcalc(X=X, W=W, k=k)
               #P <- W - W %*% X %*% stXWX %*% crossprod(X,W) # needed for se.tau2 computation below (not when using the simpler equation)

            } else {

               tau2 <- tau2.arg

            }

         }

         ### maximum-likelihood (ML), restricted maximum-likelihood (REML), and empirical Bayes (EB) estimators

         if (is.element(method[1], c("ML","REML","EB"))) {

            if (is.null(con$tau2.init)) {         # check if user specified initial value for tau2
               tau2 <- max(0, tau2, con$tau2.min) # if not, use HE estimate (or con$tau2.min) as initial estimate for tau2
            } else {
               tau2 <- con$tau2.init              # if yes, use value specified by user
            }

            while (change > con$threshold) {

               if (verbose)
                  cat(mstyle$verbose(paste(mstyle$verbose(paste("Iteration", formatC(iter, width=5, flag="-", format="f", digits=0), "tau^2 =", fmtx(tau2, digits[["var"]]), "\n")))))

               iter <- iter + 1
               old2 <- tau2
               wi   <- 1/(vi + tau2)
               if (any(tau2 + vi < 0))
                  stop(mstyle$stop("Some marginal variances are negative."))
               if (any(is.infinite(wi)))
                  stop(mstyle$stop("Division by zero when computing the inverse variance weights."))
               W     <- diag(wi, nrow=k, ncol=k)
               stXWX <- .invcalc(X=X, W=W, k=k)
               P     <- W - W %*% X %*% stXWX %*% crossprod(X,W)

               if (method[1] == "ML") {
                  PP  <- P %*% P
                  adj <- c(crossprod(Ymc,PP) %*% Ymc - sum(wi)) / sum(wi^2)
               }
               if (method[1] == "REML") {
                  PP  <- P %*% P
                  adj <- c(crossprod(Ymc,PP) %*% Ymc - .tr(P)) / .tr(PP)
               }
               if (method[1] == "EB") {
                  adj <- c(crossprod(Ymc,P) %*% Ymc * k/(k-p) - k) / sum(wi)
               }

               adj <- c(adj) * con$stepadj # apply (user-defined) step adjustment

               if (is.na(adj)) # can happen for a saturated model when fixing tau^2
                  adj <- 0

               while (tau2 + adj < con$tau2.min) # use step-halving if necessary
                  adj <- adj / 2

               tau2   <- ifelse(tau2.fix, tau2.arg, tau2 + adj)
               change <- abs(old2 - tau2)

               if (iter > con$maxiter) {
                  conv <- FALSE
                  break
               }

            }

            if (!conv) {
               if (length(method) == 1L) {
                  stop(mstyle$stop("Fisher scoring algorithm did not converge. See 'help(rma)' for possible remedies."))
               } else {
                  if (verbose)
                     warning(mstyle$warning("Fisher scoring algorithm did not converge. See 'help(rma)' for possible remedies."), call.=FALSE)
               }
            }

            ### check if ll is larger when tau^2 = 0 (only if ll0check=TRUE and only possible/sensible if allvipos and !tau2.fix)
            ### note: this doesn't catch the case where tau^2 = 0 is a local maximum

            if (conv && is.element(method[1], c("ML","REML")) && con$ll0check && allvipos && !tau2.fix) {

               wi    <- 1/vi
               W     <- diag(wi, nrow=k, ncol=k)
               stXWX <- .invcalc(X=X, W=W, k=k)
               beta  <- stXWX %*% crossprod(X,W) %*% Ymc
               RSS   <- sum(wi*(ymci - X %*% beta)^2)
               if (method[1] == "ML")
                  ll0 <- -1/2 * (k)   * log(2*base::pi) - 1/2 * sum(log(vi)) - 1/2 * RSS
               if (method[1] == "REML")
                  ll0 <- -1/2 * (k-p) * log(2*base::pi) - 1/2 * sum(log(vi)) - 1/2 * determinant(crossprod(X,W) %*% X, logarithm=TRUE)$modulus - 1/2 * RSS

               wi    <- 1/(vi + tau2)
               if (any(tau2 + vi < 0))
                  stop(mstyle$stop("Some marginal variances are negative."))
               if (any(is.infinite(wi)))
                  stop(mstyle$stop("Division by zero when computing the inverse variance weights."))
               W     <- diag(wi, nrow=k, ncol=k)
               stXWX <- .invcalc(X=X, W=W, k=k)
               beta  <- stXWX %*% crossprod(X,W) %*% Ymc
               RSS   <- sum(wi*(ymci - X %*% beta)^2)
               if (method[1] == "ML")
                  ll <- -1/2 * (k)   * log(2*base::pi) - 1/2 * sum(log(vi + tau2)) - 1/2 * RSS
               if (method[1] == "REML")
                  ll <- -1/2 * (k-p) * log(2*base::pi) - 1/2 * sum(log(vi + tau2)) - 1/2 * determinant(crossprod(X,W) %*% X, logarithm=TRUE)$modulus - 1/2 * RSS

               if (ll0 - ll > con$tol && tau2 > con$threshold) {
                  warning(mstyle$warning("Fisher scoring algorithm may have gotten stuck at a local maximum.\nSetting tau^2 = 0. Check the profile likelihood plot with profile()."), call.=FALSE)
                  tau2 <- 0
               }

            }

            ### need to run this so that wi and P are based on the final tau^2 value

            if (conv) {
               wi <- 1/(vi + tau2)
               if (any(tau2 + vi < 0))
                  stop(mstyle$stop("Some marginal variances are negative."))
               if (any(is.infinite(wi)))
                  stop(mstyle$stop("Division by zero when computing the inverse variance weights."))
               W <- diag(wi, nrow=k, ncol=k)
               stXWX <- .invcalc(X=X, W=W, k=k)
               P <- W - W %*% X %*% stXWX %*% crossprod(X,W)
            }

         }

         if (conv) {

            ### make sure that tau2 is >= con$tau2.min

            tau2 <- max(con$tau2.min, c(tau2))

            ### check if any marginal variances are negative (only possible if user has changed tau2.min)

            if (!is.na(tau2) && any(tau2 + vi < 0))
               stop(mstyle$stop("Some marginal variances are negative."))

            ### verbose output upon convergence for ML/REML/EB estimators

            if (verbose && is.element(method[1], c("ML","REML","EB"))) {
               cat(mstyle$verbose(paste("Iteration", formatC(iter, width=5, flag="-", format="f", digits=0), "tau^2 =", fmtx(tau2, digits[["var"]]), "\n")))
               cat(mstyle$verbose(paste("Fisher scoring algorithm converged after", iter, "iterations.\n")))
            }

            ### standard error of the tau^2 estimators (also when the user has fixed/specified a tau^2 value)
            ### see notes.pdf and note: .tr(P%*%P) = sum(P*t(P)) = sum(P*P) (since P is symmetric)

            if (method[1] == "HS")
               se.tau2 <- sqrt(1/sum(wi)^2 * (2*(k-p) + 4*max(tau2,0)*.tr(P) + 2*max(tau2,0)^2*sum(P*P))) # note: wi = 1/vi
            if (method[1] == "HSk")
               se.tau2 <- k/(k-p) * sqrt(1/sum(wi)^2 * (2*(k-p) + 4*max(tau2,0)*.tr(P) + 2*max(tau2,0)^2*sum(P*P)))
            if (method[1] == "HE")
               se.tau2 <- sqrt(1/(k-p)^2 * (2*sum(PV*t(PV)) + 4*max(tau2,0)*trPV + 2*max(tau2,0)^2*(k-p)))
            if (method[1] == "DL")
               se.tau2 <- sqrt(1/trP^2 * (2*(k-p) + 4*max(tau2,0)*trP + 2*max(tau2,0)^2*sum(P*P)))
            if (is.element(method[1], c("GENQ","GENQM")))
               se.tau2 <- sqrt(1/trP^2 * (2*sum(PV*t(PV)) + 4*max(tau2,0)*sum(PV*P) + 2*max(tau2,0)^2*sum(P*P)))
            if (method[1] == "SJ")
               se.tau2 <- sqrt(tau2.0^2/(k-p)^2 * (2*sum(PV*t(PV)) + 4*max(tau2,0)*sum(PV*P) + 2*max(tau2,0)^2*sum(P*P)))
            if (method[1] == "ML")
               se.tau2 <- sqrt(2/sum(wi^2)) # note: wi = 1/(vi + tau2) for ML, REML, EB, PM, PMM, and SJIT
            if (method[1] == "REML")
               se.tau2 <- sqrt(2/sum(P*P)) # based on Fisher information matrix
               #se.tau2 <- sqrt(1 / (t(Ymc) %*% P %*% P %*% P %*% Ymc - 1/2 * sum(P*P))) # based on Hessian
            if (is.element(method[1], c("EB","PM","MP","PMM","DLIT","SJIT"))) {
               wi  <- 1/(vi + tau2)
               #V  <- diag(vi, nrow=k, ncol=k)
               #PV <- P %*% V # note: this is not symmetric
               #se.tau2 <- sqrt((k/(k-p))^2 / sum(wi)^2 * (2*sum(PV*t(PV)) + 4*max(tau2,0)*sum(PV*P) + 2*max(tau2,0)^2*sum(P*P)))
               se.tau2 <- sqrt(2*k^2/(k-p) / sum(wi)^2) # these two equations are actually identical, but this one is much simpler
            }

         } else {

            method <- method[-1]

         }

      }

      if (k == 1L)
         method <- method.sav

   }

   #########################################################################

   ###### parameter estimation for the location-scale model (rma.ls)

   if (model == "rma.ls") {

      if (!is.element(method[1], c("ML","REML")))
         stop(mstyle$stop("Location-scale models can only be fitted with ML or REML estimation."))

      tau2.fix <- FALSE

      if (!is.null(tau2) && !is.na(tau2))
         warning(mstyle$warning("Argument 'tau2' ignored for location-scale models."), call.=FALSE)

      ### get optimizer arguments from control argument

      optimizer <- match.arg(con$optimizer, c("optim","nlminb","uobyqa","newuoa","bobyqa","nloptr","nlm","hjk","nmk","mads","ucminf","lbfgsb3c","subplex","BBoptim","optimParallel","constrOptim","solnp","alabama","constrOptim.nl","Nelder-Mead","BFGS","CG","L-BFGS-B","SANN","Brent","Rcgmin","Rvmmin"))
      optmethod <- match.arg(con$optmethod, c("Nelder-Mead","BFGS","CG","L-BFGS-B","SANN","Brent"))
      if (optimizer %in% c("Nelder-Mead","BFGS","CG","L-BFGS-B","SANN","Brent")) {
         optmethod <- optimizer
         optimizer <- "optim"
      }
      parallel   <- con$parallel
      cl         <- con$cl
      ncpus      <- con$ncpus
      optcontrol <- control[is.na(con.pos)] # get arguments that are control arguments for optimizer

      if (length(optcontrol) == 0L)
         optcontrol <- list()

      ### if control argument 'ncpus' is larger than 1, automatically switch to optimParallel optimizer

      if (ncpus > 1L)
         optimizer <- "optimParallel"

      ### can use optimizer="alabama" as a shortcut for optimizer="constrOptim.nl"

      if (optimizer == "alabama")
         optimizer <- "constrOptim.nl"

      ### when using an identity link, automatically set 'constrOptim' as the default optimizer (but 'solnp' by default when optbeta=TRUE)

      if (link == "identity") {
         if (optbeta) {
            if (optimizer == "nlminb") {
               optimizer <- "solnp"
            } else {
               if (!is.element(optimizer, c("solnp","nloptr","constrOptim.nl"))) {
                  optimizer <- "solnp"
                  warning(mstyle$warning(paste0("Can only use optimizers 'solnp', 'nloptr', or 'constrOptim.nl' when link='identity' and optbeta=TRUE (resetting to '", optimizer, "').")), call.=FALSE)
               }
            }
         } else {
            if (optimizer == "nlminb") {
               optimizer <- "constrOptim"
            } else {
               if (!is.element(optimizer, c("constrOptim","solnp","nloptr","constrOptim.nl"))) {
                  optimizer <- "constrOptim"
                  warning(mstyle$warning(paste0("Can only use optimizers 'constrOptim', 'solnp', 'nloptr', or 'constrOptim.nl' when link='identity' (resetting to '", optimizer, "').")), call.=FALSE)
               }
            }
         }
      }

      if (link == "log" && is.element(optimizer, c("constrOptim","constrOptim.nl")))
         stop(mstyle$stop(paste0("Cannot use '", optimizer, "' optimizer when using a log link."))) # but can use solnp and nloptr

      reml <- ifelse(method[1] == "REML", TRUE, FALSE)

      ### drop redundant predictors

      tmp <- try(lm(yi ~ Z - 1), silent=TRUE)
      if (inherits(tmp, "lm")) {
         coef.na.Z <- is.na(coef(tmp))
      } else {
         coef.na.Z <- rep(FALSE, NCOL(Z))
      }
      if (any(coef.na.Z)) {
         warning(mstyle$warning("Redundant predictors dropped from the scale model."), call.=FALSE)
         Z   <- Z[,!coef.na.Z,drop=FALSE]
         Z.f <- Z.f[,!coef.na.Z,drop=FALSE]
      }

      ### check whether intercept is included and if yes, move it to the first column (NAs already removed, so na.rm=TRUE for any() not necessary)

      is.int <- apply(Z, 2, .is.intercept)
      if (any(is.int)) {
         Z.int.incl <- TRUE
         int.indx   <- which(is.int, arr.ind=TRUE)
         Z          <- cbind(intrcpt=1,   Z[,-int.indx, drop=FALSE]) # this removes any duplicate intercepts
         Z.f        <- cbind(intrcpt=1, Z.f[,-int.indx, drop=FALSE]) # this removes any duplicate intercepts
         Z.intercept <- TRUE # set intercept appropriately so that the predict() function works
      } else {
         Z.int.incl <- FALSE
      }

      q <- NCOL(Z) # number of columns in Z (including the intercept if it is included)

      ### check whether model matrix is of full rank

      if (!.chkpd(crossprod(Z), tol=con$evtol))
         stop(mstyle$stop("Model matrix for scale part of the model not of full rank. Cannot fit model."))

      ### check whether this is an intercept-only model

      is.int <- apply(Z, 2, .is.intercept)
      if (q == 1L && is.int) {
         Z.int.only <- TRUE
      } else {
         Z.int.only <- FALSE
      }

      ### checks on alpha argument

      if (missing(alpha) || is.null(alpha) || all(is.na(alpha))) {
         alpha <- rep(NA_real_, q)
      } else {
         if (length(alpha) == 1L)
            alpha <- rep(alpha, q)
         if (length(alpha) != q)
            stop(mstyle$stop(paste0("Length of 'alpha' argument (", length(alpha), ") does not match actual number of parameters (", q, ").")))
      }

      ### checks on beta argument

      if (optbeta) {

         if (missing(beta) || is.null(beta) || all(is.na(beta))) {
            beta <- rep(NA_real_, p)
         } else {
            if (length(beta) == 1L)
               beta <- rep(beta, p)
            if (length(beta) != p)
               stop(mstyle$stop(paste0("Length of 'beta' argument (", length(beta), ") does not match actual number of parameters (", p, ").")))
         }

      }

      ### rescale Z matrix (only for models with moderators, models including a non-fixed intercept term, when not placing constraints on alpha, and when not optimizing over beta)

      if (!Z.int.only && Z.int.incl && con$scaleZ && is.na(alpha[1]) && all(is.infinite(con$alpha.min)) && all(is.infinite(con$alpha.max)) && !optbeta) {
         Zsave <- Z
         meanZ <- colMeans(Z[, 2:q, drop=FALSE])
         sdZ   <- apply(Z[, 2:q, drop=FALSE], 2, sd) # consider using colSds() from matrixStats package
         is.d  <- apply(Z, 2, .is.dummy) # is each column a dummy variable (i.e., only 0s and 1s)?
         mZ    <- rbind(c(intrcpt=1, -1*ifelse(is.d[-1], 0, meanZ/sdZ)), cbind(0, diag(ifelse(is.d[-1], 1, 1/sdZ), nrow=length(is.d)-1, ncol=length(is.d)-1)))
         imZ   <- try(suppressWarnings(solve(mZ)), silent=TRUE)
         Z[,!is.d] <- apply(Z[, !is.d, drop=FALSE], 2, scale) # rescale the non-dummy variables
         if (any(!is.na(alpha))) {
            if (inherits(imZ, "try-error"))
               stop(mstyle$stop("Unable to rescale starting values for the scale parameters."))
            alpha <- diag(imZ) * alpha
         }
      } else {
         mZ <- NULL
      }

      if (k == 1L && Z.int.only) {
         if (link == "log")
            con$alpha.init <- -10000
         if (link == "identity")
            con$alpha.init <- 0.00001
      }

      ### set/transform/check alpha.init

      if (verbose > 1)
         message(mstyle$message("Extracting/computing initial values ..."))

      if (is.null(con$alpha.init)) {

         fit <- suppressWarnings(rma.uni(yi, vi, mods=X, intercept=FALSE, method="HE", skipr2=TRUE))
         tmp <- rstandard(fit)

         if (link == "log") {

            tmp <- suppressWarnings(rma.uni(log(tmp$resid^2), 4/tmp$resid^2*tmp$se^2, mods=Z, intercept=FALSE, method="FE"))
            #tmp <- rma.uni(log(tmp$resid^2), 4/tmp$resid^2*tmp$se^2, mods=Z, intercept=FALSE, method="FE")
            #tmp <- rma.uni(log(tmp$resid^2), tmp$se^2, mods=Z, intercept=FALSE, method="FE")
            #tmp <- rma.uni(log(tmp$resid^2), 1, mods=Z, intercept=FALSE, method="FE")
            alpha.init <- coef(tmp)

         }

         if (link == "identity") {

            #tmp <- rma.uni(tmp$resid^2, 4*tmp$resid^2*tmp$se^2, mods=Z, intercept=FALSE, method="FE")
            tmp <- suppressWarnings(rma.uni(tmp$resid^2, tmp$se^2, mods=Z, intercept=FALSE, method="FE"))
            #tmp <- rma.uni(tmp$resid^2, 1, mods=Z, intercept=FALSE, method="FE")
            alpha.init <- coef(tmp)
            if (any(Z %*% alpha.init < 0))
               alpha.init <- ifelse(is.int, fit$tau2+.01, 0)
            if (any(Z %*% alpha.init < 0))
               stop(mstyle$stop("Unable to find suitable starting values for the scale parameters."))

         }

      } else {

         alpha.init <- con$alpha.init

         if (!is.null(mZ)) {
            if (inherits(imZ, "try-error"))
               stop(mstyle$stop("Unable to rescale starting values for the scale parameters."))
            alpha.init <- c(imZ %*% cbind(alpha.init))
         }

         if (link == "identity" && any(Z %*% alpha.init < 0))
            stop(mstyle$stop("Starting values for the scale parameters lead to one or more negative tau^2 values."))

         if (optbeta)
            fit <- suppressWarnings(rma.uni(yi, vi, mods=X, intercept=FALSE, method="HE", skipr2=TRUE))

      }

      if (length(alpha.init) != q)
         stop(mstyle$stop(paste0("Length of 'alpha.init' argument (", length(alpha.init), ") does not match actual number of parameters (", q, ").")))

      if (anyNA(alpha.init))
         stop(mstyle$stop("No missing values allowed in 'alpha.init'."))

      if (optbeta) {

         if (is.null(con$beta.init)) {

            beta.init <- c(fit$beta)

         } else {

            beta.init <- con$beta.init

            if (length(beta.init) != p)
               stop(mstyle$stop(paste0("Length of 'beta.init' argument (", length(beta.init), ") does not match actual number of parameters (", p, ").")))

            if (anyNA(beta.init))
               stop(mstyle$stop("No missing values allowed in 'beta.init'."))

         }

      } else {

         beta.init <- NULL

      }

      ### set potential constraints on alpha values

      if (length(con$alpha.min) == 1L)
            con$alpha.min <- rep(con$alpha.min, q)

      if (length(con$alpha.max) == 1L)
            con$alpha.max <- rep(con$alpha.max, q)

      if (length(con$alpha.min) != q)
         stop(mstyle$stop(paste0("Length of 'alpha.min' argument (", length(alpha.min), ") does not match actual number of parameters (", q, ").")))

      if (length(con$alpha.max) != q)
         stop(mstyle$stop(paste0("Length of 'alpha.max' argument (", length(alpha.max), ") does not match actual number of parameters (", q, ").")))

      if (any(xor(is.infinite(con$alpha.min),is.infinite(con$alpha.max))))
         stop(mstyle$stop("Constraints on scale coefficients must be placed on both the lower and upper bound."))

      alpha.min <- con$alpha.min
      alpha.max <- con$alpha.max

      if (link == "identity" && (any(alpha.min != -Inf) || any(alpha.max != Inf)))
         stop(mstyle$stop("Cannot use constraints on scale coefficients when using an identity link."))

      alpha.init <- pmax(alpha.init, alpha.min)
      alpha.init <- pmin(alpha.init, alpha.max)

      alpha.init <- mapply(.mapinvfun.alpha, alpha.init, alpha.min, alpha.max)

      ### estimate alpha (and beta) values

      if (verbose > 1)
         message(mstyle$message("Estimating scale parameters ...\n"))

      tmp <- .chkopt(optimizer, optcontrol, ineq=link=="identity")
      optimizer  <- tmp$optimizer
      optcontrol <- tmp$optcontrol
      par.arg    <- tmp$par.arg
      ctrl.arg   <- tmp$ctrl.arg

      ### set up default cluster when using optimParallel

      if (optimizer == "optimParallel::optimParallel") {

         parallel$cl <- NULL

         if (is.null(cl)) {

            ncpus <- as.integer(ncpus)

            if (ncpus < 1L)
               stop(mstyle$stop("Control argument 'ncpus' must be >= 1."))

            cl <- parallel::makePSOCKcluster(ncpus)
            on.exit(parallel::stopCluster(cl), add=TRUE)

         } else {

            if (!inherits(cl, "SOCKcluster"))
               stop(mstyle$stop("Specified cluster is not of class 'SOCKcluster'."))

         }

         parallel$cl <- cl

         if (is.null(parallel$forward))
            parallel$forward <- FALSE

         if (is.null(parallel$loginfo)) {
            if (verbose) {
               parallel$loginfo <- TRUE
            } else {
               parallel$loginfo <- FALSE
            }
         }

      }

      #return(list(con=con, optimizer=optimizer, optmethod=optmethod, optcontrol=optcontrol, ctrl.arg=ctrl.arg))

      if (link == "log") {

         optcall <- paste0(optimizer, "(", par.arg, "=c(beta.init, alpha.init), .ll.rma.ls, ", ifelse(optimizer=="optim", "method=optmethod, ", ""),
                          "yi=yi, vi=vi, X=X, Z=Z, reml=reml, k=k, pX=p, alpha.arg=alpha, beta.arg=beta, verbose=verbose, digits=digits,
                          REMLf=con$REMLf, link=link, mZ=mZ, alpha.min=alpha.min, alpha.max=alpha.max, alpha.transf=TRUE,
                          tau2.min=con$tau2.min, tau2.max=con$tau2.max, optbeta=optbeta", ctrl.arg, ")\n")

      }

      if (link == "identity") {

         if (optimizer == "constrOptim")
            optcall <- paste0("constrOptim(theta=c(beta.init, alpha.init), f=.ll.rma.ls, grad=NULL, ui=Z, ci=rep(0,k),
                              yi=yi, vi=vi, X=X, Z=Z, reml=reml, k=k, pX=p, alpha.arg=alpha, beta.arg=beta, verbose=verbose, digits=digits,
                              REMLf=con$REMLf, link=link, mZ=mZ, alpha.min=alpha.min, alpha.max=alpha.max, alpha.transf=TRUE,
                              tau2.min=con$tau2.min, tau2.max=con$tau2.max, optbeta=optbeta", ctrl.arg, ")\n")

         if (optimizer == "Rsolnp::solnp")
            optcall <- paste0("Rsolnp::solnp(pars=c(beta.init, alpha.init), fun=.ll.rma.ls, ineqfun=.rma.ls.ineqfun.pos, ineqLB=rep(0,k), ineqUB=rep(Inf,k),
                              yi=yi, vi=vi, X=X, Z=Z, reml=reml, k=k, pX=p, alpha.arg=alpha, beta.arg=beta, verbose=verbose, digits=digits,
                              REMLf=con$REMLf, link=link, mZ=mZ, alpha.min=alpha.min, alpha.max=alpha.max, alpha.transf=TRUE,
                              tau2.min=con$tau2.min, tau2.max=con$tau2.max, optbeta=optbeta", ctrl.arg, ")\n")

         if (optimizer == "nloptr::nloptr")
            optcall <- paste0("nloptr::nloptr(x0=c(beta.init, alpha.init), eval_f=.ll.rma.ls, eval_g_ineq=.rma.ls.ineqfun.neg,
                              yi=yi, vi=vi, X=X, Z=Z, reml=reml, k=k, pX=p, alpha.arg=alpha, beta.arg=beta, verbose=verbose, digits=digits,
                              REMLf=con$REMLf, link=link, mZ=mZ, alpha.min=alpha.min, alpha.max=alpha.max, alpha.transf=TRUE,
                              tau2.min=con$tau2.min, tau2.max=con$tau2.max, optbeta=optbeta", ctrl.arg, ")\n")

         if (optimizer == "alabama::constrOptim.nl")
            optcall <- paste0("alabama::constrOptim.nl(par=c(beta.init, alpha.init), fn=.ll.rma.ls, hin=.rma.ls.ineqfun.pos,
                              yi=yi, vi=vi, X=X, Z=Z, reml=reml, k=k, pX=p, alpha.arg=alpha, beta.arg=beta, verbose=verbose, digits=digits,
                              REMLf=con$REMLf, link=link, mZ=mZ, alpha.min=alpha.min, alpha.max=alpha.max, alpha.transf=TRUE,
                              tau2.min=con$tau2.min, tau2.max=con$tau2.max, optbeta=optbeta", ctrl.arg, ")\n")

      }

      #print(optcall)
      #return(optcall)

      if (verbose) {
         opt.res <- try(eval(str2lang(optcall)), silent=!verbose)
      } else {
         opt.res <- try(suppressWarnings(eval(str2lang(optcall))), silent=!verbose)
      }

      if (isTRUE(ddd$retopt))
         return(opt.res)

      ### convergence checks (if verbose print optimParallel log, if verbose > 2 print opt.res, and unify opt.res$par)

      opt.res$par <- .chkconv(optimizer=optimizer, opt.res=opt.res, optcontrol=optcontrol, fun="rma", verbose=verbose)

      ### back-transform in case constraints were placed on alpha values

      if (optbeta) {
         opt.res$par[-seq_len(p)] <- mapply(.mapfun.alpha, opt.res$par[-seq_len(p)], alpha.min, alpha.max)
      } else {
         opt.res$par <- mapply(.mapfun.alpha, opt.res$par, alpha.min, alpha.max)
      }

      ### replace fixed alpha (and beta) values in opt.res$par

      if (optbeta) {
         opt.res$par[seq_len(p)]  <- ifelse(is.na(beta),  opt.res$par[seq_len(p)],  beta)
         opt.res$par[-seq_len(p)] <- ifelse(is.na(alpha), opt.res$par[-seq_len(p)], alpha)
      } else {
         opt.res$par <- ifelse(is.na(alpha), opt.res$par, alpha)
      }

      ### try to compute vcov matrix for scale parameter estimates

      H <- NA_real_

      if (optbeta) {
         va <- matrix(NA_real_, nrow=p+q, ncol=p+q)
         hest <- c(is.na(beta), is.na(alpha))
      } else {
         va <- matrix(NA_real_, nrow=q, ncol=q)
         hest <- is.na(alpha)
      }

      if (any(hest) && !isTRUE(ddd$skiphes)) {

         if (verbose > 1)
            message(mstyle$message("\nComputing Hessian ..."))

         if (con$hesspack == "numDeriv")
            H <- try(numDeriv::hessian(func=.ll.rma.ls, x=opt.res$par, method.args=con$hessianCtrl, yi=yi, vi=vi, X=X,
                                       Z=Z, reml=reml, k=k, pX=p, alpha.arg=alpha, beta.arg=beta, verbose=FALSE, digits=digits,
                                       REMLf=con$REMLf, link=link, mZ=mZ, alpha.min=alpha.min, alpha.max=alpha.max, alpha.transf=FALSE,
                                       tau2.min=con$tau2.min, tau2.max=con$tau2.max, optbeta=optbeta), silent=TRUE)
         if (con$hesspack == "pracma")
            H <- try(pracma::hessian(f=.ll.rma.ls, x0=opt.res$par, yi=yi, vi=vi, X=X,
                                     Z=Z, reml=reml, k=k, pX=p, alpha.arg=alpha, beta.arg=beta, verbose=FALSE, digits=digits,
                                     REMLf=con$REMLf, link=link, mZ=mZ, alpha.min=alpha.min, alpha.max=alpha.max, alpha.transf=FALSE,
                                     tau2.min=con$tau2.min, tau2.max=con$tau2.max, optbeta=optbeta), silent=TRUE)

         if (inherits(H, "try-error")) {

            warning(mstyle$warning("Error when trying to compute the Hessian."), call.=FALSE)

         } else {

            H.hest <- H[hest, hest, drop=FALSE]
            iH.hest <- try(suppressWarnings(chol2inv(chol(H.hest))), silent=TRUE)

            if (inherits(iH.hest, "try-error") || anyNA(iH.hest) || any(is.infinite(iH.hest))) {

               warning(mstyle$warning("Error when trying to invert the Hessian."), call.=FALSE)

            } else {

               va[hest, hest] <- iH.hest

            }

         }

      }

      if (optbeta) {
         vba <- va
         vb  <- va[seq_len(p),   seq_len(p), drop=FALSE]
         va  <- va[-seq_len(p), -seq_len(p), drop=FALSE]
      }

      ### get scale (and location) parameter estimates

      alpha.arg <- alpha
      beta.arg  <- beta

      if (optbeta) {
         beta  <- cbind(opt.res$par[seq_len(p)])
         alpha <- cbind(opt.res$par[-seq_len(p)])
      } else {
         alpha <- cbind(opt.res$par)
      }

      if (any(alpha <= alpha.min + 10*.Machine$double.eps^0.25) || any(alpha >= alpha.max - 10*.Machine$double.eps^0.25))
         warning(mstyle$warning("One or more 'alpha' estimates are (almost) equal to their lower or upper bound.\nTreat results with caution (or consider adjusting 'alpha.min' and/or 'alpha.max')."), call.=FALSE)

      ### scale back alpha and va when Z matrix was rescaled

      if (!is.null(mZ)) {
         alpha <- mZ %*% alpha
         va[!hest,] <- 0
         va[,!hest] <- 0
         va <- mZ %*% va %*% t(mZ)
         va[!hest,] <- NA_real_
         va[,!hest] <- NA_real_
         Z <- Zsave
      }

      ### set/check 'att' argument

      att <- .set.btt(att, q, Z.int.incl, colnames(Z))
      m.alpha <- length(att) # number of alphas to test (m = q if all alphas are tested)

      ### ddf calculation

      if (is.element(test, c("knha","adhoc","t"))) {
         ddf.alpha <- k-q
      } else {
         ddf.alpha <- NA_integer_
      }

      ### QS calculation

      QS <- try(as.vector(t(alpha)[att] %*% chol2inv(chol(va[att,att])) %*% alpha[att]), silent=TRUE)

      if (inherits(QS, "try-error"))
         QS <- NA_real_

      se.alpha <- sqrt(diag(va))

      rownames(alpha) <- rownames(va) <- colnames(va) <- colnames(Z)

      names(se.alpha) <- NULL
      zval.alpha  <- c(alpha/se.alpha)

      if (is.element(test, c("knha","adhoc","t"))) {
         QS         <- QS / m.alpha
         QSdf       <- c(m.alpha, ddf.alpha)
         QSp        <- if (QSdf[2] > 0) pf(QS, df1=QSdf[1], df2=QSdf[2], lower.tail=FALSE) else NA_real_
         pval.alpha <- if (ddf.alpha > 0) 2*pt(abs(zval.alpha), df=ddf.alpha, lower.tail=FALSE) else rep(NA_real_,q)
         crit.alpha <- if (ddf.alpha > 0) qt(level/2, df=ddf.alpha, lower.tail=FALSE) else NA_real_
      } else {
         QSdf       <- c(m.alpha, ddf.alpha)
         QSp        <- pchisq(QS, df=QSdf[1], lower.tail=FALSE)
         pval.alpha <- 2*pnorm(abs(zval.alpha), lower.tail=FALSE)
         crit.alpha <- qnorm(level/2, lower.tail=FALSE)
      }

      ci.lb.alpha <- c(alpha - crit.alpha * se.alpha)
      ci.ub.alpha <- c(alpha + crit.alpha * se.alpha)

      if (link == "log")
         tau2 <- exp(as.vector(Z %*% alpha))
      if (link == "identity")
         tau2 <- as.vector(Z %*% alpha)

   }

   ### equal/fixed/common-effects model (note: sets tau2 to zero even when tau2 value is specified)

   if (is.element(method[1], c("FE","EE","CE")))
      tau2 <- 0

   #########################################################################

   ###### model fitting, test statistics, and confidence intervals

   if (verbose > 1)
      message(mstyle$message("\nModel fitting ..."))

   wi <- 1/(vi + tau2)
   W  <- diag(wi, nrow=k, ncol=k)
   M  <- diag(vi + tau2, nrow=k, ncol=k)

   if (weighted) {

      #########################
      ### weighted analysis ###
      #########################

      ### fit model with weighted estimation

      if (is.null(weights) || is.element(test, c("knha","adhoc"))) {

         ### if no weights are specified, use default inverse variance weights, that is, 1/vi or 1/(vi + tau2)
         ### also, even with weights, if test="knha" or "adhoc", need to run this to get RSS.knha

         ### if any vi = 0 and tau^2 is estimated to be 0 (or is set to 0 for a FE model), then get Inf for wi

         if (any(is.infinite(wi)))
            stop(mstyle$stop("Division by zero when computing the inverse variance weights."))

         ### don't recompute beta and vb when optbeta=TRUE, since these are already estimated

         if (!optbeta) {
            if (tau2.inf) {
               beta <- cbind(coef(lm(yi ~ 0 + X)))
               vb <- diag(rep(Inf,p), nrow=p, ncol=p)
            } else {
               stXWX <- .invcalc(X=X, W=W, k=k)
               beta  <- stXWX %*% crossprod(X,W) %*% Y
               vb    <- stXWX
            }
         }
         RSS.f <- sum(wi*c(yi - X %*% beta)^2)
         #P     <- W - W %*% X %*% stXWX %*% crossprod(X,W)
         #RSS.f <- crossprod(Y,P) %*% Y
         RSS.knha <- RSS.f

      }

      if (!is.null(weights)) {

         ### if weights are specified, use them (note: RSS.f is recomputed if test="knha" or "adhoc")

         A     <- diag(weights, nrow=k, ncol=k)
         stXAX <- .invcalc(X=X, W=A, k=k)
         beta  <- stXAX %*% crossprod(X,A) %*% Y
         vb    <- stXAX %*% t(X) %*% A %*% M %*% A %*% X %*% stXAX
         RSS.f <- sum(wi*c(yi - X %*% beta)^2)
         #P     <- W - W %*% X %*% stXAX %*% t(X) %*% A - A %*% X %*% stXAX %*% t(X) %*% W + A %*% X %*% stXAX %*% t(X) %*% W %*% X %*% stXAX %*% t(X) %*% A
         #RSS.f <- crossprod(Y,P) %*% Y

      }

      #return(list(beta=beta, vb=vb, se=sqrt(diag(vb)), RSS.f=RSS.f))

      ### calculate scaling factor for Knapp & Hartung method
      ### note: catch cases where RSS.knha is extremely small, which is probably due to all yi being equal
      ### then set s2w to 0 (to avoid the strange looking output we would obtain if we don't do this)

      if (is.element(test, c("knha","adhoc"))) {

         if (RSS.knha <= .Machine$double.eps) {
            s2w <- 0
         } else {
            s2w <- RSS.knha / (k-p)
         }

      }

   } else {

      ###########################
      ### unweighted analysis ###
      ###########################

      ### fit model with unweighted estimation
      ### note: 1) if user has specified weights, they are ignored
      ###       2) but if method="GENQ/GENQM", they were used to estimate tau^2

      stXX  <- .invcalc(X=X, W=diag(k), k=k)
      beta  <- stXX %*% crossprod(X,Y)
      vb    <- tcrossprod(stXX,X) %*% M %*% X %*% stXX
      RSS.f <- sum(wi*(yi - X %*% beta)^2)
      #P     <- W - W %*% X %*% tcrossprod(stXX,X) - X %*% stXX %*% crossprod(X,W) + X %*% stXX %*% crossprod(X,W) %*% X %*% tcrossprod(stXX,X)
      #RSS.f <- crossprod(Y,P) %*% Y

      ### calculate scaling factor for Knapp & Hartung method

      if (is.element(test, c("knha","adhoc"))) {

         if (any(is.infinite(wi)))
            stop(mstyle$stop("Division by zero when computing the inverse variance weights."))

         stXWX     <- .invcalc(X=X, W=W, k=k)
         beta.knha <- stXWX %*% crossprod(X,W) %*% Y
         RSS.knha  <- sum(wi*(yi - X %*% beta.knha)^2)
         #P         <- W - W %*% X %*% stXWX %*% crossprod(X,W)
         #RSS.knha  <- c(crossprod(Y,P) %*% Y)

         if (RSS.knha <= .Machine$double.eps) {
            s2w <- 0
         } else {
            s2w <- RSS.knha / (k-p)
         }

      }

   }

   if (verbose > 1)
      message(mstyle$message("Conducting tests of the fixed effects ..."))

   ### the Knapp & Hartung method as described in the literature is for random/mixed-effects models

   if (is.element(method[1], c("FE","EE","CE")) && is.element(test, c("knha","adhoc")))
      warning(mstyle$warning(paste0("Knapp and Hartung method is not meant to be used in the context of ", method[1], " models.")), call.=FALSE)

   ### Knapp & Hartung method with ad-hoc correction so that the scale factor is always >= 1

   if (test == "adhoc")
      s2w[s2w < 1] <- 1

   ### for Knapp & Hartung method, apply scaling to vb

   vb <- s2w * vb

   ### handle special case of tau2=Inf

   if (tau2.inf)
      vb <- diag(rep(Inf,p), nrow=p, ncol=p)

   ### ddf calculation

   if (is.element(test, c("knha","adhoc","t"))) {
      ddf <- .chkddd(ddd$dfs, k-p, ddd$dfs[[1]]) # would be nice to allow multiple dfs values, but tricky since some methods are set up for a single df value
   } else {
      ddf <- NA_integer_
   }

   ### QM calculation

   QM <- try(as.vector(t(beta)[btt] %*% chol2inv(chol(vb[btt,btt])) %*% beta[btt]), silent=TRUE)

   if (inherits(QM, "try-error"))
      QM <- NA_real_

   ### abbreviate some types of coefficient names

   if (.isTRUE(ddd$abbrev)) {
      tmp <- colnames(X)
      tmp <- gsub("relevel(factor(", "", tmp, fixed=TRUE)
      tmp <- gsub("\\), ref = \"[[:alnum:]]*\")", "", tmp)
      tmp <- gsub("poly(", "", tmp, fixed=TRUE)
      tmp <- gsub(", degree = [[:digit:]], raw = TRUE)", "^", tmp)
      tmp <- gsub(", degree = [[:digit:]], raw = T)", "^", tmp)
      tmp <- gsub(", degree = [[:digit:]])", "^", tmp)
      tmp <- gsub("rcs\\([[:alnum:]]*, [[:digit:]]\\)", "", tmp)
      tmp <- gsub("factor(", "", tmp, fixed=TRUE)
      tmp <- gsub("I(", "", tmp, fixed=TRUE)
      tmp <- gsub(")", "", tmp, fixed=TRUE)
      colnames(X) <- tmp
   }

   rownames(beta) <- rownames(vb) <- colnames(vb) <- colnames(X.f) <- colnames(X)

   se <- sqrt(diag(vb))
   names(se) <- NULL
   zval <- c(beta/se)

   if (is.element(test, c("knha","adhoc","t"))) {
      QM   <- QM / m
      QMdf <- c(m, ddf)
      QMp  <- if (QMdf[2] > 0) pf(QM, df1=QMdf[1], df2=QMdf[2], lower.tail=FALSE) else NA_real_
      pval <- if (ddf > 0) 2*pt(abs(zval), df=ddf, lower.tail=FALSE) else rep(NA_real_,p)
      crit <- if (ddf > 0) qt(level/2, df=ddf, lower.tail=FALSE) else NA_real_
   } else {
      QMdf <- c(m, ddf)
      QMp  <- pchisq(QM, df=QMdf[1], lower.tail=FALSE)
      pval <- 2*pnorm(abs(zval), lower.tail=FALSE)
      crit <- qnorm(level/2, lower.tail=FALSE)
   }

   ci.lb <- c(beta - crit * se)
   ci.ub <- c(beta + crit * se)

   #########################################################################

   ### heterogeneity test (Wald-type test of the extra coefficients in the saturated model)

   if (verbose > 1)
      message(mstyle$message("Conducting heterogeneity test ..."))

   if (allvipos) {

      ### heterogeneity test (always uses inverse variance method)
      # note: this is unaffected by the 'weighted' argument, since under H0, the same parameters are
      # estimated and weighted estimation provides the most efficient estimates; therefore, also any
      # arbitrary weights specified by the user are not relevant here (different from what the metan
      # command in Stata does!) see also: Chen, Z., Ng, H. K. T., & Nadarajah, S. (2014). A note on
      # Cochran test for homogeneity in one-way ANOVA and meta-analysis. Statistical Papers, 55(2),
      # 301-310. This shows that the weights used are not relevant.

      if (k > p) {

         wi    <- 1/vi
         W.FE  <- diag(wi, nrow=k, ncol=k) # note: ll.REML below involves W, so cannot overwrite W
         stXWX <- .invcalc(X=X, W=W.FE, k=k)
         P     <- W.FE - W.FE %*% X %*% stXWX %*% crossprod(X,W.FE) # need P below for calculation of I^2
         QE    <- max(0, c(crossprod(Ymc,P) %*% Ymc))
         #beta.FE <- stXWX %*% crossprod(X,W.FE) %*% Y
         #QE    <- max(0, sum(wi*(yi - X %*% beta.FE)^2))
         QEp   <- pchisq(QE, df=k-p, lower.tail=FALSE)

         ### calculation of 'typical' sampling variance

         #vt <- (k-1) / (sum(wi) - sum(wi^2)/sum(wi)) # this only applies to the RE model

         if (i2def == "1")
            vt <- (k-p) / .tr(P)
         if (i2def == "2")
            vt <- 1 / mean(wi) # harmonic mean of the vi values (see Takkouche et al., 1999)

         ### calculation of I^2 and H^2

         if (is.element(method[1], c("FE","EE","CE"))) {
            I2 <- max(0, 100 * (QE - (k-p)) / QE)
            H2 <- QE / (k-p)
         } else {
            I2 <- 100 * tau2 / (vt + tau2) # vector for location-scale models
            H2 <- tau2 / vt + 1            # vector for location-scale models
         }

      } else {

         QE  <- 0
         QEp <- 1
         I2  <- 0
         H2  <- 1
         vt  <- 0

      }

   } else {

      if (!vi0)
         warning(mstyle$warning(paste0("Cannot compute ", ifelse(int.only, "Q", "QE"), "-test, I^2, or H^2 when there are non-positive sampling variances in the data.")), call.=FALSE)

      vt <- NA_real_

   }

   #########################################################################

   ###### fit statistics

   if (verbose > 1)
      message(mstyle$message("Computing fit statistics and log-likelihood ..."))

   ### note: tau2 is not counted as a parameter when it was fixed by the user (same for fixed alpha values)
   q.est <- ifelse(model == "rma.uni", 0, sum(is.na(alpha.arg)))
   parms <- ifelse(optbeta, sum(is.na(beta.arg)), p) + ifelse(model == "rma.uni", ifelse(is.element(method[1], c("FE","EE","CE")) || tau2.fix, 0, 1), q.est)

   ll.ML    <- -1/2 * (k) * log(2*base::pi) - 1/2 * sum(log(vi + tau2)) - 1/2 * RSS.f
   ll.REML  <- -1/2 * (k-p) * log(2*base::pi) + ifelse(con$REMLf, 1/2 * determinant(crossprod(X), logarithm=TRUE)$modulus, 0) +
               -1/2 * sum(log(vi + tau2)) - 1/2 * determinant(crossprod(X,W) %*% X, logarithm=TRUE)$modulus - 1/2 * RSS.f

   if (k > p) {
      if (allvipos) {
         dev.ML <- -2 * (ll.ML - sum(dnorm(yi, mean=yi, sd=sqrt(vi), log=TRUE)))
      } else {
         dev.ML <- -2 * ll.ML
      }
   } else {
      dev.ML <- 0
   }

   AIC.ML    <- -2 * ll.ML   + 2*parms
   BIC.ML    <- -2 * ll.ML   +   parms * log(k)
   AICc.ML   <- -2 * ll.ML   + 2*parms * max(k, parms+2) / (max(k, parms+2) - parms - 1)
   dev.REML  <- -2 * (ll.REML - 0) # saturated model has ll = 0 when using the full REML likelihood
   AIC.REML  <- -2 * ll.REML + 2*parms
   BIC.REML  <- -2 * ll.REML +   parms * log(k-p)
   AICc.REML <- -2 * ll.REML + 2*parms * max(k-p, parms+2) / (max(k-p, parms+2) - parms - 1)

   fit.stats <- matrix(c(ll.ML, dev.ML, AIC.ML, BIC.ML, AICc.ML, ll.REML, dev.REML, AIC.REML, BIC.REML, AICc.REML), ncol=2, byrow=FALSE)
   dimnames(fit.stats) <- list(c("ll","dev","AIC","BIC","AICc"), c("ML","REML"))
   fit.stats <- data.frame(fit.stats)

   #########################################################################

   ### compute pseudo R^2 statistic for mixed-effects models with an intercept (only for rma.uni normal models)

   if (!int.only && int.incl && model == "rma.uni" && !isTRUE(ddd$skipr2)) {

      if (verbose > 1)
         message(mstyle$message("Computing R^2 ..."))

      if (is.element(method[1], c("FE","EE","CE"))) {

         if (identical(var(yi),0)) {
            R2 <- 0
         } else {
            if (weighted) {
               if (is.null(weights)) {
                  R2 <- max(0, 100 * summary(lm(yi ~ X, weights=wi))$adj.r.squared)
               } else {
                  R2 <- max(0, 100 * summary(lm(yi ~ X, weights=weights))$adj.r.squared)
               }
            } else {
               R2 <- max(0, 100 * summary(lm(yi ~ X))$adj.r.squared)
            }
         }

      } else {

         if (r2def %in% c("1","1v","3","3v","5","6","7","8")) {

            args <- list(yi=yi, vi=vi, weights=weights, method=method, weighted=weighted, test=test,
                         verbose=ifelse(verbose, TRUE, FALSE), control=con, digits=digits, outlist="minimal")

            if (verbose > 1) {
               res0 <- try(.do.call(rma.uni, args), silent=FALSE)
            } else {
               res0 <- try(suppressWarnings(.do.call(rma.uni, args)), silent=TRUE)
            }

            if (!inherits(res0, "try-error")) {

               tau2.RE <- res0$tau2

               if (identical(tau2.RE,0) && r2def %in% c("1","3")) {

                  R2 <- 0

               } else {

                  ll0 <- logLik(res0)
                  ll1 <- ifelse(method[1] == "ML", ll.ML, ll.REML)
                  lls <- (ifelse(method[1] == "ML", dev.ML, dev.REML) + 2*ll1) / 2

                  # based on Raudenbush (1994)
                  if (r2def == "1")
                     R2 <- (tau2.RE - tau2) / tau2.RE

                  # like Raudenbush (1994) but with total variance (including sampling variance) in the denominator
                  if (r2def == "1v")
                     R2 <- (tau2.RE - tau2) / (tau2.RE + 1/mean(1/vi))

                  # model component definition with tau^2_RE in the denominator
                  if (r2def == "3")
                     R2 <- var(c(X%*%beta)) / tau2.RE

                  # model component definition with total variance (including sampling variance) in the denominator
                  if (r2def == "3v")
                     R2 <- var(c(X%*%beta)) / (tau2.RE + 1/mean(1/vi))

                  # like McFadden's R^2
                  if (r2def == "5")
                     R2 <- 1 - ll1 / ll0

                  # like Cox & Snell R^2
                  if (r2def == "6")
                     R2 <- 1 - (exp(ll0) / exp(ll1))^(2/k)

                  # like Nagelkerke R^2
                  if (r2def == "7")
                     R2 <- (1 - (exp(ll0) / exp(ll1))^(2/k)) / (1 - exp(ll0)^(2/k))

                  # how close ME model is to the saturated model in terms of ll (same as 5 for REML)
                  if (r2def == "8")
                     R2 <- (ll1 - ll0) / (lls - ll0)

               }

            } else {

               R2 <- NA_real_

            }

         } else {

            # model component definition
            if (r2def == "2")
               R2 <- var(c(X%*%beta)) / (var(c(X%*%beta)) + tau2)

            # model component definition with total variance (including sampling variance) in the denominator
            if (r2def == "2v")
               R2 <- var(c(X%*%beta)) / (var(c(X%*%beta)) + tau2 + 1/mean(1/vi))

            # squared correlation between observed and fitted values
            if (r2def == "4")
               R2 <- cor(yi, c(X%*%beta))^2

            # squared weighted correlation between observed and fitted values
            if (r2def == "4w") {
               if (is.null(weights)) {
                  # identical to eta^2 = F * df1 / (F * df1 + df2) when test="knha"
                  R2 <- cov.wt(cbind(yi, c(X%*%beta)), cor=TRUE, wt=1/(vi+tau2))$cor[1,2]^2
               } else {
                  R2 <- cov.wt(cbind(yi, c(X%*%beta)), cor=TRUE, wt=weights)$cor[1,2]^2
               }
            }

         }

         R2 <- max(0, 100 * R2)

      }

   } else {

      R2 <- NULL

   }

   if (.isTRUE(ddd$pleasedonotreportI2thankyouverymuch)) {
      I2 <- NA
      H2 <- NA
   }

   #########################################################################

   ###### prepare output

   if (verbose > 1)
      message(mstyle$message("Preparing output ..."))

   p.eff <- p
   k.eff <- k

   if (is.null(ddd$outlist) || ddd$outlist == "nodata") {

      res <- list(b=beta, beta=beta, se=se, zval=zval, pval=pval, ci.lb=ci.lb, ci.ub=ci.ub, vb=vb,
                  tau2=tau2, se.tau2=se.tau2, tau2.fix=tau2.fix, tau2.f=tau2,
                  I2=I2, H2=H2, R2=R2, vt=vt,
                  QE=QE, QEp=QEp, QM=QM, QMdf=QMdf, QMp=QMp,
                  k=k, k.f=k.f, k.eff=k.eff, k.all=k.all, p=p, p.eff=p.eff, parms=parms,
                  int.only=int.only, int.incl=int.incl, intercept=intercept, allvipos=allvipos, coef.na=coef.na,
                  yi=yi, vi=vi, X=X, weights=weights, yi.f=yi.f, vi.f=vi.f, X.f=X.f, weights.f=weights.f, M=M,
                  outdat.f=outdat.f, ni=ni, ni.f=ni.f,
                  ids=ids, not.na=not.na, subset=subset, slab=slab, slab.null=slab.null,
                  measure=measure, method=method[1], model=model, weighted=weighted,
                  test=test, dfs=ddf, ddf=ddf, s2w=s2w, btt=btt, m=m,
                  digits=digits, level=level, control=control, verbose=verbose,
                  add=add, to=to, drop00=drop00,
                  fit.stats=fit.stats,
                  formula.yi=formula.yi, formula.mods=formula.mods,
                  version=packageVersion("metafor"), call=mf)

      if (is.null(ddd$outlist))
         res <- append(res, list(data=data), which(names(res) == "fit.stats"))

   } else {

      if (ddd$outlist == "minimal") {
         res <- list(b=beta, beta=beta, se=se, zval=zval, pval=pval, ci.lb=ci.lb, ci.ub=ci.ub, vb=vb,
                     tau2=tau2, se.tau2=se.tau2, tau2.fix=tau2.fix,
                     I2=I2, H2=H2, R2=R2,
                     QE=QE, QEp=QEp, QM=QM, QMdf=QMdf, QMp=QMp,
                     k=k, k.eff=k.eff, p=p, p.eff=p.eff, parms=parms,
                     int.only=int.only,
                     measure=measure, method=method[1], model=model,
                     test=test, dfs=ddf, ddf=ddf, btt=btt, m=m,
                     digits=digits,
                     fit.stats=fit.stats)
      } else {
         res <- eval(str2lang(paste0("list(", ddd$outlist, ")")))
      }

   }

   if (model == "rma.ls") {

      res$alpha          <- alpha
      res$va             <- va
      res$se.alpha       <- se.alpha
      res$zval.alpha     <- zval.alpha
      res$pval.alpha     <- pval.alpha
      res$ci.lb.alpha    <- ci.lb.alpha
      res$ci.ub.alpha    <- ci.ub.alpha
      res$alpha.fix      <- !is.na(alpha.arg)
      res$optbeta        <- optbeta
      if (optbeta) {
         res$vba         <- vba
         res$beta.fix    <- !is.na(beta.arg)
      }
      res$q              <- q
      res$alphas         <- q
      res$link           <- link
      res$Z              <- Z
      res$Z.f            <- Z.f
      res$tau2.f         <- rep(NA_real_, k.f)
      res$tau2.f[not.na] <- tau2
      res$att            <- att
      res$m.alpha        <- m.alpha
      res$ddf.alpha      <- ddf.alpha
      res$QS             <- QS
      res$QSdf           <- QSdf
      res$QSp            <- QSp
      res$formula.scale  <- formula.scale
      res$Z.int.incl     <- Z.int.incl
      res$Z.intercept    <- Z.int.incl
      res$Z.int.only     <- Z.int.only
      res$coef.na.Z      <- coef.na.Z
      res$H              <- H

   }

   time.end <- proc.time()
   res$time <- unname(time.end - time.start)[3]

   if (.isTRUE(ddd$time))
      .print.time(res$time)

   if (verbose || .isTRUE(ddd$time))
      cat("\n")

   if (model == "rma.ls") {
      class(res) <- c("rma.ls", "rma.uni", "rma")
   } else {
      class(res) <- c("rma.uni", "rma")
   }

   return(res)

}
