.cmicalc <-
function (mi) 
{
    cmi <- gamma(mi/2)/(sqrt(mi/2) * gamma((mi - 1)/2))
    is.na <- is.na(cmi)
    cmi[is.na] <- 1 - 3/(4 * mi[is.na] - 1)
    return(cmi)
}
.diag <-
function (x) 
{
    n <- length(x)
    y <- array(0, c(n, n))
    if ((m <- n) > 0L) 
        y[1L + 0L:(m - 1L) * (n + 1L)] <- x
    return(y)
}
.dnchg <-
function (parms, ai, bi, ci, di, X.fit, random, verbose = FALSE, 
    digits = 4, dnchgcalc, dnchgprec, control.int) 
{
    p <- ncol(X.fit)
    k <- length(ai)
    b <- parms[seq.int(p)]
    tau2 <- ifelse(random, exp(parms[p + 1]), 0)
    mu.i <- X.fit %*% cbind(b)
    lli <- rep(NA, k)
    if (!random) {
        for (i in seq.int(k)) {
            lli[i] <- log(.dnchgi(logOR = mu.i[i], ai = ai[i], 
                bi = bi[i], ci = ci[i], di = di[i], random = random, 
                dnchgcalc = dnchgcalc, dnchgprec = dnchgprec))
        }
        if (verbose) 
            cat("ll =", formatC(sum(lli), digits = digits, format = "f"), 
                " ", formatC(b, digits = digits, format = "f"), 
                "\n")
    }
    if (random) {
        for (i in seq.int(k)) {
            res <- try(integrate(.dnchgi, lower = control.int$lower, 
                upper = control.int$upper, ai = ai[i], bi = bi[i], 
                ci = ci[i], di = di[i], mu.i = mu.i[i], tau2 = tau2, 
                random = random, dnchgcalc = dnchgcalc, dnchgprec = dnchgprec, 
                rel.tol = control.int$rel.tol, subdivisions = control.int$subdivisions, 
                stop.on.error = FALSE), silent = !verbose)
            if (inherits(res, "try-error")) {
                stop(paste("Could not integrate over density of non-central hypergeometric distribution in study ", 
                  i, ".", sep = ""))
            }
            else {
                if (res$value > 0) {
                  lli[i] <- log(res$value)
                }
                else {
                  lli[i] <- -Inf
                }
            }
        }
        if (verbose) 
            cat("ll = ", formatC(sum(lli), digits = digits, format = "f"), 
                " ", formatC(tau2, digits = digits, format = "f"), 
                " ", formatC(b, digits = digits, format = "f"), 
                "\n")
    }
    return(-sum(lli))
}
.dnchgi <-
function (logOR, ai, bi, ci, di, mu.i, tau2, random, dnchgcalc, 
    dnchgprec) 
{
    k <- length(logOR)
    dnchgi <- rep(NA, k)
    logOR[logOR < log(1e-12)] <- log(1e-12)
    logOR[logOR > log(1e+12)] <- log(1e+12)
    for (i in seq.int(k)) {
        ORi <- exp(logOR[i])
        if (dnchgcalc == "dnoncenhypergeom") {
            res <- try(.dnoncenhypergeom(x = ai, n1 = ai + bi, 
                n2 = ci + di, m1 = ai + ci, psi = ORi))
        }
        else {
            res <- try(dFNCHypergeo(x = ai, m1 = ai + bi, m2 = ci + 
                di, n = ai + ci, odds = ORi, precision = dnchgprec))
        }
        if (inherits(res, "try-error")) {
            stop(paste("Could not compute density of non-central hypergeometric distribution in study ", 
                i, ".", sep = ""))
        }
        else {
            dnchgi[i] <- res
        }
    }
    if (random) 
        dnchgi <- dnchgi * dnorm(logOR, mu.i, sqrt(tau2))
    return(dnchgi)
}
.dnoncenhypergeom <-
function (x = NA, n1, n2, m1, psi) 
{
    mode.compute <- function(n1, n2, m1, psi, ll, uu) {
        a <- psi - 1
        b <- -((m1 + n1 + 2) * psi + n2 - m1)
        c <- psi * (n1 + 1) * (m1 + 1)
        q <- b + sign(b) * sqrt(b * b - 4 * a * c)
        q <- -q/2
        mode <- trunc(c/q)
        if (uu >= mode && mode >= ll) 
            return(mode)
        else return(trunc(q/a))
    }
    r.function <- function(n1, n2, m1, psi, i) {
        (n1 - i + 1) * (m1 - i + 1)/i/(n2 - m1 + i) * psi
    }
    ll <- max(0, m1 - n2)
    uu <- min(n1, m1)
    if (n1 < 0 | n2 < 0) 
        stop("n1 or n2 negative in dnoncenhypergeom().\n")
    if (m1 < 0 | m1 > (n1 + n2)) 
        stop("m1 out of range in dnoncenhypergeom().\n")
    if (psi <= 0) 
        stop("psi [odds ratio] negative in dnoncenhypergeom().\n")
    if (!is.na(x) & (x < ll | x > uu)) 
        stop("x out of bounds in dnoncenhypergeom().\n")
    if (!is.na(x) & length(x) > 1) 
        stop("x neither missing or scalar in dnoncenhypergeom().\n")
    mode <- mode.compute(n1, n2, m1, psi, ll, uu)
    pi <- array(1, uu - ll + 1)
    shift <- 1 - ll
    if (mode < uu) {
        r1 <- r.function(n1, n2, m1, psi, (mode + 1):uu)
        pi[(mode + 1 + shift):(uu + shift)] <- cumprod(r1)
    }
    if (mode > ll) {
        r1 <- 1/r.function(n1, n2, m1, psi, mode:(ll + 1))
        pi[(mode - 1 + shift):(ll + shift)] <- cumprod(r1)
    }
    pi <- pi/sum(pi)
    if (is.na(x)) {
        return(cbind(ll:uu, pi))
    }
    else {
        return(pi[x + shift])
    }
}
.genperms <-
function (k) 
{
    v <- seq.int(k)
    sub <- function(k, v) {
        if (k == 1L) {
            matrix(v, 1, k)
        }
        else {
            X <- NULL
            for (i in seq.int(k)) {
                X <- rbind(X, cbind(v[i], Recall(k - 1, v[-i])))
            }
            X
        }
    }
    return(sub(k, v[seq.int(k)]))
}
.gensigns <-
function (k) 
{
    ncols <- k
    nrows <- 2^k
    out <- matrix(NA, nrow = nrows, ncol = ncols)
    for (i in seq.int(ncols)) {
        out[, i] <- rep(c(1, -1), times = (nrows/2)/((2^i)/2), 
            each = (2^i)/2)
    }
    return(out)
}
.genuperms <-
function (x) 
{
    z <- NULL
    sub <- function(x, y) {
        len.x <- length(x)
        if (len.x == 0L) {
            return(y)
        }
        else {
            prev.num <- 0
            for (i in seq.int(len.x)) {
                num <- x[i]
                if (num > prev.num) {
                  prev.num <- num
                  z <- rbind(z, Recall(x[-i], c(y, num)))
                }
            }
            return(z)
        }
    }
    return(sub(x, y = NULL))
}
.invcalc <-
function (X, W, k) 
{
    wX <- sqrt(W) %*% X
    res.qrs <- qr.solve(wX, diag(k))
    return(tcrossprod(res.qrs))
}
.onAttach <-
function (libname, pkgname) 
{
    loadmsg <- "\nLoading 'metafor' package (version 1.8-0). For an overview \nand introduction to the package please type: help(metafor)."
    packageStartupMessage(loadmsg, domain = NULL, appendLF = TRUE)
}
.QE.func <-
function (tau2val, Y, vi, X, k, objective, verbose = FALSE, digits = 4) 
{
    wi <- 1/(vi + tau2val)
    W <- .diag(wi)
    stXWX <- .invcalc(X = X, W = W, k = k)
    P <- W - W %*% X %*% stXWX %*% crossprod(X, W)
    RSS <- crossprod(Y, P) %*% Y
    if (verbose) 
        cat("tau2 =", formatC(tau2val, digits = digits, format = "f"), 
            "\tRSS - objective =", c(RSS - objective), "\n")
    return(RSS - objective)
}
.setxlab <-
function (measure, transf.char, atransf.char, gentype) 
{
    if (gentype == 1) 
        xlab <- "Observed Outcome"
    if (gentype == 2) 
        xlab <- "Overall Estimate"
    if (measure == "RR") {
        if (transf.char == "FALSE" && atransf.char == "FALSE") {
            xlab <- "Log Relative Risk"
        }
        else {
            xlab <- "Transformed Log Relative Risk"
            if (atransf.char == "exp" || atransf.char == "transf.exp.int") 
                xlab <- "Relative Risk (log scale)"
            if (transf.char == "exp" || transf.char == "transf.exp.int") 
                xlab <- "Relative Risk"
        }
    }
    if (measure == "OR" || measure == "PETO" || measure == "D2OR") {
        if (transf.char == "FALSE" && atransf.char == "FALSE") {
            xlab <- "Log Odds Ratio"
        }
        else {
            xlab <- "Transformed Log Odds Ratio"
            if (atransf.char == "exp" || atransf.char == "transf.exp.int") 
                xlab <- "Odds Ratio (log scale)"
            if (transf.char == "exp" || transf.char == "transf.exp.int") 
                xlab <- "Odds Ratio"
        }
    }
    if (measure == "RD") {
        if (transf.char == "FALSE" && atransf.char == "FALSE") {
            xlab <- "Risk Difference"
        }
        else {
            xlab <- "Transformed Risk Difference"
        }
    }
    if (measure == "AS") {
        if (transf.char == "FALSE" && atransf.char == "FALSE") {
            xlab <- "Arcsine Transformed Risk Difference"
        }
        else {
            xlab <- "Transformed Arcsine Transformed Risk Difference"
        }
    }
    if (measure == "PHI") {
        if (transf.char == "FALSE" && atransf.char == "FALSE") {
            xlab <- "Phi Coefficient"
        }
        else {
            xlab <- "Transformed Phi Coefficient"
        }
    }
    if (measure == "YUQ") {
        if (transf.char == "FALSE" && atransf.char == "FALSE") {
            xlab <- "Yule's Q"
        }
        else {
            xlab <- "Transformed Yule's Q"
        }
    }
    if (measure == "YUY") {
        if (transf.char == "FALSE" && atransf.char == "FALSE") {
            xlab <- "Yule's Y"
        }
        else {
            xlab <- "Transformed Yule's Y"
        }
    }
    if (measure == "IRR") {
        if (transf.char == "FALSE" && atransf.char == "FALSE") {
            xlab <- "Log Incidence Rate Ratio"
        }
        else {
            xlab <- "Transformed Log Incidence Relative Risk"
            if (atransf.char == "exp" || atransf.char == "transf.exp.int") 
                xlab <- "Incidence Rate Ratio (log scale)"
            if (transf.char == "exp" || transf.char == "transf.exp.int") 
                xlab <- "Incidence Rate Ratio"
        }
    }
    if (measure == "IRD") {
        if (transf.char == "FALSE" && atransf.char == "FALSE") {
            xlab <- "Incidence Rate Difference"
        }
        else {
            xlab <- "Transformed Incidence Rate Difference"
        }
    }
    if (measure == "IRSD") {
        if (transf.char == "FALSE" && atransf.char == "FALSE") {
            xlab <- "Square-Root Transformed Incidence Rate Difference"
        }
        else {
            xlab <- "Transformed Square-Root Transformed Incidence Rate Difference"
        }
    }
    if (measure == "MD") {
        if (transf.char == "FALSE" && atransf.char == "FALSE") {
            xlab <- "Mean Difference"
        }
        else {
            xlab <- "Transformed Mean Difference"
        }
    }
    if (measure == "SMD" || measure == "SMDH" || measure == "PBIT" || 
        measure == "OR2D") {
        if (transf.char == "FALSE" && atransf.char == "FALSE") {
            xlab <- "Standardized Mean Difference"
        }
        else {
            xlab <- "Transformed Standardized Mean Difference"
        }
    }
    if (measure == "ROM") {
        if (transf.char == "FALSE" && atransf.char == "FALSE") {
            xlab <- "Log Ratio of Means"
        }
        else {
            xlab <- "Transformed Log Ratio of Means"
            if (atransf.char == "exp" || atransf.char == "transf.exp.int") 
                xlab <- "Ratio of Means (log scale)"
            if (transf.char == "exp" || transf.char == "transf.exp.int") 
                xlab <- "Ratio of Means"
        }
    }
    if (measure == "RPB") {
        if (transf.char == "FALSE" && atransf.char == "FALSE") {
            xlab <- "Point-Biserial Correlation"
        }
        else {
            xlab <- "Transformed Point-Biserial Correlation"
        }
    }
    if (measure == "COR" || measure == "UCOR" || measure == "RTET" || 
        measure == "RBIS") {
        if (transf.char == "FALSE" && atransf.char == "FALSE") {
            xlab <- "Correlation Coefficient"
        }
        else {
            xlab <- "Transformed Correlation Coefficient"
        }
    }
    if (measure == "ZCOR") {
        if (transf.char == "FALSE" && atransf.char == "FALSE") {
            xlab <- "Fisher's z Transformed Correlation Coefficient"
        }
        else {
            xlab <- "Transformed Fisher's z Transformed Correlation Coefficient"
            if (atransf.char == "transf.ztor" || atransf.char == 
                "transf.ztor.int") 
                xlab <- "Correlation Coefficient"
            if (transf.char == "transf.ztor" || transf.char == 
                "transf.ztor.int") 
                xlab <- "Correlation Coefficient"
        }
    }
    if (measure == "PR") {
        if (transf.char == "FALSE" && atransf.char == "FALSE") {
            xlab <- "Proportion"
        }
        else {
            xlab <- "Transformed Proportion"
        }
    }
    if (measure == "PLN") {
        if (transf.char == "FALSE" && atransf.char == "FALSE") {
            xlab <- "Log Proportion"
        }
        else {
            xlab <- "Transformed Log Proportion"
            if (atransf.char == "exp" || atransf.char == "transf.exp.int") 
                xlab <- "Proportion (log scale)"
            if (transf.char == "exp" || transf.char == "transf.exp.int") 
                xlab <- "Proportion"
        }
    }
    if (measure == "PLO") {
        if (transf.char == "FALSE" && atransf.char == "FALSE") {
            xlab <- "Log Odds"
        }
        else {
            xlab <- "Transformed Log Odds"
            if (atransf.char == "exp" || atransf.char == "transf.exp.int") 
                xlab <- "Proportion (logit scale)"
            if (transf.char == "exp" || transf.char == "transf.exp.int") 
                xlab <- "Proportion"
        }
    }
    if (measure == "PAS") {
        if (transf.char == "FALSE" && atransf.char == "FALSE") {
            xlab <- "Arcsine Transformed Proportion"
        }
        else {
            xlab <- "Transformed Arcsine Transformed Proportion"
            if (atransf.char == "transf.iarcsin" || atransf.char == 
                "transf.iarcsin.int") 
                xlab <- "Proportion (arcsine scale)"
            if (transf.char == "transf.iarcsin" || transf.char == 
                "transf.iarcsin.int") 
                xlab <- "Proportion"
        }
    }
    if (measure == "PFT") {
        if (transf.char == "FALSE" && atransf.char == "FALSE") {
            xlab <- "Double Arcsine Transformed Proportion"
        }
        else {
            xlab <- "Transformed Double Arcsine Transformed Proportion"
            if (atransf.char == "transf.ift.hm") 
                xlab <- "Proportion"
            if (transf.char == "transf.ift.hm") 
                xlab <- "Proportion"
        }
    }
    if (measure == "IR") {
        if (transf.char == "FALSE" && atransf.char == "FALSE") {
            xlab <- "Incidence Rate"
        }
        else {
            xlab <- "Transformed Incidence Rate"
        }
    }
    if (measure == "IRLN") {
        if (transf.char == "FALSE" && atransf.char == "FALSE") {
            xlab <- "Log Incidence Rate"
        }
        else {
            xlab <- "Transformed Log Incidence Rate"
            if (atransf.char == "exp" || atransf.char == "transf.exp.int") 
                xlab <- "Incidence Rate (log scale)"
            if (transf.char == "exp" || transf.char == "transf.exp.int") 
                xlab <- "Incidence Rate"
        }
    }
    if (measure == "IRS") {
        if (transf.char == "FALSE" && atransf.char == "FALSE") {
            xlab <- "Square-Root Transformed Incidence Rate"
        }
        else {
            xlab <- "Transformed Square-Root Transformed Incidence Rate"
            if (atransf.char == "transf.isqrt" || atransf.char == 
                "transf.isqrt.int") 
                xlab <- "Incidence Rate (square-root scale)"
            if (transf.char == "transf.isqrt" || transf.char == 
                "transf.isqrt.int") 
                xlab <- "Incidence Rate"
        }
    }
    if (measure == "IRFT") {
        if (transf.char == "FALSE" && atransf.char == "FALSE") {
            xlab <- "Freeman-Tukey Transformed Incidence Rate"
        }
        else {
            xlab <- "Transformed Freeman-Tukey Transformed Incidence Rate"
        }
    }
    if (measure == "MN") {
        if (transf.char == "FALSE" && atransf.char == "FALSE") {
            xlab <- "Mean"
        }
        else {
            xlab <- "Transformed Mean"
        }
    }
    if (measure == "MC") {
        if (transf.char == "FALSE" && atransf.char == "FALSE") {
            xlab <- "Mean Change"
        }
        else {
            xlab <- "Transformed Mean Change"
        }
    }
    if (measure == "SMCC" || measure == "SMCR") {
        if (transf.char == "FALSE" && atransf.char == "FALSE") {
            xlab <- "Standardized Mean Change"
        }
        else {
            xlab <- "Transformed Standardized Mean Change"
        }
    }
    if (measure == "ARAW") {
        if (transf.char == "FALSE" && atransf.char == "FALSE") {
            xlab <- "Coefficient alpha"
        }
        else {
            xlab <- "Transformed Coefficient alpha"
        }
    }
    if (measure == "AHW") {
        if (transf.char == "FALSE" && atransf.char == "FALSE") {
            xlab <- "Transformed Coefficient alpha"
        }
        else {
            xlab <- "Transformed Coefficient alpha"
            if (atransf.char == "transf.iahw") 
                xlab <- "Coefficient alpha"
            if (transf.char == "transf.iahw") 
                xlab <- "Coefficient alpha"
        }
    }
    if (measure == "ABT") {
        if (transf.char == "FALSE" && atransf.char == "FALSE") {
            xlab <- "Transformed Coefficient alpha"
        }
        else {
            xlab <- "Transformed Coefficient alpha"
            if (atransf.char == "transf.iabt") 
                xlab <- "Coefficient alpha"
            if (transf.char == "transf.iabt") 
                xlab <- "Coefficient alpha"
        }
    }
    return(xlab)
}
.tr <-
function (X) 
return(sum(diag(X)))
`[.escalc` <-
function (x, ...) 
{
    dat <- NextMethod("[")
    yi.name <- attr(x, "var.names")[1]
    if (!is.null(yi.name) && is.element(yi.name, names(dat))) 
        eval(parse(text = paste("attr(dat$", yi.name, ", 'measure') <- attr(x$", 
            yi.name, ", 'measure')", sep = "")))
    if (!is.null(attr(x, "var.names"))) 
        attr(dat, "var.names") <- attr(x, "var.names")
    if (!is.null(attr(x, "digits"))) 
        attr(dat, "digits") <- attr(x, "digits")
    return(dat)
}
addpoly <-
function (x, ...) 
UseMethod("addpoly")
addpoly.default <-
function (x, vi, sei, ci.lb, ci.ub, rows = -1, level = 95, digits = 2, 
    annotate = TRUE, mlab, transf = FALSE, atransf = FALSE, targs, 
    col = "black", efac = 1, cex, ...) 
{
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    if (missing(targs)) 
        targs <- NULL
    if (missing(mlab)) 
        mlab <- NULL
    if (missing(cex)) 
        cex <- NULL
    alpha <- (100 - level)/100
    yi <- x
    if (hasArg(ci.lb) && hasArg(ci.ub)) {
        if (length(ci.lb) != length(ci.ub)) 
            stop("Length of ci.lb and ci.ub do not match.")
        if (missing(vi) && missing(sei)) {
            vi <- ((ci.ub - ci.lb)/(2 * qnorm(alpha/2, lower.tail = FALSE)))^2
        }
        else {
            if (missing(vi)) 
                vi <- sei^2
        }
        if (length(ci.lb) != length(vi)) 
            stop("Length of vi (or sei) does not match length of (ci.lb, ci.ub) pairs.")
    }
    else {
        if (missing(vi)) {
            if (missing(sei)) {
                stop("Must specify either vi, sei, or (ci.lb, ci.ub) pairs.")
            }
            else {
                vi <- sei^2
                ci.lb <- yi - qnorm(alpha/2, lower.tail = FALSE) * 
                  sei
                ci.ub <- yi + qnorm(alpha/2, lower.tail = FALSE) * 
                  sei
            }
        }
        else {
            ci.lb <- yi - qnorm(alpha/2, lower.tail = FALSE) * 
                sqrt(vi)
            ci.ub <- yi + qnorm(alpha/2, lower.tail = FALSE) * 
                sqrt(vi)
        }
    }
    k <- length(yi)
    if (is.null(rows)) {
        rows <- -1:(-k)
    }
    else {
        if (length(rows) == 1L) {
            rows <- rows:(rows - k + 1)
        }
    }
    if (length(rows) != length(yi)) 
        stop("Number of outcomes does not correspond to the length of the rows argument.")
    yivi.na <- is.na(cbind(yi, vi))
    if (any(yivi.na)) {
        not.na <- apply(yivi.na, MARGIN = 1, sum) == 0L
        if (na.act == "na.omit") {
            yi <- yi[not.na]
            vi <- vi[not.na]
            ci.lb <- ci.lb[not.na]
            ci.ub <- ci.ub[not.na]
            mlab <- mlab[not.na]
            rows.new <- rows
            rows.na <- rows[!not.na]
            for (j in seq.int(length(rows.na))) {
                rows.new[rows <= rows.na[j]] <- rows.new[rows <= 
                  rows.na[j]] + 1
            }
            rows <- rows.new[not.na]
        }
        if (na.act == "na.fail") 
            stop("Missing values in results.")
    }
    k <- length(yi)
    yi.ut <- yi
    ci.lb.ut <- ci.lb
    ci.ub.ut <- ci.ub
    if (is.function(transf)) {
        if (is.null(targs)) {
            yi <- sapply(yi, transf)
            ci.lb <- sapply(ci.lb, transf)
            ci.ub <- sapply(ci.ub, transf)
        }
        else {
            yi <- sapply(yi, transf, targs)
            ci.lb <- sapply(ci.lb, transf, targs)
            ci.ub <- sapply(ci.ub, transf, targs)
        }
    }
    ci.bounds <- cbind(ci.lb, ci.ub)
    rev.order <- ifelse(ci.ub < ci.lb, TRUE, FALSE)
    rev.order[is.na(rev.order)] <- FALSE
    ci.bounds[rev.order] <- ci.bounds[rev.order, 2:1]
    ci.lb <- ci.bounds[, 1]
    ci.ub <- ci.bounds[, 2]
    par.usr <- par("usr")
    height <- par.usr[4] - par.usr[3]
    cex.adj <- min(1, 20/height)
    xlim <- par.usr[1:2]
    if (is.null(cex)) 
        cex <- par("cex") * cex.adj
    if (annotate) {
        if (is.function(atransf)) {
            if (is.null(targs)) {
                annotext <- round(cbind(sapply(yi.ut, atransf), 
                  sapply(ci.lb.ut, atransf), sapply(ci.ub.ut, 
                    atransf)), digits)
            }
            else {
                annotext <- round(cbind(sapply(yi.ut, atransf, 
                  targs), sapply(ci.lb.ut, atransf, targs), sapply(ci.ub.ut, 
                  atransf, targs)), digits)
            }
            rev.order <- ifelse(annotext[, 3] < annotext[, 2], 
                TRUE, FALSE)
            rev.order[is.na(rev.order)] <- FALSE
            annotext[rev.order, 2:3] <- annotext[rev.order, 3:2]
        }
        else {
            annotext <- round(cbind(yi, ci.lb, ci.ub), digits)
        }
        annotext <- matrix(apply(annotext, 2, format, nsmall = digits), 
            ncol = 3)
        annotext <- cbind(annotext[, 1], " [ ", annotext[, 2], 
            " , ", annotext[, 3], " ]")
        annotext <- apply(annotext, 1, paste, collapse = "")
        text(x = xlim[2], rows, labels = annotext, pos = 2, cex = cex, 
            ...)
    }
    for (i in seq.int(k)) {
        polygon(x = c(ci.lb[i], yi[i], ci.ub[i], yi[i]), y = c(rows[i], 
            rows[i] + (height/100) * cex * efac, rows[i], rows[i] - 
                (height/100) * cex * efac), col = col, ...)
        if (!is.null(mlab)) {
            text(xlim[1], rows[i], mlab[i], pos = 4, cex = cex, 
                ...)
        }
    }
}
addpoly.rma <-
function (x, row = -2, level = x$level, digits = 2, annotate = TRUE, 
    mlab, transf = FALSE, atransf = FALSE, targs, col = "black", 
    efac = 1, cex, ...) 
{
    if (!is.element("rma", class(x))) 
        stop("Argument 'x' must be an object of class \"rma\".")
    if (!x$int.only) 
        stop("Fitted model should not contain moderators.")
    if (missing(targs)) 
        targs <- NULL
    if (missing(mlab)) 
        mlab <- NULL
    if (missing(cex)) 
        cex <- NULL
    if (is.null(mlab)) 
        mlab <- ifelse((x$method == "FE"), "FE Model", "RE Model")
    addpoly(x$b, vi = x$vb, rows = row, level = level, digits = digits, 
        annotate = annotate, mlab = mlab, transf = transf, atransf = atransf, 
        col = col, targs = targs, efac = efac, cex = cex, ...)
}
AIC.rma <-
function (object, ..., k = 2) 
{
    if (!is.element("rma", class(object))) 
        stop("Argument 'object' must be an object of class \"rma\".")
    AIC(logLik(object), k = k)
}
anova.rma.uni <-
function (object, object2, digits = object$digits, ...) 
{
    if (!is.element("rma.uni", class(object))) 
        stop("Argument 'object' must be an object of class \"rma.uni\".")
    if (!is.element("rma.uni", class(object2))) 
        stop("Argument 'object2' must be an object of class \"rma.uni\".")
    m.f <- object
    m.r <- object2
    if (!(identical(c(m.f$yi), c(m.r$yi)) && identical(c(m.f$vi), 
        c(m.r$vi)))) 
        stop("Observed outcomes and/or sampling variances not equal in the full and reduced model.")
    if (m.f$method == "FE") {
        p.f <- m.f$p
    }
    else {
        p.f <- m.f$p + 1
    }
    if (m.r$method == "FE") {
        p.r <- m.r$p
    }
    else {
        p.r <- m.r$p + 1
    }
    if (p.f == p.r) 
        stop("Models have the same number of parameters. LRT not meaningful.")
    if (p.f < p.r) {
        m.f <- object2
        m.r <- object
        p.s <- p.f
        p.f <- p.r
        p.r <- p.s
    }
    if (m.f$method == "FE" && m.r$method != "FE") 
        stop("Full model uses a fixed- and reduced model uses random/mixed-effects model.")
    p.diff <- p.f - p.r
    if (m.f$method == "REML") {
        LRT <- abs(m.r$fit.stats$REML[2] - m.f$fit.stats$REML[2])
        fit.stats.f <- m.f$fit.stats$REML
        fit.stats.r <- m.r$fit.stats$REML
        if (!identical(m.f$X, m.r$X)) 
            warning("Models with different fixed effects. REML comparisons are not meaningful.")
    }
    else {
        LRT <- abs(m.r$fit.stats$ML[2] - m.f$fit.stats$ML[2])
        fit.stats.f <- m.f$fit.stats$ML
        fit.stats.r <- m.r$fit.stats$ML
    }
    pval <- pchisq(LRT, df = p.diff, lower.tail = FALSE)
    if (m.f$method == "FE" || identical(m.r$tau2, 0)) {
        VAF <- NA
    }
    else {
        VAF <- round(100 * max(0, (m.r$tau2 - m.f$tau2)/m.r$tau2), 
            2)
    }
    res <- list(fit.stats.f, fit.stats.r, p.f, p.r, LRT, pval, 
        m.f$QE, m.r$QE, m.f$tau2, m.r$tau2, VAF, m.f$method, 
        digits)
    names(res) <- c("fit.stats.f", "fit.stats.r", "p.f", "p.r", 
        "LRT", "pval", "QE.f", "QE.r", "tau2.f", "tau2.r", "VAF", 
        "method", "digits")
    class(res) <- c("anova.rma.uni")
    return(res)
}
BIC.rma <-
function (object, ...) 
{
    if (!is.element("rma", class(object))) 
        stop("Argument 'object' must be an object of class \"rma\".")
    BIC(logLik(object))
}
blup <-
function (x, ...) 
UseMethod("blup")
blup.rma.uni <-
function (x, level = x$level, digits = x$digits, transf = FALSE, 
    targs, ...) 
{
    if (!is.element("rma.uni", class(x))) 
        stop("Argument 'x' must be an object of class \"rma.uni\".")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    alpha <- (100 - level)/100
    if (x$knha || x$robust) {
        crit <- qt(alpha/2, df = x$k - x$p, lower.tail = FALSE)
    }
    else {
        crit <- qnorm(alpha/2, lower.tail = FALSE)
    }
    if (missing(targs)) 
        targs <- NULL
    pred <- rep(NA, x$k.f)
    vpred <- rep(NA, x$k.f)
    li <- x$tau2/(x$tau2 + x$vi.f)
    for (i in seq.int(x$k.f)[x$not.na]) {
        Xi <- matrix(x$X.f[i, ], nrow = 1)
        pred[i] <- li[i] * x$yi.f[i] + (1 - li[i]) * Xi %*% x$b
        vpred[i] <- li[i] * x$vi.f[i] + (1 - li[i])^2 * Xi %*% 
            tcrossprod(x$vb, Xi)
    }
    se <- sqrt(vpred)
    pi.lb <- pred - crit * se
    pi.ub <- pred + crit * se
    if (is.function(transf)) {
        if (is.null(targs)) {
            pred <- sapply(pred, transf)
            se <- rep(NA, x$k.f)
            pi.lb <- sapply(pi.lb, transf)
            pi.ub <- sapply(pi.ub, transf)
        }
        else {
            pred <- sapply(pred, transf, targs)
            se <- rep(NA, x$k.f)
            pi.lb <- sapply(pi.lb, transf, targs)
            pi.ub <- sapply(pi.ub, transf, targs)
        }
    }
    pi.bounds <- cbind(pi.lb, pi.ub)
    rev.order <- ifelse(pi.ub < pi.lb, TRUE, FALSE)
    rev.order[is.na(rev.order)] <- FALSE
    pi.bounds[rev.order] <- pi.bounds[rev.order, 2:1]
    pi.lb <- pi.bounds[, 1]
    pi.ub <- pi.bounds[, 2]
    if (na.act == "na.omit") {
        out <- list(pred = pred[x$not.na], se = se[x$not.na], 
            pi.lb = pi.lb[x$not.na], pi.ub = pi.ub[x$not.na])
        out$slab <- x$slab[x$not.na]
    }
    if (na.act == "na.exclude" || na.act == "na.pass") {
        out <- list(pred = pred, se = se, pi.lb = pi.lb, pi.ub = pi.ub)
        out$slab <- x$slab
    }
    if (na.act == "na.fail" && any(!x$not.na)) 
        stop("Missing values in results.")
    out$digits <- digits
    class(out) <- c("list.rma")
    return(out)
}
cbind.escalc <-
function (..., deparse.level = 1) 
{
    dat <- data.frame(..., check.names = FALSE)
    arguments <- list(...)
    var.name.list <- list()
    digits.list <- list()
    i <- 0
    for (arg in arguments) {
        i <- i + 1
        var.name.list[[i]] <- attr(arg, "var.names")
        digits.list[[i]] <- attr(arg, "digits")
    }
    if (length(var.name.list) > 0) 
        attr(dat, "var.names") <- var.name.list[[1]]
    if (length(digits.list) > 0) 
        attr(dat, "digits") <- digits.list[[1]]
    class(dat) <- c("escalc", "data.frame")
    return(dat)
}
coef.permutest.rma.uni <-
function (object, ...) 
{
    if (!is.element("permutest.rma.uni", class(object))) 
        stop("Argument 'object' must be an object of class \"permutest.rma.uni\".")
    x <- object
    res.table <- cbind(x$b, x$se, x$zval, x$pval, x$ci.lb, x$ci.ub)
    dimnames(res.table)[[2]] <- c("estimate", "se", "zval", "pval", 
        "ci.lb", "ci.ub")
    if (x$knha || x$robust) 
        dimnames(res.table)[[2]][3] <- c("tval")
    res.table <- data.frame(res.table)
    return(res.table)
}
coef.rma <-
function (object, ...) 
{
    if (!is.element("rma", class(object))) 
        stop("Argument 'object' must be an object of class \"rma\".")
    x <- object
    coefs <- c(x$b)
    names(coefs) <- rownames(x$b)
    return(coefs)
}
coef.summary.rma <-
function (object, ...) 
{
    if (!is.element("summary.rma", class(object))) 
        stop("Argument 'object' must be an object of class \"summary.rma\".")
    x <- object
    res.table <- cbind(x$b, x$se, x$zval, x$pval, x$ci.lb, x$ci.ub)
    dimnames(res.table)[[2]] <- c("estimate", "se", "zval", "pval", 
        "ci.lb", "ci.ub")
    if (x$knha || x$robust) 
        dimnames(res.table)[[2]][3] <- c("tval")
    res.table <- data.frame(res.table)
    return(res.table)
}
confint.rma.glmm <-
function (object, parm, level = object$level, digits = object$digits, 
    ...) 
{
    if (!is.element("rma.glmm", class(object))) 
        stop("Argument 'object' must be an object of class \"rma.glmm\".")
    stop("Method not yet implemented for objects of class \"rma.glmm\". Sorry!")
}
confint.rma.mh <-
function (object, parm, level = object$level, digits = object$digits, 
    ...) 
{
    if (!is.element("rma.mh", class(object))) 
        stop("Argument 'object' must be an object of class \"rma.mh\".")
    x <- object
    if (x$k == 1) 
        stop("Stopped because k = 1.")
    alpha <- (100 - level)/100
    if (x$measure == "OR" || x$measure == "RR" || x$measure == 
        "IRR") {
        estimate <- exp(x$b)
        ci.lb <- exp(x$b - qnorm(alpha/2, lower.tail = FALSE) * 
            x$se)
        ci.ub <- exp(x$b + qnorm(alpha/2, lower.tail = FALSE) * 
            x$se)
    }
    else {
        estimate <- x$b
        ci.lb <- x$b - qnorm(alpha/2, lower.tail = FALSE) * x$se
        ci.ub <- x$b + qnorm(alpha/2, lower.tail = FALSE) * x$se
    }
    res <- cbind(estimate, ci.lb, ci.ub)
    res <- list(fixed = res)
    rownames(res$fixed) <- ""
    res$digits <- digits
    class(res) <- c("confint.rma")
    return(res)
}
confint.rma.peto <-
function (object, parm, level = object$level, digits = object$digits, 
    ...) 
{
    if (!is.element("rma.peto", class(object))) 
        stop("Argument 'object' must be an object of class \"rma.peto\".")
    x <- object
    if (x$k == 1) 
        stop("Stopped because k = 1.")
    alpha <- (100 - level)/100
    estimate <- exp(x$b)
    ci.lb <- exp(x$b - qnorm(alpha/2, lower.tail = FALSE) * x$se)
    ci.ub <- exp(x$b + qnorm(alpha/2, lower.tail = FALSE) * x$se)
    res <- cbind(estimate, ci.lb, ci.ub)
    res <- list(fixed = res)
    rownames(res$fixed) <- ""
    res$digits <- digits
    class(res) <- c("confint.rma")
    return(res)
}
confint.rma.uni <-
function (object, parm, level = object$level, fixed = FALSE, 
    random = TRUE, digits = object$digits, verbose = FALSE, control, 
    ...) 
{
    if (!is.element("rma.uni", class(object))) 
        stop("Argument 'object' must be an object of class \"rma.uni\".")
    x <- object
    k <- x$k
    p <- x$p
    yi <- x$yi
    vi <- x$vi
    X <- x$X
    if (k == 1) 
        stop("Stopped because k = 1.")
    if (missing(control)) 
        control <- list()
    if (!fixed && !random) 
        stop("At least one of the arguments 'fixed' and 'random' must be TRUE.")
    con <- list(tol = .Machine$double.eps^0.25, maxiter = 1000, 
        tau2.min = x$control$tau2.min, tau2.max = 50, verbose = FALSE)
    con[pmatch(names(control), names(con))] <- control
    if (verbose) 
        con$verbose <- verbose
    if (random) {
        if (!x$allvipos) 
            stop("Cannot compute confidence interval for the amount of (residual)\n  heterogeneity with non-positive sampling variances in the data.")
        alpha <- (100 - level)/100
        crit.u <- qchisq(alpha/2, k - p, lower.tail = FALSE)
        crit.l <- qchisq(alpha/2, k - p, lower.tail = TRUE)
        status.lb <- 1
        status.ub <- 1
        conv <- 1
        if (.QE.func(con$tau2.min, Y = cbind(yi), vi = vi, X = X, 
            k = k, objective = 0, verbose = FALSE) < crit.l) {
            tau2.lb <- NA
            tau2.ub <- NA
        }
        else {
            if (.QE.func(con$tau2.min, Y = cbind(yi), vi = vi, 
                X = X, k = k, objective = 0, verbose = FALSE) > 
                crit.u) {
                tau2.lb <- try(uniroot(.QE.func, interval = c(con$tau2.min, 
                  con$tau2.max), tol = con$tol, maxiter = con$maxiter, 
                  Y = cbind(yi), vi = vi, X = X, k = k, objective = crit.u, 
                  verbose = con$verbose, digits = digits)$root, 
                  silent = TRUE)
                if (!is.numeric(tau2.lb)) {
                  tau2.lb <- NA
                  status.lb <- 0
                  conv <- 0
                }
            }
            else {
                tau2.lb <- con$tau2.min
            }
            tau2.ub <- try(uniroot(.QE.func, interval = c(tau2.lb, 
                con$tau2.max), tol = con$tol, maxiter = con$maxiter, 
                Y = cbind(yi), vi = vi, X = X, k = k, objective = crit.l, 
                verbose = con$verbose, digits = digits)$root, 
                silent = TRUE)
            if (!is.numeric(tau2.ub)) {
                tau2.ub <- NA
                status.ub <- 0
                conv <- 0
            }
        }
        if (status.lb == 0L) {
            warning("Error in iterative search for the lower bound.")
        }
        if (status.ub == 0L) {
            warning("Error in iterative search for the upper bound.")
        }
        if (conv == 0L) {
            stop("Try increasing tau2.max (via the 'control' argument).")
        }
        wi <- 1/vi
        W <- .diag(wi)
        stXWX <- .invcalc(X = X, W = W, k = k)
        P <- W - W %*% X %*% stXWX %*% crossprod(X, W)
        vi.avg <- (k - p)/.tr(P)
        I2.lb <- 100 * tau2.lb/(vi.avg + tau2.lb)
        I2.ub <- 100 * tau2.ub/(vi.avg + tau2.ub)
        H2.lb <- tau2.lb/vi.avg + 1
        H2.ub <- tau2.ub/vi.avg + 1
        tau2 <- c(x$tau2, tau2.lb, tau2.ub)
        tau <- sqrt(c(ifelse(x$tau2 >= 0, x$tau2, NA), ifelse(tau2.lb >= 
            0, tau2.lb, NA), ifelse(tau2.ub >= 0, tau2.ub, NA)))
        I2 <- c(x$I2, I2.lb, I2.ub)
        H2 <- c(x$H2, H2.lb, H2.ub)
        res.random <- rbind(tau2, tau, I2, H2)
        dimnames(res.random)[[1]] <- c("tau^2", "tau", "I^2(%)", 
            "H^2")
        if (x$method == "FE") 
            res.random[, 1] <- NA
        dimnames(res.random)[[2]] <- c("estimate", "ci.lb", "ci.ub")
    }
    if (fixed) {
        alpha <- (100 - level)/100
        if (x$knha || x$robust) {
            crit <- qt(alpha/2, df = k - p, lower.tail = FALSE)
        }
        else {
            crit <- qnorm(alpha/2, lower.tail = FALSE)
        }
        ci.lb <- c(x$b - crit * x$se)
        ci.ub <- c(x$b + crit * x$se)
        res.fixed <- cbind(x$b, ci.lb, ci.ub)
        dimnames(res.fixed)[[2]] <- c("estimate", "ci.lb", "ci.ub")
    }
    res <- list()
    if (fixed) 
        res$fixed <- res.fixed
    if (random) 
        res$random <- res.random
    res$digits <- digits
    res$tau2.min <- con$tau2.min
    class(res) <- c("confint.rma")
    return(res)
}
cooks.distance.rma.uni <-
function (model, ...) 
{
    if (!is.element("rma.uni", class(model))) 
        stop("Argument 'model' must be an object of class \"rma.uni\".")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    x <- model
    if (x$k == 1) 
        stop("Stopped because k = 1.")
    cook.d <- rep(NA, x$k.f)
    if (x$weighted) {
        wi <- 1/(x$vi + x$tau2)
        W <- .diag(wi)
        svb <- crossprod(x$X, W) %*% x$X/x$s2w
    }
    else {
        svb <- chol2inv(chol(x$vb))
    }
    o.warn <- getOption("warn")
    on.exit(options(warn = o.warn))
    options(warn = -1)
    for (i in seq.int(x$k.f)[x$not.na]) {
        res <- try(rma(x$yi.f[-i], x$vi.f[-i], mods = cbind(x$X.f[-i, 
            ]), method = x$method, weighted = x$weighted, intercept = FALSE, 
            knha = x$knha, control = x$control, ...), silent = TRUE)
        if (inherits(res, "try-error")) 
            next
        dfbeta <- x$b - res$b
        cook.d[i] <- (crossprod(dfbeta, svb) %*% dfbeta)
    }
    if (na.act == "na.omit") {
        out <- cook.d[x$not.na]
        names(out) <- x$slab[x$not.na]
    }
    if (na.act == "na.exclude" || na.act == "na.pass") {
        out <- cook.d
        names(out) <- x$slab
    }
    if (na.act == "na.fail" && any(!x$not.na)) 
        stop("Missing values in results.")
    return(out)
}
cumul <-
function (x, ...) 
UseMethod("cumul")
cumul.rma.mh <-
function (x, order, digits = x$digits, transf = FALSE, ...) 
{
    if (!is.element("rma.mh", class(x))) 
        stop("Argument 'x' must be an object of class \"rma.mh\".")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    if (missing(order)) 
        order <- NULL
    if (is.null(order)) 
        order <- seq.int(x$k.f)
    ai.f <- x$ai.f[order]
    bi.f <- x$bi.f[order]
    ci.f <- x$ci.f[order]
    di.f <- x$di.f[order]
    x1i.f <- x$x1i.f[order]
    x2i.f <- x$x2i.f[order]
    t1i.f <- x$t1i.f[order]
    t2i.f <- x$t2i.f[order]
    yi.f <- x$yi.f[order]
    vi.f <- x$vi.f[order]
    not.na <- x$not.na[order]
    slab <- x$slab[order]
    o.warn <- getOption("warn")
    on.exit(options(warn = o.warn))
    options(warn = -1)
    b <- rep(NA, x$k.f)
    se <- rep(NA, x$k.f)
    zval <- rep(NA, x$k.f)
    pval <- rep(NA, x$k.f)
    ci.lb <- rep(NA, x$k.f)
    ci.ub <- rep(NA, x$k.f)
    QE <- rep(NA, x$k.f)
    QEp <- rep(NA, x$k.f)
    for (i in seq.int(x$k.f)[x$not.na]) {
        if (is.element(x$measure, c("RR", "OR", "RD"))) {
            res <- try(rma.mh(ai = ai.f[seq.int(i)], bi = bi.f[seq.int(i)], 
                ci = ci.f[seq.int(i)], di = di.f[seq.int(i)], 
                measure = x$measure, add = x$add, to = x$to, 
                ...), silent = TRUE)
        }
        else {
            res <- try(rma.mh(x1i = x1i.f[seq.int(i)], x2i = x2i.f[seq.int(i)], 
                t1i = t1i.f[seq.int(i)], t2i = t2i.f[seq.int(i)], 
                measure = x$measure, add = x$add, to = x$to, 
                ...), silent = TRUE)
        }
        if (inherits(res, "try-error")) 
            next
        b[i] <- res$b
        se[i] <- res$se
        zval[i] <- res$zval
        pval[i] <- res$pval
        ci.lb[i] <- res$ci.lb
        ci.ub[i] <- res$ci.ub
        QE[i] <- res$QE
        QEp[i] <- res$QEp
    }
    alpha <- (100 - x$level)/100
    crit <- qnorm(alpha/2, lower.tail = FALSE)
    b[1] <- yi.f[1]
    se[1] <- sqrt(vi.f[1])
    zval[1] <- yi.f[1]/se[1]
    pval[1] <- 2 * pnorm(abs(zval[1]), lower.tail = FALSE)
    ci.lb[1] <- yi.f[1] - crit * se[1]
    ci.ub[1] <- yi.f[1] + crit * se[1]
    QE[1] <- 0
    QEp[1] <- 1
    if (transf) {
        if (x$measure == "OR" || x$measure == "RR" || x$measure == 
            "IRR") {
            b <- exp(b)
            se <- rep(NA, x$k.f)
            ci.lb <- exp(ci.lb)
            ci.ub <- exp(ci.ub)
        }
    }
    if (na.act == "na.omit") {
        out <- list(estimate = b[not.na], se = se[not.na], zval = zval[not.na], 
            pval = pval[not.na], ci.lb = ci.lb[not.na], ci.ub = ci.ub[not.na], 
            QE = QE[not.na], QEp = QEp[not.na])
        out$slab <- slab[not.na]
    }
    if (na.act == "na.exclude" || na.act == "na.pass") {
        out <- list(estimate = b, se = se, zval = zval, pval = pval, 
            ci.lb = ci.lb, ci.ub = ci.ub, QE = QE, QEp = QEp)
        out$slab <- slab
    }
    if (na.act == "na.fail" && any(!x$not.na)) 
        stop("Missing values in results.")
    out$digits <- digits
    out$slab.null <- x$slab.null
    out$level <- x$level
    out$measure <- x$measure
    if (x$measure == "GEN") {
        attr(out$estimate, "measure") <- "GEN"
    }
    else {
        attr(out$estimate, "measure") <- x$measure
    }
    class(out) <- c("list.rma", "cumul.rma")
    return(out)
}
cumul.rma.peto <-
function (x, order, digits = x$digits, transf = FALSE, ...) 
{
    if (!is.element("rma.peto", class(x))) 
        stop("Argument 'x' must be an object of class \"rma.peto\".")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    if (missing(order)) 
        order <- NULL
    if (is.null(order)) 
        order <- seq.int(x$k.f)
    ai.f <- x$ai.f[order]
    bi.f <- x$bi.f[order]
    ci.f <- x$ci.f[order]
    di.f <- x$di.f[order]
    yi.f <- x$yi.f[order]
    vi.f <- x$vi.f[order]
    not.na <- x$not.na[order]
    slab <- x$slab[order]
    o.warn <- getOption("warn")
    on.exit(options(warn = o.warn))
    options(warn = -1)
    b <- rep(NA, x$k.f)
    se <- rep(NA, x$k.f)
    zval <- rep(NA, x$k.f)
    pval <- rep(NA, x$k.f)
    ci.lb <- rep(NA, x$k.f)
    ci.ub <- rep(NA, x$k.f)
    QE <- rep(NA, x$k.f)
    QEp <- rep(NA, x$k.f)
    for (i in seq.int(x$k.f)[x$not.na]) {
        res <- try(rma.peto(ai = ai.f[seq.int(i)], bi = bi.f[seq.int(i)], 
            ci = ci.f[seq.int(i)], di = di.f[seq.int(i)], add = x$add, 
            to = x$to, ...), silent = TRUE)
        if (inherits(res, "try-error")) 
            next
        b[i] <- res$b
        se[i] <- res$se
        zval[i] <- res$zval
        pval[i] <- res$pval
        ci.lb[i] <- res$ci.lb
        ci.ub[i] <- res$ci.ub
        QE[i] <- res$QE
        QEp[i] <- res$QEp
    }
    alpha <- (100 - x$level)/100
    crit <- qnorm(alpha/2, lower.tail = FALSE)
    b[1] <- yi.f[1]
    se[1] <- sqrt(vi.f[1])
    zval[1] <- yi.f[1]/se[1]
    pval[1] <- 2 * pnorm(abs(zval[1]), lower.tail = FALSE)
    ci.lb[1] <- yi.f[1] - crit * se[1]
    ci.ub[1] <- yi.f[1] + crit * se[1]
    QE[1] <- 0
    QEp[1] <- 1
    if (transf) {
        b <- exp(b)
        se <- rep(NA, x$k.f)
        ci.lb <- exp(ci.lb)
        ci.ub <- exp(ci.ub)
    }
    if (na.act == "na.omit") {
        out <- list(estimate = b[not.na], se = se[not.na], zval = zval[not.na], 
            pval = pval[not.na], ci.lb = ci.lb[not.na], ci.ub = ci.ub[not.na], 
            Q = QE[not.na], Qp = QEp[not.na])
        out$slab <- slab[not.na]
    }
    if (na.act == "na.exclude" || na.act == "na.pass") {
        out <- list(estimate = b, se = se, zval = zval, pval = pval, 
            ci.lb = ci.lb, ci.ub = ci.ub, Q = QE, Qp = QEp)
        out$slab <- slab
    }
    if (na.act == "na.fail" && any(!x$not.na)) 
        stop("Missing values in results.")
    out$digits <- digits
    out$slab.null <- x$slab.null
    out$level <- x$level
    out$measure <- x$measure
    if (x$measure == "GEN") {
        attr(out$estimate, "measure") <- "GEN"
    }
    else {
        attr(out$estimate, "measure") <- x$measure
    }
    class(out) <- c("list.rma", "cumul.rma")
    return(out)
}
cumul.rma.uni <-
function (x, order, digits = x$digits, transf = FALSE, targs, 
    ...) 
{
    if (!is.element("rma.uni", class(x))) 
        stop("Argument 'x' must be an object of class \"rma.uni\".")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    if (!x$int.only) 
        stop("Method only applicable for models without moderators.")
    if (missing(order)) 
        order <- NULL
    if (missing(targs)) 
        targs <- NULL
    if (is.null(order)) 
        order <- seq.int(x$k.f)
    yi.f <- x$yi.f[order]
    vi.f <- x$vi.f[order]
    X.f <- cbind(x$X.f[order, ])
    not.na <- x$not.na[order]
    slab <- x$slab[order]
    o.warn <- getOption("warn")
    on.exit(options(warn = o.warn))
    options(warn = -1)
    b <- rep(NA, x$k.f)
    se <- rep(NA, x$k.f)
    zval <- rep(NA, x$k.f)
    pval <- rep(NA, x$k.f)
    ci.lb <- rep(NA, x$k.f)
    ci.ub <- rep(NA, x$k.f)
    QE <- rep(NA, x$k.f)
    QEp <- rep(NA, x$k.f)
    tau2 <- rep(NA, x$k.f)
    I2 <- rep(NA, x$k.f)
    H2 <- rep(NA, x$k.f)
    for (i in seq.int(x$k.f)[not.na]) {
        res <- try(rma(yi.f[seq.int(i)], vi.f[seq.int(i)], method = x$method, 
            weighted = x$weighted, intercept = TRUE, knha = x$knha, 
            control = x$control, ...), silent = TRUE)
        if (inherits(res, "try-error")) 
            next
        b[i] <- res$b
        se[i] <- res$se
        zval[i] <- res$zval
        pval[i] <- res$pval
        ci.lb[i] <- res$ci.lb
        ci.ub[i] <- res$ci.ub
        QE[i] <- res$QE
        QEp[i] <- res$QEp
        tau2[i] <- res$tau2
        I2[i] <- res$I2
        H2[i] <- res$H2
    }
    alpha <- (100 - x$level)/100
    crit <- qnorm(alpha/2, lower.tail = FALSE)
    b[1] <- yi.f[1]
    se[1] <- sqrt(vi.f[1])
    zval[1] <- yi.f[1]/se[1]
    pval[1] <- 2 * pnorm(abs(zval[1]), lower.tail = FALSE)
    ci.lb[1] <- yi.f[1] - crit * se[1]
    ci.ub[1] <- yi.f[1] + crit * se[1]
    QE[1] <- 0
    QEp[1] <- 1
    tau2[1] <- 0
    I2[1] <- 0
    H2[1] <- 1
    if (is.function(transf)) {
        if (is.null(targs)) {
            b <- sapply(b, transf)
            se <- rep(NA, x$k.f)
            ci.lb <- sapply(ci.lb, transf)
            ci.ub <- sapply(ci.ub, transf)
        }
        else {
            b <- sapply(b, transf, targs)
            se <- rep(NA, x$k.f)
            ci.lb <- sapply(ci.lb, transf, targs)
            ci.ub <- sapply(ci.ub, transf, targs)
        }
    }
    ci.bounds <- cbind(ci.lb, ci.ub)
    rev.order <- ifelse(ci.ub < ci.lb, TRUE, FALSE)
    rev.order[is.na(rev.order)] <- FALSE
    ci.bounds[rev.order] <- ci.bounds[rev.order, 2:1]
    ci.lb <- ci.bounds[, 1]
    ci.ub <- ci.bounds[, 2]
    if (na.act == "na.omit") {
        out <- list(estimate = b[not.na], se = se[not.na], zval = zval[not.na], 
            pvals = pval[not.na], ci.lb = ci.lb[not.na], ci.ub = ci.ub[not.na], 
            QE = QE[not.na], QEp = QEp[not.na], tau2 = tau2[not.na], 
            I2 = I2[not.na], H2 = H2[not.na])
        out$slab <- slab[not.na]
    }
    if (na.act == "na.exclude" || na.act == "na.pass") {
        out <- list(estimate = b, se = se, zval = zval, pvals = pval, 
            ci.lb = ci.lb, ci.ub = ci.ub, QE = QE, QEp = QEp, 
            tau2 = tau2, I2 = I2, H2 = H2)
        out$slab <- slab
    }
    if (na.act == "na.fail" && any(!x$not.na)) 
        stop("Missing values in results.")
    if (x$method == "FE") 
        out <- out[-c(9, 10, 11)]
    out$digits <- digits
    out$slab.null <- x$slab.null
    out$level <- x$level
    out$measure <- x$measure
    if (x$measure == "GEN") {
        attr(out$estimate, "measure") <- "GEN"
    }
    else {
        attr(out$estimate, "measure") <- x$measure
    }
    class(out) <- c("list.rma", "cumul.rma")
    return(out)
}
deviance.rma <-
function (object, REML, ...) 
{
    if (!is.element("rma", class(object))) 
        stop("Argument 'object' must be an object of class \"rma\".")
    if (missing(REML)) {
        if (object$method == "REML") {
            REML <- TRUE
        }
        else {
            REML <- FALSE
        }
    }
    if (REML) {
        val <- object$fit.stats$REML[2]
    }
    else {
        val <- object$fit.stats$ML[2]
    }
    return(val)
}
df.residual.rma <-
function (object, ...) 
{
    if (!is.element("rma", class(object))) 
        stop("Argument 'object' must be an object of class \"rma\".")
    df.resid <- object$k.eff - object$p.eff
    return(df.resid)
}
dfbetas.rma.uni <-
function (model, ...) 
{
    if (!is.element("rma.uni", class(model))) 
        stop("Argument 'model' must be an object of class \"rma.uni\".")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    x <- model
    if (x$k == 1) 
        stop("Stopped because k = 1.")
    tau2.del <- rep(NA, x$k.f)
    dfbetas <- matrix(NA, nrow = x$k.f, ncol = length(x$b))
    if (!x$weighted) {
        svb <- chol2inv(chol(x$vb))
        stXX <- .invcalc(X = x$X, W = diag(x$k), k = x$k)
    }
    o.warn <- getOption("warn")
    on.exit(options(warn = o.warn))
    options(warn = -1)
    for (i in seq.int(x$k.f)[x$not.na]) {
        res <- try(rma(x$yi.f[-i], x$vi.f[-i], mods = cbind(x$X.f[-i, 
            ]), method = x$method, weighted = x$weighted, intercept = FALSE, 
            knha = x$knha, control = x$control, ...), silent = TRUE)
        if (inherits(res, "try-error")) 
            next
        tau2.del[i] <- res$tau2
        dfbeta <- x$b - res$b
        if (x$weighted) {
            vb.del <- .invcalc(X = x$X, W = .diag(1/(x$vi + tau2.del[i])), 
                k = x$k)
            dfbetas[i, ] <- dfbeta/sqrt(res$s2w * diag(vb.del))
        }
        else {
            vb.del <- tcrossprod(stXX, x$X) %*% .diag(x$vi + 
                tau2.del[i]) %*% x$X %*% stXX
            dfbetas[i, ] <- dfbeta/sqrt(res$s2w * diag(vb.del))
        }
    }
    if (na.act == "na.omit") {
        out <- dfbetas[x$not.na, ]
        dimnames(out)[[1]] <- x$slab[x$not.na]
    }
    if (na.act == "na.exclude" || na.act == "na.pass") {
        out <- dfbetas
        dimnames(out)[[1]] <- x$slab
    }
    if (na.act == "na.fail" && any(!x$not.na)) 
        stop("Missing values in results.")
    dimnames(out)[[2]] <- dimnames(x$b)[[1]]
    out <- data.frame(out)
    return(out)
}
escalc <-
function (measure, formula, ...) 
{
    if (missing(measure) || class(measure) == "formula") 
        stop("Must specify an effect size or outcome measure.")
    if (missing(formula)) 
        formula <- NULL
    UseMethod("escalc", formula)
}
escalc.default <-
function (measure, formula, ai, bi, ci, di, n1i, n2i, x1i, x2i, 
    t1i, t2i, m1i, m2i, sd1i, sd2i, xi, mi, ri, ti, sdi, ni, 
    data, slab, subset, add = 1/2, to = "only0", drop00 = FALSE, 
    vtype = "LS", var.names = c("yi", "vi"), append = TRUE, replace = TRUE, 
    digits = 4, ...) 
{
    if (!is.element(measure, c("RR", "OR", "PETO", "RD", "AS", 
        "PHI", "YUQ", "YUY", "RTET", "PBIT", "OR2D", "IRR", "IRD", 
        "IRSD", "MD", "SMD", "SMDH", "ROM", "RPB", "RBIS", "D2OR", 
        "COR", "UCOR", "ZCOR", "PR", "PLN", "PLO", "PAS", "PFT", 
        "IR", "IRLN", "IRS", "IRFT", "MN", "MC", "SMCC", "SMCR", 
        "ARAW", "AHW", "ABT"))) 
        stop("Unknown 'measure' specified.")
    if (!is.element(to, c("all", "only0", "if0all", "none"))) 
        stop("Unknown 'to' argument specified.")
    if (any(!is.element(vtype, c("UB", "LS", "HO", "ST", "CS")), 
        na.rm = TRUE)) 
        stop("Unknown 'vtype' argument specified.")
    if (length(var.names) != 2) 
        stop("Argument var.names must be of length 2.")
    if (any(var.names != make.names(var.names, unique = TRUE))) {
        var.names <- make.names(var.names, unique = TRUE)
        warning(paste("Argument 'var.names' does not contain syntactically valid variable names.\n  Variable names adjusted to: var.names = c('", 
            var.names[1], "', '", var.names[2], "').", sep = ""))
    }
    if (missing(data)) 
        data <- NULL
    no.data <- is.null(data)
    if (is.null(data)) {
        data <- sys.frame(sys.parent())
    }
    else {
        if (!is.data.frame(data)) {
            data <- data.frame(data)
        }
    }
    mf <- match.call()
    mf.slab <- mf[[match("slab", names(mf))]]
    mf.subset <- mf[[match("subset", names(mf))]]
    slab <- eval(mf.slab, data, enclos = sys.frame(sys.parent()))
    subset <- eval(mf.subset, data, enclos = sys.frame(sys.parent()))
    if (is.element(measure, c("RR", "OR", "RD", "AS", "PETO", 
        "PHI", "YUQ", "YUY", "RTET", "PBIT", "OR2D"))) {
        mf.ai <- mf[[match("ai", names(mf))]]
        mf.bi <- mf[[match("bi", names(mf))]]
        mf.ci <- mf[[match("ci", names(mf))]]
        mf.di <- mf[[match("di", names(mf))]]
        mf.n1i <- mf[[match("n1i", names(mf))]]
        mf.n2i <- mf[[match("n2i", names(mf))]]
        ai <- eval(mf.ai, data, enclos = sys.frame(sys.parent()))
        bi <- eval(mf.bi, data, enclos = sys.frame(sys.parent()))
        ci <- eval(mf.ci, data, enclos = sys.frame(sys.parent()))
        di <- eval(mf.di, data, enclos = sys.frame(sys.parent()))
        n1i <- eval(mf.n1i, data, enclos = sys.frame(sys.parent()))
        n2i <- eval(mf.n2i, data, enclos = sys.frame(sys.parent()))
        if (is.null(bi)) 
            bi <- n1i - ai
        if (is.null(di)) 
            di <- n2i - ci
        if (length(ai) == 0L || length(bi) == 0L || length(ci) == 
            0L || length(di) == 0L) 
            stop("Cannot compute outcomes. Check that all of the required \n  information is specified via the appropriate arguments.")
        if (any(c(ai, bi, ci, di) < 0, na.rm = TRUE)) 
            stop("One or more counts are negative.")
        ni.u <- ai + bi + ci + di
        if (drop00) {
            id00 <- c(ai == 0L & ci == 0L) | c(bi == 0L & di == 
                0L)
            id00[is.na(id00)] <- FALSE
            ai[id00] <- NA
            bi[id00] <- NA
            ci[id00] <- NA
            di[id00] <- NA
        }
        if (to == "all") {
            ai <- ai + add
            ci <- ci + add
            bi <- bi + add
            di <- di + add
        }
        if (to == "only0") {
            id0 <- c(ai == 0L | ci == 0L | bi == 0L | di == 0L)
            id0[is.na(id0)] <- FALSE
            ai[id0] <- ai[id0] + add
            ci[id0] <- ci[id0] + add
            bi[id0] <- bi[id0] + add
            di[id0] <- di[id0] + add
        }
        if (to == "if0all") {
            id0 <- c(ai == 0L | ci == 0L | bi == 0L | di == 0L)
            id0[is.na(id0)] <- FALSE
            if (any(id0)) {
                ai <- ai + add
                ci <- ci + add
                bi <- bi + add
                di <- di + add
            }
        }
        n1i <- ai + bi
        n2i <- ci + di
        ni <- n1i + n2i
        p1i <- ai/n1i
        p2i <- ci/n2i
        if (measure == "RR") {
            yi <- log(p1i) - log(p2i)
            vi <- 1/ai - 1/n1i + 1/ci - 1/n2i
        }
        if (is.element(measure, c("OR", "OR2D"))) {
            yi <- log(p1i/(1 - p1i)) - log(p2i/(1 - p2i))
            vi <- 1/ai + 1/bi + 1/ci + 1/di
        }
        if (measure == "PETO") {
            xt <- ai + ci
            yt <- bi + di
            Oi <- ai
            Ei <- xt * n1i/ni
            Vi <- xt * yt * (n1i/ni) * (n2i/ni)/(ni - 1)
            yi <- (Oi - Ei)/Vi
            vi <- 1/Vi
        }
        if (measure == "RD") {
            yi <- p1i - p2i
            k <- length(yi)
            if (length(vtype) == 1L) 
                vtype <- rep(vtype, k)
            vi <- rep(NA, k)
            mnwp1i <- sum(ai, na.rm = TRUE)/sum(n1i, na.rm = TRUE)
            mnwp2i <- sum(ci, na.rm = TRUE)/sum(n2i, na.rm = TRUE)
            for (i in seq.int(k)) {
                if (vtype[i] == "UB") 
                  vi[i] <- p1i[i] * (1 - p1i[i])/(n1i[i] - 1) + 
                    p2i[i] * (1 - p2i[i])/(n2i[i] - 1)
                if (vtype[i] == "LS") 
                  vi[i] <- p1i[i] * (1 - p1i[i])/n1i[i] + p2i[i] * 
                    (1 - p2i[i])/n2i[i]
                if (vtype[i] == "HO") 
                  vi[i] <- mnwp1i * (1 - mnwp1i)/n1i[i] + mnwp2i * 
                    (1 - mnwp2i)/n2i[i]
            }
        }
        if (measure == "AS") {
            yi <- asin(sqrt(p1i)) - asin(sqrt(p2i))
            vi <- 1/(4 * n1i) + 1/(4 * n2i)
        }
        if (measure == "PHI") {
            yi <- (ai * di - bi * ci)/sqrt((ai + bi) * (ci + 
                di) * (ai + ci) * (bi + di))
            p1. <- (ai + bi)/ni
            p2. <- (ci + di)/ni
            p.1 <- (ai + ci)/ni
            p.2 <- (bi + di)/ni
            vi <- 1/ni * (1 - yi^2 + yi * (1 + 1/2 * yi^2) * 
                (p1. - p2.) * (p.1 - p.2)/sqrt(p1. * p2. * p.1 * 
                p.2) - 3/4 * yi^2 * ((p1. - p2.)^2/(p1. * p2.) + 
                (p.1 - p.2)^2/(p.1 * p.2)))
        }
        if (measure == "YUQ") {
            ori <- ai * di/(bi * ci)
            yi <- (ori - 1)/(ori + 1)
            vi <- 1/4 * (1 - yi^2)^2 * (1/ai + 1/bi + 1/ci + 
                1/di)
        }
        if (measure == "YUY") {
            ori <- ai * di/(bi * ci)
            yi <- (sqrt(ori) - 1)/(sqrt(ori) + 1)
            vi <- 1/16 * (1 - yi^2)^2 * (1/ai + 1/bi + 1/ci + 
                1/di)
        }
        if (measure == "RTET") {
            if (!require(polycor)) 
                stop("Please install the 'polycor' package to compute this measure.")
            warn.before <- getOption("warn")
            options(warn = 2)
            k <- length(ai)
            yi <- rep(NA, k)
            vi <- rep(NA, k)
            for (i in seq.int(k)) {
                if (is.na(ai[i]) || is.na(bi[i]) || is.na(ci[i]) || 
                  is.na(di[i])) 
                  next
                tab <- as.table(matrix(c(ai[i], bi[i], ci[i], 
                  di[i]), nrow = 2, byrow = TRUE))
                res <- try(polychor(tab, std.err = TRUE, ML = TRUE, 
                  maxcor = 0.9999), silent = TRUE)
                if (inherits(res, "try-error")) {
                  next
                }
                else {
                  yi[i] <- c(res$rho)
                  vi[i] <- c(res$var[1, 1])
                }
            }
            options(warn = warn.before)
        }
        if (measure == "PBIT") {
            z1 <- qnorm(p1i)
            z2 <- qnorm(p2i)
            yi <- z1 - z2
            vi <- 2 * pi * p1i * (1 - p1i) * exp(z1^2)/n1i + 
                2 * pi * p2i * (1 - p2i) * exp(z2^2)/n2i
        }
        if (measure == "OR2D") {
            yi <- sqrt(3)/pi * yi
            vi <- 3/pi^2 * vi
        }
    }
    if (is.element(measure, c("IRR", "IRD", "IRSD"))) {
        mf.x1i <- mf[[match("x1i", names(mf))]]
        mf.x2i <- mf[[match("x2i", names(mf))]]
        mf.t1i <- mf[[match("t1i", names(mf))]]
        mf.t2i <- mf[[match("t2i", names(mf))]]
        x1i <- eval(mf.x1i, data, enclos = sys.frame(sys.parent()))
        x2i <- eval(mf.x2i, data, enclos = sys.frame(sys.parent()))
        t1i <- eval(mf.t1i, data, enclos = sys.frame(sys.parent()))
        t2i <- eval(mf.t2i, data, enclos = sys.frame(sys.parent()))
        if (length(x1i) == 0L || length(x2i) == 0L || length(t1i) == 
            0L || length(t2i) == 0L) 
            stop("Cannot compute outcomes. Check that all of the required \n  information is specified via the appropriate arguments.")
        if (any(c(x1i, x2i) < 0, na.rm = TRUE)) 
            stop("One or more counts are negative.")
        if (any(c(t1i, t2i) < 0, na.rm = TRUE)) 
            stop("One or more person-times are negative.")
        ni.u <- t1i + t2i
        if (drop00) {
            id00 <- c(x1i == 0L & x2i == 0L)
            id00[is.na(id00)] <- FALSE
            x1i[id00] <- NA
            x2i[id00] <- NA
        }
        if (to == "all") {
            x1i <- x1i + add
            x2i <- x2i + add
        }
        if (to == "only0") {
            id0 <- c(x1i == 0L | x2i == 0L)
            id0[is.na(id0)] <- FALSE
            x1i[id0] <- x1i[id0] + add
            x2i[id0] <- x2i[id0] + add
        }
        if (to == "if0all") {
            id0 <- c(x1i == 0L | x2i == 0L)
            id0[is.na(id0)] <- FALSE
            if (any(id0)) {
                x1i <- x1i + add
                x2i <- x2i + add
            }
        }
        ir1i <- x1i/t1i
        ir2i <- x2i/t2i
        if (measure == "IRR") {
            yi <- log(ir1i) - log(ir2i)
            vi <- 1/x1i + 1/x2i
        }
        if (measure == "IRD") {
            yi <- ir1i - ir2i
            vi <- ir1i/t1i + ir2i/t2i
        }
        if (measure == "IRSD") {
            yi <- sqrt(ir1i) - sqrt(ir2i)
            vi <- 1/(4 * t1i) + 1/(4 * t2i)
        }
    }
    if (is.element(measure, c("MD", "SMD", "SMDH", "ROM", "RPB", 
        "RBIS", "D2OR"))) {
        mf.m1i <- mf[[match("m1i", names(mf))]]
        mf.m2i <- mf[[match("m2i", names(mf))]]
        mf.sd1i <- mf[[match("sd1i", names(mf))]]
        mf.sd2i <- mf[[match("sd2i", names(mf))]]
        mf.n1i <- mf[[match("n1i", names(mf))]]
        mf.n2i <- mf[[match("n2i", names(mf))]]
        m1i <- eval(mf.m1i, data, enclos = sys.frame(sys.parent()))
        m2i <- eval(mf.m2i, data, enclos = sys.frame(sys.parent()))
        sd1i <- eval(mf.sd1i, data, enclos = sys.frame(sys.parent()))
        sd2i <- eval(mf.sd2i, data, enclos = sys.frame(sys.parent()))
        n1i <- eval(mf.n1i, data, enclos = sys.frame(sys.parent()))
        n2i <- eval(mf.n2i, data, enclos = sys.frame(sys.parent()))
        if (length(m1i) == 0L || length(m2i) == 0L || length(sd1i) == 
            0L || length(sd2i) == 0L || length(n1i) == 0L || 
            length(n2i) == 0L) 
            stop("Cannot compute outcomes. Check that all of the required \n  information is specified via the appropriate arguments.")
        if (any(c(sd1i, sd2i) < 0, na.rm = TRUE)) 
            stop("One or more standard deviations are negative.")
        ni.u <- n1i + n2i
        if (measure == "MD") {
            yi <- m1i - m2i
            vi <- sd1i^2/n1i + sd2i^2/n2i
        }
        if (measure == "SMD") {
            ni <- n1i + n2i
            mi <- ni - 2
            warn.before <- getOption("warn")
            options(warn = -1)
            cmi <- .cmicalc(mi)
            options(warn = warn.before)
            spi <- sqrt(((n1i - 1) * sd1i^2 + (n2i - 1) * sd2i^2)/mi)
            yi <- cmi * (m1i - m2i)/spi
            k <- length(yi)
            if (length(vtype) == 1L) 
                vtype <- rep(vtype, k)
            vi <- rep(NA, k)
            mnwyi <- sum(ni * yi, na.rm = TRUE)/sum(ni, na.rm = TRUE)
            for (i in seq.int(k)) {
                if (vtype[i] == "UB") 
                  vi[i] <- (1/n1i[i] + 1/n2i[i]) + (1 - (mi[i] - 
                    2)/(mi[i] * cmi[i]^2)) * yi[i]^2
                if (vtype[i] == "LS") 
                  vi[i] <- (1/n1i[i] + 1/n2i[i]) + yi[i]^2/(2 * 
                    (n1i[i] + n2i[i]))
                if (vtype[i] == "HO") 
                  vi[i] <- (1/n1i[i] + 1/n2i[i]) + mnwyi^2/(2 * 
                    (n1i[i] + n2i[i]))
            }
        }
        if (measure == "SMDH") {
            mi <- n1i + n2i - 2
            warn.before <- getOption("warn")
            options(warn = -1)
            cmi <- .cmicalc(mi)
            options(warn = warn.before)
            si <- sqrt((sd1i^2 + sd2i^2)/2)
            yi <- cmi * (m1i - m2i)/si
            vi <- yi^2 * (sd1i^4/(n1i - 1) + sd2i^4/(n2i - 1))/(2 * 
                (sd1i^2 + sd2i^2)^2) + (sd1i^2/(n1i - 1) + sd2i^2/(n2i - 
                1))/((sd1i^2 + sd2i^2)/2)
            vi <- vi * cmi^2
        }
        if (measure == "ROM") {
            yi <- log(m1i/m2i)
            vi <- sd1i^2/(n1i * m1i^2) + sd2i^2/(n2i * m2i^2)
        }
        if (is.element(measure, c("RPB", "RBIS"))) {
            ni <- n1i + n2i
            mi <- ni - 2
            spi <- sqrt(((n1i - 1) * sd1i^2 + (n2i - 1) * sd2i^2)/mi)
            di <- (m1i - m2i)/spi
            hi <- mi/n1i + mi/n2i
            yi <- di/sqrt(di^2 + hi)
            k <- length(yi)
            if (length(vtype) == 1L) 
                vtype <- rep(vtype, k)
            vi <- rep(NA, k)
            for (i in seq.int(k)) {
                if (vtype[i] == "ST" || vtype[i] == "LS") 
                  vi[i] <- hi[i]^2/(hi[i] + di[i]^2)^3 * (1/n1i[i] + 
                    1/n2i[i] + di[i]^2/(2 * ni[i]))
                if (vtype[i] == "CS") 
                  vi[i] <- (1 - yi[i]^2)^2 * (ni[i] * yi[i]^2/(4 * 
                    n1i[i] * n2i[i]) + (2 - 3 * yi[i]^2)/(2 * 
                    ni[i]))
            }
        }
        if (measure == "RBIS") {
            p1i <- n1i/ni
            p2i <- n2i/ni
            zi <- qnorm(p1i, lower.tail = FALSE)
            fzi <- dnorm(zi)
            yi <- sqrt(p1i * p2i)/fzi * yi
            yi.t <- ifelse(abs(yi) > 1, sign(yi), yi)
            vi <- 1/(ni - 1) * (p1i * p2i/fzi^2 - (3/2 + (1 - 
                p1i * zi/fzi) * (1 + p2i * zi/fzi)) * yi.t^2 + 
                yi.t^4)
        }
        if (measure == "D2OR") {
            ni <- n1i + n2i
            mi <- ni - 2
            spi <- sqrt(((n1i - 1) * sd1i^2 + (n2i - 1) * sd2i^2)/mi)
            di <- (m1i - m2i)/spi
            yi <- pi/sqrt(3) * di
            vi <- pi^2/3 * ((1/n1i + 1/n2i) + di^2/(2 * ni))
        }
    }
    if (is.element(measure, c("COR", "UCOR", "ZCOR"))) {
        mf.ri <- mf[[match("ri", names(mf))]]
        mf.ni <- mf[[match("ni", names(mf))]]
        ri <- eval(mf.ri, data, enclos = sys.frame(sys.parent()))
        ni <- eval(mf.ni, data, enclos = sys.frame(sys.parent()))
        if (length(ri) == 0L || length(ni) == 0L) 
            stop("Cannot compute outcomes. Check that all of the required \n  information is specified via the appropriate arguments.")
        if (any(abs(ri) > 1, na.rm = TRUE)) 
            stop("One or more correlations are > 1 or < -1.")
        if (any(abs(ni) <= 4, na.rm = TRUE)) 
            warning("Cannot estimate sampling variance when ni <= 4.")
        ni.u <- ni
        if (measure == "COR") {
            yi <- ri
            k <- length(yi)
            if (length(vtype) == 1L) 
                vtype <- rep(vtype, k)
            vi <- rep(NA, k)
            mnwyi <- sum(ni * yi, na.rm = TRUE)/sum(ni, na.rm = TRUE)
            for (i in seq.int(k)) {
                if (vtype[i] == "UB") 
                  vi[i] <- (ri[i] + ri[i] * (1 - ri[i]^2)/(2 * 
                    (ni[i] - 4)))^2 - 1 + (ni[i] - 3)/(ni[i] - 
                    2) * ((1 - ri[i]^2) + 2 * (1 - ri[i]^2)^2/ni[i] + 
                    8 * (1 - ri[i]^2)^3/(ni[i] * (ni[i] + 2)) + 
                    48 * (1 - ri[i]^2)^4/(ni[i] * (ni[i] + 2) * 
                      (ni[i] + 4)))
                if (vtype[i] == "LS") 
                  vi[i] <- (1 - yi[i]^2)^2/(ni[i] - 1)
                if (vtype[i] == "HO") 
                  vi[i] <- (1 - mnwyi^2)^2/(ni[i] - 1)
            }
            vi[ni <= 4] <- NA
        }
        if (measure == "UCOR") {
            yi <- ri + ri * (1 - ri^2)/(2 * (ni - 4))
            yi[ni <= 4] <- NA
            k <- length(yi)
            if (length(vtype) == 1L) 
                vtype <- rep(vtype, k)
            vi <- rep(NA, k)
            mnwyi <- sum(ni * yi, na.rm = TRUE)/sum(ni, na.rm = TRUE)
            for (i in seq.int(k)) {
                if (vtype[i] == "UB") {
                  vi[i] <- yi[i]^2 - 1 + (ni[i] - 3)/(ni[i] - 
                    2) * ((1 - ri[i]^2) + 2 * (1 - ri[i]^2)^2/ni[i] + 
                    8 * (1 - ri[i]^2)^3/(ni[i] * (ni[i] + 2)) + 
                    48 * (1 - ri[i]^2)^4/(ni[i] * (ni[i] + 2) * 
                      (ni[i] + 4)))
                }
                if (vtype[i] == "LS") {
                  vi[i] <- (1 - yi[i]^2)^2/(ni[i] - 1)
                }
                if (vtype[i] == "HO") 
                  vi[i] <- (1 - mnwyi^2)^2/(ni[i] - 1)
            }
            vi[ni <= 4] <- NA
        }
        if (measure == "ZCOR") {
            yi <- 1/2 * log((1 + ri)/(1 - ri))
            vi <- 1/(ni - 3)
            vi[ni <= 4] <- NA
        }
    }
    if (is.element(measure, c("PR", "PLN", "PLO", "PAS", "PFT"))) {
        mf.xi <- mf[[match("xi", names(mf))]]
        mf.mi <- mf[[match("mi", names(mf))]]
        mf.ni <- mf[[match("ni", names(mf))]]
        xi <- eval(mf.xi, data, enclos = sys.frame(sys.parent()))
        mi <- eval(mf.mi, data, enclos = sys.frame(sys.parent()))
        ni <- eval(mf.ni, data, enclos = sys.frame(sys.parent()))
        if (is.null(mi)) 
            mi <- ni - xi
        if (length(xi) == 0L || length(mi) == 0L) 
            stop("Cannot compute outcomes. Check that all of the required \n  information is specified via the appropriate arguments.")
        if (any(c(xi, mi) < 0, na.rm = TRUE)) 
            stop("One or more counts are negative.")
        ni.u <- xi + mi
        if (to == "all") {
            xi <- xi + add
            mi <- mi + add
        }
        if (to == "only0") {
            id0 <- c(xi == 0L | mi == 0L)
            id0[is.na(id0)] <- FALSE
            xi[id0] <- xi[id0] + add
            mi[id0] <- mi[id0] + add
        }
        if (to == "if0all") {
            id0 <- c(xi == 0L | mi == 0L)
            id0[is.na(id0)] <- FALSE
            if (any(id0)) {
                xi <- xi + add
                mi <- mi + add
            }
        }
        ni <- xi + mi
        pri <- xi/ni
        if (measure == "PR") {
            yi <- pri
            k <- length(yi)
            if (length(vtype) == 1L) 
                vtype <- rep(vtype, k)
            vi <- rep(NA, k)
            mnwpri <- sum(xi, na.rm = TRUE)/sum(ni, na.rm = TRUE)
            for (i in seq.int(k)) {
                if (vtype[i] == "UB") 
                  vi[i] <- pri[i] * (1 - pri[i])/(ni[i] - 1)
                if (vtype[i] == "LS") 
                  vi[i] <- pri[i] * (1 - pri[i])/ni[i]
                if (vtype[i] == "HO") 
                  vi[i] <- mnwpri * (1 - mnwpri)/ni[i]
            }
        }
        if (measure == "PLN") {
            yi <- log(pri)
            k <- length(yi)
            if (length(vtype) == 1L) 
                vtype <- rep(vtype, k)
            vi <- rep(NA, k)
            mnwpri <- sum(xi, na.rm = TRUE)/sum(ni, na.rm = TRUE)
            for (i in seq.int(k)) {
                if (vtype[i] == "LS") {
                  vi[i] <- 1/xi[i] - 1/ni[i]
                }
                if (vtype[i] == "HO") 
                  vi[i] <- 1/(mnwpri * ni[i]) - 1/ni[i]
            }
        }
        if (measure == "PLO") {
            yi <- log(pri/(1 - pri))
            k <- length(yi)
            if (length(vtype) == 1L) 
                vtype <- rep(vtype, k)
            vi <- rep(NA, k)
            mnwpri <- sum(xi, na.rm = TRUE)/sum(ni, na.rm = TRUE)
            for (i in seq.int(k)) {
                if (vtype[i] == "LS") 
                  vi[i] <- 1/xi[i] + 1/mi[i]
                if (vtype[i] == "HO") 
                  vi[i] <- 1/(mnwpri * ni[i]) + 1/((1 - mnwpri) * 
                    ni[i])
            }
        }
        if (measure == "PAS") {
            yi <- asin(sqrt(pri))
            vi <- 1/(4 * ni)
        }
        if (measure == "PFT") {
            yi <- 1/2 * (asin(sqrt(xi/(ni + 1))) + asin(sqrt((xi + 
                1)/(ni + 1))))
            vi <- 1/(4 * ni + 2)
        }
    }
    if (is.element(measure, c("IR", "IRLN", "IRS", "IRFT"))) {
        mf.xi <- mf[[match("xi", names(mf))]]
        mf.ti <- mf[[match("ti", names(mf))]]
        xi <- eval(mf.xi, data, enclos = sys.frame(sys.parent()))
        ti <- eval(mf.ti, data, enclos = sys.frame(sys.parent()))
        if (length(xi) == 0L || length(ti) == 0L) 
            stop("Cannot compute outcomes. Check that all of the required \n  information is specified via the appropriate arguments.")
        if (any(xi < 0, na.rm = TRUE)) 
            stop("One or more counts are negative.")
        if (any(ti < 0, na.rm = TRUE)) 
            stop("One or more person-times are negative.")
        ni.u <- ti
        if (to == "all") {
            xi <- xi + add
        }
        if (to == "only0") {
            id0 <- c(xi == 0L)
            id0[is.na(id0)] <- FALSE
            xi[id0] <- xi[id0] + add
        }
        if (to == "if0all") {
            id0 <- c(xi == 0L)
            id0[is.na(id0)] <- FALSE
            if (any(id0)) {
                xi <- xi + add
            }
        }
        iri <- xi/ti
        if (measure == "IR") {
            yi <- iri
            vi <- iri/ti
        }
        if (measure == "IRLN") {
            yi <- log(iri)
            vi <- 1/xi
        }
        if (measure == "IRS") {
            yi <- sqrt(iri)
            vi <- 1/(4 * ti)
        }
        if (measure == "IRFT") {
            yi <- 1/2 * (sqrt(iri) + sqrt(iri + 1/ti))
            vi <- 1/(4 * ti)
        }
    }
    if (is.element(measure, c("MN"))) {
        mf.mi <- mf[[match("mi", names(mf))]]
        mf.sdi <- mf[[match("sdi", names(mf))]]
        mf.ni <- mf[[match("ni", names(mf))]]
        mi <- eval(mf.mi, data, enclos = sys.frame(sys.parent()))
        sdi <- eval(mf.sdi, data, enclos = sys.frame(sys.parent()))
        ni <- eval(mf.ni, data, enclos = sys.frame(sys.parent()))
        if (length(mi) == 0L || length(sdi) == 0L || length(ni) == 
            0L) 
            stop("Cannot compute outcomes. Check that all of the required \n  information is specified via the appropriate arguments.")
        if (any(sdi < 0, na.rm = TRUE)) 
            stop("One or more standard deviations are negative.")
        ni.u <- ni
        if (measure == "MN") {
            yi <- mi
            vi <- sdi^2/ni
        }
    }
    if (is.element(measure, c("MC", "SMCC", "SMCR"))) {
        mf.m1i <- mf[[match("m1i", names(mf))]]
        mf.m2i <- mf[[match("m2i", names(mf))]]
        mf.sd1i <- mf[[match("sd1i", names(mf))]]
        mf.sd2i <- mf[[match("sd2i", names(mf))]]
        mf.ni <- mf[[match("ni", names(mf))]]
        mf.ri <- mf[[match("ri", names(mf))]]
        m1i <- eval(mf.m1i, data, enclos = sys.frame(sys.parent()))
        m2i <- eval(mf.m2i, data, enclos = sys.frame(sys.parent()))
        sd1i <- eval(mf.sd1i, data, enclos = sys.frame(sys.parent()))
        sd2i <- eval(mf.sd2i, data, enclos = sys.frame(sys.parent()))
        ni <- eval(mf.ni, data, enclos = sys.frame(sys.parent()))
        ri <- eval(mf.ri, data, enclos = sys.frame(sys.parent()))
        if (is.element(measure, c("MC", "SMCC"))) {
            if (length(m1i) == 0L || length(m2i) == 0L || length(sd1i) == 
                0L || length(sd2i) == 0L || length(ni) == 0L || 
                length(ri) == 0L) 
                stop("Cannot compute outcomes. Check that all of the required \n  information is specified via the appropriate arguments.")
            if (any(c(sd1i, sd2i) < 0, na.rm = TRUE)) 
                stop("One or more standard deviations are negative.")
        }
        else {
            if (length(m1i) == 0L || length(m2i) == 0L || length(sd1i) == 
                0L || length(ni) == 0L || length(ri) == 0L) 
                stop("Cannot compute outcomes. Check that all of the required \n  information is specified via the appropriate arguments.")
            if (any(sd1i < 0, na.rm = TRUE)) 
                stop("One or more standard deviations are negative.")
        }
        if (any(abs(ri) > 1, na.rm = TRUE)) 
            stop("One or more correlations are > 1 or < -1.")
        ni.u <- ni
        if (measure == "MC") {
            yi <- m1i - m2i
            vi <- (sd1i^2 + sd2i^2 - 2 * ri * sd1i * sd2i)/ni
        }
        if (measure == "SMCC") {
            mi <- ni - 1
            warn.before <- getOption("warn")
            options(warn = -1)
            cmi <- .cmicalc(mi)
            options(warn = warn.before)
            sddi <- sqrt(sd1i^2 + sd2i^2 - 2 * ri * sd1i * sd2i)
            yi <- cmi * (m1i - m2i)/sddi
            vi <- 1/ni + yi^2/(2 * ni)
        }
        if (measure == "SMCR") {
            mi <- ni - 1
            warn.before <- getOption("warn")
            options(warn = -1)
            cmi <- .cmicalc(mi)
            options(warn = warn.before)
            yi <- cmi * (m1i - m2i)/sd1i
            vi <- 2 * (1 - ri)/ni + yi^2/(2 * ni)
        }
    }
    if (is.element(measure, c("ARAW", "AHW", "ABT"))) {
        mf.ai <- mf[[match("ai", names(mf))]]
        mf.mi <- mf[[match("mi", names(mf))]]
        mf.ni <- mf[[match("ni", names(mf))]]
        ai <- eval(mf.ai, data, enclos = sys.frame(sys.parent()))
        mi <- eval(mf.mi, data, enclos = sys.frame(sys.parent()))
        ni <- eval(mf.ni, data, enclos = sys.frame(sys.parent()))
        if (length(ai) == 0L || length(mi) == 0L || length(ni) == 
            0L) 
            stop("Cannot compute outcomes. Check that all of the required \n  information is specified via the appropriate arguments.")
        if (any(ai > 1, na.rm = TRUE)) 
            stop("One or more alphas are > 1.")
        if (any(mi < 2, na.rm = TRUE)) 
            stop("One or more mi's are < 2.")
        ni.u <- ni
        if (measure == "ARAW") {
            yi <- ai
            vi <- 2 * mi * (1 - ai)^2/((mi - 1) * (ni - 2))
        }
        if (measure == "AHW") {
            yi <- 1 - (1 - ai)^(1/3)
            vi <- 18 * mi * (ni - 1) * (1 - ai)^(2/3)/((mi - 
                1) * (9 * ni - 11)^2)
        }
        if (measure == "ABT") {
            yi <- -log(1 - ai)
            vi <- 2 * mi/((mi - 1) * (ni - 2))
        }
    }
    if (any(is.infinite(c(yi, vi)))) {
        warning("Some yi and/or vi values equal to +-Inf. Recoded to NAs.")
        is.inf <- is.infinite(yi) | is.infinite(vi)
        yi[is.inf] <- NA
        vi[is.inf] <- NA
    }
    if (any(is.nan(c(yi, vi)))) {
        is.NaN <- is.nan(yi) | is.nan(vi)
        yi[is.NaN] <- NA
        vi[is.NaN] <- NA
    }
    vi[vi < 0] <- NA
    if (!is.null(slab)) {
        if (any(duplicated(slab))) 
            stop("Study labels must be unique.")
        if (length(slab) != length(yi)) 
            stop("Study labels not of same length as data.")
        attr(yi, "slab") <- slab
    }
    if (!is.null(subset)) {
        yi <- yi[subset]
        vi <- vi[subset]
        ni.u <- ni.u[subset]
        slab <- slab[subset]
        data <- data[subset, ]
    }
    attr(yi, "measure") <- measure
    if (!no.data && append) {
        dat <- data.frame(data)
        if (replace) {
            dat[[var.names[1]]] <- yi
            dat[[var.names[2]]] <- vi
            attr(dat[[var.names[1]]], "ni") <- ni.u
        }
        else {
            if (is.element(var.names[1], names(dat))) {
                is.na.yi <- is.na(dat[[var.names[1]]])
                dat[[var.names[1]]][is.na.yi] <- yi[is.na.yi]
                attributes(dat[[var.names[1]]])$ni[is.na.yi] <- ni.u[is.na.yi]
            }
            else {
                dat[[var.names[1]]] <- yi
                attr(dat[[var.names[1]]], "ni") <- ni.u
            }
            if (is.element(var.names[2], names(dat))) {
                is.na.vi <- is.na(dat[[var.names[2]]])
                dat[[var.names[2]]][is.na.vi] <- vi[is.na.vi]
            }
            else {
                dat[[var.names[2]]] <- vi
            }
        }
    }
    else {
        dat <- data.frame(yi, vi)
        attr(dat$yi, "ni") <- ni.u
        names(dat) <- var.names
    }
    attr(dat, "var.names") <- var.names
    attr(dat, "digits") <- digits
    class(dat) <- c("escalc", "data.frame")
    return(dat)
}
escalc.formula <-
function (measure, formula, weights, data, add = 1/2, to = "only0", 
    drop00 = FALSE, vtype = "LS", var.names = c("yi", "vi"), 
    digits = 4, ...) 
{
    if (!is.element(measure, c("RR", "OR", "PETO", "RD", "AS", 
        "PHI", "YUQ", "YUY", "RTET", "PBIT", "OR2D", "IRR", "IRD", 
        "IRSD", "MD", "SMD", "SMDH", "ROM", "RPB", "RBIS", "D2OR", 
        "COR", "UCOR", "ZCOR", "PR", "PLN", "PLO", "PAS", "PFT", 
        "IR", "IRLN", "IRS", "IRFT", "MN", "MC", "SMCC", "SMCR", 
        "ARAW", "AHW", "ABT"))) 
        stop("Unknown 'measure' specified.")
    if (is.element(measure, c("MC", "SMCC", "SMCR"))) 
        stop("Formula interface not (currently) implemented for this measure.")
    na.act <- getOption("na.action")
    options(na.action = "na.pass")
    on.exit(options(na.action = na.act))
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "weights"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    formula <- as.Formula(formula)
    if (length(formula)[2] < 2L) 
        stop("Right-hand side of formula must specify both a grouping and a study factor (i.e., ~ group | study).")
    mf$formula <- formula
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    weights <- model.weights(mf)
    lhs <- model.part(formula, data = mf, lhs = 1)
    rhs1 <- model.part(formula, data = mf, rhs = 1)
    study <- model.part(formula, data = mf, rhs = 2)
    if (length(study) != 1) 
        stop("A single study factor must be specified.")
    if (!is.factor(study[[1]])) 
        stop("Study variable must be a factor.")
    study <- study[[1]]
    if (any(is.na(study))) 
        stop("Study factor must not contain NAs.")
    if (is.element(measure, c("RR", "OR", "RD", "AS", "PETO", 
        "PHI", "YUQ", "YUY", "RTET", "PBIT", "OR2D"))) {
        if (is.null(weights)) 
            stop("Must specify the 'weights' argument.")
        if (length(lhs) != 1) 
            stop("Left-hand side of formula must be a single outcome factor.")
        outcome <- lhs[[1]]
        if (!is.factor(outcome)) 
            stop("Left-hand side of formula must be a factor.")
        if (nlevels(outcome) != 2) 
            stop("Outcome factor on left-hand side of formula should have two levels.")
        if (length(rhs1) != 1) 
            stop("A single grouping factor must be specified.")
        if (!is.factor(rhs1[[1]])) 
            stop("Grouping variable must be a factor.")
        group <- rhs1[[1]]
        if (nlevels(group) != 2) 
            stop("Grouping factor should have two levels.")
        if (any(is.na(group)) || any(is.na(outcome))) 
            stop("Grouping and outcome factors must not contain NAs.")
        ai <- weights[group == levels(group)[1] & outcome == 
            levels(outcome)[1]]
        bi <- weights[group == levels(group)[1] & outcome == 
            levels(outcome)[2]]
        ci <- weights[group == levels(group)[2] & outcome == 
            levels(outcome)[1]]
        di <- weights[group == levels(group)[2] & outcome == 
            levels(outcome)[2]]
        names(ai) <- mf$study[group == levels(group)[1] & outcome == 
            levels(outcome)[1]]
        names(bi) <- mf$study[group == levels(group)[1] & outcome == 
            levels(outcome)[2]]
        names(ci) <- mf$study[group == levels(group)[2] & outcome == 
            levels(outcome)[1]]
        names(di) <- mf$study[group == levels(group)[2] & outcome == 
            levels(outcome)[2]]
        return(escalc(measure = measure, ai = ai, bi = bi, ci = ci, 
            di = di, add = add, to = to, drop00 = drop00, vtype = vtype, 
            var.names = var.names, append = "FALSE", digits = digits))
    }
    if (is.element(measure, c("IRR", "IRD", "IRSD"))) {
        if (length(lhs) != 2) 
            stop("Left-hand side of formula must specify the number of events and the total person-time at risk (i.e., events/times ~).")
        events <- lhs[, 1]
        times <- lhs[, 2]
        if (!is.vector(events) || !is.vector(times)) 
            stop("The events and person-time at risk variables should be vectors.")
        if (length(rhs1) != 1) 
            stop("A single grouping factor must be specified.")
        if (!is.factor(rhs1[[1]])) 
            stop("Grouping variable must be a factor.")
        group <- rhs1[[1]]
        if (nlevels(group) != 2) 
            stop("Grouping factor should have two levels.")
        if (any(is.na(group))) 
            stop("Grouping factor must not contain NAs.")
        x1i <- events[group == levels(group)[1]]
        x2i <- events[group == levels(group)[2]]
        t1i <- times[group == levels(group)[1]]
        t2i <- times[group == levels(group)[2]]
        names(x1i) <- mf$study[group == levels(group)[1]]
        names(x2i) <- mf$study[group == levels(group)[2]]
        return(escalc(measure = measure, x1i = x1i, x2i = x2i, 
            t1i = t1i, t2i = t2i, add = add, to = to, drop00 = drop00, 
            vtype = vtype, var.names = var.names, append = "FALSE", 
            digits = digits))
    }
    if (is.element(measure, c("MD", "SMD", "SMDH", "ROM", "RPB", 
        "RBIS", "D2OR"))) {
        if (is.null(weights)) 
            stop("Must specify the 'weights' argument.")
        if (length(lhs) != 2) 
            stop("Left-hand side of formula must specify the means and standard devations (i.e., means/sds ~).")
        means <- lhs[, 1]
        sds <- lhs[, 2]
        if (!is.vector(means) || !is.vector(sds)) 
            stop("The mean and standard devation variables should be vectors.")
        if (length(rhs1) != 1) 
            stop("A single grouping factor must be specified.")
        if (!is.factor(rhs1[[1]])) 
            stop("Grouping variable must be a factor.")
        group <- rhs1[[1]]
        if (nlevels(group) != 2) 
            stop("Grouping factor should have two levels.")
        if (any(is.na(group))) 
            stop("Grouping factor must not contain NAs.")
        m1i <- means[group == levels(group)[1]]
        m2i <- means[group == levels(group)[2]]
        sd1i <- sds[group == levels(group)[1]]
        sd2i <- sds[group == levels(group)[2]]
        n1i <- weights[group == levels(group)[1]]
        n2i <- weights[group == levels(group)[2]]
        names(m1i) <- mf$study[group == levels(group)[1]]
        names(m2i) <- mf$study[group == levels(group)[2]]
        return(escalc(measure = measure, m1i = m1i, m2i = m2i, 
            sd1i = sd1i, sd2i = sd2i, n1i = n1i, n2i = n2i, vtype = vtype, 
            var.names = var.names, append = "FALSE", digits = digits))
    }
    if (is.element(measure, c("COR", "UCOR", "ZCOR"))) {
        if (is.null(weights)) 
            stop("Must specify the 'weights' argument.")
        if (length(lhs) != 1) 
            stop("Left-hand side of formula must specify the correlations (i.e., cors ~).")
        ri <- lhs[[1]]
        if (!is.vector(ri)) 
            stop("The variable specifying the correlation should be a vector.")
        ni <- weights
        names(ri) <- mf$study
        return(escalc(measure = measure, ri = ri, ni = ni, vtype = vtype, 
            var.names = var.names, append = "FALSE", digits = digits))
    }
    if (is.element(measure, c("PR", "PLN", "PLO", "PAS", "PFT"))) {
        if (length(lhs) != 1) 
            stop("Left-hand side of formula must be a single outcome factor.")
        outcome <- lhs[[1]]
        if (!is.factor(outcome)) 
            stop("Left-hand side of formula must be a factor.")
        if (nlevels(outcome) != 2) 
            stop("Outcome factor on left-hand side of formula should have two levels.")
        if (any(is.na(outcome))) 
            stop("Outcome factor must not contain NAs.")
        xi <- weights[outcome == levels(outcome)[1]]
        mi <- weights[outcome == levels(outcome)[2]]
        names(xi) <- mf$study[outcome == levels(outcome)[1]]
        names(mi) <- mf$study[outcome == levels(outcome)[2]]
        return(escalc(measure = measure, xi = xi, mi = mi, add = add, 
            to = to, vtype = vtype, var.names = var.names, append = "FALSE", 
            digits = digits))
    }
    if (is.element(measure, c("IR", "IRLN", "IRS", "IRFT"))) {
        if (length(lhs) != 2) 
            stop("Left-hand side of formula must specify the number of cases and the total person-time at risk (i.e., cases/times ~).")
        events <- lhs[, 1]
        times <- lhs[, 2]
        if (!is.vector(events) || !is.vector(times)) 
            stop("The events and person-time at risk variables should be vectors.")
        xi <- events
        ti <- times
        names(xi) <- mf$study
        return(escalc(measure = measure, xi = xi, ti = ti, add = add, 
            to = to, vtype = vtype, var.names = var.names, append = "FALSE", 
            digits = digits))
    }
    if (is.element(measure, c("MN"))) {
        if (is.null(weights)) 
            stop("Must specify the 'weights' argument.")
        if (length(lhs) != 2) 
            stop("Left-hand side of formula must specify the means and standard devations (i.e., means/sds ~).")
        means <- lhs[, 1]
        sds <- lhs[, 2]
        if (!is.vector(means) || !is.vector(sds)) 
            stop("The mean and standard devation variables should be vectors.")
        mi <- means
        sdi <- sds
        ni <- weights
        names(mi) <- mf$study
        return(escalc(measure = measure, mi = mi, sdi = sdi, 
            ni = ni, vtype = vtype, var.names = var.names, append = "FALSE", 
            digits = digits))
    }
    if (is.element(measure, c("MC", "SMCC", "SMCR"))) {
        stop("Formula interface (currently) not implemented for these outcome measures.")
    }
    if (is.element(measure, c("ARAW", "AHW", "ABT"))) {
        if (is.null(weights)) 
            stop("Must specify the 'weights' argument.")
        if (length(lhs) != 2) 
            stop("Left-hand side of formula must specify the alpha values and number of items (i.e., alphas/items ~).")
        alphas <- lhs[, 1]
        items <- lhs[, 2]
        if (!is.vector(alphas) || !is.vector(items)) 
            stop("The alpha and item variables should be vectors.")
        ai <- alphas
        mi <- items
        ni <- weights
        names(ai) <- mf$study
        return(escalc(measure = measure, ai = ai, mi = mi, ni = ni, 
            vtype = vtype, var.names = var.names, append = "FALSE", 
            digits = digits))
    }
}
fitstats <-
function (object, ...) 
UseMethod("fitstats")
fitstats.rma <-
function (object, REML, ...) 
{
    if (!is.element("rma", class(object))) 
        stop("Argument 'object' must be an object of class \"rma\".")
    if (missing(REML)) {
        if (object$method == "REML") {
            REML <- TRUE
        }
        else {
            REML <- FALSE
        }
    }
    if (REML) {
        out <- cbind(object$fit.stats$REML)
        dimnames(out)[[1]] <- c("logLik:", "deviance:", "AIC:", 
            "BIC:")
        dimnames(out)[[2]] <- c("REML")
    }
    else {
        out <- cbind(object$fit.stats$ML)
        dimnames(out)[[1]] <- c("logLik:", "deviance:", "AIC:", 
            "BIC:")
        dimnames(out)[[2]] <- c("ML")
    }
    return(out)
}
fitted.rma <-
function (object, ...) 
{
    if (!is.element("rma", class(object))) 
        stop("Argument 'object' must be an object of class \"rma\".")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    out <- c(object$X.f %*% object$b)
    names(out) <- object$slab
    not.na <- !is.na(out)
    if (na.act == "na.omit") 
        out <- out[not.na]
    if (na.act == "na.fail" && any(!not.na)) 
        stop("Missing values in results.")
    return(out)
}
forest <-
function (x, ...) 
UseMethod("forest")
forest.cumul.rma <-
function (x, annotate = TRUE, xlim, alim, ylim, at, steps = 5, 
    level = x$level, digits = 2, refline = 0, xlab, ilab, ilab.xpos, 
    ilab.pos, transf = FALSE, atransf = FALSE, targs, rows, efac = 1, 
    pch = 15, psize = 1, cex, cex.lab, cex.axis, ...) 
{
    na.act <- getOption("na.action")
    if (!is.element("cumul.rma", class(x))) 
        stop("Argument 'x' must be an object of class \"cumul.rma\".")
    transf.char <- deparse(substitute(transf))
    atransf.char <- deparse(substitute(atransf))
    if (transf.char != "FALSE" && atransf.char != "FALSE") 
        stop("Use either 'transf' or 'atransf' to specify a transformation (not both).")
    if (missing(targs)) 
        targs <- NULL
    if (missing(at)) 
        at <- NULL
    if (missing(ilab)) 
        ilab <- NULL
    if (missing(ilab.xpos)) 
        ilab.xpos <- NULL
    if (missing(ilab.pos)) 
        ilab.pos <- NULL
    if (missing(cex)) 
        cex <- NULL
    if (missing(cex.lab)) 
        cex.lab <- NULL
    if (missing(cex.axis)) 
        cex.axis <- NULL
    if (length(digits) == 1L) 
        digits <- c(digits, digits)
    alpha <- (100 - level)/100
    yi <- x$estimate
    if (is.null(attr(yi, "measure"))) {
        measure <- "GEN"
    }
    else {
        measure <- attr(yi, "measure")
    }
    vi <- x$se^2
    if (length(yi) != length(vi)) 
        stop("Length of yi and vi (or sei) vectors is not the same.")
    k <- length(yi)
    if (x$slab.null) {
        slab <- paste("+ Study ", x$slab)
        slab[1] <- paste("Study ", x$slab[1])
    }
    else {
        slab <- paste("+", x$slab)
        slab[1] <- paste(x$slab[1])
    }
    if (is.vector(ilab)) 
        ilab <- cbind(ilab)
    if (length(pch) == 1L) 
        pch <- rep(pch, k)
    if (length(pch) != length(yi)) 
        stop("Number of outcomes does not correspond to the length of the pch argument.")
    if (!is.null(psize)) {
        if (length(psize) == 1L) 
            psize <- rep(psize, k)
        if (length(psize) != length(yi)) 
            stop("Number of outcomes does not correspond to the length of the psize argument.")
    }
    if (missing(rows)) {
        rows <- k:1
    }
    else {
        if (length(rows) == 1L) 
            rows <- rows:(rows - k + 1)
    }
    if (length(rows) != length(yi)) 
        stop("Number of outcomes does not correspond to the length of the rows argument.")
    yi <- yi[k:1]
    vi <- vi[k:1]
    slab <- slab[k:1]
    ilab <- ilab[k:1, , drop = FALSE]
    pch <- pch[k:1]
    psize <- psize[k:1]
    rows <- rows[k:1]
    yivi.na <- is.na(cbind(yi, vi))
    if (any(yivi.na)) {
        not.na <- apply(yivi.na, MARGIN = 1, sum) == 0L
        if (na.act == "na.omit") {
            yi <- yi[not.na]
            vi <- vi[not.na]
            slab <- slab[not.na]
            ilab <- ilab[not.na, , drop = FALSE]
            pch <- pch[not.na]
            psize <- psize[not.na]
            rows.new <- rows
            rows.na <- rows[!not.na]
            for (j in seq.int(length(rows.na))) {
                rows.new[rows >= rows.na[j]] <- rows.new[rows >= 
                  rows.na[j]] - 1
            }
            rows <- rows.new[not.na]
        }
        if (na.act == "na.fail") 
            stop("Missing values in results.")
    }
    k <- length(yi)
    ci.lb <- yi - qnorm(alpha/2, lower.tail = FALSE) * sqrt(vi)
    ci.ub <- yi + qnorm(alpha/2, lower.tail = FALSE) * sqrt(vi)
    if (is.function(transf)) {
        if (is.null(targs)) {
            yi <- sapply(yi, transf)
            ci.lb <- sapply(ci.lb, transf)
            ci.ub <- sapply(ci.ub, transf)
        }
        else {
            yi <- sapply(yi, transf, targs)
            ci.lb <- sapply(ci.lb, transf, targs)
            ci.ub <- sapply(ci.ub, transf, targs)
        }
    }
    ci.bounds <- cbind(ci.lb, ci.ub)
    rev.order <- ifelse(ci.ub < ci.lb, TRUE, FALSE)
    rev.order[is.na(rev.order)] <- FALSE
    ci.bounds[rev.order] <- ci.bounds[rev.order, 2:1]
    ci.lb <- ci.bounds[, 1]
    ci.ub <- ci.bounds[, 2]
    if (is.null(psize)) {
        if (any(vi <= 0, na.rm = TRUE)) {
            psize <- rep(1, k)
        }
        else {
            wi <- 1/sqrt(vi)
            psize <- wi/sum(wi, na.rm = TRUE)
            psize <- (psize - min(psize, na.rm = TRUE))/(max(psize, 
                na.rm = TRUE) - min(psize, na.rm = TRUE))
            psize <- (psize * 1) + 0.5
            if (all(is.na(psize))) 
                psize <- rep(1, k)
        }
    }
    rng <- max(ci.ub, na.rm = TRUE) - min(ci.lb, na.rm = TRUE)
    if (annotate) {
        plot.multp.l <- 1.2
        plot.multp.r <- 1.2
        axis.multp.l <- 0.2
        axis.multp.r <- 0.2
    }
    else {
        plot.multp.l <- 1.2
        plot.multp.r <- 0.4
        axis.multp.l <- 0.2
        axis.multp.r <- 0.2
    }
    if (missing(xlim)) {
        xlim <- c(min(ci.lb, na.rm = TRUE) - rng * plot.multp.l, 
            max(ci.ub, na.rm = TRUE) + rng * plot.multp.r)
        xlim <- round(xlim, digits[2])
    }
    alim.spec <- TRUE
    if (missing(alim)) {
        if (is.null(at)) {
            alim <- range(pretty(x = c(min(ci.lb, na.rm = TRUE), 
                max(ci.ub, na.rm = TRUE)), n = steps - 1))
            alim.spec <- FALSE
        }
        else {
            alim <- range(at)
        }
    }
    alim <- sort(alim)
    xlim <- sort(xlim)
    if (xlim[1] > min(yi, na.rm = TRUE)) {
        xlim[1] <- min(yi, na.rm = TRUE)
    }
    if (xlim[2] < max(yi, na.rm = TRUE)) {
        xlim[2] <- max(yi, na.rm = TRUE)
    }
    if (alim[1] > min(yi, na.rm = TRUE)) {
        alim[1] <- min(yi, na.rm = TRUE)
    }
    if (alim[2] < max(yi, na.rm = TRUE)) {
        alim[2] <- max(yi, na.rm = TRUE)
    }
    if (alim[1] < xlim[1]) {
        xlim[1] <- alim[1]
    }
    if (alim[2] > xlim[2]) {
        xlim[2] <- alim[2]
    }
    if (missing(ylim)) {
        ylim <- c(0.5, k + 3)
    }
    else {
        ylim <- sort(ylim)
    }
    if (is.null(at)) {
        if (alim.spec) {
            at <- seq.int(from = alim[1], to = alim[2], length.out = steps)
        }
        else {
            at <- pretty(x = c(min(ci.lb, na.rm = TRUE), max(ci.ub, 
                na.rm = TRUE)), n = steps - 1)
        }
    }
    else {
        at[at < alim[1]] <- alim[1]
        at[at > alim[2]] <- alim[2]
        at <- unique(at)
    }
    at.lab <- at
    if (is.function(atransf)) {
        if (is.null(targs)) {
            at.lab <- formatC(sapply(at.lab, atransf), digits = digits[2], 
                format = "f")
        }
        else {
            at.lab <- formatC(sapply(at.lab, atransf, targs), 
                digits = digits[2], format = "f")
        }
    }
    else {
        at.lab <- formatC(at.lab, digits = digits[2], format = "f")
    }
    par.mar <- par("mar")
    par.mar.adj <- par.mar - c(0, 3, 1, 1)
    par.mar.adj[par.mar.adj < 0] <- 0
    par(mar = par.mar.adj)
    on.exit(par(mar = par.mar))
    plot(NA, NA, xlim = xlim, ylim = ylim, xlab = "", ylab = "", 
        yaxt = "n", xaxt = "n", xaxs = "i", bty = "n", ...)
    abline(h = ylim[2] - 2, ...)
    par.usr <- par("usr")
    height <- par.usr[4] - par.usr[3]
    lheight <- strheight("O")
    cex.adj <- ifelse(k * lheight > height * 0.8, height/(1.25 * 
        k * lheight), 1)
    if (is.null(cex)) {
        cex <- par("cex") * cex.adj
    }
    else {
        if (is.null(cex.lab)) 
            cex.lab <- cex
        if (is.null(cex.axis)) 
            cex.axis <- cex
    }
    if (is.null(cex.lab)) 
        cex.lab <- par("cex") * cex.adj
    if (is.null(cex.axis)) 
        cex.axis <- par("cex") * cex.adj
    axis(side = 1, at = at, labels = at.lab, cex.axis = cex.axis, 
        ...)
    if (missing(xlab)) {
        xlab <- .setxlab(measure, transf.char, atransf.char, 
            gentype = 2)
    }
    else {
        mtext(xlab, side = 1, at = min(at) + (max(at) - min(at))/2, 
            line = par("mgp")[1] - 0.5, cex = cex.lab, ...)
    }
    for (i in seq.int(k)) {
        if (is.na(yi[i]) || is.na(vi)[i]) 
            next
        segments(max(ci.lb[i], alim[1]), rows[i], min(ci.ub[i], 
            alim[2]), rows[i], ...)
        if (ci.lb[i] >= alim[1]) {
            segments(ci.lb[i], rows[i] - (height/150) * cex * 
                efac, ci.lb[i], rows[i] + (height/150) * cex * 
                efac, ...)
        }
        else {
            polygon(x = c(alim[1], alim[1] + (1.4/100) * cex * 
                (xlim[2] - xlim[1]), alim[1] + (1.4/100) * cex * 
                (xlim[2] - xlim[1]), alim[1]), y = c(rows[i], 
                rows[i] + (height/150) * cex * efac, rows[i] - 
                  (height/150) * cex * efac, rows[i]), col = "black", 
                ...)
        }
        if (ci.ub[i] <= alim[2]) {
            segments(ci.ub[i], rows[i] - (height/150) * cex * 
                efac, ci.ub[i], rows[i] + (height/150) * cex * 
                efac, ...)
        }
        else {
            polygon(x = c(alim[2], alim[2] - (1.4/100) * cex * 
                (xlim[2] - xlim[1]), alim[2] - (1.4/100) * cex * 
                (xlim[2] - xlim[1]), alim[2]), y = c(rows[i], 
                rows[i] + (height/150) * cex * efac, rows[i] - 
                  (height/150) * cex * efac, rows[i]), col = "black", 
                ...)
        }
    }
    text(xlim[1], rows, slab, pos = 4, cex = cex, ...)
    if (!is.null(ilab)) {
        if (is.null(ilab.xpos)) 
            stop("Must specify 'ilab.xpos' argument when adding information with 'ilab'.")
        if (length(ilab.xpos) != NCOL(ilab)) 
            stop("Number of 'ilab' columns does not match length of 'ilab.xpos' argument.")
        for (l in seq.int(NCOL(ilab))) {
            text(ilab.xpos[l], rows, ilab[, l], pos = ilab.pos[l], 
                cex = cex, ...)
        }
    }
    if (is.numeric(refline)) 
        segments(refline, ylim[1] - 5, refline, ylim[2] - 2, 
            lty = "dotted", ...)
    if (annotate) {
        if (is.function(atransf)) {
            if (is.null(targs)) {
                annotext <- round(cbind(sapply(yi, atransf), 
                  sapply(ci.lb, atransf), sapply(ci.ub, atransf)), 
                  digits[1])
            }
            else {
                annotext <- round(cbind(sapply(yi, atransf, targs), 
                  sapply(ci.lb, atransf, targs), sapply(ci.ub, 
                    atransf, targs)), digits[1])
            }
            rev.order <- ifelse(annotext[, 3] < annotext[, 2], 
                TRUE, FALSE)
            rev.order[is.na(rev.order)] <- FALSE
            annotext[rev.order, 2:3] <- annotext[rev.order, 3:2]
        }
        else {
            annotext <- round(cbind(yi, ci.lb, ci.ub), digits[1])
        }
        annotext <- matrix(apply(annotext, 2, format, nsmall = digits[1]), 
            ncol = 3)
        annotext <- cbind(annotext[, 1], " [ ", annotext[, 2], 
            " , ", annotext[, 3], " ]")
        annotext <- apply(annotext, 1, paste, collapse = "")
        text(x = xlim[2], rows, labels = annotext, pos = 2, cex = cex, 
            ...)
    }
    points(yi, rows, pch = pch, cex = cex * psize, ...)
    invisible()
}
forest.default <-
function (x, vi, sei, ci.lb, ci.ub, annotate = TRUE, showweight = FALSE, 
    xlim, alim, ylim, at, steps = 5, level = 95, digits = 2, 
    refline = 0, xlab, slab, ilab, ilab.xpos, ilab.pos, subset, 
    transf = FALSE, atransf = FALSE, targs, rows, efac = 1, pch = 15, 
    psize, cex, cex.lab, cex.axis, ...) 
{
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    transf.char <- deparse(substitute(transf))
    atransf.char <- deparse(substitute(atransf))
    if (transf.char != "FALSE" && atransf.char != "FALSE") 
        stop("Use either 'transf' or 'atransf' to specify a transformation (not both).")
    if (missing(targs)) 
        targs <- NULL
    if (missing(at)) 
        at <- NULL
    if (missing(ilab)) 
        ilab <- NULL
    if (missing(ilab.xpos)) 
        ilab.xpos <- NULL
    if (missing(ilab.pos)) 
        ilab.pos <- NULL
    if (missing(subset)) 
        subset <- NULL
    if (missing(psize)) 
        psize <- NULL
    if (missing(cex)) 
        cex <- NULL
    if (missing(cex.lab)) 
        cex.lab <- NULL
    if (missing(cex.axis)) 
        cex.axis <- NULL
    if (length(digits) == 1L) 
        digits <- c(digits, digits)
    alpha <- (100 - level)/100
    yi <- x
    if (is.null(attr(yi, "measure"))) {
        measure <- "GEN"
    }
    else {
        measure <- attr(yi, "measure")
    }
    if (hasArg(ci.lb) && hasArg(ci.ub)) {
        if (length(ci.lb) != length(ci.ub)) 
            stop("Length of ci.lb and ci.ub do not match.")
        if (missing(vi) && missing(sei)) {
            vi <- ((ci.ub - ci.lb)/(2 * qnorm(alpha/2, lower.tail = FALSE)))^2
        }
        else {
            if (missing(vi)) 
                vi <- sei^2
        }
        if (length(ci.lb) != length(vi)) 
            stop("Length of vi (or sei) does not match length of (ci.lb, ci.ub) pairs.")
    }
    else {
        if (missing(vi)) {
            if (missing(sei)) {
                stop("Must specify either vi, sei, or (ci.lb, ci.ub) pairs.")
            }
            else {
                vi <- sei^2
                ci.lb <- yi - qnorm(alpha/2, lower.tail = FALSE) * 
                  sei
                ci.ub <- yi + qnorm(alpha/2, lower.tail = FALSE) * 
                  sei
            }
        }
        else {
            ci.lb <- yi - qnorm(alpha/2, lower.tail = FALSE) * 
                sqrt(vi)
            ci.ub <- yi + qnorm(alpha/2, lower.tail = FALSE) * 
                sqrt(vi)
        }
    }
    if (length(yi) != length(vi)) 
        stop("Length of yi does not match the length of vi, sei, or the (ci.lb, ci.ub) pairs.")
    k <- length(yi)
    if (missing(slab)) {
        if (!is.null(attr(yi, "slab"))) {
            slab <- attr(yi, "slab")
        }
        else {
            slab <- paste("Study ", seq.int(k))
        }
    }
    if (length(yi) != length(slab)) 
        stop("Number of outcomes does not correspond to the length of the slab argument.")
    if (is.vector(ilab)) 
        ilab <- cbind(ilab)
    if (length(pch) == 1L) 
        pch <- rep(pch, k)
    if (length(pch) != length(yi)) 
        stop("Number of outcomes does not correspond to the length of the pch argument.")
    if (!is.null(psize)) {
        if (length(psize) == 1L) 
            psize <- rep(psize, k)
        if (length(psize) != length(yi)) 
            stop("Number of outcomes does not correspond to the length of the psize argument.")
    }
    if (!is.null(subset)) {
        yi <- yi[subset]
        vi <- vi[subset]
        ci.lb <- ci.lb[subset]
        ci.ub <- ci.ub[subset]
        slab <- slab[subset]
        ilab <- ilab[subset, , drop = FALSE]
        pch <- pch[subset]
        psize <- psize[subset]
    }
    k <- length(yi)
    if (missing(rows)) {
        rows <- k:1
    }
    else {
        if (length(rows) == 1L) 
            rows <- rows:(rows - k + 1)
    }
    if (length(rows) != length(yi)) 
        stop("Number of outcomes does not correspond to the length of the rows argument.")
    yi <- yi[k:1]
    vi <- vi[k:1]
    ci.lb <- ci.lb[k:1]
    ci.ub <- ci.ub[k:1]
    slab <- slab[k:1]
    ilab <- ilab[k:1, , drop = FALSE]
    pch <- pch[k:1]
    psize <- psize[k:1]
    rows <- rows[k:1]
    yivi.na <- is.na(cbind(yi, vi))
    if (any(yivi.na)) {
        not.na <- apply(yivi.na, MARGIN = 1, sum) == 0L
        if (na.act == "na.omit") {
            yi <- yi[not.na]
            vi <- vi[not.na]
            ci.lb <- ci.lb[not.na]
            ci.ub <- ci.ub[not.na]
            slab <- slab[not.na]
            ilab <- ilab[not.na, , drop = FALSE]
            pch <- pch[not.na]
            psize <- psize[not.na]
            rows.new <- rows
            rows.na <- rows[!not.na]
            for (j in seq.int(length(rows.na))) {
                rows.new[rows >= rows.na[j]] <- rows.new[rows >= 
                  rows.na[j]] - 1
            }
            rows <- rows.new[not.na]
        }
        if (na.act == "na.fail") 
            stop("Missing values in results.")
    }
    k <- length(yi)
    if (is.function(transf)) {
        if (is.null(targs)) {
            yi <- sapply(yi, transf)
            ci.lb <- sapply(ci.lb, transf)
            ci.ub <- sapply(ci.ub, transf)
        }
        else {
            yi <- sapply(yi, transf, targs)
            ci.lb <- sapply(ci.lb, transf, targs)
            ci.ub <- sapply(ci.ub, transf, targs)
        }
    }
    ci.bounds <- cbind(ci.lb, ci.ub)
    rev.order <- ifelse(ci.ub < ci.lb, TRUE, FALSE)
    rev.order[is.na(rev.order)] <- FALSE
    ci.bounds[rev.order] <- ci.bounds[rev.order, 2:1]
    ci.lb <- ci.bounds[, 1]
    ci.ub <- ci.bounds[, 2]
    if (showweight) {
        weights <- 1/vi
        weights <- 100 * weights/sum(weights, na.rm = TRUE)
    }
    if (is.null(psize)) {
        if (any(vi <= 0, na.rm = TRUE)) {
            psize <- rep(1, k)
        }
        else {
            wi <- 1/sqrt(vi)
            psize <- wi/sum(wi, na.rm = TRUE)
            psize <- (psize - min(psize, na.rm = TRUE))/(max(psize, 
                na.rm = TRUE) - min(psize, na.rm = TRUE))
            psize <- (psize * 1) + 0.5
            if (all(is.na(psize))) 
                psize <- rep(1, k)
        }
    }
    rng <- max(ci.ub, na.rm = TRUE) - min(ci.lb, na.rm = TRUE)
    if (annotate) {
        if (showweight) {
            plot.multp.l <- 2
            plot.multp.r <- 2
            axis.multp.l <- 0.2
            axis.multp.r <- 0.2
        }
        else {
            plot.multp.l <- 1.2
            plot.multp.r <- 1.2
            axis.multp.l <- 0.2
            axis.multp.r <- 0.2
        }
    }
    else {
        plot.multp.l <- 1.2
        plot.multp.r <- 0.4
        axis.multp.l <- 0.2
        axis.multp.r <- 0.2
    }
    if (missing(xlim)) {
        xlim <- c(min(ci.lb, na.rm = TRUE) - rng * plot.multp.l, 
            max(ci.ub, na.rm = TRUE) + rng * plot.multp.r)
        xlim <- round(xlim, digits[2])
    }
    alim.spec <- TRUE
    if (missing(alim)) {
        if (is.null(at)) {
            alim <- range(pretty(x = c(min(ci.lb, na.rm = TRUE), 
                max(ci.ub, na.rm = TRUE)), n = steps - 1))
            alim.spec <- FALSE
        }
        else {
            alim <- range(at)
        }
    }
    alim <- sort(alim)
    xlim <- sort(xlim)
    if (xlim[1] > min(yi, na.rm = TRUE)) {
        xlim[1] <- min(yi, na.rm = TRUE)
    }
    if (xlim[2] < max(yi, na.rm = TRUE)) {
        xlim[2] <- max(yi, na.rm = TRUE)
    }
    if (alim[1] > min(yi, na.rm = TRUE)) {
        alim[1] <- min(yi, na.rm = TRUE)
    }
    if (alim[2] < max(yi, na.rm = TRUE)) {
        alim[2] <- max(yi, na.rm = TRUE)
    }
    if (alim[1] < xlim[1]) {
        xlim[1] <- alim[1]
    }
    if (alim[2] > xlim[2]) {
        xlim[2] <- alim[2]
    }
    if (missing(ylim)) {
        ylim <- c(0.5, k + 3)
    }
    else {
        ylim <- sort(ylim)
    }
    if (is.null(at)) {
        if (alim.spec) {
            at <- seq.int(from = alim[1], to = alim[2], length.out = steps)
        }
        else {
            at <- pretty(x = c(min(ci.lb, na.rm = TRUE), max(ci.ub, 
                na.rm = TRUE)), n = steps - 1)
        }
    }
    else {
        at[at < alim[1]] <- alim[1]
        at[at > alim[2]] <- alim[2]
        at <- unique(at)
    }
    at.lab <- at
    if (is.function(atransf)) {
        if (is.null(targs)) {
            at.lab <- formatC(sapply(at.lab, atransf), digits = digits[2], 
                format = "f")
        }
        else {
            at.lab <- formatC(sapply(at.lab, atransf, targs), 
                digits = digits[2], format = "f")
        }
    }
    else {
        at.lab <- formatC(at.lab, digits = digits[2], format = "f")
    }
    par.mar <- par("mar")
    par.mar.adj <- par.mar - c(0, 3, 1, 1)
    par.mar.adj[par.mar.adj < 0] <- 0
    par(mar = par.mar.adj)
    on.exit(par(mar = par.mar))
    plot(NA, NA, xlim = xlim, ylim = ylim, xlab = "", ylab = "", 
        yaxt = "n", xaxt = "n", xaxs = "i", bty = "n", ...)
    abline(h = ylim[2] - 2, ...)
    par.usr <- par("usr")
    height <- par.usr[4] - par.usr[3]
    lheight <- strheight("O")
    cex.adj <- ifelse(k * lheight > height * 0.8, height/(1.25 * 
        k * lheight), 1)
    if (is.null(cex)) {
        cex <- par("cex") * cex.adj
    }
    else {
        if (is.null(cex.lab)) 
            cex.lab <- cex
        if (is.null(cex.axis)) 
            cex.axis <- cex
    }
    if (is.null(cex.lab)) 
        cex.lab <- par("cex") * cex.adj
    if (is.null(cex.axis)) 
        cex.axis <- par("cex") * cex.adj
    axis(side = 1, at = at, labels = at.lab, cex.axis = cex.axis, 
        ...)
    if (missing(xlab)) {
        xlab <- .setxlab(measure, transf.char, atransf.char, 
            gentype = 1)
    }
    else {
        mtext(xlab, side = 1, at = min(at) + (max(at) - min(at))/2, 
            line = par("mgp")[1] - 0.5, cex = cex.lab, ...)
    }
    for (i in seq.int(k)) {
        if (is.na(yi[i]) || is.na(ci.lb)[i] || is.na(ci.ub)[i]) 
            next
        segments(max(ci.lb[i], alim[1]), rows[i], min(ci.ub[i], 
            alim[2]), rows[i], ...)
        if (ci.lb[i] >= alim[1]) {
            segments(ci.lb[i], rows[i] - (height/150) * cex * 
                efac, ci.lb[i], rows[i] + (height/150) * cex * 
                efac, ...)
        }
        else {
            polygon(x = c(alim[1], alim[1] + (1.4/100) * cex * 
                (xlim[2] - xlim[1]), alim[1] + (1.4/100) * cex * 
                (xlim[2] - xlim[1]), alim[1]), y = c(rows[i], 
                rows[i] + (height/150) * cex * efac, rows[i] - 
                  (height/150) * cex * efac, rows[i]), col = "black", 
                ...)
        }
        if (ci.ub[i] <= alim[2]) {
            segments(ci.ub[i], rows[i] - (height/150) * cex * 
                efac, ci.ub[i], rows[i] + (height/150) * cex * 
                efac, ...)
        }
        else {
            polygon(x = c(alim[2], alim[2] - (1.4/100) * cex * 
                (xlim[2] - xlim[1]), alim[2] - (1.4/100) * cex * 
                (xlim[2] - xlim[1]), alim[2]), y = c(rows[i], 
                rows[i] + (height/150) * cex * efac, rows[i] - 
                  (height/150) * cex * efac, rows[i]), col = "black", 
                ...)
        }
    }
    text(xlim[1], rows, slab, pos = 4, cex = cex, ...)
    if (!is.null(ilab)) {
        if (is.null(ilab.xpos)) 
            stop("Must specify 'ilab.xpos' argument when adding information with 'ilab'.")
        if (length(ilab.xpos) != NCOL(ilab)) 
            stop("Number of 'ilab' columns does not match length of 'ilab.xpos' argument.")
        for (l in seq.int(NCOL(ilab))) {
            text(ilab.xpos[l], rows, ilab[, l], pos = ilab.pos[l], 
                cex = cex, ...)
        }
    }
    if (is.numeric(refline)) 
        segments(refline, ylim[1] - 5, refline, ylim[2] - 2, 
            lty = "dotted", ...)
    if (annotate) {
        if (is.function(atransf)) {
            if (is.null(targs)) {
                annotext <- round(cbind(sapply(yi, atransf), 
                  sapply(ci.lb, atransf), sapply(ci.ub, atransf)), 
                  digits[1])
            }
            else {
                annotext <- round(cbind(sapply(yi, atransf, targs), 
                  sapply(ci.lb, atransf, targs), sapply(ci.ub, 
                    atransf, targs)), digits[1])
            }
            rev.order <- ifelse(annotext[, 3] < annotext[, 2], 
                TRUE, FALSE)
            rev.order[is.na(rev.order)] <- FALSE
            annotext[rev.order, 2:3] <- annotext[rev.order, 3:2]
        }
        else {
            annotext <- round(cbind(yi, ci.lb, ci.ub), digits[1])
        }
        if (showweight) {
            annotext <- cbind(round(weights, digits[1]), annotext)
            annotext <- matrix(apply(annotext, 2, format, nsmall = digits[1]), 
                ncol = 4)
            annotext <- cbind(annotext[, 1], "%    ", annotext[, 
                2], " [ ", annotext[, 3], " , ", annotext[, 4], 
                " ]")
        }
        else {
            annotext <- matrix(apply(annotext, 2, format, nsmall = digits[1]), 
                ncol = 3)
            annotext <- cbind(annotext[, 1], " [ ", annotext[, 
                2], " , ", annotext[, 3], " ]")
        }
        annotext <- apply(annotext, 1, paste, collapse = "")
        text(x = xlim[2], rows, labels = annotext, pos = 2, cex = cex, 
            ...)
    }
    points(yi, rows, pch = pch, cex = cex * psize, ...)
    invisible()
}
forest.rma <-
function (x, annotate = TRUE, addfit = TRUE, addcred = FALSE, 
    showweight = FALSE, xlim, alim, ylim, at, steps = 5, level = x$level, 
    digits = 2, refline = 0, xlab, slab, mlab, ilab, ilab.xpos, 
    ilab.pos, order, transf = FALSE, atransf = FALSE, targs, 
    rows, efac = 1, pch = 15, psize, col = "darkgray", border = "darkgray", 
    cex, cex.lab, cex.axis, ...) 
{
    if (!is.element("rma", class(x))) 
        stop("Argument 'x' must be an object of class \"rma\".")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    transf.char <- deparse(substitute(transf))
    atransf.char <- deparse(substitute(atransf))
    if (transf.char != "FALSE" && atransf.char != "FALSE") 
        stop("Use either 'transf' or 'atransf' to specify a transformation (not both).")
    if (missing(targs)) 
        targs <- NULL
    if (missing(at)) 
        at <- NULL
    if (missing(ilab)) 
        ilab <- NULL
    if (missing(ilab.xpos)) 
        ilab.xpos <- NULL
    if (missing(ilab.pos)) 
        ilab.pos <- NULL
    if (missing(order)) 
        order <- NULL
    if (missing(psize)) 
        psize <- NULL
    if (missing(cex)) 
        cex <- NULL
    if (missing(cex.lab)) 
        cex.lab <- NULL
    if (missing(cex.axis)) 
        cex.axis <- NULL
    measure <- x$measure
    if (length(digits) == 1L) 
        digits <- c(digits, digits)
    alpha <- (100 - level)/100
    yi <- x$yi.f
    vi <- x$vi.f
    X <- x$X.f
    k <- length(yi)
    if (missing(slab)) {
        if (x$slab.null) {
            slab <- paste("Study ", x$slab)
        }
        else {
            slab <- x$slab
        }
    }
    if (length(yi) != length(slab)) 
        stop("Number of outcomes does not correspond to the length of the slab argument.")
    if (is.vector(ilab)) 
        ilab <- cbind(ilab)
    if (length(pch) == 1L) 
        pch <- rep(pch, k)
    if (length(pch) != length(yi)) 
        stop("Number of outcomes does not correspond to the length of the pch argument.")
    options(na.action = "na.pass")
    if (x$int.only) {
        pred <- fitted(x)
        pred.ci.lb <- rep(NA, k)
        pred.ci.ub <- rep(NA, k)
    }
    else {
        temp <- predict(x, level = level)
        pred <- temp$pred
        pred.ci.lb <- temp$ci.lb
        pred.ci.ub <- temp$ci.ub
    }
    weights <- weights(x)
    options(na.action = na.act)
    if (!is.null(psize)) {
        if (length(psize) == 1L) 
            psize <- rep(psize, k)
        if (length(psize) != length(yi)) 
            stop("Number of outcomes does not correspond to the length of the psize argument.")
    }
    if (!is.null(order)) {
        if (is.character(order)) {
            if (length(order) != 1) 
                stop("Incorrect length of order argument.")
            if (order == "obs") 
                sort.vec <- order(yi)
            if (order == "fit") 
                sort.vec <- order(pred)
            if (order == "prec") 
                sort.vec <- order(vi, yi)
            if (order == "resid") 
                sort.vec <- order(yi - pred, yi)
            if (order == "rstandard") 
                sort.vec <- order(rstandard(x)$z, yi)
            if (order == "abs.resid") 
                sort.vec <- order(abs(yi - pred), yi)
            if (order == "abs.rstandard") 
                sort.vec <- order(abs(rstandard(x)$z), yi)
        }
        else {
            sort.vec <- order
        }
        yi <- yi[sort.vec]
        vi <- vi[sort.vec]
        X <- X[sort.vec, , drop = FALSE]
        slab <- slab[sort.vec]
        ilab <- ilab[sort.vec, , drop = FALSE]
        pred <- pred[sort.vec]
        pred.ci.lb <- pred.ci.lb[sort.vec]
        pred.ci.ub <- pred.ci.ub[sort.vec]
        weights <- weights[sort.vec]
        pch <- pch[sort.vec]
        psize <- psize[sort.vec]
    }
    k <- length(yi)
    if (missing(rows)) {
        rows <- k:1
    }
    else {
        if (length(rows) == 1L) {
            rows <- rows:(rows - k + 1)
        }
    }
    if (length(rows) != length(yi)) 
        stop("Number of outcomes does not correspond to the length of the rows argument.")
    yi <- yi[k:1]
    vi <- vi[k:1]
    X <- X[k:1, , drop = FALSE]
    slab <- slab[k:1]
    ilab <- ilab[k:1, , drop = FALSE]
    pred <- pred[k:1]
    pred.ci.lb <- pred.ci.lb[k:1]
    pred.ci.ub <- pred.ci.ub[k:1]
    weights <- weights[k:1]
    pch <- pch[k:1]
    psize <- psize[k:1]
    rows <- rows[k:1]
    yiviX.na <- is.na(cbind(yi, vi, X))
    if (any(yiviX.na)) {
        not.na <- apply(yiviX.na, MARGIN = 1, sum) == 0L
        if (na.act == "na.omit") {
            yi <- yi[not.na]
            vi <- vi[not.na]
            X <- X[not.na, , drop = FALSE]
            slab <- slab[not.na]
            ilab <- ilab[not.na, , drop = FALSE]
            pred <- pred[not.na]
            pred.ci.lb <- pred.ci.lb[not.na]
            pred.ci.ub <- pred.ci.ub[not.na]
            weights <- weights[not.na]
            pch <- pch[not.na]
            psize <- psize[not.na]
            rows.new <- rows
            rows.na <- rows[!not.na]
            for (j in seq.int(length(rows.na))) {
                rows.new[rows >= rows.na[j]] <- rows.new[rows >= 
                  rows.na[j]] - 1
            }
            rows <- rows.new[not.na]
        }
        if (na.act == "na.fail") 
            stop("Missing values in results.")
    }
    k <- length(yi)
    ci.lb <- yi - qnorm(alpha/2, lower.tail = FALSE) * sqrt(vi)
    ci.ub <- yi + qnorm(alpha/2, lower.tail = FALSE) * sqrt(vi)
    if (is.function(transf)) {
        if (is.null(targs)) {
            yi <- sapply(yi, transf)
            ci.lb <- sapply(ci.lb, transf)
            ci.ub <- sapply(ci.ub, transf)
            pred <- sapply(pred, transf)
            pred.ci.lb <- sapply(pred.ci.lb, transf)
            pred.ci.ub <- sapply(pred.ci.ub, transf)
        }
        else {
            yi <- sapply(yi, transf, targs)
            ci.lb <- sapply(ci.lb, transf, targs)
            ci.ub <- sapply(ci.ub, transf, targs)
            pred <- sapply(pred, transf, targs)
            pred.ci.lb <- sapply(pred.ci.lb, transf, targs)
            pred.ci.ub <- sapply(pred.ci.ub, transf, targs)
        }
    }
    ci.bounds <- cbind(ci.lb, ci.ub)
    rev.order <- ifelse(ci.ub < ci.lb, TRUE, FALSE)
    rev.order[is.na(rev.order)] <- FALSE
    ci.bounds[rev.order] <- ci.bounds[rev.order, 2:1]
    ci.lb <- ci.bounds[, 1]
    ci.ub <- ci.bounds[, 2]
    pred.ci.bounds <- cbind(pred.ci.lb, pred.ci.ub)
    rev.order <- ifelse(pred.ci.ub < pred.ci.lb, TRUE, FALSE)
    rev.order[is.na(rev.order)] <- FALSE
    pred.ci.bounds[rev.order] <- pred.ci.bounds[rev.order, 2:1]
    pred.ci.lb <- pred.ci.bounds[, 1]
    pred.ci.ub <- pred.ci.bounds[, 2]
    if (is.null(psize)) {
        if (any(vi <= 0, na.rm = TRUE)) {
            psize <- rep(1, k)
        }
        else {
            wi <- 1/sqrt(vi)
            psize <- wi/sum(wi, na.rm = TRUE)
            psize <- (psize - min(psize, na.rm = TRUE))/(max(psize, 
                na.rm = TRUE) - min(psize, na.rm = TRUE))
            psize <- (psize * 1) + 0.5
            if (all(is.na(psize))) 
                psize <- rep(1, k)
        }
    }
    rng <- max(ci.ub, na.rm = TRUE) - min(ci.lb, na.rm = TRUE)
    if (annotate) {
        if (showweight) {
            plot.multp.l <- 2
            plot.multp.r <- 2
            axis.multp.l <- 0.2
            axis.multp.r <- 0.2
        }
        else {
            plot.multp.l <- 1.2
            plot.multp.r <- 1.2
            axis.multp.l <- 0.2
            axis.multp.r <- 0.2
        }
    }
    else {
        plot.multp.l <- 1.2
        plot.multp.r <- 0.4
        axis.multp.l <- 0.2
        axis.multp.r <- 0.2
    }
    if (missing(xlim)) {
        xlim <- c(min(ci.lb, na.rm = TRUE) - rng * plot.multp.l, 
            max(ci.ub, na.rm = TRUE) + rng * plot.multp.r)
        xlim <- round(xlim, digits[2])
    }
    alim.spec <- TRUE
    if (missing(alim)) {
        if (is.null(at)) {
            alim <- range(pretty(x = c(min(ci.lb, na.rm = TRUE), 
                max(ci.ub, na.rm = TRUE)), n = steps - 1))
            alim.spec <- FALSE
        }
        else {
            alim <- range(at)
        }
    }
    alim <- sort(alim)
    xlim <- sort(xlim)
    if (xlim[1] > min(yi, na.rm = TRUE)) {
        xlim[1] <- min(yi, na.rm = TRUE)
    }
    if (xlim[2] < max(yi, na.rm = TRUE)) {
        xlim[2] <- max(yi, na.rm = TRUE)
    }
    if (alim[1] > min(yi, na.rm = TRUE)) {
        alim[1] <- min(yi, na.rm = TRUE)
    }
    if (alim[2] < max(yi, na.rm = TRUE)) {
        alim[2] <- max(yi, na.rm = TRUE)
    }
    if (alim[1] < xlim[1]) {
        xlim[1] <- alim[1]
    }
    if (alim[2] > xlim[2]) {
        xlim[2] <- alim[2]
    }
    if (missing(ylim)) {
        if (x$int.only && addfit) {
            ylim <- c(-1.5, k + 3)
        }
        else {
            ylim <- c(0.5, k + 3)
        }
    }
    else {
        ylim <- sort(ylim)
    }
    if (is.null(at)) {
        if (alim.spec) {
            at <- seq.int(alim[1], alim[2], length.out = steps)
        }
        else {
            at <- pretty(x = c(min(ci.lb, na.rm = TRUE), max(ci.ub, 
                na.rm = TRUE)), n = steps - 1)
        }
    }
    else {
        at[at < alim[1]] <- alim[1]
        at[at > alim[2]] <- alim[2]
        at <- unique(at)
    }
    at.lab <- at
    if (is.function(atransf)) {
        if (is.null(targs)) {
            at.lab <- formatC(sapply(at.lab, atransf), digits = digits[2], 
                format = "f")
        }
        else {
            at.lab <- formatC(sapply(at.lab, atransf, targs), 
                digits = digits[2], format = "f")
        }
    }
    else {
        at.lab <- formatC(at.lab, digits = digits[2], format = "f")
    }
    par.mar <- par("mar")
    par.mar.adj <- par.mar - c(0, 3, 1, 1)
    par.mar.adj[par.mar.adj < 0] <- 0
    par(mar = par.mar.adj)
    on.exit(par(mar = par.mar))
    plot(NA, NA, xlim = xlim, ylim = ylim, xlab = "", ylab = "", 
        yaxt = "n", xaxt = "n", xaxs = "i", bty = "n", ...)
    abline(h = ylim[2] - 2, ...)
    par.usr <- par("usr")
    height <- par.usr[4] - par.usr[3]
    lheight <- strheight("O")
    cex.adj <- ifelse(k * lheight > height * 0.8, height/(1.25 * 
        k * lheight), 1)
    if (is.null(cex)) {
        cex <- par("cex") * cex.adj
    }
    else {
        if (is.null(cex.lab)) 
            cex.lab <- cex
        if (is.null(cex.axis)) 
            cex.axis <- cex
    }
    if (is.null(cex.lab)) 
        cex.lab <- par("cex") * cex.adj
    if (is.null(cex.axis)) 
        cex.axis <- par("cex") * cex.adj
    if (addfit && !x$int.only) {
        for (i in seq.int(k)) {
            if (is.na(pred[i])) 
                next
            if ((pred.ci.lb[i] > alim[1]) && (pred.ci.ub[i] < 
                alim[2])) 
                polygon(x = c(pred.ci.lb[i], pred[i], pred.ci.ub[i], 
                  pred[i]), y = c(rows[i], rows[i] + (height/100) * 
                  cex * efac, rows[i], rows[i] - (height/100) * 
                  cex * efac), col = col, border = border, ...)
        }
    }
    if (addfit && x$int.only) {
        temp <- predict(x, level = level)
        b <- temp$pred
        b.ci.lb <- temp$ci.lb
        b.ci.ub <- temp$ci.ub
        b.cr.lb <- temp$cr.lb
        b.cr.ub <- temp$cr.ub
        if (is.function(transf)) {
            if (is.null(targs)) {
                b <- sapply(b, transf)
                b.ci.lb <- sapply(b.ci.lb, transf)
                b.ci.ub <- sapply(b.ci.ub, transf)
                b.cr.lb <- sapply(b.cr.lb, transf)
                b.cr.ub <- sapply(b.cr.ub, transf)
            }
            else {
                b <- sapply(b, transf, targs)
                b.ci.lb <- sapply(b.ci.lb, transf, targs)
                b.ci.ub <- sapply(b.ci.ub, transf, targs)
                b.cr.lb <- sapply(b.cr.lb, transf, targs)
                b.cr.ub <- sapply(b.cr.ub, transf, targs)
            }
        }
        b.ci.bounds <- cbind(b.ci.lb, b.ci.ub)
        rev.order <- ifelse(b.ci.ub < b.ci.lb, TRUE, FALSE)
        rev.order[is.na(rev.order)] <- FALSE
        b.ci.bounds[rev.order] <- b.ci.bounds[rev.order, 2:1]
        b.ci.lb <- b.ci.bounds[, 1]
        b.ci.ub <- b.ci.bounds[, 2]
        b.cr.bounds <- cbind(b.cr.lb, b.cr.ub)
        rev.order <- ifelse(b.cr.ub < b.cr.lb, TRUE, FALSE)
        rev.order[is.na(rev.order)] <- FALSE
        b.cr.bounds[rev.order] <- b.cr.bounds[rev.order, 2:1]
        b.cr.lb <- b.cr.bounds[, 1]
        b.cr.ub <- b.cr.bounds[, 2]
        if (x$method != "FE" && addcred) {
            segments(max(b.cr.lb, alim[1]), -1, min(b.cr.ub, 
                alim[2]), -1, lty = "dotted", col = col, ...)
            segments(b.cr.lb, -1 - (height/150) * cex * efac, 
                b.cr.lb, -1 + (height/150) * cex * efac, col = col, 
                ...)
            segments(b.cr.ub, -1 - (height/150) * cex * efac, 
                b.cr.ub, -1 + (height/150) * cex * efac, col = col, 
                ...)
        }
        polygon(x = c(b.ci.lb, b, b.ci.ub, b), y = c(-1, -1 + 
            (height/100) * cex * efac, -1, -1 - (height/100) * 
            cex * efac), col = "black", ...)
        if (missing(mlab)) 
            mlab <- ifelse((x$method == "FE"), "FE Model", "RE Model")
        text(xlim[1], -1, mlab, pos = 4, cex = cex, ...)
    }
    axis(side = 1, at = at, labels = at.lab, cex.axis = cex.axis, 
        ...)
    if (missing(xlab)) {
        xlab <- .setxlab(measure, transf.char, atransf.char, 
            gentype = 1)
    }
    else {
        mtext(xlab, side = 1, at = min(at) + (max(at) - min(at))/2, 
            line = par("mgp")[1] - 0.5, cex = cex.lab, ...)
    }
    for (i in seq.int(k)) {
        if (is.na(yi[i]) || is.na(vi)[i]) 
            next
        segments(max(ci.lb[i], alim[1]), rows[i], min(ci.ub[i], 
            alim[2]), rows[i], ...)
        if (ci.lb[i] >= alim[1]) {
            segments(ci.lb[i], rows[i] - (height/150) * cex * 
                efac, ci.lb[i], rows[i] + (height/150) * cex * 
                efac, ...)
        }
        else {
            polygon(x = c(alim[1], alim[1] + (1.4/100) * cex * 
                (xlim[2] - xlim[1]), alim[1] + (1.4/100) * cex * 
                (xlim[2] - xlim[1]), alim[1]), y = c(rows[i], 
                rows[i] + (height/150) * cex * efac, rows[i] - 
                  (height/150) * cex * efac, rows[i]), col = "black", 
                ...)
        }
        if (ci.ub[i] <= alim[2]) {
            segments(ci.ub[i], rows[i] - (height/150) * cex * 
                efac, ci.ub[i], rows[i] + (height/150) * cex * 
                efac, ...)
        }
        else {
            polygon(x = c(alim[2], alim[2] - (1.4/100) * cex * 
                (xlim[2] - xlim[1]), alim[2] - (1.4/100) * cex * 
                (xlim[2] - xlim[1]), alim[2]), y = c(rows[i], 
                rows[i] + (height/150) * cex * efac, rows[i] - 
                  (height/150) * cex * efac, rows[i]), col = "black", 
                ...)
        }
    }
    text(xlim[1], rows, slab, pos = 4, cex = cex, ...)
    if (!is.null(ilab)) {
        if (is.null(ilab.xpos)) 
            stop("Must specify 'ilab.xpos' argument when adding information with 'ilab'.")
        if (length(ilab.xpos) != NCOL(ilab)) 
            stop("Number of 'ilab' columns does not match length of 'ilab.xpos' argument.")
        for (l in seq.int(NCOL(ilab))) {
            text(ilab.xpos[l], rows, ilab[, l], pos = ilab.pos[l], 
                cex = cex, ...)
        }
    }
    if (is.numeric(refline)) 
        segments(refline, ylim[1] - 5, refline, ylim[2] - 2, 
            lty = "dotted", ...)
    if (annotate) {
        if (is.function(atransf)) {
            if (is.null(targs)) {
                if (addfit && x$int.only) {
                  annotext <- round(cbind(sapply(c(yi, b), atransf), 
                    sapply(c(ci.lb, b.ci.lb), atransf), sapply(c(ci.ub, 
                      b.ci.ub), atransf)), digits[1])
                }
                else {
                  annotext <- round(cbind(sapply(yi, atransf), 
                    sapply(ci.lb, atransf), sapply(ci.ub, atransf)), 
                    digits[1])
                }
            }
            else {
                if (addfit && x$int.only) {
                  annotext <- round(cbind(sapply(c(yi, b), atransf, 
                    targs), sapply(c(ci.lb, b.ci.lb), atransf, 
                    targs), sapply(c(ci.ub, b.ci.ub), atransf, 
                    targs)), digits[1])
                }
                else {
                  annotext <- round(cbind(sapply(yi, atransf, 
                    targs), sapply(ci.lb, atransf, targs), sapply(ci.ub, 
                    atransf, targs)), digits[1])
                }
            }
            rev.order <- ifelse(annotext[, 3] < annotext[, 2], 
                TRUE, FALSE)
            rev.order[is.na(rev.order)] <- FALSE
            annotext[rev.order, 2:3] <- annotext[rev.order, 3:2]
        }
        else {
            if (addfit && x$int.only) {
                annotext <- round(cbind(c(yi, b), c(ci.lb, b.ci.lb), 
                  c(ci.ub, b.ci.ub)), digits[1])
            }
            else {
                annotext <- round(cbind(yi, ci.lb, ci.ub), digits[1])
            }
        }
        if (showweight) {
            if (addfit && x$int.only) {
                annotext <- cbind(round(c(weights, 100), digits[1]), 
                  annotext)
            }
            else {
                annotext <- cbind(round(weights, digits[1]), 
                  annotext)
            }
            annotext <- matrix(apply(annotext, 2, format, nsmall = digits[1]), 
                ncol = 4)
            annotext <- cbind(annotext[, 1], "%    ", annotext[, 
                2], " [ ", annotext[, 3], " , ", annotext[, 4], 
                " ]")
        }
        else {
            annotext <- matrix(apply(annotext, 2, format, nsmall = digits[1]), 
                ncol = 3)
            annotext <- cbind(annotext[, 1], " [ ", annotext[, 
                2], " , ", annotext[, 3], " ]")
        }
        annotext <- apply(annotext, 1, paste, collapse = "")
        if (addfit && x$int.only) {
            text(x = xlim[2], c(rows, -1), labels = annotext, 
                pos = 2, cex = cex, ...)
        }
        else {
            text(x = xlim[2], rows, labels = annotext, pos = 2, 
                cex = cex, ...)
        }
    }
    points(yi, rows, pch = pch, cex = cex * psize, ...)
    if (x$int.only && addfit) 
        abline(h = 0, ...)
    invisible()
}
fsn <-
function (yi, vi, sei, data, type = "Rosenthal", alpha = 0.05, 
    target, subset, digits = 4) 
{
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    type <- match.arg(type, c("Rosenthal", "Orwin", "Rosenberg"))
    if (missing(target)) 
        target <- NULL
    if (missing(data)) 
        data <- NULL
    if (is.null(data)) {
        data <- sys.frame(sys.parent())
    }
    else {
        if (!is.data.frame(data)) {
            data <- data.frame(data)
        }
    }
    mf <- match.call()
    mf.yi <- mf[[match("yi", names(mf))]]
    mf.vi <- mf[[match("vi", names(mf))]]
    mf.sei <- mf[[match("sei", names(mf))]]
    mf.subset <- mf[[match("subset", names(mf))]]
    yi <- eval(mf.yi, data, enclos = sys.frame(sys.parent()))
    vi <- eval(mf.vi, data, enclos = sys.frame(sys.parent()))
    sei <- eval(mf.sei, data, enclos = sys.frame(sys.parent()))
    subset <- eval(mf.subset, data, enclos = sys.frame(sys.parent()))
    if (is.null(vi)) 
        vi <- sei^2
    if (length(yi) != length(vi)) 
        stop("Length of yi and vi (or sei) vectors is not the same.")
    if (!is.null(subset)) {
        yi <- yi[subset]
        vi <- vi[subset]
    }
    yivi.na <- is.na(cbind(yi, vi))
    if (any(yivi.na)) {
        not.na <- apply(yivi.na, MARGIN = 1, sum) == 0L
        if (na.act == "na.omit" || na.act == "na.exclude" || 
            na.act == "na.pass") {
            yi <- yi[not.na]
            vi <- vi[not.na]
        }
        if (na.act == "na.fail") 
            stop("Missing values in results.")
    }
    if (type == "Rosenthal") {
        k <- length(yi)
        zi <- yi/sqrt(vi)
        z.avg <- abs(sum(zi)/sqrt(k))
        pval <- pnorm(z.avg, lower.tail = FALSE)
        fsnum <- ceiling(max(0, k * (z.avg/qnorm(alpha, lower.tail = FALSE))^2 - 
            k))
        meanes <- NA
        target <- NA
    }
    if (type == "Orwin") {
        k <- length(yi)
        meanes <- mean(yi)
        if (is.null(target)) {
            target <- meanes/2
        }
        if (identical(target, 0)) {
            fsnum <- Inf
        }
        else {
            fsnum <- ceiling(max(0, k * (meanes - target)/target))
        }
        pval <- NA
    }
    if (type == "Rosenberg") {
        k <- length(yi)
        wi <- 1/vi
        meanes <- sum(wi * yi)/sum(wi)
        zval <- meanes/sqrt(1/sum(wi))
        w.p <- (sum(wi * yi)/qnorm(alpha/2, lower.tail = FALSE))^2 - 
            sum(wi)
        pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)
        fsnum <- ceiling(max(0, k * w.p/sum(wi)))
        target <- NA
    }
    res <- list(type, fsnum, alpha, pval, meanes, target, digits)
    names(res) <- c("type", "fsnum", "alpha", "pval", "meanes", 
        "target", "digits")
    class(res) <- c("fsn")
    return(res)
}
funnel <-
function (x, ...) 
UseMethod("funnel")
funnel.rma <-
function (x, yaxis = "sei", xlim, ylim, xlab, ylab, steps = 5, 
    at, atransf = FALSE, targs, digits, level = x$level, addtau2 = FALSE, 
    type = "rstandard", back = "lightgray", shade = "white", 
    hlines = "white", refline, pch = 19, pch.fill = 21, ci.res = 1000, 
    ...) 
{
    if (!is.element("rma", class(x))) 
        stop("Argument 'x' must be an object of class \"rma\".")
    yaxis <- match.arg(yaxis, c("sei", "vi", "seinv", "vinv", 
        "ni", "ninv", "sqrtni", "sqrtninv", "lni"))
    type <- match.arg(type, c("rstandard", "rstudent"))
    atransf.char <- deparse(substitute(atransf))
    if (missing(ylab)) {
        if (yaxis == "sei") 
            ylab <- "Standard Error"
        if (yaxis == "vi") 
            ylab <- "Variance"
        if (yaxis == "seinv") 
            ylab <- "Inverse Standard Error"
        if (yaxis == "vinv") 
            ylab <- "Inverse Variance"
        if (yaxis == "ni") 
            ylab <- "Sample Size"
        if (yaxis == "ninv") 
            ylab <- "Inverse Sample Size"
        if (yaxis == "sqrtni") 
            ylab <- "Square-Root Sample Size"
        if (yaxis == "sqrtninv") 
            ylab <- "Inverse Square-Root Sample Size"
        if (yaxis == "lni") 
            ylab <- "Log Sample Size"
    }
    if (missing(at)) 
        at <- NULL
    if (missing(targs)) 
        targs <- NULL
    if (yaxis == "ni" || yaxis == "ninv" || yaxis == "sqrtni" || 
        yaxis == "sqrtninv" || yaxis == "lni") {
        if (is.null(x$ni)) 
            stop("No sample size information stored in model object.")
        if (any(is.na(x$ni))) 
            warning("Sample size information stored in model object \n  contains NAs. Not all studies will be plotted.")
    }
    if (missing(digits)) {
        if (yaxis == "sei") 
            digits <- c(2, 3)
        if (yaxis == "vi") 
            digits <- c(2, 3)
        if (yaxis == "seinv") 
            digits <- c(2, 3)
        if (yaxis == "vinv") 
            digits <- c(2, 3)
        if (yaxis == "ni") 
            digits <- c(2, 0)
        if (yaxis == "ninv") 
            digits <- c(2, 3)
        if (yaxis == "sqrtni") 
            digits <- c(2, 3)
        if (yaxis == "sqrtninv") 
            digits <- c(2, 3)
        if (yaxis == "lni") 
            digits <- c(2, 3)
    }
    else {
        if (length(digits) == 1L) 
            digits <- c(digits, digits)
    }
    if (x$int.only) {
        if (missing(refline)) 
            refline <- x$b
        tau2 <- ifelse(addtau2, x$tau2, 0)
        yi <- x$yi
        vi <- x$vi
        ni <- x$ni
        sei <- sqrt(vi)
        if (missing(xlab)) 
            xlab <- .setxlab(x$measure, transf.char = "FALSE", 
                atransf.char, gentype = 1)
    }
    else {
        if (missing(refline)) 
            refline <- 0
        tau2 <- 0
        na.act <- getOption("na.action")
        options(na.action = "na.pass")
        if (type == "rstandard") {
            res <- rstandard(x)
        }
        else {
            res <- rstudent(x)
        }
        options(na.action = na.act)
        not.na <- !is.na(res$resid)
        yi <- res$resid[not.na]
        sei <- res$se[not.na]
        ni <- x$ni.f[not.na]
        vi <- sei^2
        if (missing(xlab)) 
            xlab <- "Residual Value"
    }
    if (missing(ylim)) {
        if (yaxis == "sei") 
            ylim <- c(max(sei), 0)
        if (yaxis == "vi") 
            ylim <- c(max(vi), 0)
        if (yaxis == "seinv") 
            ylim <- c(min(1/sei), max(1/sei))
        if (yaxis == "vinv") 
            ylim <- c(min(1/vi), max(1/vi))
        if (yaxis == "ni") 
            ylim <- c(min(ni, na.rm = TRUE), max(ni, na.rm = TRUE))
        if (yaxis == "ninv") 
            ylim <- c(max(1/ni, na.rm = TRUE), min(1/ni, na.rm = TRUE))
        if (yaxis == "sqrtni") 
            ylim <- c(min(sqrt(ni), na.rm = TRUE), max(sqrt(ni), 
                na.rm = TRUE))
        if (yaxis == "sqrtninv") 
            ylim <- c(max(1/sqrt(ni), na.rm = TRUE), min(1/sqrt(ni), 
                na.rm = TRUE))
        if (yaxis == "lni") 
            ylim <- c(min(log(ni), na.rm = TRUE), max(log(ni), 
                na.rm = TRUE))
    }
    else {
        if (yaxis == "sei" || yaxis == "vi" || yaxis == "ninv" || 
            yaxis == "sqrtninv") 
            ylim <- c(max(ylim), min(ylim))
        if (yaxis == "seinv" || yaxis == "vinv" || yaxis == "ni" || 
            yaxis == "sqrtni" || yaxis == "lni") 
            ylim <- c(min(ylim), max(ylim))
        if (yaxis == "sei" || yaxis == "vi" || yaxis == "ni" || 
            yaxis == "ninv" || yaxis == "sqrtni" || yaxis == 
            "sqrtninv" || yaxis == "lni") {
            if (ylim[1] < 0 || ylim[2] < 0) 
                stop("Both limits for the y axis must be >= 0.")
        }
        if (yaxis == "seinv" || yaxis == "vinv") {
            if (ylim[1] <= 0 || ylim[2] <= 0) 
                stop("Both limits for the y axis must be > 0.")
        }
    }
    if (yaxis == "sei" || yaxis == "vi" || yaxis == "seinv" || 
        yaxis == "vinv") {
        alpha <- (100 - level)/100
        alpha.min <- min(alpha)
        avals <- length(alpha)
        if (yaxis == "sei") {
            x.lb.bot <- refline - qnorm(alpha.min/2, lower.tail = FALSE) * 
                sqrt(ylim[1]^2 + tau2)
            x.ub.bot <- refline + qnorm(alpha.min/2, lower.tail = FALSE) * 
                sqrt(ylim[1]^2 + tau2)
        }
        if (yaxis == "vi") {
            x.lb.bot <- refline - qnorm(alpha.min/2, lower.tail = FALSE) * 
                sqrt(ylim[1] + tau2)
            x.ub.bot <- refline + qnorm(alpha.min/2, lower.tail = FALSE) * 
                sqrt(ylim[1] + tau2)
        }
        if (yaxis == "seinv") {
            x.lb.bot <- refline - qnorm(alpha.min/2, lower.tail = FALSE) * 
                sqrt(1/ylim[1]^2 + tau2)
            x.ub.bot <- refline + qnorm(alpha.min/2, lower.tail = FALSE) * 
                sqrt(1/ylim[1]^2 + tau2)
        }
        if (yaxis == "vinv") {
            x.lb.bot <- refline - qnorm(alpha.min/2, lower.tail = FALSE) * 
                sqrt(1/ylim[1] + tau2)
            x.ub.bot <- refline + qnorm(alpha.min/2, lower.tail = FALSE) * 
                sqrt(1/ylim[1] + tau2)
        }
        if (missing(xlim)) {
            xlim <- c(min(x.lb.bot, min(yi)), max(x.ub.bot, max(yi)))
            rxlim <- xlim[2] - xlim[1]
            xlim[1] <- xlim[1] - (rxlim * 0.1)
            xlim[2] <- xlim[2] + (rxlim * 0.1)
        }
        else {
            xlim <- sort(xlim)
        }
    }
    if (yaxis == "ni" || yaxis == "ninv" || yaxis == "sqrtni" || 
        yaxis == "sqrtninv" || yaxis == "lni") {
        if (missing(xlim)) {
            xlim <- c(min(yi), max(yi))
            rxlim <- xlim[2] - xlim[1]
            xlim[1] <- xlim[1] - (rxlim * 0.1)
            xlim[2] <- xlim[2] + (rxlim * 0.1)
        }
        else {
            xlim <- sort(xlim)
        }
    }
    if (!is.null(at)) {
        xlim[1] <- min(c(xlim[1], at), na.rm = TRUE)
        xlim[2] <- max(c(xlim[2], at), na.rm = TRUE)
    }
    plot(NA, NA, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, 
        xaxt = "n", yaxt = "n", bty = "n", ...)
    par.usr <- par("usr")
    rect(par.usr[1], par.usr[3], par.usr[2], par.usr[4], col = back, 
        border = NA, ...)
    axis(side = 2, at = seq.int(ylim[1], ylim[2], length.out = steps), 
        labels = formatC(seq.int(ylim[1], ylim[2], length.out = steps), 
            digits = digits[2], format = "f"), ...)
    abline(h = seq.int(ylim[1], ylim[2], length.out = steps), 
        col = hlines, ...)
    if (yaxis == "sei" || yaxis == "vi" || yaxis == "seinv" || 
        yaxis == "vinv") {
        if (yaxis == "sei") {
            rylim <- ylim[1] - ylim[2]
            ylim[1] <- ylim[1] + (rylim * 0.1)
            ylim[2] <- max(0, ylim[2] - (rylim * 0.1))
        }
        if (yaxis == "vi") {
            rylim <- ylim[1] - ylim[2]
            ylim[1] <- ylim[1] + (rylim * 0.1)
            ylim[2] <- max(0, ylim[2] - (rylim * 0.1))
        }
        if (yaxis == "seinv") {
            rylim <- ylim[2] - ylim[1]
            ylim[2] <- ylim[2] + (rylim * 0.1)
        }
        if (yaxis == "vinv") {
            rylim <- ylim[2] - ylim[1]
            ylim[2] <- ylim[2] + (rylim * 0.1)
        }
        yi.vals <- seq.int(ylim[1], ylim[2], length.out = ci.res)
        if (yaxis == "sei") 
            vi.vals <- yi.vals^2
        if (yaxis == "vi") 
            vi.vals <- yi.vals
        if (yaxis == "seinv") 
            vi.vals <- 1/yi.vals^2
        if (yaxis == "vinv") 
            vi.vals <- 1/yi.vals
        for (m in avals:1) {
            ci.left <- refline - qnorm(alpha[m]/2, lower.tail = FALSE) * 
                sqrt(vi.vals + tau2)
            ci.right <- refline + qnorm(alpha[m]/2, lower.tail = FALSE) * 
                sqrt(vi.vals + tau2)
            polygon(c(ci.left, ci.right[ci.res:1]), c(yi.vals, 
                yi.vals[ci.res:1]), border = NA, col = shade[m], 
                ...)
            lines(ci.left, yi.vals, lty = "dotted", ...)
            lines(ci.right, yi.vals, lty = "dotted", ...)
        }
    }
    if (yaxis == "sei" || yaxis == "vi" || yaxis == "seinv" || 
        yaxis == "vinv") 
        segments(refline, ylim[1], refline, ylim[2], ...)
    if (yaxis == "ni" || yaxis == "ninv" || yaxis == "sqrtni" || 
        yaxis == "sqrtninv" || yaxis == "lni") 
        abline(v = refline, ...)
    if (yaxis == "sei") 
        points(yi, sei, pch = pch, ...)
    if (yaxis == "vi") 
        points(yi, vi, pch = pch, ...)
    if (yaxis == "seinv") 
        points(yi, 1/sei, pch = pch, ...)
    if (yaxis == "vinv") 
        points(yi, 1/vi, pch = pch, ...)
    if (yaxis == "ni") 
        points(yi, ni, pch = pch, ...)
    if (yaxis == "ninv") 
        points(yi, 1/ni, pch = pch, ...)
    if (yaxis == "sqrtni") 
        points(yi, sqrt(ni), pch = pch, ...)
    if (yaxis == "sqrtninv") 
        points(yi, 1/sqrt(ni), pch = pch, ...)
    if (yaxis == "lni") 
        points(yi, log(ni), pch = pch, ...)
    if (is.element("rma.uni.trimfill", class(x))) {
        if (yaxis == "sei") 
            points(yi[x$fill == 1], (sei)[x$fill == 1], pch = pch.fill, 
                col = "black", bg = "white", ...)
        if (yaxis == "vi") 
            points(yi[x$fill == 1], (sei^2)[x$fill == 1], pch = pch.fill, 
                col = "black", bg = "white", ...)
        if (yaxis == "seinv") 
            points(yi[x$fill == 1], (1/sei)[x$fill == 1], pch = pch.fill, 
                col = "black", bg = "white", ...)
        if (yaxis == "vinv") 
            points(yi[x$fill == 1], (1/vi)[x$fill == 1], pch = pch.fill, 
                col = "black", bg = "white", ...)
        if (yaxis == "ni") 
            points(yi[x$fill == 1], (ni)[x$fill == 1], pch = pch.fill, 
                col = "black", bg = "white", ...)
        if (yaxis == "ninv") 
            points(yi[x$fill == 1], (1/ni)[x$fill == 1], pch = pch.fill, 
                col = "black", bg = "white", ...)
        if (yaxis == "sqrtni") 
            points(yi[x$fill == 1], sqrt(ni)[x$fill == 1], pch = pch.fill, 
                col = "black", bg = "white", ...)
        if (yaxis == "sqrtninv") 
            points(yi[x$fill == 1], sqrt(1/ni)[x$fill == 1], 
                pch = pch.fill, col = "black", bg = "white", 
                ...)
        if (yaxis == "lni") 
            points(yi[x$fill == 1], (log(ni))[x$fill == 1], pch = pch.fill, 
                col = "black", bg = "white", ...)
    }
    box(bty = "l")
    if (is.null(at)) {
        at <- axTicks(side = 1)
    }
    else {
        at <- at[at > par("usr")[1]]
        at <- at[at < par("usr")[2]]
    }
    at.lab <- at
    if (is.function(atransf)) {
        if (is.null(targs)) {
            at.lab <- formatC(sapply(at.lab, atransf), digits = digits[1], 
                format = "f")
        }
        else {
            at.lab <- formatC(sapply(at.lab, atransf, targs), 
                digits = digits[1], format = "f")
        }
    }
    else {
        at.lab <- formatC(at.lab, digits = digits[1], format = "f")
    }
    axis(side = 1, at = at, labels = at.lab, ...)
    invisible()
}
galbraith <-
function (x, ...) 
UseMethod("radial")
hatvalues.rma.uni <-
function (model, ...) 
{
    if (!is.element("rma.uni", class(model))) 
        stop("Argument 'model' must be an object of class \"rma.uni\".")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    x <- model
    if (x$weighted) {
        wi <- 1/(x$vi + x$tau2)
        W <- .diag(wi)
        stXWX <- .invcalc(X = x$X, W = W, k = x$k)
        H <- x$X %*% stXWX %*% crossprod(x$X, W)
    }
    else {
        stXX <- .invcalc(X = x$X, W = diag(x$k), k = x$k)
        H <- x$X %*% tcrossprod(stXX, x$X)
    }
    hii <- rep(NA, x$k.f)
    hii[x$not.na] <- diag(H)
    hii[hii > 1 - 10 * .Machine$double.eps] <- 1
    names(hii) <- x$slab
    if (na.act == "na.omit") 
        hii <- hii[x$not.na]
    if (na.act == "na.fail" && any(!x$not.na)) 
        stop("Missing values in results.")
    return(hii)
}
influence.rma.uni <-
function (model, digits = model$digits, ...) 
{
    if (!is.element("rma.uni", class(model))) 
        stop("Argument 'model' must be an object of class \"rma.uni\".")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    x <- model
    if (x$k == 1) 
        stop("Stopped because k = 1.")
    tau2.del <- rep(NA, x$k.f)
    delpred <- rep(NA, x$k.f)
    vdelpred <- rep(NA, x$k.f)
    QE.del <- rep(NA, x$k.f)
    dffits <- rep(NA, x$k.f)
    dfbetas <- matrix(NA, nrow = x$k.f, ncol = length(x$b))
    cook.d <- rep(NA, x$k.f)
    cov.r <- rep(NA, x$k.f)
    weight <- rep(NA, x$k.f)
    det.full <- det(x$vb)
    pred.full <- x$X.f %*% x$b
    if (x$weighted) {
        wi <- 1/(x$vi + x$tau2)
        W <- .diag(wi)
        svb <- crossprod(x$X, W) %*% x$X/x$s2w
    }
    else {
        svb <- chol2inv(chol(x$vb))
        stXX <- .invcalc(X = x$X, W = diag(x$k), k = x$k)
        H <- x$X %*% stXX %*% t(x$X)
    }
    options(na.action = "na.pass")
    hii <- hatvalues(x)
    options(na.action = na.act)
    o.warn <- getOption("warn")
    on.exit(options(warn = o.warn))
    options(warn = -1)
    for (i in seq.int(x$k.f)[x$not.na]) {
        res <- try(rma(x$yi.f[-i], x$vi.f[-i], mods = cbind(x$X.f[-i, 
            ]), method = x$method, weighted = x$weighted, intercept = FALSE, 
            knha = x$knha, control = x$control, ...), silent = TRUE)
        if (inherits(res, "try-error")) 
            next
        tau2.del[i] <- res$tau2
        Xi <- matrix(x$X.f[i, ], nrow = 1)
        delpred[i] <- Xi %*% res$b
        vdelpred[i] <- Xi %*% tcrossprod(res$vb, Xi)
        QE.del[i] <- res$QE
        if (x$weighted) {
            dffits[i] <- (pred.full[i] - delpred[i])/sqrt(res$s2w * 
                hii[i] * (tau2.del[i] + x$vi.f[i]))
        }
        else {
            dffits[i] <- (pred.full[i] - delpred[i])/(sqrt(res$s2w * 
                diag(H %*% .diag(tau2.del[i] + x$vi) %*% t(H))))[i - 
                x$k.f + sum(x$not.na)]
        }
        dfbeta <- x$b - res$b
        if (x$weighted) {
            vb.del <- .invcalc(X = x$X, W = .diag(1/(x$vi + tau2.del[i])), 
                k = x$k)
            dfbetas[i, ] <- dfbeta/sqrt(res$s2w * diag(vb.del))
        }
        else {
            vb.del <- tcrossprod(stXX, x$X) %*% .diag(x$vi + 
                tau2.del[i]) %*% x$X %*% stXX
            dfbetas[i, ] <- dfbeta/sqrt(res$s2w * diag(vb.del))
        }
        cook.d[i] <- (crossprod(dfbeta, svb) %*% dfbeta)
        cov.r[i] <- det(res$vb)/det.full
    }
    delresid <- x$yi.f - delpred
    sedelresid <- sqrt(x$vi.f + vdelpred + tau2.del)
    standelres <- delresid/sedelresid
    if (x$weighted) {
        weight[x$not.na] <- wi/sum(wi) * 100
    }
    else {
        weight[x$not.na] <- 1/x$k * 100
    }
    inf <- cbind(standelres, dffits, cook.d, cov.r, tau2.del, 
        QE.del, hii, weight)
    dfb <- cbind(dfbetas)
    inf <- data.frame(inf)
    dfb <- data.frame(dfb)
    is.infl <- abs(inf$dffits) > 3 * sqrt(x$p/(x$k - x$p)) | 
        pchisq(inf$cook.d, df = x$p) > 0.5 | inf$hii > 3 * x$p/x$k | 
        apply(abs(dfb) > 1, 1, any)
    out <- list(inf = inf, dfb = dfb, tau2 = x$tau2, QE = x$QE, 
        ids = x$ids, not.na = x$not.na, is.infl = is.infl, k = x$k, 
        p = x$p, digits = digits)
    dimnames(out$inf)[[1]] <- x$slab
    dimnames(out$dfb)[[1]] <- x$slab
    dimnames(out$dfb)[[2]] <- dimnames(x$b)[[1]]
    dimnames(out$inf)[[2]] <- c("rstudent", "dffits", "cook.d", 
        "cov.r", "tau2.del", "QE.del", "hat", "weight")
    class(out) <- "infl.rma.uni"
    return(out)
}
labbe <-
function (x, ...) 
UseMethod("labbe")
labbe.rma <-
function (x, xlim, ylim, xlab, ylab, add = x$add, to = x$to, 
    transf = FALSE, targs, pch = 21, psize, bg = "gray", ...) 
{
    if (!is.element("rma", class(x))) 
        stop("Argument 'x' must be an object of class \"rma\".")
    if (!x$int.only) 
        stop("L'Abbe plot only applicable for models without moderators.")
    if (!is.element(x$measure, c("RR", "OR", "RD", "AS", "IRR", 
        "IRD", "IRSD"))) 
        stop("Argument 'measure' must be one of the following: 'RR','OR','RD','AS','IRR','IRD','IRSD'.")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    if (length(add) == 2) 
        add <- add[1]
    if (length(to) == 2) 
        to <- to[1]
    if (!is.element(to, c("all", "only0", "if0all", "none"))) 
        stop("Unknown 'to' argument specified.")
    transf.char <- deparse(substitute(transf))
    if (missing(targs)) 
        targs <- NULL
    if (missing(psize)) 
        psize <- NULL
    k <- x$k.f
    if (length(pch) == 1L) 
        pch <- rep(pch, k)
    if (length(pch) != k) 
        stop("Number of tables does not correspond to the length of the pch argument.")
    if (!is.null(psize)) {
        if (length(psize) == 1L) 
            psize <- rep(psize, k)
        if (length(psize) != k) 
            stop("Number of tables does not correspond to the length of the psize argument.")
    }
    if (length(bg) == 1L) 
        bg <- rep(bg, k)
    if (length(bg) != k) 
        stop("Number of tables does not correspond to the length of the bg argument.")
    ai <- x$ai.f
    bi <- x$bi.f
    ci <- x$ci.f
    di <- x$di.f
    x1i <- x$x1i.f
    x2i <- x$x2i.f
    t1i <- x$t1i.f
    t2i <- x$t2i.f
    yi.is.na <- is.na(x$yi.f)
    ai[yi.is.na] <- NA
    bi[yi.is.na] <- NA
    ci[yi.is.na] <- NA
    di[yi.is.na] <- NA
    x1i[yi.is.na] <- NA
    x2i[yi.is.na] <- NA
    t1i[yi.is.na] <- NA
    t2i[yi.is.na] <- NA
    options(na.action = "na.pass")
    if (x$measure == "RR") {
        dat.t <- escalc(measure = "PLN", xi = ai, mi = bi, add = add, 
            to = to)
        dat.c <- escalc(measure = "PLN", xi = ci, mi = di, add = add, 
            to = to)
    }
    if (x$measure == "OR") {
        dat.t <- escalc(measure = "PLO", xi = ai, mi = bi, add = add, 
            to = to)
        dat.c <- escalc(measure = "PLO", xi = ci, mi = di, add = add, 
            to = to)
    }
    if (x$measure == "RD") {
        dat.t <- escalc(measure = "PR", xi = ai, mi = bi, add = add, 
            to = to)
        dat.c <- escalc(measure = "PR", xi = ci, mi = di, add = add, 
            to = to)
    }
    if (x$measure == "AS") {
        dat.t <- escalc(measure = "PAS", xi = ai, mi = bi, add = add, 
            to = to)
        dat.c <- escalc(measure = "PAS", xi = ci, mi = di, add = add, 
            to = to)
    }
    if (x$measure == "IRR") {
        dat.t <- escalc(measure = "IRLN", xi = x1i, ti = t1i, 
            add = add, to = to)
        dat.c <- escalc(measure = "IRLN", xi = x2i, ti = t2i, 
            add = add, to = to)
    }
    if (x$measure == "IRD") {
        dat.t <- escalc(measure = "IR", xi = x1i, ti = t1i, add = add, 
            to = to)
        dat.c <- escalc(measure = "IR", xi = x2i, ti = t2i, add = add, 
            to = to)
    }
    if (x$measure == "IRSD") {
        dat.t <- escalc(measure = "IRS", xi = x1i, ti = t1i, 
            add = add, to = to)
        dat.c <- escalc(measure = "IRS", xi = x2i, ti = t2i, 
            add = add, to = to)
    }
    options(na.action = na.act)
    dat.t.dat.c.na <- is.na(cbind(dat.t, dat.c))
    if (any(dat.t.dat.c.na)) {
        not.na <- apply(dat.t.dat.c.na, MARGIN = 1, sum) == 0L
        dat.t <- dat.t[not.na, ]
        dat.c <- dat.c[not.na, ]
    }
    if (length(dat.t$yi) == 0 || length(dat.c$yi) == 0) 
        stop("No information in object to compute arm-level outcomes.")
    if (is.null(psize)) {
        vi <- dat.t$vi + dat.c$vi
        wi <- 1/sqrt(vi)
        psize <- 0.5 + 3 * (wi - min(wi))/(max(wi) - min(wi))
    }
    min.yi <- min(c(dat.t$yi, dat.c$yi))
    max.yi <- max(c(dat.t$yi, dat.c$yi))
    rng.yi <- max.yi - min.yi
    if (x$measure == "RD") 
        c.vals <- seq(ifelse(x$b > 0, 0, -x$b), ifelse(x$b > 
            0, 1 - x$b, 1), length.out = 1000)
    if (x$measure == "RR") 
        c.vals <- seq(min.yi - rng.yi, ifelse(x$b > 0, 0 - x$b, 
            0), length.out = 1000)
    if (x$measure == "OR") 
        c.vals <- seq(min.yi - rng.yi, max.yi + rng.yi, length.out = 1000)
    if (x$measure == "AS") 
        c.vals <- seq(ifelse(x$b > 0, 0, -x$b), ifelse(x$b > 
            0, asin(sqrt(1)) - x$b, asin(sqrt(1))), length.out = 1000)
    if (x$measure == "IRR") 
        c.vals <- seq(min.yi - rng.yi, ifelse(x$b > 0, 0 - x$b, 
            0), length.out = 1000)
    if (x$measure == "IRD") 
        c.vals <- seq(ifelse(x$b > 0, 0, -x$b), ifelse(x$b > 
            0, 1 - x$b, 1), length.out = 1000)
    if (x$measure == "IRSD") 
        c.vals <- seq(ifelse(x$b > 0, 0, -x$b), ifelse(x$b > 
            0, 1 - x$b, 1), length.out = 1000)
    t.vals <- x$b + 1 * c.vals
    if (is.function(transf)) {
        if (is.null(targs)) {
            dat.t$yi <- sapply(dat.t$yi, transf)
            dat.c$yi <- sapply(dat.c$yi, transf)
            c.vals <- sapply(c.vals, transf)
            t.vals <- sapply(t.vals, transf)
        }
        else {
            dat.t$yi <- sapply(dat.t$yi, transf, targs)
            dat.c$yi <- sapply(dat.c$yi, transf, targs)
            c.vals <- sapply(c.vals, transf, targs)
            t.vals <- sapply(t.vals, transf, targs)
        }
    }
    min.yi <- min(c(dat.t$yi, dat.c$yi))
    max.yi <- max(c(dat.t$yi, dat.c$yi))
    if (missing(xlim)) 
        xlim <- c(min.yi, max.yi)
    if (missing(ylim)) 
        ylim <- c(min.yi, max.yi)
    if (T) {
        order.vec <- order(psize, decreasing = TRUE)
        dat.t$yi <- dat.t$yi[order.vec]
        dat.c$yi <- dat.c$yi[order.vec]
        psize <- psize[order.vec]
        pch <- pch[order.vec]
    }
    if (missing(xlab)) {
        xlab <- "Observed Outcome"
        if (x$measure == "OR") {
            if (transf.char == "FALSE") {
                xlab <- "Log Odds"
            }
            else {
                xlab <- "Transformed Log Odds"
                if (transf.char == "exp" || transf.char == "transf.exp.int") 
                  xlab <- "Odds"
                if (transf.char == "transf.ilogit" || transf.char == 
                  "transf.ilogit.int") 
                  xlab <- "Proportion"
            }
        }
        if (x$measure == "RR") {
            if (transf.char == "FALSE") {
                xlab <- "Log Proportion"
            }
            else {
                xlab <- "Transformed Log Proportion"
                if (transf.char == "exp" || transf.char == "transf.exp.int") 
                  xlab <- "Proportion"
            }
        }
        if (x$measure == "RD") {
            if (transf.char == "FALSE") {
                xlab <- "Proportion"
            }
            else {
                xlab <- "Transformed Proportion"
            }
        }
        if (x$measure == "AS") {
            if (transf.char == "FALSE") {
                xlab <- "Arcsine Transformed Proportion"
            }
            else {
                xlab <- "Transformed Arcsine Transformed Risk Difference"
                if (transf.char == "transf.iarcsin" || transf.char == 
                  "transf.iarcsin.int") 
                  xlab <- "Proportion"
            }
        }
        if (x$measure == "IRR") {
            if (transf.char == "FALSE") {
                xlab <- "Log Incidence Rate"
            }
            else {
                xlab <- "Transformed Log Incidence Rate"
                if (transf.char == "exp" || transf.char == "transf.exp.int") 
                  xlab <- "Incidence Rate"
            }
        }
        if (x$measure == "IRD") {
            if (transf.char == "FALSE") {
                xlab <- "Incidence Rate"
            }
            else {
                xlab <- "Transformed Incidence Rate"
            }
        }
        if (x$measure == "IRSD") {
            if (transf.char == "FALSE") {
                xlab <- "Square-Root Transformed Incidence Rate"
            }
            else {
                xlab <- "Transformed Square-Root Transformed Incidence Rate"
                if (transf.char == "transf.isqrt" || transf.char == 
                  "transf.isqrt.int") 
                  xlab <- "Incidence Rate"
            }
        }
    }
    if (missing(ylab)) {
        ylab <- "Observed Outcome"
        if (x$measure == "OR") {
            if (transf.char == "FALSE") {
                ylab <- "Log Odds"
            }
            else {
                ylab <- "Transformed Log Odds"
                if (transf.char == "exp" || transf.char == "transf.exp.int") 
                  ylab <- "Odds"
                if (transf.char == "transf.ilogit" || transf.char == 
                  "transf.ilogit.int") 
                  ylab <- "Proportion"
            }
        }
        if (x$measure == "RR") {
            if (transf.char == "FALSE") {
                ylab <- "Log Proportion"
            }
            else {
                ylab <- "Transformed Log Proportion"
                if (transf.char == "exp" || transf.char == "transf.exp.int") 
                  ylab <- "Proportion"
            }
        }
        if (x$measure == "RD") {
            if (transf.char == "FALSE") {
                ylab <- "Proportion"
            }
            else {
                ylab <- "Transformed Proportion"
            }
        }
        if (x$measure == "AS") {
            if (transf.char == "FALSE") {
                ylab <- "Arcsine Transformed Proportion"
            }
            else {
                ylab <- "Transformed Arcsine Transformed Risk Difference"
                if (transf.char == "transf.iarcsin" || transf.char == 
                  "transf.iarcsin.int") 
                  ylab <- "Proportion"
            }
        }
        if (x$measure == "IRR") {
            if (transf.char == "FALSE") {
                ylab <- "Log Incidence Rate"
            }
            else {
                ylab <- "Transformed Log Incidence Rate"
                if (transf.char == "exp" || transf.char == "transf.exp.int") 
                  ylab <- "Incidence Rate"
            }
        }
        if (x$measure == "IRD") {
            if (transf.char == "FALSE") {
                ylab <- "Incidence Rate"
            }
            else {
                ylab <- "Transformed Incidence Rate"
            }
        }
        if (x$measure == "IRSD") {
            if (transf.char == "FALSE") {
                ylab <- "Square-Root Transformed Incidence Rate"
            }
            else {
                ylab <- "Transformed Square-Root Transformed Incidence Rate"
                if (transf.char == "transf.isqrt" || transf.char == 
                  "transf.isqrt.int") 
                  ylab <- "Incidence Rate"
            }
        }
    }
    plot(NA, NA, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, 
        cex = psize, pch = pch, bg = bg, ...)
    abline(a = 0, b = 1, ...)
    lines(c.vals, t.vals, lty = "dashed", ...)
    points(dat.c$yi, dat.t$yi, cex = psize, pch = pch, bg = bg, 
        ...)
    invisible()
}
leave1out <-
function (x, ...) 
UseMethod("leave1out")
leave1out.rma.mh <-
function (x, digits = x$digits, transf = FALSE, ...) 
{
    if (!is.element("rma.mh", class(x))) 
        stop("Argument 'x' must be an object of class \"rma.mh\".")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    if (x$k == 1) 
        stop("Stopped because k = 1.")
    o.warn <- getOption("warn")
    on.exit(options(warn = o.warn))
    options(warn = -1)
    b <- rep(NA, x$k.f)
    se <- rep(NA, x$k.f)
    zval <- rep(NA, x$k.f)
    pval <- rep(NA, x$k.f)
    ci.lb <- rep(NA, x$k.f)
    ci.ub <- rep(NA, x$k.f)
    QE <- rep(NA, x$k.f)
    QEp <- rep(NA, x$k.f)
    for (i in seq.int(x$k.f)[x$not.na]) {
        if (is.element(x$measure, c("RR", "OR", "RD"))) {
            res <- try(rma.mh(ai = x$ai.f[-i], bi = x$bi.f[-i], 
                ci = x$ci.f[-i], di = x$di.f[-i], measure = x$measure, 
                add = x$add, to = x$to, ...), silent = TRUE)
        }
        else {
            res <- try(rma.mh(x1i = x$x1i.f[-i], x2i = x$x2i.f[-i], 
                t1i = x$t1i.f[-i], t2i = x$t2i.f[-i], measure = x$measure, 
                add = x$add, to = x$to, ...), silent = TRUE)
        }
        if (inherits(res, "try-error")) 
            next
        b[i] <- res$b
        se[i] <- res$se
        zval[i] <- res$zval
        pval[i] <- res$pval
        ci.lb[i] <- res$ci.lb
        ci.ub[i] <- res$ci.ub
        QE[i] <- res$QE
        QEp[i] <- res$QEp
    }
    if (transf) {
        if (x$measure == "OR" || x$measure == "RR" || x$measure == 
            "IRR") {
            b <- exp(b)
            se <- rep(NA, x$k.f)
            ci.lb <- exp(ci.lb)
            ci.ub <- exp(ci.ub)
        }
    }
    if (na.act == "na.omit") {
        out <- list(estimate = b[x$not.na], se = se[x$not.na], 
            zval = zval[x$not.na], pval = pval[x$not.na], ci.lb = ci.lb[x$not.na], 
            ci.ub = ci.ub[x$not.na], Q = QE[x$not.na], Qp = QEp[x$not.na])
        out$slab <- x$slab[x$not.na]
    }
    if (na.act == "na.exclude" || na.act == "na.pass") {
        out <- list(estimate = b, se = se, zval = zval, pval = pval, 
            ci.lb = ci.lb, ci.ub = ci.ub, Q = QE, Qp = QEp)
        out$slab <- x$slab
    }
    if (na.act == "na.fail" && any(!x$not.na)) 
        stop("Missing values in results.")
    out$digits <- digits
    class(out) <- c("list.rma")
    return(out)
}
leave1out.rma.peto <-
function (x, digits = x$digits, transf = FALSE, ...) 
{
    if (!is.element("rma.peto", class(x))) 
        stop("Argument 'x' must be an object of class \"rma.peto\".")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    if (x$k == 1) 
        stop("Stopped because k = 1.")
    o.warn <- getOption("warn")
    on.exit(options(warn = o.warn))
    options(warn = -1)
    b <- rep(NA, x$k.f)
    se <- rep(NA, x$k.f)
    zval <- rep(NA, x$k.f)
    pval <- rep(NA, x$k.f)
    ci.lb <- rep(NA, x$k.f)
    ci.ub <- rep(NA, x$k.f)
    QE <- rep(NA, x$k.f)
    QEp <- rep(NA, x$k.f)
    for (i in seq.int(x$k.f)[x$not.na]) {
        res <- try(rma.peto(ai = x$ai.f[-i], bi = x$bi.f[-i], 
            ci = x$ci.f[-i], di = x$di.f[-i], add = x$add, to = x$to, 
            ...), silent = TRUE)
        if (inherits(res, "try-error")) 
            next
        b[i] <- res$b
        se[i] <- res$se
        zval[i] <- res$zval
        pval[i] <- res$pval
        ci.lb[i] <- res$ci.lb
        ci.ub[i] <- res$ci.ub
        QE[i] <- res$QE
        QEp[i] <- res$QEp
    }
    if (transf) {
        b <- exp(b)
        se <- rep(NA, x$k.f)
        ci.lb <- exp(ci.lb)
        ci.ub <- exp(ci.ub)
    }
    if (na.act == "na.omit") {
        out <- list(estimate = b[x$not.na], se = se[x$not.na], 
            zval = zval[x$not.na], pval = pval[x$not.na], ci.lb = ci.lb[x$not.na], 
            ci.ub = ci.ub[x$not.na], Q = QE[x$not.na], Qp = QEp[x$not.na])
        out$slab <- x$slab[x$not.na]
    }
    if (na.act == "na.exclude" || na.act == "na.pass") {
        out <- list(estimate = b, se = se, zval = zval, pval = pval, 
            ci.lb = ci.lb, ci.ub = ci.ub, Q = QE, Qp = QEp)
        out$slab <- x$slab
    }
    if (na.act == "na.fail" && any(!x$not.na)) 
        stop("Missing values in results.")
    out$digits <- digits
    class(out) <- c("list.rma")
    return(out)
}
leave1out.rma.uni <-
function (x, digits = x$digits, transf = FALSE, targs, ...) 
{
    if (!is.element("rma.uni", class(x))) 
        stop("Argument 'x' must be an object of class \"rma.uni\".")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    if (!x$int.only) 
        stop("Method only applicable for models without moderators.")
    if (x$k == 1) 
        stop("Stopped because k = 1.")
    if (missing(targs)) 
        targs <- NULL
    o.warn <- getOption("warn")
    on.exit(options(warn = o.warn))
    options(warn = -1)
    b <- rep(NA, x$k.f)
    se <- rep(NA, x$k.f)
    zval <- rep(NA, x$k.f)
    pval <- rep(NA, x$k.f)
    ci.lb <- rep(NA, x$k.f)
    ci.ub <- rep(NA, x$k.f)
    QE <- rep(NA, x$k.f)
    QEp <- rep(NA, x$k.f)
    tau2 <- rep(NA, x$k.f)
    I2 <- rep(NA, x$k.f)
    H2 <- rep(NA, x$k.f)
    for (i in seq.int(x$k.f)[x$not.na]) {
        res <- try(rma(x$yi.f[-i], x$vi.f[-i], method = x$method, 
            weighted = x$weighted, intercept = TRUE, knha = x$knha, 
            control = x$control, ...), silent = TRUE)
        if (inherits(res, "try-error")) 
            next
        b[i] <- res$b
        se[i] <- res$se
        zval[i] <- res$zval
        pval[i] <- res$pval
        ci.lb[i] <- res$ci.lb
        ci.ub[i] <- res$ci.ub
        QE[i] <- res$QE
        QEp[i] <- res$QEp
        tau2[i] <- res$tau2
        I2[i] <- res$I2
        H2[i] <- res$H2
    }
    if (is.function(transf)) {
        if (is.null(targs)) {
            b <- sapply(b, transf)
            se <- rep(NA, x$k.f)
            ci.lb <- sapply(ci.lb, transf)
            ci.ub <- sapply(ci.ub, transf)
        }
        else {
            b <- sapply(b, transf, targs)
            se <- rep(NA, x$k.f)
            ci.lb <- sapply(ci.lb, transf, targs)
            ci.ub <- sapply(ci.ub, transf, targs)
        }
    }
    ci.bounds <- cbind(ci.lb, ci.ub)
    rev.order <- ifelse(ci.ub < ci.lb, TRUE, FALSE)
    rev.order[is.na(rev.order)] <- FALSE
    ci.bounds[rev.order] <- ci.bounds[rev.order, 2:1]
    ci.lb <- ci.bounds[, 1]
    ci.ub <- ci.bounds[, 2]
    if (na.act == "na.omit") {
        out <- list(estimate = b[x$not.na], se = se[x$not.na], 
            zval = zval[x$not.na], pval = pval[x$not.na], ci.lb = ci.lb[x$not.na], 
            ci.ub = ci.ub[x$not.na], Q = QE[x$not.na], Qp = QEp[x$not.na], 
            tau2 = tau2[x$not.na], I2 = I2[x$not.na], H2 = H2[x$not.na])
        out$slab <- x$slab[x$not.na]
    }
    if (na.act == "na.exclude" || na.act == "na.pass") {
        out <- list(estimate = b, se = se, zval = zval, pval = pval, 
            ci.lb = ci.lb, ci.ub = ci.ub, Q = QE, Qp = QEp, tau2 = tau2, 
            I2 = I2, H2 = H2)
        out$slab <- x$slab
    }
    if (na.act == "na.fail" && any(!x$not.na)) 
        stop("Missing values in results.")
    if (x$method == "FE") 
        out <- out[-c(9, 10, 11)]
    out$digits <- digits
    class(out) <- c("list.rma")
    return(out)
}
llplot <-
function (measure = "OR", ai, bi, ci, di, n1i, n2i, data, subset, 
    drop00 = TRUE, xvals = 1000, xlim, ylim, xlab, ylab, scale = TRUE, 
    lty, lwd, col, level = 99.99, refline = 0, ...) 
{
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    if (!require(BiasedUrn)) 
        stop("Please install the 'BiasedUrn' package to use this function.")
    if (missing(xlab)) 
        xlab <- "Log Odds Ratio"
    if (missing(ylab)) {
        if (scale) {
            ylab <- "Scaled Log Likelihood"
        }
        else {
            ylab <- "Log Likelihood"
        }
    }
    alpha <- (100 - level)/100
    if (missing(data)) 
        data <- NULL
    if (is.null(data)) {
        data <- sys.frame(sys.parent())
    }
    else {
        if (!is.data.frame(data)) {
            data <- data.frame(data)
        }
    }
    mf <- match.call()
    mf.subset <- mf[[match("subset", names(mf))]]
    mf.lty <- mf[[match("lty", names(mf))]]
    mf.lwd <- mf[[match("lwd", names(mf))]]
    mf.col <- mf[[match("col", names(mf))]]
    subset <- eval(mf.subset, data, enclos = sys.frame(sys.parent()))
    lty <- eval(mf.lty, data, enclos = sys.frame(sys.parent()))
    lwd <- eval(mf.lwd, data, enclos = sys.frame(sys.parent()))
    col <- eval(mf.col, data, enclos = sys.frame(sys.parent()))
    mf.ai <- mf[[match("ai", names(mf))]]
    mf.bi <- mf[[match("bi", names(mf))]]
    mf.ci <- mf[[match("ci", names(mf))]]
    mf.di <- mf[[match("di", names(mf))]]
    mf.n1i <- mf[[match("n1i", names(mf))]]
    mf.n2i <- mf[[match("n2i", names(mf))]]
    ai <- eval(mf.ai, data, enclos = sys.frame(sys.parent()))
    bi <- eval(mf.bi, data, enclos = sys.frame(sys.parent()))
    ci <- eval(mf.ci, data, enclos = sys.frame(sys.parent()))
    di <- eval(mf.di, data, enclos = sys.frame(sys.parent()))
    n1i <- eval(mf.n1i, data, enclos = sys.frame(sys.parent()))
    n2i <- eval(mf.n2i, data, enclos = sys.frame(sys.parent()))
    if (is.null(bi)) 
        bi <- n1i - ai
    if (is.null(di)) 
        di <- n2i - ci
    dat <- escalc(measure = "OR", ai = ai, bi = bi, ci = ci, 
        di = di, drop00 = drop00)
    yi <- dat$yi
    vi <- dat$vi
    k <- length(ai)
    ids <- seq.int(k)
    if (!is.null(lty)) {
        if (length(lty) == 1L) {
            lty <- rep(lty, k)
        }
        else {
            if (length(lty) != k) 
                stop("Length of 'lty' argument does not match data.")
        }
    }
    if (!is.null(lwd)) {
        if (length(lwd) == 1L) {
            lwd <- rep(lwd, k)
        }
        else {
            if (length(lwd) != k) 
                stop("Length of 'lwd' argument does not match data.")
        }
    }
    if (!is.null(col)) {
        if (length(col) == 1L) {
            col <- rep(col, k)
        }
        else {
            if (length(col) != k) 
                stop("Length of 'col' argument does not match data.")
        }
    }
    id0 <- c(ai == 0L | bi == 0L | ci == 0L | di == 0L)
    id0[is.na(id0)] <- FALSE
    id00 <- c(ai == 0L & ci == 0L) | c(bi == 0L & di == 0L)
    id00[is.na(id00)] <- FALSE
    if (drop00) {
        ai[id00] <- NA
        bi[id00] <- NA
        ci[id00] <- NA
        di[id00] <- NA
    }
    if (!is.null(subset)) {
        yi <- yi[subset]
        vi <- vi[subset]
        ai <- ai[subset]
        bi <- bi[subset]
        ci <- ci[subset]
        di <- di[subset]
        ids <- ids[subset]
        lty <- lty[subset]
        lwd <- lwd[subset]
        col <- col[subset]
        id0 <- id0[subset]
        id00 <- id00[subset]
        k <- length(ai)
    }
    aibicidi.na <- is.na(cbind(ai, bi, ci, di))
    if (any(aibicidi.na)) {
        not.na <- apply(aibicidi.na, MARGIN = 1, sum) == 0L
        if (na.act == "na.omit" || na.act == "na.exclude" || 
            na.act == "na.pass") {
            yi <- yi[not.na]
            vi <- vi[not.na]
            ai <- ai[not.na]
            bi <- bi[not.na]
            ci <- ci[not.na]
            di <- di[not.na]
            ids <- ids[not.na]
            lty <- lty[not.na]
            lwd <- lwd[not.na]
            col <- col[not.na]
            id0 <- id0[not.na]
            id00 <- id00[not.na]
            k <- length(ai)
            warning("Studies with NAs omitted from plotting.")
        }
        if (na.act == "na.fail") 
            stop("Missing values in studies.")
    }
    else {
        not.na <- rep(TRUE, k)
    }
    if (k < 1) 
        stop("Processing terminated since k = 0.")
    if (is.null(lty)) 
        lty <- ifelse(id0 | id00, ifelse(id00, "dotted", "dashed"), 
            "solid")
    if (is.null(lwd)) 
        lwd <- seq(from = 0.2, to = 4, length = k)[rank(1/vi)]
    if (is.null(col)) 
        col <- paste("gray", round(seq(from = 0, to = 80, length = k))[rank(vi)], 
            sep = "")
    ci.lb <- yi - qnorm(alpha/2, lower.tail = FALSE) * sqrt(vi)
    ci.ub <- yi + qnorm(alpha/2, lower.tail = FALSE) * sqrt(vi)
    if (missing(xlim)) {
        xlim <- c(min(ci.lb, na.rm = TRUE), max(ci.ub, na.rm = TRUE))
    }
    else {
        xlim <- sort(xlim)
    }
    logORs <- seq(from = xlim[1], to = xlim[2], length.out = xvals)
    lls <- matrix(NA, nrow = k, ncol = xvals)
    out <- matrix(TRUE, nrow = k, ncol = xvals)
    for (i in 1:k) {
        for (j in 1:xvals) {
            lls[i, j] <- .dnchgi(logORs[j], ai = ai[i], bi = bi[i], 
                ci = ci[i], di = di[i], random = FALSE, dnchgcalc = "dFNCHypergeo", 
                dnchgprec = 1e-10)
            if (logORs[j] >= ci.lb[i] & logORs[j] <= ci.ub[i]) 
                out[i, j] <- FALSE
        }
    }
    if (scale) {
        trapezoid <- function(x, y) sum(diff(x) * (y[-1] + y[-length(y)]))/2
        lls.sum <- rep(NA, k)
        for (i in 1:k) {
            lls.sum[i] <- trapezoid(logORs[!is.na(lls[i, ])], 
                lls[i, !is.na(lls[i, ])])
        }
        lls <- apply(lls, 2, "/", lls.sum)
    }
    lls[out] <- NA
    if (missing(ylim)) {
        ylim <- c(0, max(lls, na.rm = TRUE))
    }
    else {
        ylim <- sort(ylim)
    }
    plot(NA, NA, xlim = c(xlim[1], xlim[2]), ylim = ylim, xlab = xlab, 
        ylab = ylab, ...)
    if (is.numeric(refline)) 
        abline(v = refline, lty = "solid", lwd = 2, ...)
    for (i in (1:k)[order(1/vi)]) {
        lines(logORs, lls[i, ], lty = lty[i], lwd = lwd[i], col = col[i], 
            ...)
    }
    invisible(lls)
}
logLik.rma <-
function (object, REML, ...) 
{
    if (!is.element("rma", class(object))) 
        stop("Argument 'object' must be an object of class \"rma\".")
    if (missing(REML)) {
        if (object$method == "REML") {
            REML <- TRUE
        }
        else {
            REML <- FALSE
        }
    }
    if (REML) {
        val <- object$fit.stats$REML[1]
    }
    else {
        val <- object$fit.stats$ML[1]
    }
    attr(val, "nall") <- object$k.eff
    attr(val, "nobs") <- object$k.eff - ifelse(REML, 1, 0) * 
        object$p.eff
    attr(val, "df") <- object$parms
    class(val) <- "logLik"
    return(val)
}
metafor.news <-
function () 
{
    news(package = "metafor")
}
nobs.rma <-
function (object, ...) 
{
    if (!is.element("rma", class(object))) 
        stop("Argument 'object' must be an object of class \"rma\".")
    n.obs <- object$k.eff - ifelse(object$method == "REML", 1, 
        0) * object$p.eff
    return(n.obs)
}
permutest <-
function (x, ...) 
UseMethod("permutest")
permutest.rma.uni <-
function (x, exact = FALSE, iter = 1000, progbar = TRUE, retpermdist = FALSE, 
    digits = x$digits, ...) 
{
    if (!is.element("rma.uni", class(x))) 
        stop("Argument 'x' must be an object of class \"rma.uni\".")
    if (x$int.only) {
        exact.iter <- 2^x$k
    }
    else {
        X <- as.data.frame(x$X)[do.call(order, as.data.frame(x$X)), 
            ]
        indices <- cumsum(c(TRUE, !duplicated(X)[-1]))
        indices <- rep(cumsum(rle(indices)$length) - (rle(indices)$length - 
            1), rle(indices)$length)
        ind.table <- table(indices)
        exact.iter <- round(prod((max(ind.table) + 1):x$k)/prod(factorial(ind.table[-which(ind.table == 
            max(ind.table))[1]])))
    }
    if (exact || (exact.iter <= iter)) {
        exact <- TRUE
        iter <- exact.iter
        if (iter == Inf) 
            stop("Too many iterations required for exact permutation test.\n")
        if (progbar) 
            cat("Running ", iter, " iterations for exact permutation test.\n", 
                sep = "")
    }
    else {
        if (progbar) 
            cat("Running ", iter, " iterations for approximate permutation test.\n", 
                sep = "")
    }
    if (progbar) 
        pbar <- txtProgressBar(min = 0, max = iter, style = 3)
    if (x$int.only) {
        zval.perm <- rep(NA, iter)
        QM.perm <- rep(NA, iter)
        if (exact) {
            signmat <- .gensigns(x$k)
            for (i in seq.int(iter)) {
                res <- try(rma(signmat[i, ] * x$yi, x$vi, method = x$method, 
                  weighted = x$weighted, intercept = TRUE, knha = x$knha, 
                  control = x$control, btt = 1, ...), silent = FALSE)
                if (inherits(res, "try-error")) 
                  next
                zval.perm[i] <- res$zval
                QM.perm[i] <- res$QM
                if (progbar) 
                  setTxtProgressBar(pbar, i)
            }
        }
        else {
            i <- 1
            while (i <= iter) {
                signs <- 2 * rbinom(x$k, 1, 0.5) - 1
                res <- try(rma(signs * x$yi, x$vi, method = x$method, 
                  weighted = x$weighted, intercept = TRUE, knha = x$knha, 
                  control = x$control, btt = 1, ...), silent = FALSE)
                if (inherits(res, "try-error")) 
                  next
                zval.perm[i] <- res$zval
                QM.perm[i] <- res$QM
                i <- i + 1
                if (progbar) 
                  setTxtProgressBar(pbar, i)
            }
        }
        if (x$zval > 0) {
            pval <- 2 * mean(zval.perm >= x$zval, na.rm = TRUE)
        }
        else {
            pval <- 2 * mean(zval.perm <= x$zval, na.rm = TRUE)
        }
        pval[pval > 1] <- 1
        QMp <- mean(QM.perm >= x$QM, na.rm = TRUE)
    }
    else {
        zval.perm <- matrix(NA, nrow = iter, ncol = x$p)
        QM.perm <- rep(NA, iter)
        if (exact) {
            permmat <- .genuperms(indices)
            for (i in seq.int(iter)) {
                res <- try(rma(x$yi, x$vi, mods = cbind(X[permmat[i, 
                  ], ]), method = x$method, weighted = x$weighted, 
                  intercept = FALSE, knha = x$knha, control = x$control, 
                  btt = x$btt, ...), silent = FALSE)
                if (inherits(res, "try-error")) 
                  next
                zval.perm[i, ] <- res$zval
                QM.perm[i] <- res$QM
                if (progbar) 
                  setTxtProgressBar(pbar, i)
            }
        }
        else {
            i <- 1
            while (i <= iter) {
                res <- try(rma(x$yi, x$vi, mods = cbind(X[sample(x$k), 
                  ]), method = x$method, weighted = x$weighted, 
                  intercept = FALSE, knha = x$knha, control = x$control, 
                  btt = x$btt, ...), silent = FALSE)
                if (inherits(res, "try-error")) 
                  next
                zval.perm[i, ] <- res$zval
                QM.perm[i] <- res$QM
                i <- i + 1
                if (progbar) 
                  setTxtProgressBar(pbar, i)
            }
        }
        pval <- rep(NA, x$p)
        for (j in seq.int(x$p)) {
            if (x$zval[j] > 0) {
                pval[j] <- 2 * mean(zval.perm[, j] >= x$zval[j], 
                  na.rm = TRUE)
            }
            else {
                pval[j] <- 2 * mean(zval.perm[, j] <= x$zval[j], 
                  na.rm = TRUE)
            }
        }
        pval[pval > 1] <- 1
        QMp <- mean(QM.perm >= x$QM, na.rm = TRUE)
    }
    if (progbar) 
        close(pbar)
    out <- list(pval = pval, QMp = QMp, b = x$b, se = x$se, zval = x$zval, 
        ci.lb = x$ci.lb, ci.ub = x$ci.ub, QM = x$QM, k = x$k, 
        p = x$p, btt = x$btt, m = x$m, knha = x$knha, robust = x$robust, 
        int.only = x$int.only, digits = digits)
    if (retpermdist) {
        out$QM.perm <- QM.perm
        out$zval.perm <- data.frame(zval.perm)
        names(out$zval.perm) <- colnames(x$X)
    }
    class(out) <- "permutest.rma.uni"
    return(out)
}
plot.infl.rma.uni <-
function (x, plotinf = TRUE, plotdfb = FALSE, dfbnew = FALSE, 
    logcov = TRUE, layout, slab.style = 1, las = 0, pch = 21, 
    bg = "black", bg.infl = "red", col.na = "lightgray", ...) 
{
    if (class(x) != "infl.rma.uni") 
        stop("Argument 'x' must be an object of class \"infl.rma.uni\".")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    any.na <- is.na(cbind(x$inf, x$dfb))
    if (any(any.na)) {
        if (na.act == "na.fail") 
            stop("Missing values in results.")
    }
    if (is.logical(plotinf)) {
        if (plotinf) {
            which.inf <- 1:8
        }
    }
    else {
        which.inf <- plotinf
        which.inf <- which.inf[(which.inf >= 1) & (which.inf <= 
            8)]
        which.inf <- unique(round(which.inf))
        if (length(which.inf) == 0L) 
            stop("Incorrect specification of 'plotinf' argument.")
        plotinf <- TRUE
    }
    if (is.logical(plotdfb)) {
        if (plotdfb) {
            which.dfb <- seq.int(x$p)
        }
    }
    else {
        which.dfb <- plotdfb
        which.dfb <- which.dfb[(which.dfb >= 1) & (which.dfb <= 
            x$p)]
        which.dfb <- unique(round(which.dfb))
        if (length(which.dfb) == 0L) 
            stop("Incorrect specification of 'plotdfb' argument.")
        plotdfb <- TRUE
    }
    if (!plotinf & !plotdfb) 
        stop("At least one of the arguments 'plotinf' or 'plotdfb' argument must be TRUE.")
    if (!plotinf & dfbnew) 
        dfbnew <- FALSE
    par.mar <- par("mar")
    par.mar.adj <- par.mar - c(2, 2, 2, 1)
    par.mar.adj[par.mar.adj < 1] <- 1
    par(mar = par.mar.adj)
    on.exit(par(mar = par.mar))
    lplot <- function(..., minlength, strict) {
        plot(...)
    }
    lpoints <- function(..., minlength, strict) {
        points(...)
    }
    llines <- function(..., minlength, strict) {
        lines(...)
    }
    laxis <- function(..., minlength, strict) {
        axis(...)
    }
    labline <- function(..., minlength, strict) {
        abline(...)
    }
    ids <- switch(slab.style, `1` = x$ids, `2` = dimnames(x$inf)[[1]], 
        `3` = abbreviate(dimnames(x$inf)[[1]], ...))
    if (plotinf) {
        par.mfrow <- par("mfrow")
        on.exit(par(mfrow = par.mfrow), add = TRUE)
        if (missing(layout)) {
            if (length(which.inf) == 2) 
                par(mfrow = c(2, 1))
            if (length(which.inf) == 3) 
                par(mfrow = c(3, 1))
            if (length(which.inf) == 4) 
                par(mfrow = c(2, 2))
            if (length(which.inf) == 5) 
                par(mfrow = c(5, 1))
            if (length(which.inf) == 6) 
                par(mfrow = c(3, 2))
            if (length(which.inf) == 7) 
                par(mfrow = c(7, 1))
            if (length(which.inf) == 8) 
                par(mfrow = c(4, 2))
        }
        else {
            layout <- layout[(layout >= 1)]
            layout <- round(layout)
            if (length(layout) != 2L) 
                stop("Incorrect specification of 'layout' argument.")
            par(mfrow = layout)
        }
        for (i in seq.int(length(which.inf))) {
            if (which.inf[i] == 1) {
                zi <- x$inf$rstudent
                not.na <- !is.na(zi)
                if (na.act == "na.omit") {
                  zi <- zi[not.na]
                  len.ids <- length(x$ids) - sum(!not.na)
                  ids.infl <- x$is.infl[not.na]
                  lab.ids <- ids[not.na]
                }
                if (na.act == "na.exclude" || na.act == "na.pass") {
                  len.ids <- length(x$ids)
                  ids.infl <- x$is.infl
                  lab.ids <- ids
                }
                zi.min <- min(zi, -2, na.rm = TRUE)
                zi.max <- max(zi, 2, na.rm = TRUE)
                lplot(NA, NA, xlim = c(1, len.ids), ylim = c(zi.min, 
                  zi.max), xaxt = "n", main = "rstudent", xlab = "", 
                  ylab = "", las = las, ...)
                laxis(side = 1, at = seq.int(len.ids), labels = lab.ids, 
                  xlab = "", las = las, ...)
                labline(h = 0, lty = "dashed", ...)
                labline(h = c(qnorm(0.025), qnorm(0.975)), lty = "dotted", 
                  ...)
                if (na.act == "na.exclude" || na.act == "na.pass") 
                  llines(seq.int(len.ids)[not.na], zi[not.na], 
                    col = col.na, ...)
                llines(seq.int(len.ids), zi, ...)
                lpoints(seq.int(len.ids), zi, pch = pch, bg = bg, 
                  ...)
                lpoints(seq.int(len.ids)[ids.infl], zi[ids.infl], 
                  bg = bg.infl, pch = pch, ...)
            }
            if (which.inf[i] == 2) {
                zi <- x$inf$dffits
                not.na <- !is.na(zi)
                if (na.act == "na.omit") {
                  zi <- zi[not.na]
                  len.ids <- length(x$ids) - sum(!not.na)
                  ids.infl <- x$is.infl[not.na]
                  lab.ids <- ids[not.na]
                }
                if (na.act == "na.exclude" || na.act == "na.pass") {
                  len.ids <- length(x$ids)
                  ids.infl <- x$is.infl
                  lab.ids <- ids
                }
                zi.min <- min(zi, na.rm = TRUE)
                zi.max <- max(zi, na.rm = TRUE)
                lplot(NA, NA, xlim = c(1, len.ids), ylim = c(zi.min, 
                  zi.max), xaxt = "n", main = "dffits", xlab = "", 
                  ylab = "", las = las, ...)
                laxis(side = 1, at = seq.int(len.ids), labels = lab.ids, 
                  xlab = "", las = las, ...)
                labline(h = 0, lty = "dashed", ...)
                labline(h = 3 * sqrt(x$p/(x$k - x$p)), lty = "dotted", 
                  ...)
                labline(h = -3 * sqrt(x$p/(x$k - x$p)), lty = "dotted", 
                  ...)
                if (na.act == "na.exclude" || na.act == "na.pass") 
                  llines(seq.int(len.ids)[not.na], zi[not.na], 
                    col = col.na, ...)
                llines(seq.int(len.ids), zi, ...)
                lpoints(seq.int(len.ids), zi, pch = pch, bg = bg, 
                  ...)
                lpoints(seq.int(len.ids)[ids.infl], zi[ids.infl], 
                  bg = bg.infl, pch = pch, ...)
            }
            if (which.inf[i] == 3) {
                zi <- x$inf$cook.d
                not.na <- !is.na(zi)
                if (na.act == "na.omit") {
                  zi <- zi[not.na]
                  len.ids <- length(x$ids) - sum(!not.na)
                  ids.infl <- x$is.infl[not.na]
                  lab.ids <- ids[not.na]
                }
                if (na.act == "na.exclude" || na.act == "na.pass") {
                  len.ids <- length(x$ids)
                  ids.infl <- x$is.infl
                  lab.ids <- ids
                }
                zi.min <- min(zi, na.rm = TRUE)
                zi.max <- max(zi, na.rm = TRUE)
                lplot(NA, NA, xlim = c(1, len.ids), ylim = c(zi.min, 
                  zi.max), xaxt = "n", main = "cook.d", xlab = "", 
                  ylab = "", las = las, ...)
                laxis(side = 1, at = seq.int(len.ids), labels = lab.ids, 
                  xlab = "", las = las, ...)
                labline(h = qchisq(0.5, df = x$p), lty = "dotted", 
                  ...)
                if (na.act == "na.exclude" || na.act == "na.pass") 
                  llines(seq.int(len.ids)[not.na], zi[not.na], 
                    col = col.na, ...)
                llines(seq.int(len.ids), zi, ...)
                lpoints(seq.int(len.ids), zi, pch = pch, bg = bg, 
                  ...)
                lpoints(seq.int(len.ids)[ids.infl], zi[ids.infl], 
                  bg = bg.infl, pch = pch, ...)
            }
            if (which.inf[i] == 4) {
                zi <- x$inf$cov.r
                not.na <- !is.na(zi)
                if (na.act == "na.omit") {
                  zi <- zi[not.na]
                  len.ids <- length(x$ids) - sum(!not.na)
                  ids.infl <- x$is.infl[not.na]
                  lab.ids <- ids[not.na]
                }
                if (na.act == "na.exclude" || na.act == "na.pass") {
                  len.ids <- length(x$ids)
                  ids.infl <- x$is.infl
                  lab.ids <- ids
                }
                zi.min <- min(zi, na.rm = TRUE)
                zi.max <- max(zi, na.rm = TRUE)
                if (logcov) {
                  lplot(NA, NA, xlim = c(1, len.ids), ylim = c(zi.min, 
                    zi.max), xaxt = "n", main = "cov.r", xlab = "", 
                    ylab = "", las = las, log = "y", ...)
                }
                else {
                  lplot(NA, NA, xlim = c(1, len.ids), ylim = c(zi.min, 
                    zi.max), xaxt = "n", main = "cov.r", xlab = "", 
                    ylab = "", las = las, ...)
                }
                laxis(side = 1, at = seq.int(len.ids), labels = lab.ids, 
                  xlab = "", las = las, ...)
                labline(h = 1, lty = "dashed", ...)
                if (na.act == "na.exclude" || na.act == "na.pass") 
                  llines(seq.int(len.ids)[not.na], zi[not.na], 
                    col = col.na, ...)
                llines(seq.int(len.ids), zi, ...)
                lpoints(seq.int(len.ids), zi, pch = pch, bg = bg, 
                  ...)
                lpoints(seq.int(len.ids)[ids.infl], zi[ids.infl], 
                  bg = bg.infl, pch = pch, ...)
            }
            if (which.inf[i] == 5) {
                zi <- x$inf$tau2.del
                not.na <- !is.na(zi)
                if (na.act == "na.omit") {
                  zi <- zi[not.na]
                  len.ids <- length(x$ids) - sum(!not.na)
                  ids.infl <- x$is.infl[not.na]
                  lab.ids <- ids[not.na]
                }
                if (na.act == "na.exclude" || na.act == "na.pass") {
                  len.ids <- length(x$ids)
                  ids.infl <- x$is.infl
                  lab.ids <- ids
                }
                zi.min <- min(zi, na.rm = TRUE)
                zi.max <- max(zi, na.rm = TRUE)
                lplot(NA, NA, xlim = c(1, len.ids), ylim = c(zi.min, 
                  zi.max), xaxt = "n", main = "tau2.del", xlab = "", 
                  ylab = "", las = las, ...)
                laxis(side = 1, at = seq.int(len.ids), labels = lab.ids, 
                  xlab = "", las = las, ...)
                labline(h = x$tau2, lty = "dashed", ...)
                if (na.act == "na.exclude" || na.act == "na.pass") 
                  llines(seq.int(len.ids)[not.na], zi[not.na], 
                    col = col.na, ...)
                llines(seq.int(len.ids), zi, ...)
                lpoints(seq.int(len.ids), zi, pch = pch, bg = bg, 
                  ...)
                lpoints(seq.int(len.ids)[ids.infl], zi[ids.infl], 
                  bg = bg.infl, pch = pch, ...)
            }
            if (which.inf[i] == 6) {
                zi <- x$inf$QE.del
                not.na <- !is.na(zi)
                if (na.act == "na.omit") {
                  zi <- zi[not.na]
                  len.ids <- length(x$ids) - sum(!not.na)
                  ids.infl <- x$is.infl[not.na]
                  lab.ids <- ids[not.na]
                }
                if (na.act == "na.exclude" || na.act == "na.pass") {
                  len.ids <- length(x$ids)
                  ids.infl <- x$is.infl
                  lab.ids <- ids
                }
                zi.min <- min(zi, na.rm = TRUE)
                zi.max <- max(zi, na.rm = TRUE)
                lplot(NA, NA, xlim = c(1, len.ids), ylim = c(zi.min, 
                  zi.max), xaxt = "n", main = "QE.del", xlab = "", 
                  ylab = "", las = las, ...)
                laxis(side = 1, at = seq.int(len.ids), labels = lab.ids, 
                  xlab = "", las = las, ...)
                labline(h = x$QE, lty = "dashed", ...)
                if (na.act == "na.exclude" || na.act == "na.pass") 
                  llines(seq.int(len.ids)[not.na], zi[not.na], 
                    col = col.na, ...)
                llines(seq.int(len.ids), zi, ...)
                lpoints(seq.int(len.ids), zi, pch = pch, bg = bg, 
                  ...)
                lpoints(seq.int(len.ids)[ids.infl], zi[ids.infl], 
                  bg = bg.infl, pch = pch, ...)
            }
            if (which.inf[i] == 7) {
                zi <- x$inf$hat
                not.na <- !is.na(zi)
                if (na.act == "na.omit") {
                  zi <- zi[not.na]
                  len.ids <- length(x$ids) - sum(!not.na)
                  ids.infl <- x$is.infl[not.na]
                  lab.ids <- ids[not.na]
                }
                if (na.act == "na.exclude" || na.act == "na.pass") {
                  len.ids <- length(x$ids)
                  ids.infl <- x$is.infl
                  lab.ids <- ids
                }
                zi.min <- 0
                zi.max <- max(zi, na.rm = TRUE)
                lplot(NA, NA, xlim = c(1, len.ids), ylim = c(zi.min, 
                  zi.max), xaxt = "n", main = "hat", xlab = "", 
                  ylab = "", las = las, ...)
                laxis(side = 1, at = seq.int(len.ids), labels = lab.ids, 
                  xlab = "", las = las, ...)
                labline(h = x$p/x$k, lty = "dashed", ...)
                labline(h = 3 * x$p/x$k, lty = "dotted", ...)
                if (na.act == "na.exclude" || na.act == "na.pass") 
                  llines(seq.int(len.ids)[not.na], zi[not.na], 
                    col = col.na, ...)
                llines(seq.int(len.ids), zi, ...)
                lpoints(seq.int(len.ids), zi, pch = pch, bg = bg, 
                  ...)
                lpoints(seq.int(len.ids)[ids.infl], zi[ids.infl], 
                  bg = bg.infl, pch = pch, ...)
            }
            if (which.inf[i] == 8) {
                zi <- x$inf$weight
                not.na <- !is.na(zi)
                if (na.act == "na.omit") {
                  zi <- zi[not.na]
                  len.ids <- length(x$ids) - sum(!not.na)
                  ids.infl <- x$is.infl[not.na]
                  lab.ids <- ids[not.na]
                }
                if (na.act == "na.exclude" || na.act == "na.pass") {
                  len.ids <- length(x$ids)
                  ids.infl <- x$is.infl
                  lab.ids <- ids
                }
                zi.min <- 0
                zi.max <- max(zi, na.rm = TRUE)
                lplot(NA, NA, xlim = c(1, len.ids), ylim = c(zi.min, 
                  zi.max), xaxt = "n", main = "weight", xlab = "", 
                  ylab = "", las = las, ...)
                laxis(side = 1, at = seq.int(len.ids), labels = lab.ids, 
                  xlab = "", las = las, ...)
                labline(h = 100/x$k, lty = "dashed", ...)
                if (na.act == "na.exclude" || na.act == "na.pass") 
                  llines(seq.int(len.ids)[not.na], zi[not.na], 
                    col = col.na, ...)
                llines(seq.int(len.ids), zi, ...)
                lpoints(seq.int(len.ids), zi, pch = pch, bg = bg, 
                  ...)
                lpoints(seq.int(len.ids)[ids.infl], zi[ids.infl], 
                  bg = bg.infl, pch = pch, ...)
            }
        }
    }
    if (plotdfb) {
        if (dfbnew) {
            dev.new()
            par.mar <- par("mar")
            par.mar.adj <- par.mar - c(2, 2, 2, 1)
            par.mar.adj[par.mar.adj < 1] <- 1
            par(mar = par.mar.adj)
            on.exit(par(mar = par.mar), add = TRUE)
        }
        else {
            if (plotinf) {
                par.ask <- par("ask")
                par(ask = TRUE)
                on.exit(par(ask = par.ask), add = TRUE)
            }
        }
        par(mfrow = c(length(which.dfb), 1))
        for (i in seq.int(length(which.dfb))) {
            zi <- x$dfb[, which.dfb[i]]
            not.na <- !is.na(zi)
            if (na.act == "na.omit") {
                zi <- zi[not.na]
                len.ids <- length(x$ids) - sum(!not.na)
                ids.infl <- x$is.infl[not.na]
                lab.ids <- ids[not.na]
            }
            if (na.act == "na.exclude" || na.act == "na.pass") {
                len.ids <- length(x$ids)
                ids.infl <- x$is.infl
                lab.ids <- ids
            }
            lplot(NA, NA, xlim = c(1, len.ids), ylim = range(zi, 
                na.rm = TRUE), xaxt = "n", main = paste("dfb: ", 
                dimnames(x$dfb)[[2]][which.dfb[i]]), xlab = "", 
                ylab = "", las = las, ...)
            laxis(side = 1, at = seq.int(len.ids), labels = lab.ids, 
                xlab = "", las = las, ...)
            labline(h = 0, lty = "dashed", ...)
            labline(h = 1, lty = "dotted", ...)
            labline(h = -1, lty = "dotted", ...)
            if (na.act == "na.exclude" || na.act == "na.pass") 
                llines(seq.int(len.ids)[not.na], zi[not.na], 
                  col = col.na, ...)
            llines(seq.int(len.ids), zi, ...)
            lpoints(seq.int(len.ids), zi, pch = pch, bg = bg, 
                ...)
            lpoints(seq.int(len.ids)[ids.infl], zi[ids.infl], 
                bg = bg.infl, pch = pch, ...)
        }
    }
    invisible()
}
plot.rma.glmm <-
function (x, qqplot = FALSE, ...) 
{
    if (!is.element("rma.glmm", class(x))) 
        stop("Argument 'x' must be an object of class \"rma.glmm\".")
    stop("Method not yet implemented for objects of class \"rma.glmm\". Sorry!")
}
plot.rma.mh <-
function (x, qqplot = FALSE, ...) 
{
    if (!is.element("rma.mh", class(x))) 
        stop("Argument 'x' must be an object of class \"rma.mh\".")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    par.mfrow <- par("mfrow")
    par(mfrow = c(2, 2))
    on.exit(par(mfrow = par.mfrow))
    forest(x, ...)
    title("Forest Plot", ...)
    funnel(x, ...)
    title("Funnel Plot", ...)
    radial(x, ...)
    title("Radial Plot", ...)
    if (qqplot) {
        qqnorm(x, ...)
    }
    else {
        options(na.action = "na.pass")
        z <- rstandard(x)$z
        options(na.action = na.act)
        not.na <- !is.na(z)
        if (na.act == "na.omit") {
            z <- z[not.na]
            ids <- x$ids[not.na]
            not.na <- not.na[not.na]
        }
        if (na.act == "na.exclude" || na.act == "na.pass") 
            ids <- x$ids
        k <- length(z)
        plot(NA, NA, xlim = c(1, k), ylim = c(min(z, -2, na.rm = TRUE), 
            max(z, 2, na.rm = TRUE)), xaxt = "n", xlab = "Study", 
            ylab = "", bty = "l", ...)
        lines(seq.int(k)[not.na], z[not.na], col = "lightgray", 
            ...)
        lines(seq.int(k), z, ...)
        points(seq.int(k), z, pch = 21, bg = "black", ...)
        axis(side = 1, at = seq.int(k), labels = ids, ...)
        abline(h = 0, lty = "dashed", ...)
        abline(h = c(qnorm(0.025), qnorm(0.975)), lty = "dotted", 
            ...)
        title("Standardized Residuals", ...)
    }
    invisible()
}
plot.rma.peto <-
function (x, qqplot = FALSE, ...) 
{
    if (!is.element("rma.peto", class(x))) 
        stop("Argument 'x' must be an object of class \"rma.peto\".")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    par.mfrow <- par("mfrow")
    par(mfrow = c(2, 2))
    on.exit(par(mfrow = par.mfrow))
    forest(x, ...)
    title("Forest Plot", ...)
    funnel(x, ...)
    title("Funnel Plot", ...)
    radial(x, ...)
    title("Radial Plot", ...)
    if (qqplot) {
        qqnorm(x, ...)
    }
    else {
        options(na.action = "na.pass")
        z <- rstandard(x)$z
        options(na.action = na.act)
        not.na <- !is.na(z)
        if (na.act == "na.omit") {
            z <- z[not.na]
            ids <- x$ids[not.na]
            not.na <- not.na[not.na]
        }
        if (na.act == "na.exclude" || na.act == "na.pass") 
            ids <- x$ids
        k <- length(z)
        plot(NA, NA, xlim = c(1, k), ylim = c(min(z, -2, na.rm = TRUE), 
            max(z, 2, na.rm = TRUE)), xaxt = "n", xlab = "Study", 
            ylab = "", bty = "l", ...)
        lines(seq.int(k)[not.na], z[not.na], col = "lightgray", 
            ...)
        lines(seq.int(k), z, ...)
        points(seq.int(k), z, pch = 21, bg = "black", ...)
        axis(side = 1, at = seq.int(k), labels = ids, ...)
        abline(h = 0, lty = "dashed", ...)
        abline(h = c(qnorm(0.025), qnorm(0.975)), lty = "dotted", 
            ...)
        title("Standardized Residuals", ...)
    }
    invisible()
}
plot.rma.uni <-
function (x, qqplot = FALSE, ...) 
{
    if (!is.element("rma.uni", class(x))) 
        stop("Argument 'x' must be an object of class \"rma.uni\".")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    par.mfrow <- par("mfrow")
    par(mfrow = c(2, 2))
    on.exit(par(mfrow = par.mfrow))
    if (x$int.only) {
        forest(x, ...)
        title("Forest Plot", ...)
        funnel(x, ...)
        title("Funnel Plot", ...)
        radial(x, ...)
        title("Radial Plot", ...)
        if (qqplot) {
            qqnorm(x, ...)
        }
        else {
            options(na.action = "na.pass")
            z <- rstandard(x)$z
            options(na.action = na.act)
            not.na <- !is.na(z)
            if (na.act == "na.omit") {
                z <- z[not.na]
                ids <- x$ids[not.na]
                not.na <- not.na[not.na]
            }
            if (na.act == "na.exclude" || na.act == "na.pass") 
                ids <- x$ids
            k <- length(z)
            plot(NA, NA, xlim = c(1, k), ylim = c(min(z, -2, 
                na.rm = TRUE), max(z, 2, na.rm = TRUE)), xaxt = "n", 
                xlab = "Study", ylab = "", bty = "l", ...)
            lines(seq.int(k)[not.na], z[not.na], col = "lightgray", 
                ...)
            lines(seq.int(k), z, ...)
            points(seq.int(k), z, pch = 21, bg = "black", ...)
            axis(side = 1, at = seq.int(k), labels = ids, ...)
            abline(h = 0, lty = "dashed", ...)
            abline(h = c(qnorm(0.025), qnorm(0.975)), lty = "dotted", 
                ...)
            title("Standardized Residuals", ...)
        }
    }
    else {
        forest(x, ...)
        title("Forest Plot", ...)
        funnel(x, ...)
        title("Residual Funnel Plot", ...)
        options(na.action = "na.pass")
        z <- rstandard(x)$z
        pred <- fitted(x)
        options(na.action = na.act)
        plot(pred, z, ylim = c(min(z, -2, na.rm = TRUE), max(z, 
            2, na.rm = TRUE)), pch = 19, bty = "l", xlab = "Fitted Value", 
            ylab = "Standardized Residual", ...)
        abline(h = 0, lty = "dashed", ...)
        abline(h = c(qnorm(0.025), qnorm(0.975)), lty = "dotted", 
            ...)
        title("Fitted vs. Standardized Residuals", ...)
        if (qqplot) {
            qqnorm(x, ...)
        }
        else {
            options(na.action = "na.pass")
            z <- rstandard(x)$z
            options(na.action = na.act)
            not.na <- !is.na(z)
            if (na.act == "na.omit") {
                z <- z[not.na]
                ids <- x$ids[not.na]
                not.na <- not.na[not.na]
            }
            if (na.act == "na.exclude" || na.act == "na.pass") 
                ids <- x$ids
            k <- length(z)
            plot(NA, NA, xlim = c(1, k), ylim = c(min(z, -2, 
                na.rm = TRUE), max(z, 2, na.rm = TRUE)), xaxt = "n", 
                xlab = "Study", ylab = "", bty = "l", ...)
            lines(seq.int(k)[not.na], z[not.na], col = "lightgray", 
                ...)
            lines(seq.int(k), z, ...)
            points(seq.int(k), z, pch = 21, bg = "black", ...)
            axis(side = 1, at = seq.int(k), labels = ids, ...)
            abline(h = 0, lty = "dashed", ...)
            abline(h = c(qnorm(0.025), qnorm(0.975)), lty = "dotted", 
                ...)
            title("Standardized Residuals", ...)
        }
    }
    invisible()
}
predict.rma <-
function (object, newmods, addx = FALSE, level = object$level, 
    digits = object$digits, transf = FALSE, targs, ...) 
{
    if (!is.element("rma", class(object))) 
        stop("Argument 'object' must be an object of class \"rma\".")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    x <- object
    alpha <- (100 - level)/100
    if (missing(newmods)) 
        newmods <- NULL
    if (missing(targs)) 
        targs <- NULL
    if (x$knha || x$robust) {
        crit <- qt(alpha/2, df = x$k - x$p, lower.tail = FALSE)
    }
    else {
        crit <- qnorm(alpha/2, lower.tail = FALSE)
    }
    if (x$int.only && !is.null(newmods)) 
        stop("Cannot specify new moderator values for models without moderators.")
    if (is.null(newmods)) {
        if (x$int.only) {
            k.new <- 1
            X.new <- cbind(1)
        }
        else {
            k.new <- x$k.f
            X.new <- x$X.f
        }
    }
    else {
        if ((!x$intercept && x$p == 1L) || (x$intercept && x$p == 
            2L)) {
            k.new <- length(newmods)
            X.new <- cbind(c(newmods))
        }
        else {
            if (is.vector(newmods) || nrow(newmods) == 1L) {
                k.new <- 1
                X.new <- rbind(newmods)
            }
            else {
                k.new <- NROW(newmods)
                X.new <- cbind(newmods)
            }
        }
        if (x$intercept) {
            X.new <- cbind(intrcpt = rep(1, k.new), X.new)
        }
    }
    pred <- rep(NA, k.new)
    vpred <- rep(NA, k.new)
    for (i in seq.int(k.new)) {
        Xi.new <- matrix(X.new[i, ], nrow = 1)
        pred[i] <- Xi.new %*% x$b
        vpred[i] <- Xi.new %*% tcrossprod(x$vb, Xi.new)
    }
    se <- sqrt(vpred)
    ci.lb <- pred - crit * se
    ci.ub <- pred + crit * se
    cr.lb <- pred - crit * sqrt(vpred + x$tau2)
    cr.ub <- pred + crit * sqrt(vpred + x$tau2)
    if (is.function(transf)) {
        if (is.null(targs)) {
            pred <- sapply(pred, transf)
            se <- rep(NA, k.new)
            ci.lb <- sapply(ci.lb, transf)
            ci.ub <- sapply(ci.ub, transf)
            cr.lb <- sapply(cr.lb, transf)
            cr.ub <- sapply(cr.ub, transf)
        }
        else {
            pred <- sapply(pred, transf, targs)
            se <- rep(NA, k.new)
            ci.lb <- sapply(ci.lb, transf, targs)
            ci.ub <- sapply(ci.ub, transf, targs)
            cr.lb <- sapply(cr.lb, transf, targs)
            cr.ub <- sapply(cr.ub, transf, targs)
        }
    }
    ci.bounds <- cbind(ci.lb, ci.ub)
    rev.order <- ifelse(ci.ub < ci.lb, TRUE, FALSE)
    rev.order[is.na(rev.order)] <- FALSE
    ci.bounds[rev.order] <- ci.bounds[rev.order, 2:1]
    ci.lb <- ci.bounds[, 1]
    ci.ub <- ci.bounds[, 2]
    cr.bounds <- cbind(cr.lb, cr.ub)
    rev.order <- ifelse(cr.ub < cr.lb, TRUE, FALSE)
    rev.order[is.na(rev.order)] <- FALSE
    cr.bounds[rev.order] <- cr.bounds[rev.order, 2:1]
    cr.lb <- cr.bounds[, 1]
    cr.ub <- cr.bounds[, 2]
    if (is.null(newmods) && !x$int.only) {
        slab <- x$slab
    }
    else {
        slab <- seq.int(k.new)
    }
    if (x$int.only) 
        slab <- ""
    if (na.act == "na.omit") {
        not.na <- !is.na(pred)
        out <- list(pred = pred[not.na], se = se[not.na], ci.lb = ci.lb[not.na], 
            ci.ub = ci.ub[not.na], cr.lb = cr.lb[not.na], cr.ub = cr.ub[not.na])
        if (addx) 
            out$X <- matrix(X.new[not.na, ], ncol = x$p)
        out$slab <- slab[not.na]
    }
    if (na.act == "na.exclude" || na.act == "na.pass") {
        out <- list(pred = pred, se = se, ci.lb = ci.lb, ci.ub = ci.ub, 
            cr.lb = cr.lb, cr.ub = cr.ub)
        if (addx) 
            out$X <- matrix(X.new, ncol = x$p)
        out$slab <- slab
    }
    if (addx) 
        colnames(out$X) <- colnames(x$X)
    if (na.act == "na.fail" && any(!x$not.na)) 
        stop("Missing values in results.")
    if (x$method == "FE") {
        out$cr.lb <- NULL
        out$cr.ub <- NULL
    }
    out$digits <- digits
    out$method <- x$method
    class(out) <- c("predict.rma")
    return(out)
}
print.anova.rma.uni <-
function (x, digits = x$digits, ...) 
{
    if (class(x) != "anova.rma.uni") 
        stop("Argument 'x' must be an object of class \"anova.rma.uni\".")
    cutoff <- paste(c(".", rep(0, digits - 1), 1), collapse = "")
    ncutoff <- as.numeric(cutoff)
    res.table <- rbind(c(x$p.f, x$fit.stats.f[3], x$fit.stats.f[4], 
        x$fit.stats.f[1], NA, NA, x$QE.f, x$tau2.f, NA), c(x$p.r, 
        x$fit.stats.r[3], x$fit.stats.r[4], x$fit.stats.r[1], 
        x$LRT, x$pval, x$QE.r, x$tau2.r, NA))
    res.table[, seq.int(2, 9)] <- formatC(res.table[, seq.int(2, 
        9)], digits = digits, format = "f")
    colnames(res.table) <- c("df", "AIC", "BIC", "logLik", "LRT", 
        "pval", "QE", "tau^2", "VAF")
    rownames(res.table) <- c("Full", "Reduced")
    pval <- x$pval
    if (pval > ncutoff) {
        res.table[2, 6] <- formatC(pval, digits = digits, format = "f")
    }
    else {
        res.table[2, 6] <- paste("<", cutoff, sep = "", collapse = "")
    }
    res.table[1, c(5, 6)] <- ""
    res.table[1, 9] <- ""
    res.table[2, 9] <- paste(x$VAF, "%", sep = "")
    if (x$method == "FE") {
        res.table <- res.table[, seq.int(7)]
    }
    print(res.table, quote = FALSE, right = TRUE)
    invisible()
}
print.confint.rma <-
function (x, digits = x$digits, ...) 
{
    if (!is.element("confint.rma", class(x))) 
        stop("Argument 'x' must be an object of class \"confint.rma\".")
    cat("\n")
    if (names(x)[1] == "fixed") {
        res.fixed <- formatC(x$fixed, digits = digits, format = "f")
        print(res.fixed, quote = FALSE, right = TRUE)
    }
    if (names(x)[1] == "random" | names(x)[2] == "random") {
        if (names(x)[1] == "fixed") 
            cat("\n")
        res.random <- formatC(x$random, digits = digits, format = "f")
        print(res.random, quote = FALSE, right = TRUE)
        if (is.na(x$random[1, 2]) && is.na(x$random[1, 3])) 
            message("\nThe upper and lower CI bounds for tau^2 both fall below ", 
                x$tau2.min, ".\nThe CIs given above are therefore equal to the null set.", 
                sep = "")
    }
    cat("\n")
    invisible()
}
print.escalc <-
function (x, digits, ...) 
{
    if (!is.element("escalc", class(x))) 
        stop("Argument 'x' must be an object of class \"escalc\".")
    attr(x, "class") <- NULL
    if (is.null(attr(x, "var.names"))) {
        x <- data.frame(x)
        print(x, ...)
    }
    else {
        if (missing(digits)) 
            digits <- attr(x, "digits")
        if (is.null(digits)) 
            digits <- 4
        yi.pos <- which(names(x) == attr(x, "var.names")[1])
        vi.pos <- which(names(x) == attr(x, "var.names")[2])
        x <- data.frame(x)
        if (length(yi.pos) > 0) 
            x[yi.pos] <- apply(x[yi.pos], 2, formatC, digits = digits, 
                format = "f")
        if (length(vi.pos) > 0) 
            x[vi.pos] <- apply(x[vi.pos], 2, formatC, digits = digits, 
                format = "f")
        print(x, quote = FALSE, right = TRUE)
    }
}
print.fsn <-
function (x, digits = x$digits, ...) 
{
    if (class(x) != "fsn") 
        stop("Argument 'x' must be an object of class \"fsn\".")
    cutoff <- paste(c(".", rep(0, digits - 1), 1), collapse = "")
    ncutoff <- as.numeric(cutoff)
    cat("\n")
    cat("Fail-safe N Calculation Using the", x$type, "Approach", 
        "\n\n")
    if (x$type == "Rosenthal") {
        pval <- x$pval
        if (pval > ncutoff) {
            pval <- formatC(pval, digits = digits, format = "f")
        }
        else {
            pval <- paste("<", cutoff, sep = "", collapse = "")
        }
        cat("Observed Significance Level:", formatC(pval, digits = digits, 
            format = "f"), "\n")
        cat("Target Significance Level:  ", x$alpha, "\n\n")
        cat("Fail-safe N:", x$fsnum, "\n\n")
    }
    if (x$type == "Orwin") {
        cat("Average Effect Size:", formatC(x$meanes, digits = digits, 
            format = "f"), "\n")
        cat("Target Effect Size: ", formatC(x$target, digits = digits, 
            format = "f"), "\n\n")
        cat("Fail-safe N:", x$fsnum, "\n\n")
    }
    if (x$type == "Rosenberg") {
        pval <- x$pval
        if (pval > ncutoff) {
            pval <- formatC(pval, digits = digits, format = "f")
        }
        else {
            pval <- paste("<", cutoff, sep = "", collapse = "")
        }
        cat("Average Effect Size:        ", formatC(x$meanes, 
            digits = digits, format = "f"), "\n")
        cat("Observed Significance Level:", formatC(pval, digits = digits, 
            format = "f"), "\n")
        cat("Target Significance Level:  ", x$alpha, "\n\n")
        cat("Fail-safe N:", x$fsnum, "\n\n")
    }
    invisible()
}
print.infl.rma.uni <-
function (x, digits = x$digits, ...) 
{
    if (class(x) != "infl.rma.uni") 
        stop("Argument 'x' must be an object of class \"infl.rma.uni\".")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    inf <- round(x$inf, digits)
    dfb <- round(x$dfb, digits)
    inf$inf <- ifelse(!is.na(x$is.infl) & x$is.infl, "*", "")
    any.na <- is.na(cbind(inf, dfb))
    if (any(any.na)) {
        if (na.act == "na.omit") {
            inf <- inf[x$not.na, ]
            dfb <- dfb[x$not.na, ]
            out <- list(inf = inf, dfb = dfb)
        }
        if (na.act == "na.exclude" || na.act == "na.pass") {
            out <- list(inf = inf, dfb = dfb)
        }
        if (na.act == "na.fail") 
            stop("Missing values in results.")
    }
    else {
        out <- list(inf = inf, dfb = dfb)
    }
    print(out)
}
print.list.rma <-
function (x, digits = x$digits, ...) 
{
    if (!is.element("list.rma", class(x))) 
        stop("Argument 'x' must be an object of class \"list.rma\".")
    force(digits)
    attr(x, "class") <- NULL
    out <- x[seq.int(which(names(x) == "slab") - 1)]
    out <- data.frame(out, row.names = x$slab)
    out <- apply(out, 2, formatC, digits = digits, format = "f")
    print(out, quote = FALSE, right = TRUE)
}
print.permutest.rma.uni <-
function (x, digits = x$digits, signif.stars = getOption("show.signif.stars"), 
    signif.legend = signif.stars, ...) 
{
    if (class(x) != "permutest.rma.uni") 
        stop("Argument 'x' must be an object of class \"permutest.rma.uni\".")
    cutoff <- paste(c(".", rep(0, digits - 1), 1), collapse = "")
    ncutoff <- as.numeric(cutoff)
    QMp <- x$QMp
    if (QMp > ncutoff) {
        QMp <- paste("=", formatC(QMp, digits = digits, format = "f"))
    }
    else {
        QMp <- paste("< ", cutoff, sep = "", collapse = "")
    }
    cat("\n")
    if (!x$int.only) {
        cat("Test of Moderators (coefficient(s) ", paste(x$btt, 
            collapse = ","), "): \n", sep = "")
        if (x$knha || x$robust) {
            cat("F(df1 = ", x$m, ", df2 = ", x$k - x$p, ") = ", 
                formatC(x$QM, digits = digits, format = "f"), 
                ", p-val* ", QMp, "\n\n", sep = "")
        }
        else {
            cat("QM(df = ", x$m, ") = ", formatC(x$QM, digits = digits, 
                format = "f"), ", p-val* ", QMp, "\n\n", sep = "")
        }
    }
    res.table <- cbind(x$b, x$se, x$zval, x$pval, x$ci.lb, x$ci.ub)
    dimnames(res.table)[[2]] <- c("estimate", "se", "zval", "pval*", 
        "ci.lb", "ci.ub")
    if (x$knha || x$robust) 
        dimnames(res.table)[[2]][3] <- c("tval")
    signif <- symnum(x$pval, corr = FALSE, na = FALSE, cutpoints = c(0, 
        0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", 
        "*", ".", " "))
    if (signif.stars) {
        res.table <- cbind(formatC(res.table, digits = digits, 
            format = "f"), signif)
        dimnames(res.table)[[2]][7] <- ""
    }
    else {
        res.table <- formatC(res.table, digits = digits, format = "f")
    }
    res.table[x$pval > ncutoff, 4] <- formatC(x$pval[x$pval > 
        ncutoff], digits = digits, format = "f")
    res.table[x$pval < ncutoff, 4] <- paste("<", cutoff, sep = "", 
        collapse = "")
    cat("Model Results:")
    cat("\n\n")
    print(res.table, quote = FALSE, right = TRUE, print.gap = 2)
    cat("\n")
    if (signif.legend) 
        cat("---\nSignif. codes: ", attr(signif, "legend"), "\n\n")
    invisible()
}
print.predict.rma <-
function (x, digits = x$digits, ...) 
{
    if (!is.element("predict.rma", class(x))) 
        stop("Argument 'x' must be an object of class \"predict.rma\".")
    force(digits)
    attr(x, "class") <- NULL
    out <- x[seq.int(which(names(x) == "slab") - 1)]
    out <- data.frame(out, row.names = x$slab)
    if (NROW(out) == 0L) 
        stop("All predicted values are NA.", call. = FALSE)
    if (x$method == "FE") {
        out[, seq.int(4)] <- apply(out[, seq.int(4)], 2, formatC, 
            digits = digits, format = "f")
    }
    else {
        out[, seq.int(6)] <- apply(out[, seq.int(6)], 2, formatC, 
            digits = digits, format = "f")
    }
    print(out, quote = FALSE, right = TRUE)
}
print.ranktest.rma <-
function (x, digits = x$digits, ...) 
{
    if (class(x) != "ranktest.rma") 
        stop("Argument 'x' must be an object of class \"ranktest.rma\".")
    cutoff <- paste(c(".", rep(0, digits - 1), 1), collapse = "")
    ncutoff <- as.numeric(cutoff)
    pval <- x$pval
    if (pval > ncutoff) {
        pval <- paste("=", formatC(pval, digits = digits, format = "f"))
    }
    else {
        pval <- paste("< ", cutoff, sep = "", collapse = "")
    }
    cat("\n")
    cat("Rank Correlation Test for Funnel Plot Asymmetry\n\n")
    cat("Kendall's tau = ", formatC(x$tau, digits = digits, format = "f"), 
        ", p ", pval, "\n\n", sep = "")
    invisible()
}
print.regtest.rma <-
function (x, digits = x$digits, ret.fit = x$ret.fit, ...) 
{
    if (class(x) != "regtest.rma") 
        stop("Argument 'x' must be an object of class \"regtest.rma\".")
    cutoff <- paste(c(".", rep(0, digits - 1), 1), collapse = "")
    ncutoff <- as.numeric(cutoff)
    pval <- x$pval
    if (pval > ncutoff) {
        pval <- paste("=", formatC(pval, digits = digits, format = "f"))
    }
    else {
        pval <- paste("< ", cutoff, sep = "", collapse = "")
    }
    cat("\n")
    cat("Regression Test for Funnel Plot Asymmetry\n\n")
    if (x$model == "lm") {
        cat("model:     weighted regression with multiplicative dispersion\n")
    }
    else {
        cat("model:    ", ifelse(x$method == "FE", "fixed-effects", 
            "mixed-effects"), "meta-regression model\n")
    }
    if (x$predictor == "sei") 
        cat("predictor: standard error\n")
    if (x$predictor == "vi") 
        cat("predictor: sampling variance\n")
    if (x$predictor == "ni") 
        cat("predictor: sample size\n")
    if (x$predictor == "ninv") 
        cat("predictor: inverse of the sample size\n")
    if (ret.fit) {
        print(x$fit)
    }
    else {
        cat("\n")
    }
    if (is.na(x$dfs)) {
        cat("test for funnel plot asymmetry: z = ", formatC(x$zval, 
            digits = digits, format = "f"), ", p ", pval, "\n\n", 
            sep = "")
    }
    else {
        cat("test for funnel plot asymmetry: t = ", formatC(x$zval, 
            digits = digits, format = "f"), ", df = ", x$dfs, 
            ", p ", pval, "\n\n", sep = "")
    }
    invisible()
}
print.rma.glmm <-
function (x, digits = x$digits, showfit = FALSE, signif.stars = getOption("show.signif.stars"), 
    signif.legend = signif.stars, ...) 
{
    if (!is.element("rma.glmm", class(x))) 
        stop("Argument 'x' must be an object of class \"rma.glmm\".")
    cutoff <- paste(c(".", rep(0, digits - 1), 1), collapse = "")
    ncutoff <- as.numeric(cutoff)
    cat("\n")
    if (x$method == "FE") {
        if (x$int.only) {
            cat("Fixed-Effects Model (k = ", x$k, ")", sep = "")
        }
        else {
            cat("Fixed-Effects with Moderators Model (k = ", 
                x$k, ")", sep = "")
        }
    }
    else {
        if (x$int.only) {
            cat("Random-Effects Model (k = ", x$k, "; ", sep = "")
        }
        else {
            cat("Mixed-Effects Model (k = ", x$k, "; ", sep = "")
        }
        cat("tau^2 estimator: ", x$method, ")", sep = "")
    }
    if (is.element(x$measure, c("OR", "IRR"))) {
        cat("\n")
        if (x$model == "UM.FS") 
            cat("Model Type: Unconditional Model with Fixed Study Effects")
        if (x$model == "UM.RS") 
            cat("Model Type: Unconditional Model with Random Study Effects")
        if (x$model == "CM.AL") 
            cat("Model Type: Conditional Model with Approximate Likelihood")
        if (x$model == "CM.EL") 
            cat("Model Type: Conditional Model with Exact Likelihood")
    }
    if (showfit) {
        cat("\n")
        fs <- c(formatC(round(x$fit.stats$ML, digits = digits), 
            digits = digits, format = "f"))
        names(fs) <- c("logLik", "deviance", "AIC", "BIC")
        cat("\n")
        print(fs, quote = FALSE, print.gap = 2)
        cat("\n")
    }
    else {
        cat("\n\n")
    }
    if (x$method != "FE") {
        if (x$int.only) {
            cat("tau^2 (estimated amount of total heterogeneity): ", 
                formatC(x$tau2, digits = ifelse(abs(x$tau2) <= 
                  .Machine$double.eps * 10, 0, digits), format = "f"), 
                ifelse(is.na(x$se.tau2), "", paste(" (SE = ", 
                  formatC(x$se.tau2, digits = digits, format = "f"), 
                  ")", sep = "")), "\n", sep = "")
            cat("tau (square root of estimated tau^2 value):      ", 
                ifelse(x$tau2 >= 0, formatC(sqrt(x$tau2), digits = ifelse(x$tau2 <= 
                  .Machine$double.eps * 10, 0, digits), format = "f"), 
                  NA), "\n\n", sep = "")
            cat("I^2 (total heterogeneity / total variability):   ", 
                ifelse(is.na(x$I2), NA, formatC(x$I2, digits = 2, 
                  format = "f")), "%", "\n", sep = "")
            cat("H^2 (total variability / sampling variability):  ", 
                ifelse(is.na(x$H2), NA, formatC(x$H2, digits = 2, 
                  format = "f")), sep = "")
        }
        else {
            cat("tau^2 (estimated amount of residual heterogeneity):     ", 
                formatC(x$tau2, digits = ifelse(abs(x$tau2) <= 
                  .Machine$double.eps * 10, 0, digits), format = "f"), 
                ifelse(is.na(x$se.tau2), "", paste(" (SE = ", 
                  formatC(x$se.tau2, digits = digits, format = "f"), 
                  ")", sep = "")), "\n", sep = "")
            cat("tau (square root of estimated tau^2 value):             ", 
                ifelse(x$tau2 >= 0, formatC(sqrt(x$tau2), digits = ifelse(x$tau2 <= 
                  .Machine$double.eps * 10, 0, digits), format = "f"), 
                  NA), "\n\n", sep = "")
            cat("I^2 (residual heterogeneity / unaccounted variability): ", 
                ifelse(is.na(x$I2), NA, formatC(x$I2, digits = 2, 
                  format = "f")), "%", "\n", sep = "")
            cat("H^2 (unaccounted variability / sampling variability):   ", 
                ifelse(is.na(x$H2), NA, formatC(x$H2, digits = 2, 
                  format = "f")), sep = "")
        }
        cat("\n\n")
    }
    if (!is.na(x$sigma2)) {
        cat("sigma^2 (estimated amount of study level variability): ", 
            formatC(x$sigma2, digits = ifelse(abs(x$sigma2) <= 
                .Machine$double.eps * 10, 0, digits), format = "f"), 
            "\n", sep = "")
        cat("sigma (square root of estimated sigma^2 value):        ", 
            ifelse(x$sigma2 >= 0, formatC(sqrt(x$sigma2), digits = ifelse(x$sigma2 <= 
                .Machine$double.eps * 10, 0, digits), format = "f"), 
                NA), "\n\n", sep = "")
    }
    if (!is.na(x$QE.Wld) && !is.na(x$QE.LRT)) {
        QEp.Wld <- x$QEp.Wld
        QEp.LRT <- x$QEp.LRT
        if (QEp.Wld > ncutoff) {
            QEp.Wld <- paste("=", formatC(QEp.Wld, digits = digits, 
                format = "f"))
        }
        else {
            QEp.Wld <- paste("< ", cutoff, sep = "", collapse = "")
        }
        if (QEp.LRT > ncutoff) {
            QEp.LRT <- paste("=", formatC(QEp.LRT, digits = digits, 
                format = "f"))
        }
        else {
            QEp.LRT <- paste("< ", cutoff, sep = "", collapse = "")
        }
        QE.Wld <- formatC(round(x$QE.Wld, digits = digits), digits = digits, 
            format = "f")
        QE.LRT <- formatC(round(x$QE.LRT, digits = digits), digits = digits, 
            format = "f")
        if (nchar(QE.Wld) > nchar(QE.LRT)) 
            QE.LRT <- paste(paste(rep(" ", nchar(QE.Wld) - nchar(QE.LRT)), 
                collapse = ""), QE.LRT, sep = "")
        if (nchar(QE.LRT) > nchar(QE.Wld)) 
            QE.Wld <- paste(paste(rep(" ", nchar(QE.LRT) - nchar(QE.Wld)), 
                collapse = ""), QE.Wld, sep = "")
        if (x$int.only) {
            cat("Tests for Heterogeneity: \n")
            cat("Wld(df = ", x$QE.df, ") = ", QE.Wld, ", p-val ", 
                QEp.Wld, "\n", sep = "")
            cat("LRT(df = ", x$QE.df, ") = ", QE.LRT, ", p-val ", 
                QEp.LRT, "\n\n", sep = "")
        }
        else {
            cat("Tests for Residual Heterogeneity: \n")
            cat("Wld(df = ", x$QE.df, ") = ", QE.Wld, ", p-val ", 
                QEp.Wld, "\n", sep = "")
            cat("LRT(df = ", x$QE.df, ") = ", QE.LRT, ", p-val ", 
                QEp.LRT, "\n\n", sep = "")
        }
    }
    QMp <- x$QMp
    if (QMp > ncutoff) {
        QMp <- paste("=", formatC(QMp, digits = digits, format = "f"))
    }
    else {
        QMp <- paste("< ", cutoff, sep = "", collapse = "")
    }
    if (x$p > 1) {
        cat("Test of Moderators (coefficient(s) ", paste(x$btt, 
            collapse = ","), "): \n", sep = "")
        if (x$knha || x$robust) {
            cat("F(df1 = ", x$m, ", df2 = ", x$k - x$p, ") = ", 
                formatC(x$QM, digits = digits, format = "f"), 
                ", p-val ", QMp, "\n\n", sep = "")
        }
        else {
            cat("QM(df = ", x$m, ") = ", formatC(x$QM, digits = digits, 
                format = "f"), ", p-val ", QMp, "\n\n", sep = "")
        }
    }
    if (x$int.only) {
        res.table <- c(x$b, x$se, x$zval, x$pval, x$ci.lb, x$ci.ub)
        names(res.table) <- c("estimate", "se", "zval", "pval", 
            "ci.lb", "ci.ub")
        if (x$knha || x$robust) 
            names(res.table)[3] <- c("tval")
        res.table <- formatC(res.table, digits = digits, format = "f")
        signif <- symnum(x$pval, corr = FALSE, na = FALSE, cutpoints = c(0, 
            0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", 
            "*", ".", " "))
        if (signif.stars) {
            res.table <- c(formatC(res.table, digits = digits, 
                format = "f"), signif)
            names(res.table)[7] <- ""
        }
        else {
            res.table <- formatC(res.table, digits = digits, 
                format = "f")
        }
        res.table[4][x$pval > ncutoff] <- formatC(x$pval[x$pval > 
            ncutoff], digits = digits, format = "f")
        res.table[4][x$pval < ncutoff] <- paste("<", cutoff, 
            sep = "", collapse = "")
    }
    else {
        res.table <- cbind(x$b, x$se, x$zval, x$pval, x$ci.lb, 
            x$ci.ub)
        dimnames(res.table)[[2]] <- c("estimate", "se", "zval", 
            "pval", "ci.lb", "ci.ub")
        if (x$knha || x$robust) 
            dimnames(res.table)[[2]][3] <- c("tval")
        signif <- symnum(x$pval, corr = FALSE, na = FALSE, cutpoints = c(0, 
            0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", 
            "*", ".", " "))
        if (signif.stars) {
            res.table <- cbind(formatC(res.table, digits = digits, 
                format = "f"), signif)
            dimnames(res.table)[[2]][7] <- ""
        }
        else {
            res.table <- formatC(res.table, digits = digits, 
                format = "f")
        }
        res.table[x$pval > ncutoff, 4] <- formatC(x$pval[x$pval > 
            ncutoff], digits = digits, format = "f")
        res.table[x$pval < ncutoff, 4] <- paste("<", cutoff, 
            sep = "", collapse = "")
    }
    cat("Model Results:")
    cat("\n\n")
    if (x$int.only) {
        print(res.table, quote = FALSE, right = TRUE)
    }
    else {
        print(res.table, quote = FALSE, right = TRUE, print.gap = 2)
    }
    cat("\n")
    if (signif.legend) 
        cat("---\nSignif. codes: ", attr(signif, "legend"), "\n\n")
    invisible()
}
print.rma.mh <-
function (x, digits = x$digits, showfit = FALSE, ...) 
{
    if (!is.element("rma.mh", class(x))) 
        stop("Argument 'x' must be an object of class \"rma.mh\".")
    cutoff <- paste(c(".", rep(0, digits - 1), 1), collapse = "")
    ncutoff <- as.numeric(cutoff)
    cat("\n")
    cat("Fixed-Effects Model (k = ", x$k, ")", sep = "")
    if (showfit) {
        cat("\n")
        fs <- c(formatC(x$fit.stats$ML, digits = digits, format = "f"))
        names(fs) <- c("logLik", "deviance", "AIC", "BIC")
        cat("\n")
        print(fs, quote = FALSE, print.gap = 2)
        cat("\n")
    }
    else {
        cat("\n\n")
    }
    if (!is.na(x$QE)) {
        QEp <- x$QEp
        if (QEp > ncutoff) {
            QEp <- paste("=", formatC(QEp, digits = digits, format = "f"))
        }
        else {
            QEp <- paste("< ", cutoff, sep = "", collapse = "")
        }
        cat("Test for Heterogeneity: \n")
        cat("Q(df = ", x$k.yi - 1, ") = ", formatC(x$QE, digits = digits, 
            format = "f"), ", p-val ", QEp, sep = "")
    }
    if (x$measure == "OR" || x$measure == "RR" || x$measure == 
        "IRR") {
        res.table <- c(x$b, x$se, x$zval, x$pval, x$ci.lb, x$ci.ub)
        res.table.exp <- c(exp(x$b), exp(x$ci.lb), exp(x$ci.ub))
        if (!is.na(x$b)) {
            res.table <- formatC(res.table, digits = digits, 
                format = "f")
            res.table[4][x$pval > ncutoff] <- formatC(x$pval[x$pval > 
                ncutoff], digits = digits, format = "f")
            res.table[4][x$pval < ncutoff] <- paste("<", cutoff, 
                sep = "", collapse = "")
        }
        if (!is.na(x$b)) {
            res.table.exp <- formatC(res.table.exp, digits = digits, 
                format = "f")
        }
        names(res.table) <- c("estimate", "se", "zval", "pval", 
            "ci.lb", "ci.ub")
        names(res.table.exp) <- c("estimate", "ci.lb", "ci.ub")
        cat("\n\n")
        cat("Model Results (log scale):")
        cat("\n\n")
        print(res.table, quote = FALSE, right = TRUE)
        cat("\n")
        cat("Model Results (", x$measure, " scale):", sep = "")
        cat("\n\n")
        print(res.table.exp, quote = FALSE, right = TRUE)
        cat("\n")
        if (x$measure == "OR") {
            if (is.na(x$CMH)) {
                cat("Cochran-Mantel-Haenszel Test:     test value not computable for these data \n", 
                  sep = "")
            }
            else {
                pval <- x$CMHp
                if (pval > ncutoff) {
                  pval <- paste("=", formatC(pval, digits = digits, 
                    format = "f"))
                }
                else {
                  pval <- paste("< ", cutoff, sep = "", collapse = "")
                }
                cat("Cochran-Mantel-Haenszel Test:     CMH = ", 
                  formatC(x$CMH, digits, format = "f"), ", df = 1,", 
                  paste(rep(" ", nchar(x$k.pos) - 1, collapse = "")), 
                  " p-val ", pval, "\n", sep = "")
            }
            if (is.na(x$TAp)) {
                cat("Tarone's Test for Heterogeneity:  test value not computable for these data \n\n", 
                  sep = "")
            }
            else {
                pval <- x$TAp
                if (pval > ncutoff) {
                  pval <- paste("=", formatC(pval, digits = digits, 
                    format = "f"))
                }
                else {
                  pval <- paste("< ", cutoff, sep = "", collapse = "")
                }
                cat("Tarone's Test for Heterogeneity:  X^2 = ", 
                  formatC(x$TA, digits, format = "f"), ", df = ", 
                  x$k.pos - 1, ", p-val ", pval, "\n\n", sep = "")
            }
        }
    }
    else {
        res.table <- c(x$b, x$se, x$zval, x$pval, x$ci.lb, x$ci.ub)
        if (!is.na(x$b)) {
            res.table <- formatC(res.table, digits = digits, 
                format = "f")
            res.table[4][x$pval > ncutoff] <- formatC(x$pval[x$pval > 
                ncutoff], digits = digits, format = "f")
            res.table[4][x$pval < ncutoff] <- paste("<", cutoff, 
                sep = "", collapse = "")
        }
        names(res.table) <- c("estimate", "se", "zval", "pval", 
            "ci.lb", "ci.ub")
        cat("\n\n")
        cat("Model Results:")
        cat("\n\n")
        print(res.table, quote = FALSE, right = TRUE)
    }
    invisible()
}
print.rma.peto <-
function (x, digits = x$digits, showfit = FALSE, ...) 
{
    if (!is.element("rma.peto", class(x))) 
        stop("Argument 'x' must be an object of class \"rma.peto\".")
    cutoff <- paste(c(".", rep(0, digits - 1), 1), collapse = "")
    ncutoff <- as.numeric(cutoff)
    cat("\n")
    cat("Fixed-Effects Model (k = ", x$k, ")", sep = "")
    if (showfit) {
        cat("\n")
        fs <- c(formatC(x$fit.stats$ML, digits = digits, format = "f"))
        names(fs) <- c("logLik", "deviance", "AIC", "BIC")
        cat("\n")
        print(fs, quote = FALSE, print.gap = 2)
        cat("\n")
    }
    else {
        cat("\n\n")
    }
    if (!is.na(x$QE)) {
        QEp <- x$QEp
        if (QEp > ncutoff) {
            QEp <- paste("=", formatC(QEp, digits = digits, format = "f"))
        }
        else {
            QEp <- paste("< ", cutoff, sep = "", collapse = "")
        }
        cat("Test for Heterogeneity: \n")
        cat("Q(df = ", x$k.pos - 1, ") = ", formatC(x$QE, digits = digits, 
            format = "f"), ", p-val ", QEp, sep = "")
    }
    res.table <- c(x$b, x$se, x$zval, x$pval, x$ci.lb, x$ci.ub)
    res.table.exp <- c(exp(x$b), exp(x$ci.lb), exp(x$ci.ub))
    if (!is.na(x$b)) {
        res.table <- formatC(res.table, digits = digits, format = "f")
        res.table[4][x$pval > ncutoff] <- formatC(x$pval[x$pval > 
            ncutoff], digits = digits, format = "f")
        res.table[4][x$pval < ncutoff] <- paste("<", cutoff, 
            sep = "", collapse = "")
    }
    if (!is.na(x$b)) {
        res.table.exp <- formatC(res.table.exp, digits = digits, 
            format = "f")
    }
    names(res.table) <- c("estimate", "se", "zval", "pval", "ci.lb", 
        "ci.ub")
    names(res.table.exp) <- c("estimate", "ci.lb", "ci.ub")
    cat("\n\n")
    cat("Model Results (log scale):")
    cat("\n\n")
    print(res.table, quote = FALSE, right = TRUE)
    cat("\n")
    cat("Model Results (OR scale):", sep = "")
    cat("\n\n")
    print(res.table.exp, quote = FALSE, right = TRUE)
    cat("\n")
    invisible()
}
print.rma.uni <-
function (x, digits = x$digits, showfit = FALSE, signif.stars = getOption("show.signif.stars"), 
    signif.legend = signif.stars, ...) 
{
    if (!is.element("rma.uni", class(x))) 
        stop("Argument 'x' must be an object of class \"rma.uni\".")
    if (is.element("rma.uni.trimfill", class(x))) {
        cat("\n")
        if (x$k0 <= 0) {
            cat("Estimated number of missing studies on the", 
                x$side, "side is zero.\n")
        }
        cat("Estimated number of missing studies on the", x$side, 
            "side:", x$k0, "\n")
    }
    cutoff <- paste(c(".", rep(0, digits - 1), 1), collapse = "")
    ncutoff <- as.numeric(cutoff)
    cat("\n")
    if (x$method == "FE") {
        if (x$int.only) {
            cat("Fixed-Effects Model (k = ", x$k, ")", sep = "")
        }
        else {
            cat("Fixed-Effects with Moderators Model (k = ", 
                x$k, ")", sep = "")
        }
    }
    else {
        if (x$int.only) {
            cat("Random-Effects Model (k = ", x$k, "; ", sep = "")
        }
        else {
            cat("Mixed-Effects Model (k = ", x$k, "; ", sep = "")
        }
        cat("tau^2 estimator: ", x$method, ")", sep = "")
    }
    if (showfit) {
        cat("\n")
        if (x$method == "REML") {
            fs <- c(formatC(round(x$fit.stats$REML, digits = digits), 
                digits = digits, format = "f"))
            names(fs) <- c("logLik", "deviance", "AIC", "BIC")
        }
        else {
            fs <- c(formatC(round(x$fit.stats$ML, digits = digits), 
                digits = digits, format = "f"))
            names(fs) <- c("logLik", "deviance", "AIC", "BIC")
        }
        cat("\n")
        print(fs, quote = FALSE, print.gap = 2)
        cat("\n")
    }
    else {
        cat("\n\n")
    }
    if (x$method != "FE") {
        if (x$int.only) {
            cat("tau^2 (estimated amount of total heterogeneity): ", 
                formatC(x$tau2, digits = ifelse(abs(x$tau2) <= 
                  .Machine$double.eps * 10, 0, digits), format = "f"), 
                ifelse(is.na(x$se.tau2), "", paste(" (SE = ", 
                  formatC(x$se.tau2, digits = digits, format = "f"), 
                  ")", sep = "")), "\n", sep = "")
            cat("tau (square root of estimated tau^2 value):      ", 
                ifelse(x$tau2 >= 0, formatC(sqrt(x$tau2), digits = ifelse(x$tau2 <= 
                  .Machine$double.eps * 10, 0, digits), format = "f"), 
                  NA), "\n", sep = "")
            cat("I^2 (total heterogeneity / total variability):   ", 
                ifelse(is.na(x$I2), NA, formatC(x$I2, digits = 2, 
                  format = "f")), "%", "\n", sep = "")
            cat("H^2 (total variability / sampling variability):  ", 
                ifelse(is.na(x$H2), NA, formatC(x$H2, digits = 2, 
                  format = "f")), sep = "")
        }
        else {
            cat("tau^2 (estimated amount of residual heterogeneity):     ", 
                formatC(x$tau2, digits = ifelse(abs(x$tau2) <= 
                  .Machine$double.eps * 10, 0, digits), format = "f"), 
                ifelse(is.na(x$se.tau2), "", paste(" (SE = ", 
                  formatC(x$se.tau2, digits = digits, format = "f"), 
                  ")", sep = "")), "\n", sep = "")
            cat("tau (square root of estimated tau^2 value):             ", 
                ifelse(x$tau2 >= 0, formatC(sqrt(x$tau2), digits = ifelse(x$tau2 <= 
                  .Machine$double.eps * 10, 0, digits), format = "f"), 
                  NA), "\n", sep = "")
            cat("I^2 (residual heterogeneity / unaccounted variability): ", 
                ifelse(is.na(x$I2), NA, formatC(x$I2, digits = 2, 
                  format = "f")), "%", "\n", sep = "")
            cat("H^2 (unaccounted variability / sampling variability):   ", 
                ifelse(is.na(x$H2), NA, formatC(x$H2, digits = 2, 
                  format = "f")), sep = "")
        }
        cat("\n\n")
    }
    if (!is.na(x$QE)) {
        QEp <- x$QEp
        if (QEp > ncutoff) {
            QEp <- paste("=", formatC(QEp, digits = digits, format = "f"))
        }
        else {
            QEp <- paste("< ", cutoff, sep = "", collapse = "")
        }
        if (x$int.only) {
            cat("Test for Heterogeneity: \n")
            cat("Q(df = ", x$k - x$p, ") = ", formatC(x$QE, digits = digits, 
                format = "f"), ", p-val ", QEp, "\n\n", sep = "")
        }
        else {
            cat("Test for Residual Heterogeneity: \n")
            cat("QE(df = ", x$k - x$p, ") = ", formatC(x$QE, 
                digits = digits, format = "f"), ", p-val ", QEp, 
                "\n\n", sep = "")
        }
    }
    QMp <- x$QMp
    if (QMp > ncutoff) {
        QMp <- paste("=", formatC(QMp, digits = digits, format = "f"))
    }
    else {
        QMp <- paste("< ", cutoff, sep = "", collapse = "")
    }
    if (x$p > 1) {
        cat("Test of Moderators (coefficient(s) ", paste(x$btt, 
            collapse = ","), "): \n", sep = "")
        if (x$knha || x$robust) {
            cat("F(df1 = ", x$m, ", df2 = ", x$k - x$p, ") = ", 
                formatC(x$QM, digits = digits, format = "f"), 
                ", p-val ", QMp, "\n\n", sep = "")
        }
        else {
            cat("QM(df = ", x$m, ") = ", formatC(x$QM, digits = digits, 
                format = "f"), ", p-val ", QMp, "\n\n", sep = "")
        }
    }
    if (x$int.only) {
        res.table <- c(x$b, x$se, x$zval, x$pval, x$ci.lb, x$ci.ub)
        names(res.table) <- c("estimate", "se", "zval", "pval", 
            "ci.lb", "ci.ub")
        if (x$knha || x$robust) 
            names(res.table)[3] <- c("tval")
        res.table <- formatC(res.table, digits = digits, format = "f")
        signif <- symnum(x$pval, corr = FALSE, na = FALSE, cutpoints = c(0, 
            0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", 
            "*", ".", " "))
        if (signif.stars) {
            res.table <- c(formatC(res.table, digits = digits, 
                format = "f"), signif)
            names(res.table)[7] <- ""
        }
        else {
            res.table <- formatC(res.table, digits = digits, 
                format = "f")
        }
        res.table[4][x$pval > ncutoff] <- formatC(x$pval[x$pval > 
            ncutoff], digits = digits, format = "f")
        res.table[4][x$pval < ncutoff] <- paste("<", cutoff, 
            sep = "", collapse = "")
    }
    else {
        res.table <- cbind(x$b, x$se, x$zval, x$pval, x$ci.lb, 
            x$ci.ub)
        dimnames(res.table)[[2]] <- c("estimate", "se", "zval", 
            "pval", "ci.lb", "ci.ub")
        if (x$knha || x$robust) 
            dimnames(res.table)[[2]][3] <- c("tval")
        signif <- symnum(x$pval, corr = FALSE, na = FALSE, cutpoints = c(0, 
            0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", 
            "*", ".", " "))
        if (signif.stars) {
            res.table <- cbind(formatC(res.table, digits = digits, 
                format = "f"), signif)
            dimnames(res.table)[[2]][7] <- ""
        }
        else {
            res.table <- formatC(res.table, digits = digits, 
                format = "f")
        }
        res.table[x$pval > ncutoff, 4] <- formatC(x$pval[x$pval > 
            ncutoff], digits = digits, format = "f")
        res.table[x$pval < ncutoff, 4] <- paste("<", cutoff, 
            sep = "", collapse = "")
    }
    cat("Model Results:")
    cat("\n\n")
    if (x$int.only) {
        print(res.table, quote = FALSE, right = TRUE)
    }
    else {
        print(res.table, quote = FALSE, right = TRUE, print.gap = 2)
    }
    cat("\n")
    if (signif.legend) 
        cat("---\nSignif. codes: ", attr(signif, "legend"), "\n\n")
    invisible()
}
print.summary.escalc <-
function (x, digits, ...) 
{
    if (!is.element("summary.escalc", class(x))) 
        stop("Argument 'x' must be an object of class \"summary.escalc\".")
    attr(x, "class") <- NULL
    if (is.null(attr(x, "var.names")) || is.null(attr(x, "out.names"))) {
        x <- data.frame(x)
        print(x, ...)
    }
    else {
        if (missing(digits)) 
            digits <- attr(x, "digits")
        if (is.null(digits)) 
            digits <- 4
        yi.pos <- which(names(x) == attr(x, "var.names")[1])
        vi.pos <- which(names(x) == attr(x, "var.names")[2])
        sei.pos <- which(names(x) == attr(x, "out.names")[1])
        zi.pos <- which(names(x) == attr(x, "out.names")[2])
        ci.lb.pos <- which(names(x) == attr(x, "out.names")[3])
        ci.ub.pos <- which(names(x) == attr(x, "out.names")[4])
        x <- data.frame(x)
        if (length(yi.pos) > 0) 
            x[yi.pos] <- apply(x[yi.pos], 2, formatC, digits = digits, 
                format = "f")
        if (length(vi.pos) > 0) 
            x[vi.pos] <- apply(x[vi.pos], 2, formatC, digits = digits, 
                format = "f")
        if (length(sei.pos) > 0) 
            x[sei.pos] <- apply(x[sei.pos], 2, formatC, digits = digits, 
                format = "f")
        if (length(zi.pos) > 0) 
            x[zi.pos] <- apply(x[zi.pos], 2, formatC, digits = digits, 
                format = "f")
        if (length(ci.lb.pos) > 0) 
            x[ci.lb.pos] <- apply(x[ci.lb.pos], 2, formatC, digits = digits, 
                format = "f")
        if (length(ci.ub.pos) > 0) 
            x[ci.ub.pos] <- apply(x[ci.ub.pos], 2, formatC, digits = digits, 
                format = "f")
        print(x, quote = FALSE, right = TRUE)
    }
}
print.summary.rma <-
function (x, digits = x$digits, showfit = TRUE, signif.stars = getOption("show.signif.stars"), 
    signif.legend = signif.stars, ...) 
{
    if (!is.element("summary.rma", class(x))) 
        stop("Argument 'x' must be an object of class \"summary.rma\".")
    class(x) <- class(x)[-1]
    print(x, digits = digits, showfit = showfit, signif.stars = signif.stars, 
        signif.legend = signif.legend, ...)
    invisible()
}
qqnorm.rma.glmm <-
function (y, ...) 
{
    if (!is.element("rma.glmm", class(y))) 
        stop("Argument 'y' must be an object of class \"rma.glmm\".")
    stop("Method not yet implemented for objects of class \"rma.glmm\". Sorry!")
}
qqnorm.rma.mh <-
function (y, type = "rstandard", pch = 19, label = FALSE, offset = 0.3, 
    ...) 
{
    if (!is.element("rma.mh", class(y))) 
        stop("Argument 'y' must be an object of class \"rma.mh\".")
    x <- y
    type <- match.arg(type, c("rstandard", "rstudent"))
    if (x$k == 1) 
        stop("Stopped because k = 1.")
    if (length(label) != 1) 
        stop("Argument 'label' should be of length 1.")
    if (type == "rstandard") {
        res <- rstandard(x)
        not.na <- !is.na(res$z)
        zi <- res$z[not.na]
        slab <- res$slab[not.na]
        pos <- order(zi)
        slab <- slab[pos]
    }
    else {
        res <- rstudent(x)
        not.na <- !is.na(res$z)
        zi <- res$z[not.na]
        slab <- res$slab[not.na]
        pos <- order(zi)
        slab <- slab[pos]
    }
    sav <- qqnorm(zi, pch = pch, bty = "l", ...)
    abline(a = 0, b = 1, lty = "solid", ...)
    if (is.numeric(label)) {
        label <- round(label)
        if (label < 1 | label > x$k) 
            stop("Out of range value for 'label' argument.")
        pos.x <- sav$x[pos]
        pos.y <- sav$y[pos]
        dev <- abs(pos.x - pos.y)
        for (i in seq.int(x$k)) {
            if (sum(dev > dev[i]) < label) 
                text(pos.x[i], pos.y[i], slab[i], pos = ifelse(pos.x[i] >= 
                  0, 2, 4), offset = offset, ...)
        }
    }
    else {
        if (is.logical(label)) 
            label <- ifelse(label, "all", "none")
        pos.x <- sav$x[pos]
        pos.y <- sav$y[pos]
        if (label != "none") {
            for (i in seq.int(x$k)) {
                text(pos.x[i], pos.y[i], slab[i], pos = ifelse(pos.x[i] >= 
                  0, 2, 4), offset = offset, ...)
            }
        }
    }
    invisible(sav)
}
qqnorm.rma.peto <-
function (y, type = "rstandard", pch = 19, label = FALSE, offset = 0.3, 
    ...) 
{
    if (!is.element("rma.peto", class(y))) 
        stop("Argument 'y' must be an y of class \"rma.peto\".")
    x <- y
    type <- match.arg(type, c("rstandard", "rstudent"))
    if (x$k == 1) 
        stop("Stopped because k = 1.")
    if (length(label) != 1) 
        stop("Argument 'label' should be of length 1.")
    if (type == "rstandard") {
        res <- rstandard(x)
        not.na <- !is.na(res$z)
        zi <- res$z[not.na]
        slab <- res$slab[not.na]
        pos <- order(zi)
        slab <- slab[pos]
    }
    else {
        res <- rstudent(x)
        not.na <- !is.na(res$z)
        zi <- res$z[not.na]
        slab <- res$slab[not.na]
        pos <- order(zi)
        slab <- slab[pos]
    }
    sav <- qqnorm(zi, pch = pch, bty = "l", ...)
    abline(a = 0, b = 1, lty = "solid", ...)
    if (is.numeric(label)) {
        label <- round(label)
        if (label < 1 | label > x$k) 
            stop("Out of range value for 'label' argument.")
        pos.x <- sav$x[pos]
        pos.y <- sav$y[pos]
        dev <- abs(pos.x - pos.y)
        for (i in seq.int(x$k)) {
            if (sum(dev > dev[i]) < label) 
                text(pos.x[i], pos.y[i], slab[i], pos = ifelse(pos.x[i] >= 
                  0, 2, 4), offset = offset, ...)
        }
    }
    else {
        if (is.logical(label)) 
            label <- ifelse(label, "all", "none")
        pos.x <- sav$x[pos]
        pos.y <- sav$y[pos]
        if (label != "none") {
            for (i in seq.int(x$k)) {
                text(pos.x[i], pos.y[i], slab[i], pos = ifelse(pos.x[i] >= 
                  0, 2, 4), offset = offset, ...)
            }
        }
    }
    invisible(sav)
}
qqnorm.rma.uni <-
function (y, type = "rstandard", pch = 19, envelope = TRUE, level = y$level, 
    bonferroni = FALSE, reps = 1000, smooth = TRUE, bass = 0, 
    label = FALSE, offset = 0.3, ...) 
{
    if (!is.element("rma.uni", class(y))) 
        stop("Argument 'y' must be an y of class \"rma.uni\".")
    x <- y
    type <- match.arg(type, c("rstandard", "rstudent"))
    if (x$k == 1) 
        stop("Stopped because k = 1.")
    draw.envelope <- envelope
    if (label == "out" & !envelope) {
        envelope <- TRUE
        draw.envelope <- FALSE
    }
    if (length(label) != 1) 
        stop("Argument 'label' should be of length 1.")
    if (type == "rstandard") {
        res <- rstandard(x)
        not.na <- !is.na(res$z)
        zi <- res$z[not.na]
        slab <- res$slab[not.na]
        pos <- order(zi)
        slab <- slab[pos]
    }
    else {
        res <- rstudent(x)
        not.na <- !is.na(res$z)
        zi <- res$z[not.na]
        slab <- res$slab[not.na]
        pos <- order(zi)
        slab <- slab[pos]
    }
    sav <- qqnorm(zi, pch = pch, bty = "l", ...)
    abline(a = 0, b = 1, lty = "solid", ...)
    if (envelope) {
        alpha <- (100 - level)/100
        dat <- matrix(rnorm(x$k * reps), nrow = x$k, ncol = reps)
        if (x$weighted) {
            wi <- 1/(x$vi + x$tau2)
            W <- .diag(wi)
            stXWX <- .invcalc(X = x$X, W = W, k = x$k)
            H <- x$X %*% stXWX %*% crossprod(x$X, W)
        }
        else {
            stXX <- .invcalc(X = x$X, W = diag(x$k), k = x$k)
            H <- x$X %*% tcrossprod(stXX, x$X)
        }
        ImH <- diag(x$k) - H
        ei <- ImH %*% dat
        ei <- apply(ei, 2, sort)
        if (bonferroni) {
            lb <- apply(ei, 1, quantile, (alpha/2)/x$k)
            ub <- apply(ei, 1, quantile, 1 - (alpha/2)/x$k)
        }
        else {
            lb <- apply(ei, 1, quantile, (alpha/2))
            ub <- apply(ei, 1, quantile, 1 - (alpha/2))
        }
        temp.lb <- qqnorm(lb, plot.it = FALSE)
        if (smooth) 
            temp.lb <- supsmu(temp.lb$x, temp.lb$y, bass = bass)
        if (draw.envelope) 
            lines(temp.lb$x, temp.lb$y, lty = "dotted", ...)
        temp.ub <- qqnorm(ub, plot.it = FALSE)
        if (smooth) 
            temp.ub <- supsmu(temp.ub$x, temp.ub$y, bass = bass)
        if (draw.envelope) 
            lines(temp.ub$x, temp.ub$y, lty = "dotted", ...)
    }
    if (is.numeric(label)) {
        label <- round(label)
        if (label < 1 | label > x$k) 
            stop("Out of range value for 'label' argument.")
        pos.x <- sav$x[pos]
        pos.y <- sav$y[pos]
        dev <- abs(pos.x - pos.y)
        for (i in seq.int(x$k)) {
            if (sum(dev > dev[i]) < label) 
                text(pos.x[i], pos.y[i], slab[i], pos = ifelse(pos.x[i] >= 
                  0, 2, 4), offset = offset, ...)
        }
    }
    else {
        if (is.logical(label)) 
            label <- ifelse(label, "out", "none")
        pos.x <- sav$x[pos]
        pos.y <- sav$y[pos]
        if (label != "none") {
            for (i in seq.int(x$k)) {
                if (label == "all") {
                  text(pos.x[i], pos.y[i], slab[i], pos = ifelse(pos.x[i] >= 
                    0, 2, 4), offset = offset, ...)
                }
                else {
                  if (pos.y[i] < temp.lb$y[i] || pos.y[i] > temp.ub$y[i]) 
                    text(pos.x[i], pos.y[i], slab[i], pos = ifelse(pos.x[i] >= 
                      0, 2, 4), offset = offset, ...)
                }
            }
        }
    }
    invisible(sav)
}
radial <-
function (x, ...) 
UseMethod("radial")
radial.rma <-
function (x, center = FALSE, xlim = NULL, zlim, xlab, zlab, atz, 
    aty, steps = 7, level = x$level, digits = 2, back = "lightgray", 
    transf = FALSE, targs, pch = 19, arc.res = 100, cex, ...) 
{
    if (!is.element("rma", class(x))) 
        stop("Argument 'x' must be an object of class \"rma\".")
    if (missing(targs)) 
        targs <- NULL
    if (missing(atz)) 
        atz <- NULL
    if (missing(aty)) 
        aty <- NULL
    if (x$int.only) {
        yi <- x$yi
        yi.c <- yi
        vi <- x$vi
        b <- c(x$b)
        ci.lb <- x$ci.lb
        ci.ub <- x$ci.ub
        tau2 <- c(x$tau2)
        if (is.null(aty)) {
            atyis <- range(yi)
        }
        else {
            atyis <- range(aty)
            aty.c <- aty
        }
    }
    else {
        stop("Radial plots only applicable for models without moderators.")
    }
    if (center) {
        yi <- yi - x$b
        b <- 0
        ci.lb <- ci.lb - x$b
        ci.ub <- ci.ub - x$b
        atyis <- atyis - x$b
        if (!is.null(aty)) 
            aty <- aty - x$b
    }
    alpha <- (100 - level)/100
    zcrit <- qnorm(alpha/2, lower.tail = FALSE)
    zi <- yi/sqrt(vi + tau2)
    xi <- 1/sqrt(vi + tau2)
    if (missing(xlim)) {
        xlims <- c(0, (1.3 * max(xi)))
    }
    else {
        xlims <- sort(xlim)
    }
    ci.xpos <- xlims[2] + 0.12 * (xlims[2] - xlims[1])
    ya.xpos <- xlims[2] + 0.14 * (xlims[2] - xlims[1])
    xaxismax <- xlims[2]
    if (missing(zlim)) {
        zlims <- c(min(-5, 1.1 * min(zi), 1.1 * ci.lb * ci.xpos, 
            1.1 * min(atyis) * ya.xpos, 1.1 * min(yi) * ya.xpos, 
            -1.1 * zcrit + xaxismax * b), max(5, 1.1 * max(zi), 
            1.1 * ci.ub * ci.xpos, 1.1 * max(atyis) * ya.xpos, 
            1.1 * max(yi) * ya.xpos, 1.1 * zcrit + xaxismax * 
                b))
    }
    else {
        zlims <- sort(zlim)
    }
    par.mar <- par("mar")
    par.mar.adj <- par.mar - c(0, -3, 0, -5)
    par.mar.adj[par.mar.adj < 1] <- 1
    par(mar = par.mar.adj)
    on.exit(par(mar = par.mar))
    if (missing(xlab)) {
        if (x$method == "FE") {
            xlab <- expression(x[i] == 1/sqrt(v[i]), ...)
        }
        else {
            xlab <- expression(x[i] == 1/sqrt(v[i] + tau^2), 
                ...)
        }
    }
    par.pty <- par("pty")
    par(pty = "s")
    on.exit(par(pty = par.pty), add = TRUE)
    plot(NA, NA, ylim = zlims, xlim = xlims, bty = "n", xaxt = "n", 
        yaxt = "n", xlab = xlab, ylab = "", xaxs = "i", yaxs = "i", 
        ...)
    if (missing(cex)) 
        cex <- par("cex")
    polygon(c(0, xaxismax, xaxismax, 0), c(zcrit, zcrit + xaxismax * 
        b, -zcrit + xaxismax * b, -zcrit), border = NA, col = back, 
        ...)
    segments(0, 0, xaxismax, xaxismax * b, lty = "solid", ...)
    segments(0, -zcrit, xaxismax, -zcrit + xaxismax * b, lty = "dotted", 
        ...)
    segments(0, zcrit, xaxismax, zcrit + xaxismax * b, lty = "dotted", 
        ...)
    axis(side = 1, ...)
    if (is.null(atz)) {
        axis(side = 2, at = seq(-4, 4, length = 9), labels = NA, 
            las = 1, tcl = par("tcl")/2, ...)
        axis(side = 2, at = seq(-2, 2, length = 3), las = 1, 
            ...)
    }
    else {
        axis(side = 2, at = atz, labels = atz, las = 1, ...)
    }
    if (missing(zlab)) {
        if (center) {
            if (x$method == "FE") {
                mtext(expression(z[i] == frac(y[i] - hat(theta), 
                  sqrt(v[i]))), side = 2, line = par.mar.adj[2] - 
                  1, at = 0, adj = 0, las = 1, cex = cex, ...)
            }
            else {
                mtext(expression(z[i] == frac(y[i] - hat(mu), 
                  sqrt(v[i] + tau^2))), side = 2, line = par.mar.adj[2] - 
                  1, adj = 0, at = 0, las = 1, cex = cex, ...)
            }
        }
        else {
            if (x$method == "FE") {
                mtext(expression(z[i] == frac(y[i], sqrt(v[i]))), 
                  side = 2, line = par.mar.adj[2] - 2, at = 0, 
                  adj = 0, las = 1, cex = cex, ...)
            }
            else {
                mtext(expression(z[i] == frac(y[i], sqrt(v[i] + 
                  tau^2))), side = 2, line = par.mar.adj[2] - 
                  1, at = 0, adj = 0, las = 1, cex = cex, ...)
            }
        }
    }
    else {
        mtext(zlab, side = 2, line = par.mar.adj[2] - 4, at = 0, 
            cex = cex, ...)
    }
    par.xpd <- par("xpd")
    par(xpd = TRUE)
    par.usr <- par("usr")
    asp.rat <- (par.usr[4] - par.usr[3])/(par.usr[2] - par.usr[1])
    if (length(arc.res) == 1L) 
        arc.res <- c(arc.res, arc.res/4)
    if (is.null(aty)) {
        atyis <- seq(min(yi), max(yi), length = arc.res[1])
    }
    else {
        atyis <- seq(min(aty), max(aty), length = arc.res[1])
    }
    len <- ya.xpos
    xis <- rep(NA, length(atyis))
    zis <- rep(NA, length(atyis))
    for (i in seq.int(length(atyis))) {
        xis[i] <- sqrt(len^2/(1 + (atyis[i]/asp.rat)^2))
        zis[i] <- xis[i] * atyis[i]
    }
    valid <- zis > zlims[1] & zis < zlims[2]
    lines(xis[valid], zis[valid], ...)
    if (is.null(aty)) {
        atyis <- seq(min(yi), max(yi), length = steps)
    }
    else {
        atyis <- aty
    }
    len.l <- ya.xpos
    len.u <- ya.xpos + 0.015 * (xlims[2] - xlims[1])
    xis.l <- rep(NA, length(atyis))
    zis.l <- rep(NA, length(atyis))
    xis.u <- rep(NA, length(atyis))
    zis.u <- rep(NA, length(atyis))
    for (i in seq.int(length(atyis))) {
        xis.l[i] <- sqrt(len.l^2/(1 + (atyis[i]/asp.rat)^2))
        zis.l[i] <- xis.l[i] * atyis[i]
        xis.u[i] <- sqrt(len.u^2/(1 + (atyis[i]/asp.rat)^2))
        zis.u[i] <- xis.u[i] * atyis[i]
    }
    valid <- zis.l > zlims[1] & zis.u > zlims[1] & zis.l < zlims[2] & 
        zis.u < zlims[2]
    if (any(valid)) 
        segments(xis.l[valid], zis.l[valid], xis.u[valid], (xis.u * 
            atyis)[valid], ...)
    if (is.null(aty)) {
        atyis <- seq(min(yi), max(yi), length = steps)
        atyis.lab <- seq(min(yi.c), max(yi.c), length = steps)
    }
    else {
        atyis <- aty
        atyis.lab <- aty.c
    }
    len <- ya.xpos + 0.02 * (xlims[2] - xlims[1])
    xis <- rep(NA, length(atyis))
    zis <- rep(NA, length(atyis))
    for (i in seq.int(length(atyis))) {
        xis[i] <- sqrt(len^2/(1 + (atyis[i]/asp.rat)^2))
        zis[i] <- xis[i] * atyis[i]
    }
    if (is.function(transf)) {
        if (is.null(targs)) {
            atyis.lab <- sapply(atyis.lab, transf)
        }
        else {
            atyis.lab <- sapply(atyis.lab, transf, targs)
        }
    }
    valid <- zis > zlims[1] & zis < zlims[2]
    if (any(valid)) 
        text(xis[valid], zis[valid], formatC(atyis.lab[valid], 
            digits = digits, format = "f"), pos = 4, cex = cex, 
            ...)
    atyis <- seq(ci.lb, ci.ub, length = arc.res[2])
    len <- ci.xpos
    xis <- rep(NA, length(atyis))
    zis <- rep(NA, length(atyis))
    for (i in seq.int(length(atyis))) {
        xis[i] <- sqrt(len^2/(1 + (atyis[i]/asp.rat)^2))
        zis[i] <- xis[i] * atyis[i]
    }
    valid <- zis > zlims[1] & zis < zlims[2]
    if (any(valid)) 
        lines(xis[valid], zis[valid], ...)
    atyis <- c(ci.lb, b, ci.ub)
    len.l <- ci.xpos - 0.007 * (xlims[2] - xlims[1])
    len.u <- ci.xpos + 0.007 * (xlims[2] - xlims[1])
    xis.l <- rep(NA, 3)
    zis.l <- rep(NA, 3)
    xis.u <- rep(NA, 3)
    zis.u <- rep(NA, 3)
    for (i in seq.int(length(atyis))) {
        xis.l[i] <- sqrt(len.l^2/(1 + (atyis[i]/asp.rat)^2))
        zis.l[i] <- xis.l[i] * atyis[i]
        xis.u[i] <- sqrt(len.u^2/(1 + (atyis[i]/asp.rat)^2))
        zis.u[i] <- xis.u[i] * atyis[i]
    }
    valid <- zis.l > zlims[1] & zis.u > zlims[1] & zis.l < zlims[2] & 
        zis.u < zlims[2]
    if (any(valid)) 
        segments(xis.l[valid], zis.l[valid], xis.u[valid], (xis.u * 
            atyis)[valid], ...)
    par(xpd = par.xpd)
    points(xi, zi, pch = pch, cex = cex, ...)
    invisible()
}
ranktest <-
function (x, ...) 
UseMethod("ranktest")
ranktest.rma <-
function (x, ...) 
{
    if (!is.element("rma", class(x))) 
        stop("Argument 'x' must be an object of class \"rma\".")
    yi <- x$yi
    vi <- x$vi
    res <- rma(yi, vi, method = "FE")
    b <- res$b
    vb <- res$vb
    vi.star <- vi - c(vb)
    yi.star <- (yi - c(b))/sqrt(vi.star)
    res <- cor.test(yi.star, vi, method = "kendall", exact = TRUE)
    pval <- res$p.value
    tau <- c(res$estimate)
    res <- list(tau, pval, x$digits)
    names(res) <- c("tau", "pval", "digits")
    class(res) <- c("ranktest.rma")
    return(res)
}
regtest <-
function (x, ...) 
UseMethod("regtest")
regtest.rma <-
function (x, model = "rma", predictor = "sei", ret.fit = FALSE, 
    ...) 
{
    if (!is.element("rma", class(x))) 
        stop("Argument 'x' must be an object of class \"rma\".")
    if (is.element("rma.glmm", class(x))) 
        stop("Method not yet implemented for objects of class \"rma.glmm\". Sorry!")
    model <- match.arg(model, c("lm", "rma"))
    predictor <- match.arg(predictor, c("sei", "vi", "ni", "ninv"))
    yi <- x$yi
    vi <- x$vi
    ni <- x$ni
    X <- x$X
    p <- x$p
    if (predictor == "sei") 
        X <- cbind(X, sei = sqrt(vi))
    if (predictor == "vi") 
        X <- cbind(X, vi = vi)
    if (predictor == "ni" || predictor == "ninv") {
        if (is.null(ni)) {
            stop("No sample size information stored in model object.")
        }
        else {
            ni <- c(scale(ni))
            if (predictor == "ni") {
                X <- cbind(X, ni = ni)
            }
            else {
                X <- cbind(X, ninv = 1/ni)
            }
        }
    }
    if (model == "rma") {
        fit <- rma(yi, vi, mods = X, method = x$method, weighted = x$weighted, 
            intercept = FALSE, knha = x$knha, control = x$control, 
            ...)
        zval <- fit$zval[p + 1]
        pval <- fit$pval[p + 1]
        dfs <- ifelse(x$knha || x$robust, fit$k - fit$p, NA)
    }
    else {
        if (NCOL(X) >= x$k) 
            stop("Too few observations to carry out the test.")
        if (x$int.incl) {
            X <- X[, -1, drop = FALSE]
            fit <- lm(yi ~ X, weights = 1/vi)
        }
        else {
            fit <- lm(yi ~ X - 1, weights = 1/vi)
        }
        fit <- summary(fit)
        zval <- coef(fit)[p + 1, 3]
        pval <- coef(fit)[p + 1, 4]
        dfs <- x$k - x$p - 1
    }
    res <- list(model, predictor, zval, pval, dfs, x$method, 
        x$digits, ret.fit, fit)
    names(res) <- c("model", "predictor", "zval", "pval", "dfs", 
        "method", "digits", "ret.fit", "fit")
    class(res) <- c("regtest.rma")
    return(res)
}
residuals.rma <-
function (object, ...) 
{
    if (!is.element("rma", class(object))) 
        stop("Argument 'object' must be an object of class \"rma\".")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    x <- object
    out <- c(x$yi.f - x$X.f %*% x$b)
    out[abs(out) < 100 * .Machine$double.eps] <- 0
    names(out) <- x$slab
    not.na <- !is.na(out)
    if (na.act == "na.omit") 
        out <- out[not.na]
    if (na.act == "na.fail" && any(!not.na)) 
        stop("Missing values in results.")
    return(out)
}
rma <-
function (yi, vi, sei, weights, ai, bi, ci, di, n1i, n2i, x1i, 
    x2i, t1i, t2i, m1i, m2i, sd1i, sd2i, xi, mi, ri, ti, sdi, 
    ni, mods, measure = "GEN", intercept = TRUE, data, slab, 
    subset, add = 1/2, to = "only0", drop00 = FALSE, vtype = "LS", 
    method = "REML", weighted = TRUE, knha = FALSE, level = 95, 
    digits = 4, btt, tau2, verbose = FALSE, control) 
{
    if (!is.element(measure, c("GEN", "RR", "OR", "PETO", "RD", 
        "AS", "PHI", "YUQ", "YUY", "RTET", "PBIT", "OR2D", "IRR", 
        "IRD", "IRSD", "MD", "SMD", "SMDH", "ROM", "RPB", "RBIS", 
        "D2OR", "COR", "UCOR", "ZCOR", "PR", "PLN", "PLO", "PAS", 
        "PFT", "IR", "IRLN", "IRS", "IRFT", "MN", "MC", "SMCC", 
        "SMCR", "ARAW", "AHW", "ABT"))) 
        stop("Unknown 'measure' specified.")
    if (!is.element(method, c("FE", "HS", "HE", "DL", "SJ", "ML", 
        "REML", "EB", "DLit", "SJit", "PM"))) 
        stop("Unknown 'method' specified.")
    if (length(add) > 1) 
        add <- add[1]
    if (length(to) > 1) 
        to <- to[1]
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    if (missing(btt)) 
        btt <- NULL
    if (missing(tau2)) 
        tau2 <- NULL
    if (missing(control)) 
        control <- list()
    if (missing(data)) 
        data <- NULL
    if (is.null(data)) {
        data <- sys.frame(sys.parent())
    }
    else {
        if (!is.data.frame(data)) {
            data <- data.frame(data)
        }
    }
    mf <- match.call()
    mf.slab <- mf[[match("slab", names(mf))]]
    mf.subset <- mf[[match("subset", names(mf))]]
    mf.mods <- mf[[match("mods", names(mf))]]
    slab <- eval(mf.slab, data, enclos = sys.frame(sys.parent()))
    subset <- eval(mf.subset, data, enclos = sys.frame(sys.parent()))
    mods <- eval(mf.mods, data, enclos = sys.frame(sys.parent()))
    ai <- bi <- ci <- di <- x1i <- x2i <- t1i <- t2i <- NA
    mf.yi <- mf[[match("yi", names(mf))]]
    yi <- eval(mf.yi, data, enclos = sys.frame(sys.parent()))
    if (!is.null(yi)) {
        yi.is.formula <- ifelse(class(yi) == "formula", TRUE, 
            FALSE)
        if (yi.is.formula) {
            options(na.action = "na.pass")
            mods <- model.matrix(yi, data = data)
            yi <- model.response(model.frame(yi, data = data))
            options(na.action = na.act)
            names(yi) <- NULL
            intercept <- FALSE
        }
        if (!is.null(attr(yi, "measure"))) 
            measure <- attr(yi, "measure")
        attr(yi, "measure") <- measure
        mf.vi <- mf[[match("vi", names(mf))]]
        mf.sei <- mf[[match("sei", names(mf))]]
        mf.weights <- mf[[match("weights", names(mf))]]
        mf.ni <- mf[[match("ni", names(mf))]]
        vi <- eval(mf.vi, data, enclos = sys.frame(sys.parent()))
        sei <- eval(mf.sei, data, enclos = sys.frame(sys.parent()))
        weights <- eval(mf.weights, data, enclos = sys.frame(sys.parent()))
        ni <- eval(mf.ni, data, enclos = sys.frame(sys.parent()))
        if (is.null(vi)) {
            if (is.null(sei)) {
                if (is.null(weights)) {
                  stop("Need to specify vi, sei, or weights argument.")
                }
                else {
                  vi <- 1/weights
                }
            }
            else {
                vi <- sei^2
            }
        }
        if (is.null(ni) && !is.null(attr(yi, "ni"))) 
            ni <- attr(yi, "ni")
        if (length(yi) != length(vi)) 
            stop("Length of yi and vi (or sei) vectors is not the same.")
        if (!is.null(ni) && (length(yi) != length(ni))) 
            stop("Length of yi and ni vectors is not the same.")
    }
    else {
        if (is.element(measure, c("RR", "OR", "PETO", "RD", "AS", 
            "PHI", "YUQ", "YUY", "RTET", "PBIT", "OR2D"))) {
            mf.ai <- mf[[match("ai", names(mf))]]
            mf.bi <- mf[[match("bi", names(mf))]]
            mf.ci <- mf[[match("ci", names(mf))]]
            mf.di <- mf[[match("di", names(mf))]]
            mf.n1i <- mf[[match("n1i", names(mf))]]
            mf.n2i <- mf[[match("n2i", names(mf))]]
            ai <- eval(mf.ai, data, enclos = sys.frame(sys.parent()))
            bi <- eval(mf.bi, data, enclos = sys.frame(sys.parent()))
            ci <- eval(mf.ci, data, enclos = sys.frame(sys.parent()))
            di <- eval(mf.di, data, enclos = sys.frame(sys.parent()))
            n1i <- eval(mf.n1i, data, enclos = sys.frame(sys.parent()))
            n2i <- eval(mf.n2i, data, enclos = sys.frame(sys.parent()))
            if (is.null(bi)) 
                bi <- n1i - ai
            if (is.null(di)) 
                di <- n2i - ci
            dat <- escalc(measure, ai = ai, bi = bi, ci = ci, 
                di = di, add = add, to = to, drop00 = drop00, 
                vtype = vtype)
        }
        if (is.element(measure, c("IRR", "IRD", "IRSD"))) {
            mf.x1i <- mf[[match("x1i", names(mf))]]
            mf.x2i <- mf[[match("x2i", names(mf))]]
            mf.t1i <- mf[[match("t1i", names(mf))]]
            mf.t2i <- mf[[match("t2i", names(mf))]]
            x1i <- eval(mf.x1i, data, enclos = sys.frame(sys.parent()))
            x2i <- eval(mf.x2i, data, enclos = sys.frame(sys.parent()))
            t1i <- eval(mf.t1i, data, enclos = sys.frame(sys.parent()))
            t2i <- eval(mf.t2i, data, enclos = sys.frame(sys.parent()))
            dat <- escalc(measure, x1i = x1i, x2i = x2i, t1i = t1i, 
                t2i = t2i, add = add, to = to, drop00 = drop00, 
                vtype = vtype)
        }
        if (is.element(measure, c("MD", "SMD", "SMDH", "ROM", 
            "RPB", "RBIS", "D2OR"))) {
            mf.m1i <- mf[[match("m1i", names(mf))]]
            mf.m2i <- mf[[match("m2i", names(mf))]]
            mf.sd1i <- mf[[match("sd1i", names(mf))]]
            mf.sd2i <- mf[[match("sd2i", names(mf))]]
            mf.n1i <- mf[[match("n1i", names(mf))]]
            mf.n2i <- mf[[match("n2i", names(mf))]]
            m1i <- eval(mf.m1i, data, enclos = sys.frame(sys.parent()))
            m2i <- eval(mf.m2i, data, enclos = sys.frame(sys.parent()))
            sd1i <- eval(mf.sd1i, data, enclos = sys.frame(sys.parent()))
            sd2i <- eval(mf.sd2i, data, enclos = sys.frame(sys.parent()))
            n1i <- eval(mf.n1i, data, enclos = sys.frame(sys.parent()))
            n2i <- eval(mf.n2i, data, enclos = sys.frame(sys.parent()))
            dat <- escalc(measure, m1i = m1i, m2i = m2i, sd1i = sd1i, 
                sd2i = sd2i, n1i = n1i, n2i = n2i, vtype = vtype)
        }
        if (is.element(measure, c("COR", "UCOR", "ZCOR"))) {
            mf.ri <- mf[[match("ri", names(mf))]]
            mf.ni <- mf[[match("ni", names(mf))]]
            ri <- eval(mf.ri, data, enclos = sys.frame(sys.parent()))
            ni <- eval(mf.ni, data, enclos = sys.frame(sys.parent()))
            dat <- escalc(measure, ri = ri, ni = ni, vtype = vtype)
        }
        if (is.element(measure, c("PR", "PLN", "PLO", "PAS", 
            "PFT"))) {
            mf.xi <- mf[[match("xi", names(mf))]]
            mf.mi <- mf[[match("mi", names(mf))]]
            mf.ni <- mf[[match("ni", names(mf))]]
            xi <- eval(mf.xi, data, enclos = sys.frame(sys.parent()))
            mi <- eval(mf.mi, data, enclos = sys.frame(sys.parent()))
            ni <- eval(mf.ni, data, enclos = sys.frame(sys.parent()))
            if (is.null(mi)) 
                mi <- ni - xi
            dat <- escalc(measure, xi = xi, mi = mi, add = add, 
                to = to, vtype = vtype)
        }
        if (is.element(measure, c("IR", "IRLN", "IRS", "IRFT"))) {
            mf.xi <- mf[[match("xi", names(mf))]]
            mf.ti <- mf[[match("ti", names(mf))]]
            xi <- eval(mf.xi, data, enclos = sys.frame(sys.parent()))
            ti <- eval(mf.ti, data, enclos = sys.frame(sys.parent()))
            dat <- escalc(measure, xi = xi, ti = ti, add = add, 
                to = to, vtype = vtype)
        }
        if (is.element(measure, c("MN"))) {
            mf.mi <- mf[[match("mi", names(mf))]]
            mf.sdi <- mf[[match("sdi", names(mf))]]
            mf.ni <- mf[[match("ni", names(mf))]]
            mi <- eval(mf.mi, data, enclos = sys.frame(sys.parent()))
            sdi <- eval(mf.sdi, data, enclos = sys.frame(sys.parent()))
            ni <- eval(mf.ni, data, enclos = sys.frame(sys.parent()))
            dat <- escalc(measure, mi = mi, sdi = sdi, ni = ni, 
                vtype = vtype)
        }
        if (is.element(measure, c("MC", "SMCC", "SMCR"))) {
            mf.m1i <- mf[[match("m1i", names(mf))]]
            mf.m2i <- mf[[match("m2i", names(mf))]]
            mf.sd1i <- mf[[match("sd1i", names(mf))]]
            mf.sd2i <- mf[[match("sd2i", names(mf))]]
            mf.ri <- mf[[match("ri", names(mf))]]
            mf.ni <- mf[[match("ni", names(mf))]]
            m1i <- eval(mf.m1i, data, enclos = sys.frame(sys.parent()))
            m2i <- eval(mf.m2i, data, enclos = sys.frame(sys.parent()))
            sd1i <- eval(mf.sd1i, data, enclos = sys.frame(sys.parent()))
            sd2i <- eval(mf.sd2i, data, enclos = sys.frame(sys.parent()))
            ri <- eval(mf.ri, data, enclos = sys.frame(sys.parent()))
            ni <- eval(mf.ni, data, enclos = sys.frame(sys.parent()))
            dat <- escalc(measure, m1i = m1i, m2i = m2i, sd1i = sd1i, 
                sd2i = sd2i, ri = ri, ni = ni, vtype = vtype)
        }
        if (is.element(measure, c("ARAW", "AHW", "ABT"))) {
            mf.ai <- mf[[match("ai", names(mf))]]
            mf.mi <- mf[[match("mi", names(mf))]]
            mf.ni <- mf[[match("ni", names(mf))]]
            ai <- eval(mf.ai, data, enclos = sys.frame(sys.parent()))
            mi <- eval(mf.mi, data, enclos = sys.frame(sys.parent()))
            ni <- eval(mf.ni, data, enclos = sys.frame(sys.parent()))
            dat <- escalc(measure, ai = ai, mi = mi, ni = ni, 
                vtype = vtype)
        }
        if (is.element(measure, c("GEN"))) 
            stop("Specify the desired outcome measure via the 'measure' argument.")
        yi <- dat$yi
        vi <- dat$vi
        ni <- attr(yi, "ni")
    }
    is.formula <- ifelse(class(mods) == "formula", TRUE, FALSE)
    if (is.formula) {
        options(na.action = "na.pass")
        mods <- model.matrix(mods, data = data)
        attr(mods, "assign") <- NULL
        options(na.action = na.act)
        intercept <- FALSE
    }
    if (is.vector(mods)) 
        mods <- cbind(mods)
    if (is.data.frame(mods)) 
        mods <- as.matrix(mods)
    if (!is.null(mods) && (nrow(mods) != length(yi))) 
        stop("Number of rows of the design matrix does not match length of yi.")
    k <- length(yi)
    ids <- seq.int(k)
    if (is.null(slab)) {
        if (!is.null(attr(yi, "slab"))) {
            slab.null <- FALSE
            slab <- attr(yi, "slab")
        }
        else {
            slab.null <- TRUE
            slab <- seq.int(k)
        }
    }
    else {
        if (any(duplicated(slab))) 
            stop("Study labels must be unique.")
        if (length(slab) != length(yi)) 
            stop("Study labels not of same length as data.")
        slab.null <- FALSE
    }
    if (!is.null(subset)) {
        yi <- yi[subset]
        vi <- vi[subset]
        ni <- ni[subset]
        ai <- ai[subset]
        bi <- bi[subset]
        ci <- ci[subset]
        di <- di[subset]
        x1i <- x1i[subset]
        x2i <- x2i[subset]
        t1i <- t1i[subset]
        t2i <- t2i[subset]
        mods <- mods[subset, , drop = FALSE]
        slab <- slab[subset]
        ids <- ids[subset]
        k <- length(yi)
    }
    ai.f <- ai
    bi.f <- bi
    ci.f <- ci
    di.f <- di
    x1i.f <- x1i
    x2i.f <- x2i
    t1i.f <- t1i
    t2i.f <- t2i
    yi.f <- yi
    vi.f <- vi
    ni.f <- ni
    mods.f <- mods
    k.f <- k
    YVM.na <- is.na(cbind(yi, vi, mods))
    if (any(YVM.na)) {
        not.na <- apply(YVM.na, MARGIN = 1, sum) == 0L
        if (na.act == "na.omit" || na.act == "na.exclude" || 
            na.act == "na.pass") {
            yi <- yi[not.na]
            vi <- vi[not.na]
            ni <- ni[not.na]
            mods <- mods[not.na, , drop = FALSE]
            k <- length(yi)
            warning("Studies with NAs omitted from model fitting.")
        }
        if (na.act == "na.fail") 
            stop("Missing values in data.")
    }
    else {
        not.na <- rep(TRUE, k)
    }
    if (k < 1) 
        stop("Processing terminated since k = 0.")
    if (k == 1) {
        method <- "FE"
        knha <- FALSE
    }
    if (any(vi <= 0)) {
        allvipos <- FALSE
        warning("There are outcomes with non-positive sampling variances.")
        if (any(vi < 0)) {
            vi[vi <= 0] <- 0
            warning("Negative sampling variances constrained to zero.")
        }
    }
    else {
        allvipos <- TRUE
    }
    if (is.null(mods) && !intercept) {
        warning("Must either include an intercept (intercept=TRUE) and/or moderators in model.\n  Coerced intercept into the model.")
        intercept <- TRUE
    }
    if (intercept) {
        X <- cbind(intrcpt = rep(1, k), mods)
        X.f <- cbind(intrcpt = rep(1, k.f), mods.f)
    }
    else {
        X <- mods
        X.f <- mods.f
    }
    is.int <- apply(X == 1L, 2, sum) == k
    if (any(is.int)) {
        int.incl <- TRUE
        int.indx <- which(is.int, arr.ind = TRUE)
        X <- cbind(intrcpt = 1, X[, -int.indx, drop = FALSE])
        X.f <- cbind(intrcpt = 1, X.f[, -int.indx, drop = FALSE])
        if (is.formula) 
            intercept <- TRUE
    }
    else {
        int.incl <- FALSE
    }
    p <- NCOL(X)
    if ((p == 1L) && (all(sapply(X, identical, 1)))) {
        int.only <- TRUE
    }
    else {
        int.only <- FALSE
    }
    evstXX <- eigen(t(X) %*% X, symmetric = TRUE)$values
    if (any(evstXX <= .Machine$double.eps)) 
        stop("Design matrix not of full rank. Cannot fit model.")
    if (method == "FE") {
        if (p > k) {
            stop("Number of parameters to be estimated is larger than the number of observations.")
        }
    }
    else {
        if (is.numeric(tau2)) {
            if (p > k) {
                stop("Number of parameters to be estimated is larger than the number of observations.")
            }
        }
        else {
            if ((p + 1) > k) {
                stop("Number of parameters to be estimated is larger than the number of observations.")
            }
        }
    }
    if (is.null(btt)) {
        if (p > 1) {
            if (int.incl) {
                btt <- seq.int(2, p)
            }
            else {
                btt <- seq.int(p)
            }
        }
        else {
            btt <- 1
        }
    }
    else {
        btt <- btt[(btt >= 1) & (btt <= p)]
        btt <- unique(round(btt))
        if (length(intersect(btt, seq.int(p))) == 0L) {
            stop("Non-existent coefficients specified with 'btt'.")
        }
    }
    bntt <- setdiff(seq.int(p), btt)
    m <- length(btt)
    con <- list(tau2.init = NULL, tau2.min = 0, tau2.max = 50, 
        threshold = 10^-5, maxiter = 100, stepadj = 1, REML2 = FALSE, 
        verbose = FALSE)
    con[pmatch(names(control), names(con))] <- control
    if (verbose) 
        con$verbose <- verbose
    iter <- 0
    se.tau2 <- I2 <- H2 <- QE <- QEp <- NA
    s2w <- 1
    alpha <- (100 - level)/100
    robust <- FALSE
    Y <- as.matrix(yi)
    if (!is.numeric(tau2)) {
        if (method == "HS") {
            if (!allvipos) 
                stop("HS estimator cannot be used with non-positive sampling variances.")
            wi <- 1/vi
            W <- .diag(wi)
            stXWX <- .invcalc(X = X, W = W, k = k)
            P <- W - W %*% X %*% stXWX %*% crossprod(X, W)
            RSS <- crossprod(Y, P) %*% Y
            tau2 <- (RSS - k)/sum(wi)
            se.tau2 <- sqrt(1/sum(wi)^2 * (2 * (k - p) + 4 * 
                max(tau2, 0) * .tr(P) + 2 * max(tau2, 0)^2 * 
                .tr(P %*% P)))
        }
        if (is.element(method, c("HE", "ML", "REML", "EB"))) {
            stXX <- .invcalc(X = X, W = diag(k), k = k)
            P <- diag(k) - X %*% tcrossprod(stXX, X)
            RSS <- crossprod(Y, P) %*% Y
            trPV <- .tr(P %*% .diag(vi))
            tau2 <- (RSS - trPV)/(k - p)
            se.tau2 <- sqrt(1/(k - p)^2 * (2 * .tr(P %*% .diag(vi) %*% 
                P %*% .diag(vi)) + 4 * max(tau2, 0) * trPV + 
                2 * max(tau2, 0)^2 * (k - p)))
        }
        if (method == "DL") {
            if (!allvipos) 
                stop("DL estimator cannot be used with non-positive sampling variances.")
            wi <- 1/vi
            W <- .diag(wi)
            stXWX <- .invcalc(X = X, W = W, k = k)
            P <- W - W %*% X %*% stXWX %*% crossprod(X, W)
            RSS <- crossprod(Y, P) %*% Y
            trP <- .tr(P)
            tau2 <- (RSS - (k - p))/trP
            se.tau2 <- sqrt(1/trP^2 * (2 * (k - p) + 4 * max(tau2, 
                0) * trP + 2 * max(tau2, 0)^2 * .tr(P %*% P)))
        }
        if (method == "SJ") {
            if (is.null(con$tau2.init)) {
                tau2.0 <- var(yi) * (k - 1)/k
            }
            else {
                tau2.0 <- con$tau2.init
            }
            wi <- 1/(vi + tau2.0)
            W <- .diag(wi)
            stXWX <- .invcalc(X = X, W = W, k = k)
            P <- W - W %*% X %*% stXWX %*% crossprod(X, W)
            RSS <- crossprod(Y, P) %*% Y
            tau2 <- tau2.0 * RSS/(k - p)
            se.tau2 <- sqrt(tau2.0^2/(k - p)^2 * (2 * .tr(P %*% 
                .diag(vi) %*% P %*% .diag(vi)) + 4 * max(tau2, 
                0) * .tr(P %*% .diag(vi) %*% P) + 2 * max(tau2, 
                0)^2 * .tr(P %*% P)))
        }
        if (method == "DLit") {
            conv <- 1
            change <- con$threshold + 1
            if (is.null(con$tau2.init)) {
                tau2 <- 0
            }
            else {
                tau2 <- con$tau2.init
            }
            while (change > con$threshold) {
                if (con$verbose) 
                  cat("Iteration", iter, "\ttau^2 =", round(tau2, 
                    digits), "\n")
                iter <- iter + 1
                tau2.old <- tau2
                wi <- 1/(vi + tau2)
                W <- .diag(wi)
                stXWX <- .invcalc(X = X, W = W, k = k)
                P <- W - W %*% X %*% stXWX %*% crossprod(X, W)
                RSS <- crossprod(Y, P) %*% Y
                trP <- .tr(P)
                tau2 <- (RSS - (k - p))/trP
                tau2[tau2 < con$tau2.min] <- con$tau2.min
                change <- abs(tau2.old - tau2)
                if (iter > con$maxiter) {
                  conv <- 0
                  break
                }
            }
            if (conv == 0L) 
                stop("Algorithm did not converge.")
            se.tau2 <- sqrt(1/trP^2 * (2 * (k - p) + 4 * max(tau2, 
                0) * trP + 2 * max(tau2, 0)^2 * .tr(P %*% P)))
        }
        if (method == "SJit") {
            conv <- 1
            change <- con$threshold + 1
            if (is.null(con$tau2.init)) {
                tau2 <- var(yi) * (k - 1)/k
            }
            else {
                tau2 <- con$tau2.init
            }
            while (change > con$threshold) {
                if (con$verbose) 
                  cat("Iteration", iter, "\ttau^2 =", round(tau2, 
                    digits), "\n")
                iter <- iter + 1
                tau2.old <- tau2
                wi <- 1/(vi + tau2)
                W <- .diag(wi)
                stXWX <- .invcalc(X = X, W = W, k = k)
                P <- W - W %*% X %*% stXWX %*% crossprod(X, W)
                RSS <- crossprod(Y, P) %*% Y
                tau2 <- tau2 * RSS/(k - p)
                change <- abs(tau2.old - tau2)
                if (iter > con$maxiter) {
                  conv <- 0
                  break
                }
            }
            if (conv == 0L) 
                stop("Algorithm did not converge.")
        }
        if (method == "PM") {
            if (!allvipos) 
                stop("PM estimator cannot be used with non-positive sampling variances.")
            if (.QE.func(con$tau2.min, Y = Y, vi = vi, X = X, 
                k = k, objective = k - p) < 0) {
                tau2 <- con$tau2.min
            }
            else {
                tau2 <- try(uniroot(.QE.func, interval = c(con$tau2.min, 
                  con$tau2.max), tol = con$threshold, Y = Y, 
                  vi = vi, X = X, k = k, objective = k - p)$root, 
                  silent = TRUE)
                if (!is.numeric(tau2)) {
                  stop("Error in iterative search for tau2. Try increasing tau2.max or switch to another 'method'.")
                }
            }
        }
        if (is.element(method, c("ML", "REML", "EB"))) {
            conv <- 1
            change <- con$threshold + 1
            if (is.null(con$tau2.init)) {
                tau2 <- max(0, tau2)
            }
            else {
                tau2 <- con$tau2.init
            }
            while (change > con$threshold) {
                if (con$verbose) 
                  cat("Iteration", iter, "\ttau^2 =", round(tau2, 
                    digits), "\n")
                iter <- iter + 1
                tau2.old <- tau2
                wi <- 1/(vi + tau2)
                W <- .diag(wi)
                stXWX <- .invcalc(X = X, W = W, k = k)
                P <- W - W %*% X %*% stXWX %*% crossprod(X, W)
                if (method == "ML") {
                  PP <- P %*% P
                  adj <- (crossprod(Y, PP) %*% Y - sum(wi))/sum(wi^2)
                }
                if (method == "REML") {
                  PP <- P %*% P
                  adj <- (crossprod(Y, PP) %*% Y - .tr(P))/.tr(PP)
                }
                if (method == "EB") {
                  adj <- (crossprod(Y, P) %*% Y * k/(k - p) - 
                    k)/sum(wi)
                }
                adj <- adj * con$stepadj
                while (tau2 + adj < con$tau2.min) {
                  adj <- adj/2
                }
                tau2 <- tau2 + adj
                change <- abs(tau2.old - tau2)
                if (iter > con$maxiter) {
                  conv <- 0
                  break
                }
            }
            if (conv == 0L) 
                stop("Fisher scoring algorithm did not converge. Try increasing the number of iterations,\n  adjusting the threshold, or use a different estimator for tau^2.")
            if (method == "ML") 
                se.tau2 <- sqrt(2/sum(wi^2))
            if (method == "REML") 
                se.tau2 <- sqrt(2/.tr(PP))
            if (method == "EB") 
                se.tau2 <- sqrt((k/(k - p))^2/sum(wi)^2 * (2 * 
                  .tr(P %*% .diag(vi) %*% P %*% .diag(vi)) + 
                  4 * max(tau2, 0) * .tr(P %*% .diag(vi) %*% 
                    P) + 2 * max(tau2, 0)^2 * .tr(P %*% P)))
        }
        tau2 <- max(con$tau2.min, c(tau2))
        if (con$verbose && is.element(method, c("ML", "REML", 
            "EB"))) 
            cat("Fisher scoring algorithm converged after", iter, 
                "iterations.", "\n")
    }
    else {
        if (method == "HS") {
            wi <- 1/vi
            W <- .diag(wi)
            stXWX <- .invcalc(X = X, W = W, k = k)
            P <- W - W %*% X %*% stXWX %*% crossprod(X, W)
            se.tau2 <- sqrt(1/sum(wi)^2 * (2 * (k - p) + 4 * 
                max(tau2, 0) * .tr(P) + 2 * max(tau2, 0)^2 * 
                .tr(P %*% P)))
        }
        if (method == "HE") {
            stXX <- .invcalc(X = X, W = diag(k), k = k)
            P <- diag(k) - X %*% tcrossprod(stXX, X)
            RSS <- crossprod(Y, P) %*% Y
            trPV <- .tr(P %*% .diag(vi))
            se.tau2 <- sqrt(1/(k - p)^2 * (2 * .tr(P %*% .diag(vi) %*% 
                P %*% .diag(vi)) + 4 * max(tau2, 0) * trPV + 
                2 * max(tau2, 0)^2 * (k - p)))
        }
        if (method == "DL") {
            wi <- 1/vi
            W <- .diag(wi)
            stXWX <- .invcalc(X = X, W = W, k = k)
            P <- W - W %*% X %*% stXWX %*% crossprod(X, W)
            trP <- .tr(P)
            se.tau2 <- sqrt(1/trP^2 * (2 * (k - p) + 4 * max(tau2, 
                0) * trP + 2 * max(tau2, 0)^2 * .tr(P %*% P)))
        }
        if (method == "SJ") {
            if (is.null(con$tau2.init)) {
                tau2.0 <- var(yi) * (k - 1)/k
            }
            else {
                tau2.0 <- con$tau2.init
            }
            wi <- 1/(vi + tau2.0)
            W <- .diag(wi)
            stXWX <- .invcalc(X = X, W = W, k = k)
            P <- W - W %*% X %*% stXWX %*% crossprod(X, W)
            se.tau2 <- sqrt(tau2.0^2/(k - p)^2 * (2 * .tr(P %*% 
                .diag(vi) %*% P %*% .diag(vi)) + 4 * max(tau2, 
                0) * .tr(P %*% .diag(vi) %*% P) + 2 * max(tau2, 
                0)^2 * .tr(P %*% P)))
        }
        if (method == "ML") {
            wi <- 1/(vi + tau2)
            se.tau2 <- sqrt(2/sum(wi^2))
        }
        if (method == "REML") {
            wi <- 1/(vi + tau2)
            W <- .diag(wi)
            stXWX <- .invcalc(X = X, W = W, k = k)
            P <- W - W %*% X %*% stXWX %*% crossprod(X, W)
            se.tau2 <- sqrt(2/.tr(P %*% P))
        }
        if (method == "EB") {
            wi <- 1/(vi + tau2)
            W <- .diag(wi)
            stXWX <- .invcalc(X = X, W = W, k = k)
            P <- W - W %*% X %*% stXWX %*% crossprod(X, W)
            se.tau2 <- sqrt((k/(k - p))^2/sum(wi)^2 * (2 * .tr(P %*% 
                .diag(vi) %*% P %*% .diag(vi)) + 4 * tau2 * .tr(P %*% 
                .diag(vi) %*% P) + 2 * tau2^2 * .tr(P %*% P)))
        }
    }
    if (method == "FE") {
        tau2 <- 0
        if (!allvipos && weighted) 
            stop("Weighted estimation cannot be used with a fixed-effects\n  model when there are non-positive sampling variances.")
    }
    if (allvipos) {
        wi <- 1/vi
        W <- .diag(wi)
        stXWX <- .invcalc(X = X, W = W, k = k)
        P <- W - W %*% X %*% stXWX %*% crossprod(X, W)
        QE <- max(0, c(crossprod(Y, P) %*% Y))
        QEp <- ifelse(k - p >= 1, pchisq(QE, df = k - p, lower.tail = FALSE), 
            1)
        vi.avg <- (k - p)/.tr(P)
        I2 <- 100 * tau2/(vi.avg + tau2)
        H2 <- tau2/vi.avg + 1
    }
    wi <- 1/(vi + tau2)
    W <- .diag(wi)
    if (weighted) {
        stXWX <- .invcalc(X = X, W = W, k = k)
        b <- stXWX %*% crossprod(X, W) %*% Y
        vb <- stXWX
        P <- W - W %*% X %*% stXWX %*% crossprod(X, W)
        RSS.f <- crossprod(Y, P) %*% Y
        if (robust) {
            resid <- c(Y - X %*% b)
            vb <- vb %*% t(X) %*% W %*% diag(resid^2) %*% W %*% 
                X %*% vb
            vb <- vb * k/(k - p)
        }
        if (knha) {
            if (RSS.f <= .Machine$double.eps) {
                s2w <- 1
            }
            else {
                s2w <- c(RSS.f)/(k - p)
            }
            vb <- s2w * vb
            if (method == "FE") 
                warning("Knapp & Hartung (2003) method is not meant to be used in the context of fixed-effects models.")
        }
        if (robust) {
            QM <- c(t(b)[btt] %*% chol2inv(chol(vb[btt, btt])) %*% 
                b[btt])
        }
        else {
            if (length(bntt) == 0L) {
                QM <- c(sum(wi * yi^2) - RSS.f)/s2w
            }
            else {
                Xr <- X[, bntt, drop = FALSE]
                stXWX <- .invcalc(X = Xr, W = W, k = k)
                P <- W - W %*% Xr %*% stXWX %*% crossprod(Xr, 
                  W)
                RSS.r <- crossprod(Y, P) %*% Y
                QM <- c(RSS.r - RSS.f)/s2w
            }
        }
    }
    else {
        stXX <- .invcalc(X = X, W = diag(k), k = k)
        b <- stXX %*% crossprod(X, Y)
        vb <- tcrossprod(stXX, X) %*% .diag(vi + tau2) %*% X %*% 
            stXX
        P <- W - W %*% X %*% tcrossprod(stXX, X) - X %*% stXX %*% 
            crossprod(X, W) + X %*% stXX %*% crossprod(X, W) %*% 
            X %*% tcrossprod(stXX, X)
        RSS.f <- crossprod(Y, P) %*% Y
        if (robust) {
            resid <- c(Y - X %*% b)
            vb <- stXX %*% t(X) %*% diag(resid^2) %*% X %*% stXX
            vb <- vb * k/(k - p)
        }
        if (knha) {
            stXWX <- .invcalc(X = X, W = W, k = k)
            P <- W - W %*% X %*% stXWX %*% crossprod(X, W)
            RSS.k <- c(crossprod(Y, P) %*% Y)
            if (RSS.k <= .Machine$double.eps) {
                s2w <- 1
            }
            else {
                s2w <- c(RSS.k)/(k - p)
            }
            vb <- s2w * vb
            if (method == "FE") 
                warning("Knapp & Hartung (2003) method is not meant to be used in the context of fixed-effects models.")
        }
        QM <- c(t(b)[btt] %*% chol2inv(chol(vb[btt, btt])) %*% 
            b[btt])
    }
    se <- sqrt(diag(vb))
    names(se) <- NULL
    zval <- c(b/se)
    if (knha || robust) {
        QM <- QM/m
        QMp <- pf(QM, df1 = m, df2 = k - p, lower.tail = FALSE)
        pval <- 2 * pt(abs(zval), df = k - p, lower.tail = FALSE)
        crit <- qt(alpha/2, df = k - p, lower.tail = FALSE)
    }
    else {
        QMp <- pchisq(QM, df = m, lower.tail = FALSE)
        pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)
        crit <- qnorm(alpha/2, lower.tail = FALSE)
    }
    ci.lb <- c(b - crit * se)
    ci.ub <- c(b + crit * se)
    parms <- p + ifelse(method == "FE", 0, 1)
    ll.ML <- -1/2 * (k) * log(2 * base::pi) - 1/2 * sum(log(vi + 
        tau2)) - 1/2 * RSS.f
    ll.REML <- -1/2 * (k - p) * log(2 * base::pi) + ifelse(con$REML2, 
        0, 1/2 * determinant(crossprod(X, X), logarithm = TRUE)$modulus) - 
        1/2 * sum(log(vi + tau2)) - 1/2 * determinant(crossprod(X, 
        W) %*% X, logarithm = TRUE)$modulus - 1/2 * RSS.f
    if (k - p > 0L) {
        dev.ML <- -2 * (ll.ML - sum(dnorm(yi, mean = yi, sd = sqrt(vi), 
            log = TRUE)))
    }
    else {
        dev.ML <- 0
    }
    AIC.ML <- -2 * ll.ML + 2 * parms
    BIC.ML <- -2 * ll.ML + parms * log(k)
    dev.REML <- -2 * (ll.REML - 0)
    AIC.REML <- -2 * ll.REML + 2 * parms
    BIC.REML <- -2 * ll.REML + parms * log(k - p)
    fit.stats <- matrix(c(ll.ML, dev.ML, AIC.ML, BIC.ML, ll.REML, 
        dev.REML, AIC.REML, BIC.REML), ncol = 2, byrow = FALSE)
    dimnames(fit.stats) <- list(c("ll", "dev", "AIC", "BIC"), 
        c("ML", "REML"))
    fit.stats <- data.frame(fit.stats)
    p.eff <- p
    k.eff <- k
    res <- list(b, se, zval, pval, ci.lb, ci.ub, vb, tau2, se.tau2, 
        k, k.f, k.eff, p, p.eff, parms, m, QE, QEp, QM, QMp, 
        I2, H2, int.only, int.incl, allvipos, yi, vi, X, yi.f, 
        vi.f, X.f, ai.f, bi.f, ci.f, di.f, x1i.f, x2i.f, t1i.f, 
        t2i.f, ni, ni.f, ids, not.na, slab, slab.null, measure, 
        method, weighted, knha, robust, s2w, btt, intercept, 
        digits, level, con, add, to, drop00, fit.stats)
    names(res) <- c("b", "se", "zval", "pval", "ci.lb", "ci.ub", 
        "vb", "tau2", "se.tau2", "k", "k.f", "k.eff", "p", "p.eff", 
        "parms", "m", "QE", "QEp", "QM", "QMp", "I2", "H2", "int.only", 
        "int.incl", "allvipos", "yi", "vi", "X", "yi.f", "vi.f", 
        "X.f", "ai.f", "bi.f", "ci.f", "di.f", "x1i.f", "x2i.f", 
        "t1i.f", "t2i.f", "ni", "ni.f", "ids", "not.na", "slab", 
        "slab.null", "measure", "method", "weighted", "knha", 
        "robust", "s2w", "btt", "intercept", "digits", "level", 
        "control", "add", "to", "drop00", "fit.stats")
    class(res) <- c("rma.uni", "rma")
    return(res)
}
rma.glmm <-
function (ai, bi, ci, di, n1i, n2i, x1i, x2i, t1i, t2i, xi, mi, 
    ti, ni, mods, measure, intercept = TRUE, data, slab, subset, 
    add = 1/2, to = "only0", drop00 = TRUE, vtype = "LS", model = "UM.FS", 
    method = "ML", tdist = FALSE, level = 95, digits = 4, btt, 
    nAGQ = 100, verbose = FALSE, control) 
{
    if (missing(measure)) 
        stop("Need to specify 'measure' argument.")
    if (!is.element(measure, c("OR", "IRR", "PLO", "IRLN"))) 
        stop("Unknown 'measure' specified.")
    if (!is.element(method, c("FE", "ML"))) 
        stop("Unknown 'method' specified.")
    if (length(add) > 1) 
        add <- add[1]
    if (length(to) > 1) 
        to <- to[1]
    if (!is.element(model, c("UM.FS", "UM.RS", "CM.EL", "CM.AL"))) 
        stop("Unknown 'model' specified.")
    if (model == "CM.AL" & measure == "IRR") 
        model <- "CM.EL"
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    if (missing(btt)) 
        btt <- NULL
    if (missing(control)) 
        control <- list()
    knha <- tdist
    if (missing(data)) 
        data <- NULL
    if (is.null(data)) {
        data <- sys.frame(sys.parent())
    }
    else {
        if (!is.data.frame(data)) {
            data <- data.frame(data)
        }
    }
    mf <- match.call()
    mf.slab <- mf[[match("slab", names(mf))]]
    mf.subset <- mf[[match("subset", names(mf))]]
    mf.mods <- mf[[match("mods", names(mf))]]
    slab <- eval(mf.slab, data, enclos = sys.frame(sys.parent()))
    subset <- eval(mf.subset, data, enclos = sys.frame(sys.parent()))
    mods <- eval(mf.mods, data, enclos = sys.frame(sys.parent()))
    ai <- bi <- ci <- di <- x1i <- x2i <- t1i <- t2i <- xi <- mi <- ti <- NA
    if (is.element(measure, c("OR"))) {
        mf.ai <- mf[[match("ai", names(mf))]]
        mf.bi <- mf[[match("bi", names(mf))]]
        mf.ci <- mf[[match("ci", names(mf))]]
        mf.di <- mf[[match("di", names(mf))]]
        mf.n1i <- mf[[match("n1i", names(mf))]]
        mf.n2i <- mf[[match("n2i", names(mf))]]
        ai <- eval(mf.ai, data, enclos = sys.frame(sys.parent()))
        bi <- eval(mf.bi, data, enclos = sys.frame(sys.parent()))
        ci <- eval(mf.ci, data, enclos = sys.frame(sys.parent()))
        di <- eval(mf.di, data, enclos = sys.frame(sys.parent()))
        n1i <- eval(mf.n1i, data, enclos = sys.frame(sys.parent()))
        n2i <- eval(mf.n2i, data, enclos = sys.frame(sys.parent()))
        if (is.null(bi)) 
            bi <- n1i - ai
        if (is.null(di)) 
            di <- n2i - ci
        dat <- escalc(measure, ai = ai, bi = bi, ci = ci, di = di, 
            add = add, to = to, drop00 = drop00, vtype = vtype)
    }
    if (is.element(measure, c("IRR"))) {
        mf.x1i <- mf[[match("x1i", names(mf))]]
        mf.x2i <- mf[[match("x2i", names(mf))]]
        mf.t1i <- mf[[match("t1i", names(mf))]]
        mf.t2i <- mf[[match("t2i", names(mf))]]
        x1i <- eval(mf.x1i, data, enclos = sys.frame(sys.parent()))
        x2i <- eval(mf.x2i, data, enclos = sys.frame(sys.parent()))
        t1i <- eval(mf.t1i, data, enclos = sys.frame(sys.parent()))
        t2i <- eval(mf.t2i, data, enclos = sys.frame(sys.parent()))
        dat <- escalc(measure, x1i = x1i, x2i = x2i, t1i = t1i, 
            t2i = t2i, add = add, to = to, drop00 = drop00, vtype = vtype)
    }
    if (is.element(measure, c("PLO"))) {
        mf.xi <- mf[[match("xi", names(mf))]]
        mf.mi <- mf[[match("mi", names(mf))]]
        mf.ni <- mf[[match("ni", names(mf))]]
        xi <- eval(mf.xi, data, enclos = sys.frame(sys.parent()))
        mi <- eval(mf.mi, data, enclos = sys.frame(sys.parent()))
        ni <- eval(mf.ni, data, enclos = sys.frame(sys.parent()))
        if (is.null(mi)) 
            mi <- ni - xi
        dat <- escalc(measure, xi = xi, mi = mi, add = add, to = to, 
            vtype = vtype)
    }
    if (is.element(measure, c("IRLN"))) {
        mf.xi <- mf[[match("xi", names(mf))]]
        mf.ti <- mf[[match("ti", names(mf))]]
        xi <- eval(mf.xi, data, enclos = sys.frame(sys.parent()))
        ti <- eval(mf.ti, data, enclos = sys.frame(sys.parent()))
        dat <- escalc(measure, xi = xi, ti = ti, add = add, to = to, 
            vtype = vtype)
    }
    yi <- dat$yi
    vi <- dat$vi
    ni <- attr(yi, "ni")
    is.formula <- ifelse(class(mods) == "formula", TRUE, FALSE)
    if (is.formula) {
        options(na.action = "na.pass")
        mods <- model.matrix(mods, data = data)
        attr(mods, "assign") <- NULL
        options(na.action = na.act)
        intercept <- FALSE
    }
    if (is.vector(mods)) 
        mods <- cbind(mods)
    if (is.data.frame(mods)) 
        mods <- as.matrix(mods)
    if (!is.null(mods) && (nrow(mods) != length(yi))) 
        stop("Number of rows of the design matrix does not match length of yi.")
    k <- length(yi)
    ids <- seq.int(k)
    if (is.null(slab)) {
        if (!is.null(attr(yi, "slab"))) {
            slab.null <- FALSE
            slab <- attr(yi, "slab")
        }
        else {
            slab.null <- TRUE
            slab <- seq.int(k)
        }
    }
    else {
        if (any(duplicated(slab))) 
            stop("Study labels must be unique.")
        if (length(slab) != length(yi)) 
            stop("Study labels not of same length as data.")
        slab.null <- FALSE
    }
    if (!is.null(subset)) {
        yi <- yi[subset]
        vi <- vi[subset]
        ni <- ni[subset]
        ai <- ai[subset]
        bi <- bi[subset]
        ci <- ci[subset]
        di <- di[subset]
        x1i <- x1i[subset]
        x2i <- x2i[subset]
        t1i <- t1i[subset]
        t2i <- t2i[subset]
        xi <- xi[subset]
        mi <- mi[subset]
        ti <- ti[subset]
        mods <- mods[subset, , drop = FALSE]
        slab <- slab[subset]
        ids <- ids[subset]
        k <- length(yi)
    }
    if (measure == "OR") {
        if (drop00) {
            id00 <- c(ai == 0L & ci == 0L) | c(bi == 0L & di == 
                0L)
            id00[is.na(id00)] <- FALSE
            ai[id00] <- NA
            bi[id00] <- NA
            ci[id00] <- NA
            di[id00] <- NA
        }
    }
    if (measure == "IRR") {
        if (drop00) {
            id00 <- c(x1i == 0L & x2i == 0L)
            id00[is.na(id00)] <- FALSE
            x1i[id00] <- NA
            x2i[id00] <- NA
        }
    }
    ai.f <- ai
    bi.f <- bi
    ci.f <- ci
    di.f <- di
    x1i.f <- x1i
    x2i.f <- x2i
    t1i.f <- t1i
    t2i.f <- t2i
    xi.f <- xi
    mi.f <- mi
    ti.f <- ti
    yi.f <- yi
    vi.f <- vi
    ni.f <- ni
    mods.f <- mods
    k.f <- k
    if (is.element(measure, c("OR"))) {
        aibicidimods.na <- is.na(cbind(ai, bi, ci, di, mods))
        if (any(aibicidimods.na)) {
            not.na <- apply(aibicidimods.na, MARGIN = 1, sum) == 
                0L
            if (na.act == "na.omit" || na.act == "na.exclude" || 
                na.act == "na.pass") {
                ai <- ai[not.na]
                bi <- bi[not.na]
                ci <- ci[not.na]
                di <- di[not.na]
                mods <- mods[not.na, , drop = FALSE]
                k <- length(ai)
                warning("Studies with NAs omitted from model fitting.")
            }
            if (na.act == "na.fail") 
                stop("Missing values in studies.")
        }
        else {
            not.na <- rep(TRUE, k)
        }
    }
    if (is.element(measure, c("IRR"))) {
        x1ix2it1it2imods.na <- is.na(cbind(x1i, x2i, t1i, t2i, 
            mods))
        if (any(x1ix2it1it2imods.na)) {
            not.na <- apply(x1ix2it1it2imods.na, MARGIN = 1, 
                sum) == 0L
            if (na.act == "na.omit" || na.act == "na.exclude" || 
                na.act == "na.pass") {
                x1i <- x1i[not.na]
                x2i <- x2i[not.na]
                t1i <- t1i[not.na]
                t2i <- t2i[not.na]
                mods <- mods[not.na, , drop = FALSE]
                k <- length(ai)
                warning("Studies with NAs omitted from model fitting.")
            }
            if (na.act == "na.fail") 
                stop("Missing values in studies.")
        }
        else {
            not.na <- rep(TRUE, k)
        }
    }
    if (is.element(measure, c("PLO"))) {
        ximimods.na <- is.na(cbind(xi, mi, mods))
        if (any(ximimods.na)) {
            not.na <- apply(ximimods.na, MARGIN = 1, sum) == 
                0L
            if (na.act == "na.omit" || na.act == "na.exclude" || 
                na.act == "na.pass") {
                xi <- xi[not.na]
                mi <- mi[not.na]
                mods <- mods[not.na, , drop = FALSE]
                k <- length(ai)
                warning("Studies with NAs omitted from model fitting.")
            }
            if (na.act == "na.fail") 
                stop("Missing values in studies.")
        }
        else {
            not.na <- rep(TRUE, k)
        }
    }
    if (is.element(measure, c("IRLN"))) {
        xitimods.na <- is.na(cbind(xi, ti, mods))
        if (any(xitimods.na)) {
            not.na <- apply(xitimods.na, MARGIN = 1, sum) == 
                0L
            if (na.act == "na.omit" || na.act == "na.exclude" || 
                na.act == "na.pass") {
                xi <- xi[not.na]
                ti <- ti[not.na]
                mods <- mods[not.na, , drop = FALSE]
                k <- length(ai)
                warning("Studies with NAs omitted from model fitting.")
            }
            if (na.act == "na.fail") 
                stop("Missing values in studies.")
        }
        else {
            not.na <- rep(TRUE, k)
        }
    }
    if (k < 1) 
        stop("Processing terminated since k = 0.")
    if (k < 2) 
        stop("Need at least k=2 studies to fit models.")
    yivi.na <- is.na(cbind(yi, vi, mods.f))
    if (any(yivi.na)) {
        not.na.yivi <- apply(yivi.na, MARGIN = 1, sum) == 0L
        if (na.act == "na.omit" || na.act == "na.exclude" || 
            na.act == "na.pass") {
            yi <- yi[not.na.yivi]
            ni <- ni[not.na.yivi]
            vi <- vi[not.na.yivi]
            warning("Some yi/vi values are NA.")
        }
        if (na.act == "na.fail") 
            stop("Missing yi/vi values.")
    }
    else {
        not.na.yivi <- rep(TRUE, k)
    }
    k.yi <- length(yi)
    vi.I2 <- vi.f
    mods.I2 <- mods.f
    YVM.na <- is.na(cbind(yi.f, vi.f, mods.f))
    if (any(YVM.na)) {
        not.na.I2 <- apply(YVM.na, MARGIN = 1, sum) == 0L
        if (na.act == "na.omit" || na.act == "na.exclude" || 
            na.act == "na.pass") {
            vi.I2 <- vi.f[not.na.I2]
            mods.I2 <- mods.f[not.na.I2, , drop = FALSE]
        }
        if (na.act == "na.fail") 
            stop("Missing values in data.")
    }
    k.I2 <- length(vi.I2)
    if (is.null(mods) && !intercept) {
        warning("Must either include an intercept (intercept=TRUE) and/or moderators in model.\n  Coerced intercept into the model.")
        intercept <- TRUE
    }
    if (intercept) {
        X <- cbind(intrcpt = rep(1, k), mods)
        X.f <- cbind(intrcpt = rep(1, k.f), mods.f)
        X.I2 <- cbind(intrcpt = rep(1, k.I2), mods.I2)
    }
    else {
        X <- mods
        X.f <- mods.f
        X.I2 <- mods.I2
    }
    is.int <- apply(X == 1L, 2, sum) == k
    if (any(is.int)) {
        int.incl <- TRUE
        int.indx <- which(is.int, arr.ind = TRUE)
        X <- cbind(intrcpt = 1, X[, -int.indx, drop = FALSE])
        X.f <- cbind(intrcpt = 1, X.f[, -int.indx, drop = FALSE])
        X.I2 <- cbind(intrcpt = 1, X.I2[, -int.indx, drop = FALSE])
        if (is.formula) 
            intercept <- TRUE
    }
    else {
        int.incl <- FALSE
    }
    p <- NCOL(X)
    if ((p == 1L) && (all(sapply(X, identical, 1)))) {
        int.only <- TRUE
    }
    else {
        int.only <- FALSE
    }
    evstXX <- eigen(t(X) %*% X, symmetric = TRUE)$values
    if (any(evstXX <= .Machine$double.eps)) 
        stop("Design matrix not of full rank. Cannot fit model.")
    if (method == "FE" && p > k) 
        stop("Number of parameters to be estimated is larger than the number of observations.")
    if (method != "FE" && (p + 1) > k) 
        stop("Number of parameters to be estimated is larger than the number of observations.")
    if (is.null(btt)) {
        if (p > 1) {
            if (int.incl) {
                btt <- seq.int(2, p)
            }
            else {
                btt <- seq.int(p)
            }
        }
        else {
            btt <- 1
        }
    }
    else {
        btt <- btt[(btt >= 1) & (btt <= p)]
        btt <- unique(round(btt))
        if (length(intersect(btt, seq.int(p))) == 0L) {
            stop("Non-existent coefficients specified with 'btt'.")
        }
    }
    m <- length(btt)
    con <- list(verbose = FALSE, epsilon = 1e-08, maxit = 25, 
        maxIter = 300, maxFN = 900, routine = "optim", method = "BFGS", 
        reltol = 1e-08, REPORT = 1, rel.tol = 1e-08, subdivisions = 100, 
        lower = -Inf, upper = Inf, dnchgcalc = "dFNCHypergeo", 
        dnchgprec = 1e-10)
    con[pmatch(names(control), names(con))] <- control
    control.glm <- list(epsilon = con$epsilon, maxit = con$maxit, 
        trace = con$verbose)
    control.lmer <- list(msVerbose = con$verbose, maxIter = con$maxIter, 
        maxFN = con$maxFN)
    control.optim <- list(trace = ifelse(con$verbose, 1, 0), 
        reltol = con$reltol, REPORT = con$REPORT)
    control.int <- list(rel.tol = con$rel.tol, subdivisions = con$subdivisions, 
        lower = con$lower, upper = con$upper)
    if (verbose) 
        con$verbose <- verbose
    if (!is.element(con$routine, c("optim", "clogit", "clogistic"))) 
        stop("Unknown routine specified.")
    if (con$dnchgcalc != "dnoncenhypergeom" && con$dnchgcalc != 
        "dFNCHypergeo") 
        stop("Unknown dnchgcalc method specified.")
    if (is.element(measure, c("OR", "IRR"))) {
        if ((model == "UM.FS" && method == "ML") || (model == 
            "UM.RS") || (model == "CM.AL" && method == "ML") || 
            (model == "CM.EL" && method == "ML")) {
            if (!require(lme4)) 
                stop("Please install the 'lme4' package to fit this model.")
        }
    }
    if (measure == "OR" && model == "CM.EL") {
        if (con$routine == "optim") {
            if (con$dnchgcalc == "dFNCHypergeo") {
                if (!require(BiasedUrn)) 
                  stop("Please install the 'BiasedUrn' package to fit this model.")
            }
        }
        if (con$routine == "clogit") {
            if (!require(survival)) 
                stop("Please install the 'survival' package to fit this model.")
        }
        if (con$routine == "clogistic") {
            if (!require(Epi)) 
                stop("Please install the 'Epi' package to fit this model.")
        }
    }
    if (is.element(measure, c("PLO", "IRLN"))) {
        if (method == "ML") {
            if (!require(lme4)) 
                stop("Please install the 'lme4' package to fit this model.")
        }
    }
    silent <- !con$verbose
    se.tau2 <- I2 <- H2 <- QE <- QEp <- NA
    alpha <- (100 - level)/100
    robust <- FALSE
    if (is.element(measure, c("OR", "IRR"))) {
        if (is.element(model, c("UM.FS", "UM.RS"))) {
            if (measure == "OR") {
                dat.grp <- cbind(xi = c(rbind(ai, ci)), mi = c(rbind(bi, 
                  di)))
                dat.fam <- binomial
                dat.off <- NULL
            }
            if (measure == "IRR") {
                dat.grp <- cbind(xi = c(rbind(x1i, x2i)))
                dat.fam <- poisson
                dat.off <- log(c(rbind(t1i, t2i)))
            }
            group1 <- rep(c(1, 0), times = k)
            group2 <- rep(c(0, 1), times = k)
            group12 <- rep(c(1/2, -1/2), times = k)
            study <- factor(rep(seq.int(k), each = 2))
            intrcpt <- cbind(rep(1, 2 * k))
            X.fit <- X[rep(seq(k), each = 2), , drop = FALSE]
            X.fit <- cbind(group1 * X.fit[, , drop = FALSE])
            row.names(X.fit) <- seq.int(2 * k)
            if (model == "UM.FS") {
                if (con$verbose) 
                  message("Fitting FE model.")
                res.FE <- try(glm(dat.grp ~ -1 + X.fit + study, 
                  offset = dat.off, family = dat.fam, control = control.glm), 
                  silent = silent)
                if (inherits(res.FE, "try-error")) 
                  stop("Cannot fit FE model.")
                if (con$verbose) 
                  message("Fitting saturated model.")
                X.QE <- model.matrix(~-1 + X.fit + study + study:group1)
                res.QE <- try(glm(dat.grp ~ -1 + X.QE, offset = dat.off, 
                  family = dat.fam, control = control.glm), silent = silent)
                if (inherits(res.QE, "try-error")) 
                  stop("Cannot fit saturated model.")
                ll.FE <- c(logLik(res.FE))
                ll.QE <- c(logLik(res.QE))
                b2.QE <- cbind(na.omit(coef(res.QE)[-seq.int(k + 
                  p)]))
                vb2.QE <- vcov(res.QE)[-seq.int(k + p), -seq.int(k + 
                  p), drop = FALSE]
                if (method == "ML") {
                  if (con$verbose) 
                    message("Fitting ML model.")
                  res.ML <- try(lmer(dat.grp ~ -1 + X.fit + study + 
                    (group12 - 1 | study), offset = dat.off, 
                    family = dat.fam, nAGQ = nAGQ, verbose = con$verbose, 
                    control = control.lmer), silent = silent)
                  if (inherits(res.ML, "try-error")) 
                    stop("Cannot fit ML model.")
                  ll.ML <- ll.QE - 1/2 * lme4::deviance(res.ML)
                }
                if (method == "FE") {
                  b <- cbind(coef(res.FE)[seq.int(p)])
                  vb <- vcov(res.FE)[seq.int(p), seq.int(p), 
                    drop = FALSE]
                  tau2 <- 0
                  sigma2 <- NA
                  parms <- p + k
                  p.eff <- p + k
                  k.eff <- 2 * k
                }
                if (method == "ML") {
                  b <- cbind(fixef(res.ML)[seq.int(p)])
                  vb <- as.matrix(lme4::vcov(res.ML))[seq.int(p), 
                    seq.int(p), drop = FALSE]
                  tau2 <- VarCorr(res.ML)[[1]][1]
                  sigma2 <- NA
                  parms <- p + k + 1
                  p.eff <- p + k
                  k.eff <- 2 * k
                }
            }
            if (model == "UM.RS") {
                if (con$verbose) 
                  message("Fitting FE model.")
                res.FE <- try(lmer(dat.grp ~ -1 + X.fit + intrcpt + 
                  (1 | study), offset = dat.off, family = dat.fam, 
                  nAGQ = nAGQ, verbose = con$verbose, control = control.lmer), 
                  silent = silent)
                if (inherits(res.FE, "try-error")) 
                  stop("Cannot fit FE model.")
                if (con$verbose) 
                  message("Fitting saturated model.")
                X.QE <- model.matrix(~-1 + X.fit + intrcpt + 
                  study:group1)
                res.QE <- try(glm(dat.grp ~ -1 + X.QE, offset = dat.off, 
                  family = dat.fam, control = control.glm), silent = silent)
                X.QE <- X.QE[, !is.na(coef(res.QE)), drop = FALSE]
                res.QE <- try(lmer(dat.grp ~ -1 + X.QE + (1 | 
                  study), offset = dat.off, family = dat.fam, 
                  start = c(sqrt(VarCorr(res.FE)[[1]][1])), nAGQ = nAGQ, 
                  verbose = con$verbose, control = control.lmer), 
                  silent = silent)
                if (inherits(res.QE, "try-error")) 
                  stop("Cannot fit saturated model.")
                ll.FE <- c(lme4::logLik(res.FE))
                ll.QE <- c(lme4::logLik(res.QE))
                b2.QE <- cbind(fixef(res.QE)[-seq.int(p + 1)])
                vb2.QE <- as.matrix(lme4::vcov(res.QE))[-seq.int(p + 
                  1), -seq.int(p + 1), drop = FALSE]
                if (method == "ML") {
                  if (con$verbose) 
                    message("Fitting ML model.")
                  res.ML <- try(lmer(dat.grp ~ -1 + X.fit + intrcpt + 
                    (1 | study) + (group12 - 1 | study), offset = dat.off, 
                    family = dat.fam, nAGQ = nAGQ, verbose = con$verbose, 
                    control = control.lmer), silent = silent)
                  if (inherits(res.ML, "try-error")) 
                    stop("Cannot fit ML model.")
                  ll.ML <- c(lme4::logLik(res.ML))
                }
                if (method == "FE") {
                  b <- cbind(fixef(res.FE)[seq.int(p)])
                  vb <- as.matrix(lme4::vcov(res.FE))[seq.int(p), 
                    seq.int(p), drop = FALSE]
                  tau2 <- 0
                  sigma2 <- VarCorr(res.FE)[[1]][1]
                  parms <- p + 1 + 1
                  p.eff <- p + 1
                  k.eff <- 2 * k
                }
                if (method == "ML") {
                  b <- cbind(fixef(res.ML)[seq.int(p)])
                  vb <- as.matrix(lme4::vcov(res.ML))[seq.int(p), 
                    seq.int(p), drop = FALSE]
                  tau2 <- VarCorr(res.ML)[[2]][1]
                  sigma2 <- VarCorr(res.ML)[[1]][1]
                  parms <- p + 1 + 2
                  p.eff <- p + 1
                  k.eff <- 2 * k
                }
            }
        }
        if ((measure == "IRR" && model == "CM.EL") || (measure == 
            "OR" && model == "CM.AL") || (measure == "OR" && 
            model == "CM.EL")) {
            if (measure == "OR") {
                dat.grp <- cbind(xi = ai, mi = ci)
                dat.off <- log((ai + bi)/(ci + di))
            }
            if (measure == "IRR") {
                dat.grp <- cbind(xi = x1i, mi = x2i)
                dat.off <- log(t1i/t2i)
            }
            study <- factor(seq.int(k))
            X.fit <- X
            if (con$verbose) 
                message("Fitting FE model.")
            res.FE <- try(glm(dat.grp ~ -1 + X.fit, offset = dat.off, 
                family = binomial, control = control.glm), silent = silent)
            if (inherits(res.FE, "try-error")) 
                stop("Cannot fit FE model.")
            if (con$verbose) 
                message("Fitting saturated model.")
            X.QE <- model.matrix(~-1 + X.fit + study)
            res.QE <- try(glm(dat.grp ~ -1 + X.QE, offset = dat.off, 
                family = binomial, control = control.glm), silent = silent)
            if (inherits(res.QE, "try-error")) 
                stop("Cannot fit saturated model.")
            ll.FE <- c(logLik(res.FE))
            ll.QE <- c(logLik(res.QE))
            b2.QE <- cbind(na.omit(coef(res.QE)[-seq.int(p)]))
            vb2.QE <- vcov(res.QE)[-seq.int(p), -seq.int(p), 
                drop = FALSE]
            if (method == "ML") {
                if (con$verbose) 
                  message("Fitting ML model.")
                if (con$verbose) {
                  res.ML <- try(lmer(dat.grp ~ -1 + X.fit + (1 | 
                    study), offset = dat.off, family = binomial, 
                    nAGQ = nAGQ, verbose = con$verbose, control = control.lmer), 
                    silent = silent)
                }
                else {
                  res.ML <- suppressMessages(try(lmer(dat.grp ~ 
                    -1 + X.fit + (1 | study), offset = dat.off, 
                    family = binomial, nAGQ = nAGQ, verbose = con$verbose, 
                    control = control.lmer), silent = silent))
                }
                if (inherits(res.ML, "try-error")) 
                  stop("Cannot fit ML model.")
                ll.ML <- ll.QE - 1/2 * lme4::deviance(res.ML)
            }
            if (method == "FE") {
                b <- cbind(coef(res.FE)[seq.int(p)])
                vb <- vcov(res.FE)[seq.int(p), seq.int(p), drop = FALSE]
                tau2 <- 0
                sigma2 <- NA
                parms <- p
                p.eff <- p
                k.eff <- k
            }
            if (method == "ML") {
                b <- cbind(fixef(res.ML)[seq.int(p)])
                vb <- as.matrix(lme4::vcov(res.ML))[seq.int(p), 
                  seq.int(p), drop = FALSE]
                tau2 <- VarCorr(res.ML)[[1]][1]
                sigma2 <- NA
                parms <- p + 1
                p.eff <- p
                k.eff <- k
            }
        }
        if (measure == "OR" && model == "CM.EL") {
            if (con$verbose) 
                message("Fitting FE model.")
            if (con$routine == "optim") {
                res.FE <- try(optim(par = c(coef(res.FE)[seq.int(p)], 
                  0), .dnchg, method = con$method, hessian = TRUE, 
                  ai = ai, bi = bi, ci = ci, di = di, X.fit = X.fit, 
                  random = FALSE, verbose = con$verbose, digits = digits, 
                  dnchgcalc = con$dnchgcalc, dnchgprec = con$dnchgprec, 
                  control = control.optim), silent = silent)
                if (inherits(res.FE, "try-error") || res.FE$convergence != 
                  0) 
                  stop("Cannot fit FE model.")
                if (con$verbose) 
                  message("Fitting saturated model.")
                X.QE <- X.QE[, !is.na(coef(res.QE)), drop = FALSE]
                res.QE <- try(optim(par = c(c(na.omit(coef(res.QE))), 
                  0), .dnchg, method = con$method, hessian = TRUE, 
                  ai = ai, bi = bi, ci = ci, di = di, X.fit = X.QE, 
                  random = FALSE, verbose = con$verbose, digits = digits, 
                  dnchgcalc = con$dnchgcalc, dnchgprec = con$dnchgprec, 
                  control = control.optim), silent = silent)
                if (inherits(res.QE, "try-error") || res.QE$convergence != 
                  0) 
                  stop("Cannot fit saturated model.")
                ll.FE <- -1 * res.FE$value
                ll.QE <- -1 * res.QE$value
                b2.QE <- res.QE$par
                hessian <- res.QE$hessian
                p.QE <- length(b2.QE)
                b2.QE <- b2.QE[-p.QE]
                hessian <- hessian[-p.QE, -p.QE, drop = FALSE]
                p.QE <- length(b2.QE)
                is.0 <- apply(hessian == 0L, 2, sum) == p.QE
                b2.QE <- b2.QE[!is.0]
                hessian <- hessian[!is.0, !is.0, drop = FALSE]
                b2.QE <- cbind(b2.QE[-seq.int(p)])
                h.A <- hessian[seq.int(p), seq.int(p), drop = FALSE]
                h.B <- hessian[seq.int(p), -seq.int(p), drop = FALSE]
                h.C <- hessian[-seq.int(p), seq.int(p), drop = FALSE]
                h.D <- hessian[-seq.int(p), -seq.int(p), drop = FALSE]
                chol.h.A <- try(chol(h.A), silent = silent)
                if (class(chol.h.A) == "try-error") 
                  stop("Cannot invert Hessian for saturated model.")
                Ivb2.QE <- h.D - h.C %*% chol2inv(chol.h.A) %*% 
                  h.B
                QE.Wld <- c(t(b2.QE) %*% Ivb2.QE %*% b2.QE)
            }
            if (con$routine == "clogit" || con$routine == "clogistic") {
                event <- unlist(sapply(seq.int(k), function(i) c(rep(1, 
                  ai[i]), rep(0, bi[i]), rep(1, ci[i]), rep(0, 
                  di[i]))))
                group1 <- unlist(sapply(seq.int(k), function(i) c(rep(1, 
                  ai[i]), rep(1, bi[i]), rep(0, ci[i]), rep(0, 
                  di[i]))))
                study.l <- factor(rep(seq.int(k), times = ni))
                X.fit.l <- X[rep(seq.int(k), times = ni), ]
                X.fit.l <- cbind(group1 * X.fit.l)
                if (con$routine == "clogit") 
                  res.FE <- try(clogit(event ~ X.fit.l + strata(study.l), 
                    iter.max = control.glm$maxit), silent = silent)
                if (con$routine == "clogistic") 
                  res.FE <- try(clogistic(event ~ X.fit.l, strata = study.l, 
                    iter.max = control.glm$maxit, model = FALSE, 
                    x = FALSE, y = FALSE), silent = silent)
                if (inherits(res.FE, "try-error")) 
                  stop("Cannot fit FE model.")
                if (con$verbose) 
                  message("Fitting saturated model.")
                X.QE.l <- model.matrix(~-1 + X.fit.l + study.l:group1)
                X.QE.l <- X.QE.l[, !is.na(coef(res.QE)), drop = FALSE]
                X.QE <- X.QE[, !is.na(coef(res.QE)), drop = FALSE]
                if (con$routine == "clogit") 
                  res.QE <- try(clogit(event ~ X.QE.l + strata(study.l), 
                    iter.max = control.glm$maxit), silent = silent)
                if (con$routine == "clogistic") 
                  res.QE <- try(clogistic(event ~ X.QE.l, strata = study.l, 
                    iter.max = control.glm$maxit, model = FALSE, 
                    x = FALSE, y = FALSE), silent = silent)
                if (inherits(res.QE, "try-error")) 
                  stop("Cannot fit saturated model.")
                ll.FE <- -1 * .dnchg(c(cbind(coef(res.FE)), 0), 
                  ai = ai, bi = bi, ci = ci, di = di, X.fit = X.fit, 
                  random = FALSE, verbose = con$verbose, digits = digits, 
                  dnchgcalc = con$dnchgcalc, dnchgprec = con$dnchgprec)
                ll.QE <- -1 * .dnchg(c(cbind(coef(res.QE)), 0), 
                  ai = ai, bi = bi, ci = ci, di = di, X.fit = X.QE, 
                  random = FALSE, verbose = con$verbose, digits = digits, 
                  dnchgcalc = con$dnchgcalc, dnchgprec = con$dnchgprec)
                b2.QE <- cbind(coef(res.QE)[-seq.int(p)])
                vb2.QE <- vcov(res.QE)[-seq.int(p), -seq.int(p), 
                  drop = FALSE]
            }
            if (method == "ML") {
                if (con$verbose) 
                  message("Fitting ML model.")
                res.ML <- try(optim(par = c(b, log(tau2 + 0.001)), 
                  .dnchg, method = con$method, hessian = TRUE, 
                  ai = ai, bi = bi, ci = ci, di = di, X.fit = X.fit, 
                  random = TRUE, verbose = con$verbose, digits = digits, 
                  dnchgcalc = con$dnchgcalc, dnchgprec = con$dnchgprec, 
                  control.int = control.int, control = control.optim), 
                  silent = silent)
                if (inherits(res.ML, "try-error") || res.ML$convergence != 
                  0) 
                  stop("Cannot fit GLMM model.")
                ll.ML <- -1 * res.ML$value
            }
            if (method == "FE") {
                if (con$routine == "optim") {
                  b <- cbind(res.FE$par[seq.int(p)])
                  chol.h <- try(chol(res.FE$hessian[seq.int(p), 
                    seq.int(p)]), silent = silent)
                  if (class(chol.h) == "try-error") 
                    stop("Cannot invert Hessian for fixed-effects model.")
                  vb <- chol2inv(chol.h)
                }
                if (con$routine == "clogit" || con$routine == 
                  "clogistic") {
                  b <- cbind(coef(res.FE))
                  vb <- vcov(res.FE)
                }
                tau2 <- 0
                sigma2 <- NA
                parms <- p
                p.eff <- p
                k.eff <- k
            }
            if (method == "ML") {
                b <- cbind(res.ML$par[seq.int(p)])
                chol.h <- try(chol(res.ML$hessian), silent = silent)
                if (class(chol.h) == "try-error") 
                  stop("Cannot invert Hessian for random/mixed-effects model.")
                vb.f <- chol2inv(chol.h)
                vb <- vb.f[seq.int(p), seq.int(p), drop = FALSE]
                tau2 <- exp(res.ML$par[p + 1])
                sigma2 <- NA
                parms <- p + 1
                p.eff <- p
                k.eff <- k
                se.tau2 <- sqrt(vb.f[p + 1, p + 1]) * tau2
            }
        }
    }
    if (is.element(measure, c("PLO", "IRLN"))) {
        if (measure == "PLO") {
            dat.grp <- cbind(xi = xi, mi = mi)
            dat.fam <- binomial
            dat.off <- NULL
        }
        if (measure == "IRLN") {
            dat.grp <- xi
            dat.fam <- poisson
            dat.off <- log(ti)
        }
        study <- factor(seq.int(k))
        X.fit <- X
        if (con$verbose) 
            message("Fitting FE model.")
        res.FE <- try(glm(dat.grp ~ -1 + X.fit, offset = dat.off, 
            family = dat.fam, control = control.glm), silent = silent)
        if (inherits(res.FE, "try-error")) 
            stop("Cannot fit FE model.")
        if (con$verbose) 
            message("Fitting saturated model.")
        X.QE <- model.matrix(~-1 + X.fit + study)
        if (con$verbose) {
            res.QE <- try(glm(dat.grp ~ -1 + X.QE, offset = dat.off, 
                family = dat.fam, control = control.glm), silent = silent)
        }
        else {
            res.QE <- suppressWarnings(try(glm(dat.grp ~ -1 + 
                X.QE, offset = dat.off, family = dat.fam, control = control.glm), 
                silent = silent))
        }
        if (inherits(res.QE, "try-error")) 
            stop("Cannot fit saturated model.")
        ll.FE <- c(logLik(res.FE))
        ll.QE <- c(logLik(res.QE))
        b2.QE <- cbind(na.omit(coef(res.QE)[-seq.int(p)]))
        vb2.QE <- vcov(res.QE)[-seq.int(p), -seq.int(p), drop = FALSE]
        if (method == "ML") {
            if (con$verbose) 
                message("Fitting ML model.")
            if (con$verbose) {
                res.ML <- try(lmer(dat.grp ~ -1 + X.fit + (1 | 
                  study), offset = dat.off, family = dat.fam, 
                  nAGQ = nAGQ, verbose = con$verbose, control = control.lmer), 
                  silent = silent)
            }
            else {
                res.ML <- suppressMessages(try(lmer(dat.grp ~ 
                  -1 + X.fit + (1 | study), offset = dat.off, 
                  family = dat.fam, nAGQ = nAGQ, verbose = con$verbose, 
                  control = control.lmer), silent = silent))
            }
            if (inherits(res.ML, "try-error")) 
                stop("Cannot fit ML model.")
            ll.ML <- ll.QE - 1/2 * lme4::deviance(res.ML)
        }
        if (method == "FE") {
            b <- cbind(coef(res.FE)[seq.int(p)])
            vb <- vcov(res.FE)[seq.int(p), seq.int(p), drop = FALSE]
            tau2 <- 0
            sigma2 <- NA
            parms <- p
            p.eff <- p
            k.eff <- k
        }
        if (method == "ML") {
            b <- cbind(fixef(res.ML)[seq.int(p)])
            vb <- as.matrix(lme4::vcov(res.ML))[seq.int(p), seq.int(p), 
                drop = FALSE]
            tau2 <- VarCorr(res.ML)[[1]][1]
            sigma2 <- NA
            parms <- p + 1
            p.eff <- p
            k.eff <- k
        }
    }
    rownames(vb) <- colnames(vb) <- rownames(b) <- colnames(X)
    if (measure != "OR" || model != "CM.EL" || con$routine != 
        "optim") {
        if (dim(vb2.QE)[1] > 0) {
            chol.h <- try(chol(vb2.QE), silent = silent)
            if (class(chol.h) == "try-error") 
                stop("Cannot invert Hessian for saturated model.")
            QE.Wld <- c(t(b2.QE) %*% chol2inv(chol.h) %*% b2.QE)
        }
        else {
            QE.Wld <- 0
        }
    }
    QE.LRT <- -2 * (ll.FE - ll.QE)
    QE.Wld[QE.Wld < 0] <- 0
    QE.LRT[QE.LRT < 0] <- 0
    QE.df <- k - p
    if (QE.df > 0L) {
        QEp.Wld <- pchisq(QE.Wld, df = QE.df, lower.tail = FALSE)
        QEp.LRT <- pchisq(QE.LRT, df = QE.df, lower.tail = FALSE)
    }
    else {
        QEp.Wld <- 1
        QEp.LRT <- 1
    }
    wi <- 1/vi.I2
    W <- .diag(wi)
    stXWX <- .invcalc(X = X.I2, W = W, k = k.I2)
    P <- W - W %*% X.I2 %*% stXWX %*% crossprod(X.I2, W)
    vi.avg <- (k.I2 - p)/.tr(P)
    I2 <- 100 * tau2/(vi.avg + tau2)
    H2 <- tau2/vi.avg + 1
    chol.h <- try(chol(vb[btt, btt]), silent = silent)
    if (class(chol.h) == "try-error") 
        stop("Cannot invert Hessian for QM test.")
    QM <- c(t(b)[btt] %*% chol2inv(chol.h) %*% b[btt])
    se <- sqrt(diag(vb))
    names(se) <- NULL
    zval <- c(b/se)
    if (knha || robust) {
        QM <- QM/m
        QMp <- pf(QM, df1 = m, df2 = k - p, lower.tail = FALSE)
        pval <- 2 * pt(abs(zval), df = k - p, lower.tail = FALSE)
        crit <- qt(alpha/2, df = k - p, lower.tail = FALSE)
    }
    else {
        QMp <- pchisq(QM, df = m, lower.tail = FALSE)
        pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)
        crit <- qnorm(alpha/2, lower.tail = FALSE)
    }
    ci.lb <- c(b - crit * se)
    ci.ub <- c(b + crit * se)
    ll.ML <- ifelse(method == "FE", ll.FE, ll.ML)
    ll.REML <- NA
    dev.ML <- -2 * (ll.ML - ll.QE)
    AIC.ML <- -2 * ll.ML + 2 * (parms)
    BIC.ML <- -2 * ll.ML + (parms) * log(k.eff)
    dev.REML <- NA
    AIC.REML <- NA
    BIC.REML <- NA
    fit.stats <- matrix(c(ll.ML, dev.ML, AIC.ML, BIC.ML, ll.REML, 
        dev.REML, AIC.REML, BIC.REML), ncol = 2, byrow = FALSE)
    dimnames(fit.stats) <- list(c("ll", "dev", "AIC", "BIC"), 
        c("ML", "REML"))
    fit.stats <- data.frame(fit.stats)
    weighted <- TRUE
    robust <- FALSE
    res <- list(b, se, zval, pval, ci.lb, ci.ub, vb, tau2, se.tau2, 
        sigma2, k, k.f, k.yi, k.eff, p, p.eff, parms, m, QE.Wld, 
        QEp.Wld, QE.LRT, QEp.LRT, QE.df, QM, QMp, I2, H2, int.only, 
        int.incl, yi, vi, X, yi.f, vi.f, X.f, ai, bi, ci, di, 
        ai.f, bi.f, ci.f, di.f, x1i, x2i, t1i, t2i, x1i.f, x2i.f, 
        t1i.f, t2i.f, xi, mi, ti, xi.f, mi.f, ti.f, ni, ni.f, 
        ids, not.na, not.na.yivi, slab, slab.null, measure, method, 
        model, weighted, knha, robust, btt, intercept, digits, 
        level, con, add, to, drop00, fit.stats)
    names(res) <- c("b", "se", "zval", "pval", "ci.lb", "ci.ub", 
        "vb", "tau2", "se.tau2", "sigma2", "k", "k.f", "k.yi", 
        "k.eff", "p", "p.eff", "parms", "m", "QE.Wld", "QEp.Wld", 
        "QE.LRT", "QEp.LRT", "QE.df", "QM", "QMp", "I2", "H2", 
        "int.only", "int.incl", "yi", "vi", "X", "yi.f", "vi.f", 
        "X.f", "ai", "bi", "ci", "di", "ai.f", "bi.f", "ci.f", 
        "di.f", "x1i", "x2i", "t1i", "t2i", "x1i.f", "x2i.f", 
        "t1i.f", "t2i.f", "xi", "mi", "ti", "xi.f", "mi.f", "ti.f", 
        "ni", "ni.f", "ids", "not.na", "not.na.yivi", "slab", 
        "slab.null", "measure", "method", "model", "weighted", 
        "knha", "robust", "btt", "intercept", "digits", "level", 
        "control", "add", "to", "drop00", "fit.stats")
    class(res) <- c("rma.glmm", "rma")
    return(res)
}
rma.mh <-
function (ai, bi, ci, di, n1i, n2i, x1i, x2i, t1i, t2i, measure = "OR", 
    data, slab, subset, add = 1/2, to = "only0", drop00 = TRUE, 
    level = 95, digits = 4) 
{
    if (!is.element(measure, c("OR", "RR", "RD", "IRR"))) 
        stop("Mantel-Haenszel method can only be used with measures OR, RR, RD, and IRR.")
    if (length(add) == 1) 
        add <- c(add, 0)
    if (length(add) != 2) 
        stop("Argument 'add' should specify two values (see 'help(rma.mh)').")
    if (length(to) == 1) 
        to <- c(to, "none")
    if (length(to) != 2) 
        stop("Argument 'to' should specify two values (see 'help(rma.mh)').")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    if (!is.element(to[1], c("all", "only0", "if0all", "none"))) 
        stop("Unknown 'to' argument specified.")
    if (!is.element(to[2], c("all", "only0", "if0all", "none"))) 
        stop("Unknown 'to' argument specified.")
    if (missing(data)) 
        data <- NULL
    if (is.null(data)) {
        data <- sys.frame(sys.parent())
    }
    else {
        if (!is.data.frame(data)) {
            data <- data.frame(data)
        }
    }
    mf <- match.call()
    mf.slab <- mf[[match("slab", names(mf))]]
    mf.subset <- mf[[match("subset", names(mf))]]
    slab <- eval(mf.slab, data, enclos = sys.frame(sys.parent()))
    subset <- eval(mf.subset, data, enclos = sys.frame(sys.parent()))
    if (is.element(measure, c("RR", "OR", "RD"))) {
        x1i <- x2i <- t1i <- t2i <- x1i.f <- x2i.f <- t1i.f <- t2i.f <- NA
        mf.ai <- mf[[match("ai", names(mf))]]
        mf.bi <- mf[[match("bi", names(mf))]]
        mf.ci <- mf[[match("ci", names(mf))]]
        mf.di <- mf[[match("di", names(mf))]]
        mf.n1i <- mf[[match("n1i", names(mf))]]
        mf.n2i <- mf[[match("n2i", names(mf))]]
        ai <- eval(mf.ai, data, enclos = sys.frame(sys.parent()))
        bi <- eval(mf.bi, data, enclos = sys.frame(sys.parent()))
        ci <- eval(mf.ci, data, enclos = sys.frame(sys.parent()))
        di <- eval(mf.di, data, enclos = sys.frame(sys.parent()))
        n1i <- eval(mf.n1i, data, enclos = sys.frame(sys.parent()))
        n2i <- eval(mf.n2i, data, enclos = sys.frame(sys.parent()))
        if (is.null(bi)) 
            bi <- n1i - ai
        if (is.null(di)) 
            di <- n2i - ci
        ni <- ai + bi + ci + di
        k <- length(ai)
        ids <- seq.int(k)
        if (is.null(slab)) {
            slab.null <- TRUE
            slab <- seq.int(k)
        }
        else {
            if (any(duplicated(slab))) 
                stop("Study labels must be unique.")
            if (length(slab) != length(ai)) 
                stop("Study labels not of same length as data.")
            slab.null <- FALSE
        }
        if (!is.null(subset)) {
            ai <- ai[subset]
            bi <- bi[subset]
            ci <- ci[subset]
            di <- di[subset]
            ni <- ni[subset]
            slab <- slab[subset]
            ids <- ids[subset]
            k <- length(ai)
        }
        dat <- escalc(ai = ai, bi = bi, ci = ci, di = di, add = add[1], 
            to = to[1], drop00 = drop00, measure = measure)
        yi <- dat$yi
        vi <- dat$vi
        ai.f <- ai
        bi.f <- bi
        ci.f <- ci
        di.f <- di
        yi.f <- yi
        vi.f <- vi
        ni.f <- ni
        k.f <- k
        aibicidi.na <- is.na(cbind(ai, bi, ci, di))
        if (any(aibicidi.na)) {
            not.na <- apply(aibicidi.na, MARGIN = 1, sum) == 
                0L
            if (na.act == "na.omit" || na.act == "na.exclude" || 
                na.act == "na.pass") {
                ai <- ai[not.na]
                bi <- bi[not.na]
                ci <- ci[not.na]
                di <- di[not.na]
                k <- length(ai)
                warning("Tables with NAs omitted from model fitting.")
            }
            if (na.act == "na.fail") 
                stop("Missing values in tables.")
        }
        else {
            not.na <- rep(TRUE, k)
        }
        if (k < 1) 
            stop("Processing terminated since k = 0.")
        yivi.na <- is.na(cbind(yi, vi))
        if (any(yivi.na)) {
            not.na.yivi <- apply(yivi.na, MARGIN = 1, sum) == 
                0L
            if (na.act == "na.omit" || na.act == "na.exclude" || 
                na.act == "na.pass") {
                yi <- yi[not.na.yivi]
                vi <- vi[not.na.yivi]
                ni <- ni[not.na.yivi]
                warning("Some yi/vi values are NA.")
            }
            if (na.act == "na.fail") 
                stop("Missing yi/vi values.")
        }
        else {
            not.na.yivi <- rep(TRUE, k)
        }
        k.yi <- length(yi)
        if (to[2] == "all") {
            ai <- ai + add[2]
            bi <- bi + add[2]
            ci <- ci + add[2]
            di <- di + add[2]
        }
        if (to[2] == "only0") {
            id0 <- c(ai == 0L | bi == 0L | ci == 0L | di == 0L)
            ai[id0] <- ai[id0] + add[2]
            bi[id0] <- bi[id0] + add[2]
            ci[id0] <- ci[id0] + add[2]
            di[id0] <- di[id0] + add[2]
        }
        if (to[2] == "if0all") {
            id0 <- c(ai == 0L | bi == 0L | ci == 0L | di == 0L)
            if (any(id0)) {
                ai <- ai + add[2]
                bi <- bi + add[2]
                ci <- ci + add[2]
                di <- di + add[2]
            }
        }
        alpha <- (100 - level)/100
        n1i <- ai + bi
        n2i <- ci + di
        Ni <- ai + bi + ci + di
    }
    if (is.element(measure, c("IRR"))) {
        ai <- bi <- ci <- di <- ai.f <- bi.f <- ci.f <- di.f <- NA
        mf.x1i <- mf[[match("x1i", names(mf))]]
        mf.x2i <- mf[[match("x2i", names(mf))]]
        mf.t1i <- mf[[match("t1i", names(mf))]]
        mf.t2i <- mf[[match("t2i", names(mf))]]
        x1i <- eval(mf.x1i, data, enclos = sys.frame(sys.parent()))
        x2i <- eval(mf.x2i, data, enclos = sys.frame(sys.parent()))
        t1i <- eval(mf.t1i, data, enclos = sys.frame(sys.parent()))
        t2i <- eval(mf.t2i, data, enclos = sys.frame(sys.parent()))
        ni <- t1i + t2i
        k <- length(x1i)
        ids <- seq.int(k)
        if (is.null(slab)) {
            slab.null <- TRUE
            slab <- seq.int(k)
        }
        else {
            if (any(duplicated(slab))) 
                stop("Study labels must be unique.")
            if (length(slab) != length(ai)) 
                stop("Study labels not of same length as data.")
            slab.null <- FALSE
        }
        if (!is.null(subset)) {
            x1i <- x1i[subset]
            x2i <- x2i[subset]
            t1i <- t1i[subset]
            t2i <- t2i[subset]
            ni <- ni[subset]
            slab <- slab[subset]
            ids <- ids[subset]
            k <- length(x1i)
        }
        dat <- escalc(x1i = x1i, x2i = x2i, t1i = t1i, t2i = t2i, 
            add = add[1], to = to[1], drop00 = drop00, measure = measure)
        yi <- dat$yi
        vi <- dat$vi
        x1i.f <- x1i
        x2i.f <- x2i
        t1i.f <- t1i
        t2i.f <- t2i
        yi.f <- yi
        vi.f <- vi
        ni.f <- ni
        k.f <- k
        x1ix2it1it2i.na <- is.na(cbind(x1i, x2i, t1i, t2i))
        if (any(x1ix2it1it2i.na)) {
            not.na <- apply(x1ix2it1it2i.na, MARGIN = 1, sum) == 
                0L
            if (na.act == "na.omit" || na.act == "na.exclude" || 
                na.act == "na.pass") {
                x1i <- x1i[not.na]
                x2i <- x2i[not.na]
                t1i <- t1i[not.na]
                t2i <- t2i[not.na]
                k <- length(x1i)
                warning("Tables with NAs omitted from model fitting.")
            }
            if (na.act == "na.fail") 
                stop("Missing values in tables.")
        }
        else {
            not.na <- rep(TRUE, k)
        }
        if (k < 1) 
            stop("Processing terminated since k = 0.")
        yivi.na <- is.na(cbind(yi, vi))
        if (any(yivi.na)) {
            not.na.yivi <- apply(yivi.na, MARGIN = 1, sum) == 
                0L
            if (na.act == "na.omit" || na.act == "na.exclude" || 
                na.act == "na.pass") {
                yi <- yi[not.na.yivi]
                vi <- vi[not.na.yivi]
                ni <- ni[not.na.yivi]
                warning("Some yi/vi values are NA.")
            }
            if (na.act == "na.fail") 
                stop("Missing yi/vi values.")
        }
        else {
            not.na.yivi <- rep(TRUE, k)
        }
        k.yi <- length(yi)
        if (to[2] == "all") {
            x1i <- x1i + add[2]
            x2i <- x2i + add[2]
        }
        if (to[2] == "only0") {
            id0 <- c(x1i == 0L | x2i == 0L)
            x1i[id0] <- x1i[id0] + add[2]
            x2i[id0] <- x2i[id0] + add[2]
        }
        if (to[2] == "if0all") {
            id0 <- c(x1i == 0L | x2i == 0L)
            if (any(id0)) {
                x1i <- x1i + add[2]
                x2i <- x2i + add[2]
            }
        }
        alpha <- (100 - level)/100
        Ti <- t1i + t2i
    }
    CO <- COp <- CMH <- CMHp <- BD <- BDp <- TA <- TAp <- k.pos <- NA
    if (measure == "OR") {
        Pi <- ai/Ni + di/Ni
        Qi <- bi/Ni + ci/Ni
        Ri <- (ai/Ni) * di
        Si <- (bi/Ni) * ci
        R <- sum(Ri)
        S <- sum(Si)
        if (identical(R, 0) || identical(S, 0)) {
            b.exp <- NA
            b <- NA
            se <- NA
            zval <- NA
            pval <- NA
            ci.lb <- NA
            ci.ub <- NA
        }
        else {
            b.exp <- R/S
            b <- log(b.exp)
            se <- sqrt(1/2 * (sum(Pi * Ri)/R^2 + sum(Pi * Si + 
                Qi * Ri)/(R * S) + sum(Qi * Si)/S^2))
            zval <- b/se
            pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)
            ci.lb <- b - qnorm(alpha/2, lower.tail = FALSE) * 
                se
            ci.ub <- b + qnorm(alpha/2, lower.tail = FALSE) * 
                se
        }
        names(b) <- "intrcpt"
        vb <- matrix(se^2, dimnames = list("intrcpt", "intrcpt"))
        xt <- ai + ci
        yt <- bi + di
        if (identical(sum(xt), 0) || identical(sum(yt), 0)) {
            CO <- NA
            COp <- NA
            CMH <- NA
            CMHp <- NA
        }
        else {
            CO <- (abs(sum(ai - (n1i/Ni) * xt)))^2/sum((n1i/Ni) * 
                (n2i/Ni) * (xt * yt/Ni))
            COp <- pchisq(CO, df = 1, lower.tail = FALSE)
            CMH <- (abs(sum(ai - (n1i/Ni) * xt)) - 0.5)^2/sum((n1i/Ni) * 
                (n2i/Ni) * (xt * yt/(Ni - 1)))
            CMHp <- pchisq(CMH, df = 1, lower.tail = FALSE)
        }
        if (is.na(b)) {
            BD <- NA
            TA <- NA
            BDp <- NA
            TAp <- NA
            k.pos <- 0
        }
        else {
            if (identical(b.exp, 1)) {
                N11 <- (n1i/Ni) * xt
            }
            else {
                A <- b.exp * (n1i + xt) + (n2i - xt)
                B <- sqrt(A^2 - 4 * n1i * xt * b.exp * (b.exp - 
                  1))
                N11 <- (A - B)/(2 * (b.exp - 1))
            }
            pos <- (N11 > 0) & (xt > 0) & (yt > 0)
            k.pos <- sum(pos)
            N11 <- N11[pos]
            N12 <- n1i[pos] - N11
            N21 <- xt[pos] - N11
            N22 <- N11 - n1i[pos] - xt[pos] + Ni[pos]
            BD <- sum((ai[pos] - N11)^2/(1/N11 + 1/N12 + 1/N21 + 
                1/N22)^(-1))
            TA <- BD - sum(ai[pos] - N11)^2/sum((1/N11 + 1/N12 + 
                1/N21 + 1/N22)^(-1))
            if (k.pos > 1) {
                BDp <- pchisq(BD, df = k.pos - 1, lower.tail = FALSE)
                TAp <- pchisq(TA, df = k.pos - 1, lower.tail = FALSE)
            }
            else {
                BDp <- NA
                TAp <- NA
            }
        }
    }
    if (measure == "RR") {
        R <- sum(ai * (n2i/Ni))
        S <- sum(ci * (n1i/Ni))
        if (identical(sum(ai), 0) || identical(sum(ci), 0)) {
            b.exp <- NA
            b <- NA
            se <- NA
            zval <- NA
            pval <- NA
            ci.lb <- NA
            ci.ub <- NA
        }
        else {
            b.exp <- R/S
            b <- log(b.exp)
            se <- sqrt(sum(((n1i/Ni) * (n2i/Ni) * (ai + ci) - 
                (ai/Ni) * ci))/(R * S))
            zval <- b/se
            pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)
            ci.lb <- b - qnorm(alpha/2, lower.tail = FALSE) * 
                se
            ci.ub <- b + qnorm(alpha/2, lower.tail = FALSE) * 
                se
        }
        names(b) <- "intrcpt"
        vb <- matrix(se^2, dimnames = list("intrcpt", "intrcpt"))
    }
    if (measure == "RD") {
        b <- sum(ai * (n2i/Ni) - ci * (n1i/Ni))/sum(n1i * (n2i/Ni))
        se <- sqrt((b * (sum(ci * (n1i/Ni)^2 - ai * (n2i/Ni)^2 + 
            (n1i/Ni) * (n2i/Ni) * (n2i - n1i)/2)) + sum(ai * 
            (n2i - ci)/Ni + ci * (n1i - ai)/Ni)/2)/sum(n1i * 
            (n2i/Ni))^2)
        zval <- b/se
        pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)
        ci.lb <- b - qnorm(alpha/2, lower.tail = FALSE) * se
        ci.ub <- b + qnorm(alpha/2, lower.tail = FALSE) * se
        names(b) <- "intrcpt"
        vb <- matrix(se^2, dimnames = list("intrcpt", "intrcpt"))
    }
    if (measure == "IRR") {
        R <- sum(x1i * (t2i/Ti))
        S <- sum(x2i * (t1i/Ti))
        if (identical(sum(x1i), 0) || identical(sum(x2i), 0)) {
            b.exp <- NA
            b <- NA
            se <- NA
            zval <- NA
            pval <- NA
            ci.lb <- NA
            ci.ub <- NA
        }
        else {
            b.exp <- R/S
            b <- log(b.exp)
            se <- sqrt(sum((t1i/Ti) * (t2i/Ti) * (x1i + x2i))/(R * 
                S))
            zval <- b/se
            pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)
            ci.lb <- b - qnorm(alpha/2, lower.tail = FALSE) * 
                se
            ci.ub <- b + qnorm(alpha/2, lower.tail = FALSE) * 
                se
        }
        names(b) <- "intrcpt"
        vb <- matrix(se^2, dimnames = list("intrcpt", "intrcpt"))
    }
    wi <- 1/vi
    QE <- sum(wi * (yi - b)^2)
    if (k.yi - 1 >= 1) {
        QEp <- pchisq(QE, df = k.yi - 1, lower.tail = FALSE)
    }
    else {
        QEp <- 1
    }
    ll.ML <- -1/2 * (k.yi) * log(2 * base::pi) - 1/2 * sum(log(vi)) - 
        1/2 * QE
    ll.REML <- -1/2 * (k.yi - 1) * log(2 * base::pi) + 1/2 * 
        log(k.yi) - 1/2 * sum(log(vi)) - 1/2 * log(sum(wi)) - 
        1/2 * QE
    dev.ML <- -2 * (ll.ML - sum(dnorm(yi, mean = yi, sd = sqrt(vi), 
        log = TRUE)))
    AIC.ML <- -2 * ll.ML + 2
    BIC.ML <- -2 * ll.ML + log(k.yi)
    dev.REML <- -2 * (ll.REML - 0)
    AIC.REML <- -2 * ll.REML + 2
    BIC.REML <- -2 * ll.REML + log(k.yi - 1)
    fit.stats <- matrix(c(ll.ML, dev.ML, AIC.ML, BIC.ML, ll.REML, 
        dev.REML, AIC.REML, BIC.REML), ncol = 2, byrow = FALSE)
    dimnames(fit.stats) <- list(c("ll", "dev", "AIC", "BIC"), 
        c("ML", "REML"))
    fit.stats <- data.frame(fit.stats)
    parms <- 1
    p <- 1
    p.eff <- 1
    k.eff <- k
    tau2 <- 0
    X.f <- cbind(rep(1, k.f))
    intercept <- TRUE
    int.only <- TRUE
    method <- "FE"
    weighted <- TRUE
    knha <- FALSE
    robust <- FALSE
    res <- list(b, se, zval, pval, ci.lb, ci.ub, vb, tau2, k, 
        k.f, k.yi, k.pos, k.eff, p, parms, QE, QEp, CO, COp, 
        CMH, CMHp, BD, BDp, TA, TAp, int.only, yi, vi, yi.f, 
        vi.f, X.f, ai, bi, ci, di, x1i, x2i, t1i, t2i, ai.f, 
        bi.f, ci.f, di.f, x1i.f, x2i.f, t1i.f, t2i.f, ni, ni.f, 
        ids, not.na, not.na.yivi, slab, slab.null, measure, method, 
        weighted, knha, robust, intercept, digits, level, add, 
        to, drop00, fit.stats)
    names(res) <- c("b", "se", "zval", "pval", "ci.lb", "ci.ub", 
        "vb", "tau2", "k", "k.f", "k.yi", "k.pos", "k.eff", "p", 
        "parms", "QE", "QEp", "CO", "COp", "CMH", "CMHp", "BD", 
        "BDp", "TA", "TAp", "int.only", "yi", "vi", "yi.f", "vi.f", 
        "X.f", "ai", "bi", "ci", "di", "x1i", "x2i", "t1i", "t2i", 
        "ai.f", "bi.f", "ci.f", "di.f", "x1i.f", "x2i.f", "t1i.f", 
        "t2i.f", "ni", "ni.f", "ids", "not.na", "not.na.yivi", 
        "slab", "slab.null", "measure", "method", "weighted", 
        "knha", "robust", "intercept", "digits", "level", "add", 
        "to", "drop00", "fit.stats")
    class(res) <- c("rma.mh", "rma")
    return(res)
}
rma.peto <-
function (ai, bi, ci, di, n1i, n2i, data, slab, subset, add = 1/2, 
    to = "only0", drop00 = TRUE, level = 95, digits = 4) 
{
    if (length(add) == 1) 
        add <- c(add, 0)
    if (length(add) != 2) 
        stop("Argument 'add' should specify two values (see 'help(rma.peto)').")
    if (length(to) == 1) 
        to <- c(to, "none")
    if (length(to) != 2) 
        stop("Argument 'to' should specify two values (see 'help(rma.peto)').")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    if (!is.element(to[1], c("all", "only0", "if0all", "none"))) 
        stop("Unknown 'to' argument specified.")
    if (!is.element(to[2], c("all", "only0", "if0all", "none"))) 
        stop("Unknown 'to' argument specified.")
    if (missing(data)) 
        data <- NULL
    if (is.null(data)) {
        data <- sys.frame(sys.parent())
    }
    else {
        if (!is.data.frame(data)) {
            data <- data.frame(data)
        }
    }
    mf <- match.call()
    mf.slab <- mf[[match("slab", names(mf))]]
    mf.subset <- mf[[match("subset", names(mf))]]
    slab <- eval(mf.slab, data, enclos = sys.frame(sys.parent()))
    subset <- eval(mf.subset, data, enclos = sys.frame(sys.parent()))
    mf.ai <- mf[[match("ai", names(mf))]]
    mf.bi <- mf[[match("bi", names(mf))]]
    mf.ci <- mf[[match("ci", names(mf))]]
    mf.di <- mf[[match("di", names(mf))]]
    mf.n1i <- mf[[match("n1i", names(mf))]]
    mf.n2i <- mf[[match("n2i", names(mf))]]
    ai <- eval(mf.ai, data, enclos = sys.frame(sys.parent()))
    bi <- eval(mf.bi, data, enclos = sys.frame(sys.parent()))
    ci <- eval(mf.ci, data, enclos = sys.frame(sys.parent()))
    di <- eval(mf.di, data, enclos = sys.frame(sys.parent()))
    n1i <- eval(mf.n1i, data, enclos = sys.frame(sys.parent()))
    n2i <- eval(mf.n2i, data, enclos = sys.frame(sys.parent()))
    if (is.null(bi)) 
        bi <- n1i - ai
    if (is.null(di)) 
        di <- n2i - ci
    ni <- ai + bi + ci + di
    k <- length(ai)
    ids <- seq.int(k)
    if (is.null(slab)) {
        slab.null <- TRUE
        slab <- seq.int(k)
    }
    else {
        if (any(duplicated((slab)))) 
            stop("Study labels must be unique.")
        if (length(slab) != length(ai)) 
            stop("Study labels not of same length as data.")
        slab.null <- FALSE
    }
    if (!is.null(subset)) {
        ai <- ai[subset]
        bi <- bi[subset]
        ci <- ci[subset]
        di <- di[subset]
        ni <- ni[subset]
        slab <- slab[subset]
        ids <- ids[subset]
        k <- length(ai)
    }
    dat <- escalc(ai = ai, bi = bi, ci = ci, di = di, add = add[1], 
        to = to[1], drop00 = drop00, measure = "PETO")
    yi <- dat$yi
    vi <- dat$vi
    ai.f <- ai
    bi.f <- bi
    ci.f <- ci
    di.f <- di
    yi.f <- yi
    vi.f <- vi
    ni.f <- ni
    k.f <- k
    aibicidi.na <- is.na(cbind(ai, bi, ci, di))
    if (any(aibicidi.na)) {
        not.na <- apply(aibicidi.na, MARGIN = 1, sum) == 0L
        if (na.act == "na.omit" || na.act == "na.exclude" || 
            na.act == "na.pass") {
            ai <- ai[not.na]
            bi <- bi[not.na]
            ci <- ci[not.na]
            di <- di[not.na]
            k <- length(ai)
            warning("Tables with NAs omitted from model fitting.")
        }
        if (na.act == "na.fail") 
            stop("Missing values in tables.")
    }
    else {
        not.na <- rep(TRUE, k)
    }
    if (k < 1) 
        stop("Processing terminated since k = 0.")
    yivi.na <- is.na(cbind(yi, vi))
    if (any(yivi.na)) {
        not.na.yivi <- apply(yivi.na, MARGIN = 1, sum) == 0L
        if (na.act == "na.omit" || na.act == "na.exclude" || 
            na.act == "na.pass") {
            yi <- yi[not.na.yivi]
            vi <- vi[not.na.yivi]
            ni <- ni[not.na.yivi]
            warning("Some yi/vi values are NA.")
        }
        if (na.act == "na.fail") 
            stop("Missing yi/vi values.")
    }
    else {
        not.na.yivi <- rep(TRUE, k)
    }
    k.yi <- length(yi)
    if (to[2] == "all") {
        ai <- ai + add[2]
        bi <- bi + add[2]
        ci <- ci + add[2]
        di <- di + add[2]
    }
    if (to[2] == "only0") {
        id0 <- c(ai == 0L | bi == 0L | ci == 0L | di == 0L)
        ai[id0] <- ai[id0] + add[2]
        bi[id0] <- bi[id0] + add[2]
        ci[id0] <- ci[id0] + add[2]
        di[id0] <- di[id0] + add[2]
    }
    if (to[2] == "if0all") {
        id0 <- c(ai == 0L | bi == 0L | ci == 0L | di == 0L)
        if (any(id0)) {
            ai <- ai + add[2]
            bi <- bi + add[2]
            ci <- ci + add[2]
            di <- di + add[2]
        }
    }
    alpha <- (100 - level)/100
    n1i <- ai + bi
    n2i <- ci + di
    Ni <- ai + bi + ci + di
    xt <- ai + ci
    yt <- bi + di
    Ei <- xt * n1i/Ni
    Vi <- xt * yt * (n1i/Ni) * (n2i/Ni)/(Ni - 1)
    sumVi <- sum(Vi)
    if (sumVi == 0L) 
        stop("All tables have either only events or no events at all. Peto's method cannot be used.")
    b <- sum(ai - Ei)/sumVi
    se <- sqrt(1/sumVi)
    zval <- b/se
    pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)
    ci.lb <- b - qnorm(alpha/2, lower.tail = FALSE) * se
    ci.ub <- b + qnorm(alpha/2, lower.tail = FALSE) * se
    names(b) <- "intrcpt"
    vb <- matrix(se^2, dimnames = list("intrcpt", "intrcpt"))
    k.pos <- sum(Vi > 0)
    Vi[Vi == 0] <- NA
    QE <- sum((ai - Ei)^2/Vi, na.rm = TRUE) - sum(ai - Ei)^2/sum(Vi, 
        na.rm = TRUE)
    if (k.pos >= 1) {
        QEp <- pchisq(QE, df = k.yi - 1, lower.tail = FALSE)
    }
    else {
        QEp <- 1
    }
    wi <- 1/vi
    RSS <- sum(wi * (yi - b)^2)
    ll.ML <- -1/2 * (k.yi) * log(2 * base::pi) - 1/2 * sum(log(vi)) - 
        1/2 * RSS
    ll.REML <- -1/2 * (k.yi - 1) * log(2 * base::pi) + 1/2 * 
        log(k.yi) - 1/2 * sum(log(vi)) - 1/2 * log(sum(wi)) - 
        1/2 * RSS
    dev.ML <- -2 * (ll.ML - sum(dnorm(yi, mean = yi, sd = sqrt(vi), 
        log = TRUE)))
    AIC.ML <- -2 * ll.ML + 2
    BIC.ML <- -2 * ll.ML + log(k.yi)
    dev.REML <- -2 * (ll.REML - 0)
    AIC.REML <- -2 * ll.REML + 2
    BIC.REML <- -2 * ll.REML + log(k.yi - 1)
    fit.stats <- matrix(c(ll.ML, dev.ML, AIC.ML, BIC.ML, ll.REML, 
        dev.REML, AIC.REML, BIC.REML), ncol = 2, byrow = FALSE)
    dimnames(fit.stats) <- list(c("ll", "dev", "AIC", "BIC"), 
        c("ML", "REML"))
    fit.stats <- data.frame(fit.stats)
    parms <- 1
    p <- 1
    p.eff <- 1
    k.eff <- k
    tau2 <- 0
    X.f <- cbind(rep(1, k.f))
    intercept <- TRUE
    int.only <- TRUE
    measure <- "PETO"
    method <- "FE"
    weighted <- TRUE
    knha <- FALSE
    robust <- FALSE
    res <- list(b, se, zval, pval, ci.lb, ci.ub, vb, tau2, k, 
        k.f, k.yi, k.pos, k.eff, p, parms, QE, QEp, int.only, 
        yi, vi, yi.f, vi.f, X.f, ai, bi, ci, di, ai.f, bi.f, 
        ci.f, di.f, ni, ni.f, ids, not.na, not.na.yivi, slab, 
        slab.null, measure, method, weighted, knha, robust, intercept, 
        digits, level, add, to, drop00, fit.stats)
    names(res) <- c("b", "se", "zval", "pval", "ci.lb", "ci.ub", 
        "vb", "tau2", "k", "k.f", "k.yi", "k.pos", "k.eff", "p", 
        "parms", "QE", "QEp", "int.only", "yi", "vi", "yi.f", 
        "vi.f", "X.f", "ai", "bi", "ci", "di", "ai.f", "bi.f", 
        "ci.f", "di.f", "ni", "ni.f", "ids", "not.na", "not.na.yivi", 
        "slab", "slab.null", "measure", "method", "weighted", 
        "knha", "robust", "intercept", "digits", "level", "add", 
        "to", "drop00", "fit.stats")
    class(res) <- c("rma.peto", "rma")
    return(res)
}
rma.uni <-
function (yi, vi, sei, weights, ai, bi, ci, di, n1i, n2i, x1i, 
    x2i, t1i, t2i, m1i, m2i, sd1i, sd2i, xi, mi, ri, ti, sdi, 
    ni, mods, measure = "GEN", intercept = TRUE, data, slab, 
    subset, add = 1/2, to = "only0", drop00 = FALSE, vtype = "LS", 
    method = "REML", weighted = TRUE, knha = FALSE, level = 95, 
    digits = 4, btt, tau2, verbose = FALSE, control) 
{
    if (!is.element(measure, c("GEN", "RR", "OR", "PETO", "RD", 
        "AS", "PHI", "YUQ", "YUY", "RTET", "PBIT", "OR2D", "IRR", 
        "IRD", "IRSD", "MD", "SMD", "SMDH", "ROM", "RPB", "RBIS", 
        "D2OR", "COR", "UCOR", "ZCOR", "PR", "PLN", "PLO", "PAS", 
        "PFT", "IR", "IRLN", "IRS", "IRFT", "MN", "MC", "SMCC", 
        "SMCR", "ARAW", "AHW", "ABT"))) 
        stop("Unknown 'measure' specified.")
    if (!is.element(method, c("FE", "HS", "HE", "DL", "SJ", "ML", 
        "REML", "EB", "DLit", "SJit", "PM"))) 
        stop("Unknown 'method' specified.")
    if (length(add) > 1) 
        add <- add[1]
    if (length(to) > 1) 
        to <- to[1]
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    if (missing(btt)) 
        btt <- NULL
    if (missing(tau2)) 
        tau2 <- NULL
    if (missing(control)) 
        control <- list()
    if (missing(data)) 
        data <- NULL
    if (is.null(data)) {
        data <- sys.frame(sys.parent())
    }
    else {
        if (!is.data.frame(data)) {
            data <- data.frame(data)
        }
    }
    mf <- match.call()
    mf.slab <- mf[[match("slab", names(mf))]]
    mf.subset <- mf[[match("subset", names(mf))]]
    mf.mods <- mf[[match("mods", names(mf))]]
    slab <- eval(mf.slab, data, enclos = sys.frame(sys.parent()))
    subset <- eval(mf.subset, data, enclos = sys.frame(sys.parent()))
    mods <- eval(mf.mods, data, enclos = sys.frame(sys.parent()))
    ai <- bi <- ci <- di <- x1i <- x2i <- t1i <- t2i <- NA
    mf.yi <- mf[[match("yi", names(mf))]]
    yi <- eval(mf.yi, data, enclos = sys.frame(sys.parent()))
    if (!is.null(yi)) {
        yi.is.formula <- ifelse(class(yi) == "formula", TRUE, 
            FALSE)
        if (yi.is.formula) {
            options(na.action = "na.pass")
            mods <- model.matrix(yi, data = data)
            yi <- model.response(model.frame(yi, data = data))
            options(na.action = na.act)
            names(yi) <- NULL
            intercept <- FALSE
        }
        if (!is.null(attr(yi, "measure"))) 
            measure <- attr(yi, "measure")
        attr(yi, "measure") <- measure
        mf.vi <- mf[[match("vi", names(mf))]]
        mf.sei <- mf[[match("sei", names(mf))]]
        mf.weights <- mf[[match("weights", names(mf))]]
        mf.ni <- mf[[match("ni", names(mf))]]
        vi <- eval(mf.vi, data, enclos = sys.frame(sys.parent()))
        sei <- eval(mf.sei, data, enclos = sys.frame(sys.parent()))
        weights <- eval(mf.weights, data, enclos = sys.frame(sys.parent()))
        ni <- eval(mf.ni, data, enclos = sys.frame(sys.parent()))
        if (is.null(vi)) {
            if (is.null(sei)) {
                if (is.null(weights)) {
                  stop("Need to specify vi, sei, or weights argument.")
                }
                else {
                  vi <- 1/weights
                }
            }
            else {
                vi <- sei^2
            }
        }
        if (is.null(ni) && !is.null(attr(yi, "ni"))) 
            ni <- attr(yi, "ni")
        if (length(yi) != length(vi)) 
            stop("Length of yi and vi (or sei) vectors is not the same.")
        if (!is.null(ni) && (length(yi) != length(ni))) 
            stop("Length of yi and ni vectors is not the same.")
    }
    else {
        if (is.element(measure, c("RR", "OR", "PETO", "RD", "AS", 
            "PHI", "YUQ", "YUY", "RTET", "PBIT", "OR2D"))) {
            mf.ai <- mf[[match("ai", names(mf))]]
            mf.bi <- mf[[match("bi", names(mf))]]
            mf.ci <- mf[[match("ci", names(mf))]]
            mf.di <- mf[[match("di", names(mf))]]
            mf.n1i <- mf[[match("n1i", names(mf))]]
            mf.n2i <- mf[[match("n2i", names(mf))]]
            ai <- eval(mf.ai, data, enclos = sys.frame(sys.parent()))
            bi <- eval(mf.bi, data, enclos = sys.frame(sys.parent()))
            ci <- eval(mf.ci, data, enclos = sys.frame(sys.parent()))
            di <- eval(mf.di, data, enclos = sys.frame(sys.parent()))
            n1i <- eval(mf.n1i, data, enclos = sys.frame(sys.parent()))
            n2i <- eval(mf.n2i, data, enclos = sys.frame(sys.parent()))
            if (is.null(bi)) 
                bi <- n1i - ai
            if (is.null(di)) 
                di <- n2i - ci
            dat <- escalc(measure, ai = ai, bi = bi, ci = ci, 
                di = di, add = add, to = to, drop00 = drop00, 
                vtype = vtype)
        }
        if (is.element(measure, c("IRR", "IRD", "IRSD"))) {
            mf.x1i <- mf[[match("x1i", names(mf))]]
            mf.x2i <- mf[[match("x2i", names(mf))]]
            mf.t1i <- mf[[match("t1i", names(mf))]]
            mf.t2i <- mf[[match("t2i", names(mf))]]
            x1i <- eval(mf.x1i, data, enclos = sys.frame(sys.parent()))
            x2i <- eval(mf.x2i, data, enclos = sys.frame(sys.parent()))
            t1i <- eval(mf.t1i, data, enclos = sys.frame(sys.parent()))
            t2i <- eval(mf.t2i, data, enclos = sys.frame(sys.parent()))
            dat <- escalc(measure, x1i = x1i, x2i = x2i, t1i = t1i, 
                t2i = t2i, add = add, to = to, drop00 = drop00, 
                vtype = vtype)
        }
        if (is.element(measure, c("MD", "SMD", "SMDH", "ROM", 
            "RPB", "RBIS", "D2OR"))) {
            mf.m1i <- mf[[match("m1i", names(mf))]]
            mf.m2i <- mf[[match("m2i", names(mf))]]
            mf.sd1i <- mf[[match("sd1i", names(mf))]]
            mf.sd2i <- mf[[match("sd2i", names(mf))]]
            mf.n1i <- mf[[match("n1i", names(mf))]]
            mf.n2i <- mf[[match("n2i", names(mf))]]
            m1i <- eval(mf.m1i, data, enclos = sys.frame(sys.parent()))
            m2i <- eval(mf.m2i, data, enclos = sys.frame(sys.parent()))
            sd1i <- eval(mf.sd1i, data, enclos = sys.frame(sys.parent()))
            sd2i <- eval(mf.sd2i, data, enclos = sys.frame(sys.parent()))
            n1i <- eval(mf.n1i, data, enclos = sys.frame(sys.parent()))
            n2i <- eval(mf.n2i, data, enclos = sys.frame(sys.parent()))
            dat <- escalc(measure, m1i = m1i, m2i = m2i, sd1i = sd1i, 
                sd2i = sd2i, n1i = n1i, n2i = n2i, vtype = vtype)
        }
        if (is.element(measure, c("COR", "UCOR", "ZCOR"))) {
            mf.ri <- mf[[match("ri", names(mf))]]
            mf.ni <- mf[[match("ni", names(mf))]]
            ri <- eval(mf.ri, data, enclos = sys.frame(sys.parent()))
            ni <- eval(mf.ni, data, enclos = sys.frame(sys.parent()))
            dat <- escalc(measure, ri = ri, ni = ni, vtype = vtype)
        }
        if (is.element(measure, c("PR", "PLN", "PLO", "PAS", 
            "PFT"))) {
            mf.xi <- mf[[match("xi", names(mf))]]
            mf.mi <- mf[[match("mi", names(mf))]]
            mf.ni <- mf[[match("ni", names(mf))]]
            xi <- eval(mf.xi, data, enclos = sys.frame(sys.parent()))
            mi <- eval(mf.mi, data, enclos = sys.frame(sys.parent()))
            ni <- eval(mf.ni, data, enclos = sys.frame(sys.parent()))
            if (is.null(mi)) 
                mi <- ni - xi
            dat <- escalc(measure, xi = xi, mi = mi, add = add, 
                to = to, vtype = vtype)
        }
        if (is.element(measure, c("IR", "IRLN", "IRS", "IRFT"))) {
            mf.xi <- mf[[match("xi", names(mf))]]
            mf.ti <- mf[[match("ti", names(mf))]]
            xi <- eval(mf.xi, data, enclos = sys.frame(sys.parent()))
            ti <- eval(mf.ti, data, enclos = sys.frame(sys.parent()))
            dat <- escalc(measure, xi = xi, ti = ti, add = add, 
                to = to, vtype = vtype)
        }
        if (is.element(measure, c("MN"))) {
            mf.mi <- mf[[match("mi", names(mf))]]
            mf.sdi <- mf[[match("sdi", names(mf))]]
            mf.ni <- mf[[match("ni", names(mf))]]
            mi <- eval(mf.mi, data, enclos = sys.frame(sys.parent()))
            sdi <- eval(mf.sdi, data, enclos = sys.frame(sys.parent()))
            ni <- eval(mf.ni, data, enclos = sys.frame(sys.parent()))
            dat <- escalc(measure, mi = mi, sdi = sdi, ni = ni, 
                vtype = vtype)
        }
        if (is.element(measure, c("MC", "SMCC", "SMCR"))) {
            mf.m1i <- mf[[match("m1i", names(mf))]]
            mf.m2i <- mf[[match("m2i", names(mf))]]
            mf.sd1i <- mf[[match("sd1i", names(mf))]]
            mf.sd2i <- mf[[match("sd2i", names(mf))]]
            mf.ri <- mf[[match("ri", names(mf))]]
            mf.ni <- mf[[match("ni", names(mf))]]
            m1i <- eval(mf.m1i, data, enclos = sys.frame(sys.parent()))
            m2i <- eval(mf.m2i, data, enclos = sys.frame(sys.parent()))
            sd1i <- eval(mf.sd1i, data, enclos = sys.frame(sys.parent()))
            sd2i <- eval(mf.sd2i, data, enclos = sys.frame(sys.parent()))
            ri <- eval(mf.ri, data, enclos = sys.frame(sys.parent()))
            ni <- eval(mf.ni, data, enclos = sys.frame(sys.parent()))
            dat <- escalc(measure, m1i = m1i, m2i = m2i, sd1i = sd1i, 
                sd2i = sd2i, ri = ri, ni = ni, vtype = vtype)
        }
        if (is.element(measure, c("ARAW", "AHW", "ABT"))) {
            mf.ai <- mf[[match("ai", names(mf))]]
            mf.mi <- mf[[match("mi", names(mf))]]
            mf.ni <- mf[[match("ni", names(mf))]]
            ai <- eval(mf.ai, data, enclos = sys.frame(sys.parent()))
            mi <- eval(mf.mi, data, enclos = sys.frame(sys.parent()))
            ni <- eval(mf.ni, data, enclos = sys.frame(sys.parent()))
            dat <- escalc(measure, ai = ai, mi = mi, ni = ni, 
                vtype = vtype)
        }
        if (is.element(measure, c("GEN"))) 
            stop("Specify the desired outcome measure via the 'measure' argument.")
        yi <- dat$yi
        vi <- dat$vi
        ni <- attr(yi, "ni")
    }
    is.formula <- ifelse(class(mods) == "formula", TRUE, FALSE)
    if (is.formula) {
        options(na.action = "na.pass")
        mods <- model.matrix(mods, data = data)
        attr(mods, "assign") <- NULL
        options(na.action = na.act)
        intercept <- FALSE
    }
    if (is.vector(mods)) 
        mods <- cbind(mods)
    if (is.data.frame(mods)) 
        mods <- as.matrix(mods)
    if (!is.null(mods) && (nrow(mods) != length(yi))) 
        stop("Number of rows of the design matrix does not match length of yi.")
    k <- length(yi)
    ids <- seq.int(k)
    if (is.null(slab)) {
        if (!is.null(attr(yi, "slab"))) {
            slab.null <- FALSE
            slab <- attr(yi, "slab")
        }
        else {
            slab.null <- TRUE
            slab <- seq.int(k)
        }
    }
    else {
        if (any(duplicated(slab))) 
            stop("Study labels must be unique.")
        if (length(slab) != length(yi)) 
            stop("Study labels not of same length as data.")
        slab.null <- FALSE
    }
    if (!is.null(subset)) {
        yi <- yi[subset]
        vi <- vi[subset]
        ni <- ni[subset]
        ai <- ai[subset]
        bi <- bi[subset]
        ci <- ci[subset]
        di <- di[subset]
        x1i <- x1i[subset]
        x2i <- x2i[subset]
        t1i <- t1i[subset]
        t2i <- t2i[subset]
        mods <- mods[subset, , drop = FALSE]
        slab <- slab[subset]
        ids <- ids[subset]
        k <- length(yi)
    }
    ai.f <- ai
    bi.f <- bi
    ci.f <- ci
    di.f <- di
    x1i.f <- x1i
    x2i.f <- x2i
    t1i.f <- t1i
    t2i.f <- t2i
    yi.f <- yi
    vi.f <- vi
    ni.f <- ni
    mods.f <- mods
    k.f <- k
    YVM.na <- is.na(cbind(yi, vi, mods))
    if (any(YVM.na)) {
        not.na <- apply(YVM.na, MARGIN = 1, sum) == 0L
        if (na.act == "na.omit" || na.act == "na.exclude" || 
            na.act == "na.pass") {
            yi <- yi[not.na]
            vi <- vi[not.na]
            ni <- ni[not.na]
            mods <- mods[not.na, , drop = FALSE]
            k <- length(yi)
            warning("Studies with NAs omitted from model fitting.")
        }
        if (na.act == "na.fail") 
            stop("Missing values in data.")
    }
    else {
        not.na <- rep(TRUE, k)
    }
    if (k < 1) 
        stop("Processing terminated since k = 0.")
    if (k == 1) {
        method <- "FE"
        knha <- FALSE
    }
    if (any(vi <= 0)) {
        allvipos <- FALSE
        warning("There are outcomes with non-positive sampling variances.")
        if (any(vi < 0)) {
            vi[vi <= 0] <- 0
            warning("Negative sampling variances constrained to zero.")
        }
    }
    else {
        allvipos <- TRUE
    }
    if (is.null(mods) && !intercept) {
        warning("Must either include an intercept (intercept=TRUE) and/or moderators in model.\n  Coerced intercept into the model.")
        intercept <- TRUE
    }
    if (intercept) {
        X <- cbind(intrcpt = rep(1, k), mods)
        X.f <- cbind(intrcpt = rep(1, k.f), mods.f)
    }
    else {
        X <- mods
        X.f <- mods.f
    }
    is.int <- apply(X == 1L, 2, sum) == k
    if (any(is.int)) {
        int.incl <- TRUE
        int.indx <- which(is.int, arr.ind = TRUE)
        X <- cbind(intrcpt = 1, X[, -int.indx, drop = FALSE])
        X.f <- cbind(intrcpt = 1, X.f[, -int.indx, drop = FALSE])
        if (is.formula) 
            intercept <- TRUE
    }
    else {
        int.incl <- FALSE
    }
    p <- NCOL(X)
    if ((p == 1L) && (all(sapply(X, identical, 1)))) {
        int.only <- TRUE
    }
    else {
        int.only <- FALSE
    }
    evstXX <- eigen(t(X) %*% X, symmetric = TRUE)$values
    if (any(evstXX <= .Machine$double.eps)) 
        stop("Design matrix not of full rank. Cannot fit model.")
    if (method == "FE") {
        if (p > k) {
            stop("Number of parameters to be estimated is larger than the number of observations.")
        }
    }
    else {
        if (is.numeric(tau2)) {
            if (p > k) {
                stop("Number of parameters to be estimated is larger than the number of observations.")
            }
        }
        else {
            if ((p + 1) > k) {
                stop("Number of parameters to be estimated is larger than the number of observations.")
            }
        }
    }
    if (is.null(btt)) {
        if (p > 1) {
            if (int.incl) {
                btt <- seq.int(2, p)
            }
            else {
                btt <- seq.int(p)
            }
        }
        else {
            btt <- 1
        }
    }
    else {
        btt <- btt[(btt >= 1) & (btt <= p)]
        btt <- unique(round(btt))
        if (length(intersect(btt, seq.int(p))) == 0L) {
            stop("Non-existent coefficients specified with 'btt'.")
        }
    }
    bntt <- setdiff(seq.int(p), btt)
    m <- length(btt)
    con <- list(tau2.init = NULL, tau2.min = 0, tau2.max = 50, 
        threshold = 10^-5, maxiter = 100, stepadj = 1, REML2 = FALSE, 
        verbose = FALSE)
    con[pmatch(names(control), names(con))] <- control
    if (verbose) 
        con$verbose <- verbose
    iter <- 0
    se.tau2 <- I2 <- H2 <- QE <- QEp <- NA
    s2w <- 1
    alpha <- (100 - level)/100
    robust <- FALSE
    Y <- as.matrix(yi)
    if (!is.numeric(tau2)) {
        if (method == "HS") {
            if (!allvipos) 
                stop("HS estimator cannot be used with non-positive sampling variances.")
            wi <- 1/vi
            W <- .diag(wi)
            stXWX <- .invcalc(X = X, W = W, k = k)
            P <- W - W %*% X %*% stXWX %*% crossprod(X, W)
            RSS <- crossprod(Y, P) %*% Y
            tau2 <- (RSS - k)/sum(wi)
            se.tau2 <- sqrt(1/sum(wi)^2 * (2 * (k - p) + 4 * 
                max(tau2, 0) * .tr(P) + 2 * max(tau2, 0)^2 * 
                .tr(P %*% P)))
        }
        if (is.element(method, c("HE", "ML", "REML", "EB"))) {
            stXX <- .invcalc(X = X, W = diag(k), k = k)
            P <- diag(k) - X %*% tcrossprod(stXX, X)
            RSS <- crossprod(Y, P) %*% Y
            trPV <- .tr(P %*% .diag(vi))
            tau2 <- (RSS - trPV)/(k - p)
            se.tau2 <- sqrt(1/(k - p)^2 * (2 * .tr(P %*% .diag(vi) %*% 
                P %*% .diag(vi)) + 4 * max(tau2, 0) * trPV + 
                2 * max(tau2, 0)^2 * (k - p)))
        }
        if (method == "DL") {
            if (!allvipos) 
                stop("DL estimator cannot be used with non-positive sampling variances.")
            wi <- 1/vi
            W <- .diag(wi)
            stXWX <- .invcalc(X = X, W = W, k = k)
            P <- W - W %*% X %*% stXWX %*% crossprod(X, W)
            RSS <- crossprod(Y, P) %*% Y
            trP <- .tr(P)
            tau2 <- (RSS - (k - p))/trP
            se.tau2 <- sqrt(1/trP^2 * (2 * (k - p) + 4 * max(tau2, 
                0) * trP + 2 * max(tau2, 0)^2 * .tr(P %*% P)))
        }
        if (method == "SJ") {
            if (is.null(con$tau2.init)) {
                tau2.0 <- var(yi) * (k - 1)/k
            }
            else {
                tau2.0 <- con$tau2.init
            }
            wi <- 1/(vi + tau2.0)
            W <- .diag(wi)
            stXWX <- .invcalc(X = X, W = W, k = k)
            P <- W - W %*% X %*% stXWX %*% crossprod(X, W)
            RSS <- crossprod(Y, P) %*% Y
            tau2 <- tau2.0 * RSS/(k - p)
            se.tau2 <- sqrt(tau2.0^2/(k - p)^2 * (2 * .tr(P %*% 
                .diag(vi) %*% P %*% .diag(vi)) + 4 * max(tau2, 
                0) * .tr(P %*% .diag(vi) %*% P) + 2 * max(tau2, 
                0)^2 * .tr(P %*% P)))
        }
        if (method == "DLit") {
            conv <- 1
            change <- con$threshold + 1
            if (is.null(con$tau2.init)) {
                tau2 <- 0
            }
            else {
                tau2 <- con$tau2.init
            }
            while (change > con$threshold) {
                if (con$verbose) 
                  cat("Iteration", iter, "\ttau^2 =", round(tau2, 
                    digits), "\n")
                iter <- iter + 1
                tau2.old <- tau2
                wi <- 1/(vi + tau2)
                W <- .diag(wi)
                stXWX <- .invcalc(X = X, W = W, k = k)
                P <- W - W %*% X %*% stXWX %*% crossprod(X, W)
                RSS <- crossprod(Y, P) %*% Y
                trP <- .tr(P)
                tau2 <- (RSS - (k - p))/trP
                tau2[tau2 < con$tau2.min] <- con$tau2.min
                change <- abs(tau2.old - tau2)
                if (iter > con$maxiter) {
                  conv <- 0
                  break
                }
            }
            if (conv == 0L) 
                stop("Algorithm did not converge.")
            se.tau2 <- sqrt(1/trP^2 * (2 * (k - p) + 4 * max(tau2, 
                0) * trP + 2 * max(tau2, 0)^2 * .tr(P %*% P)))
        }
        if (method == "SJit") {
            conv <- 1
            change <- con$threshold + 1
            if (is.null(con$tau2.init)) {
                tau2 <- var(yi) * (k - 1)/k
            }
            else {
                tau2 <- con$tau2.init
            }
            while (change > con$threshold) {
                if (con$verbose) 
                  cat("Iteration", iter, "\ttau^2 =", round(tau2, 
                    digits), "\n")
                iter <- iter + 1
                tau2.old <- tau2
                wi <- 1/(vi + tau2)
                W <- .diag(wi)
                stXWX <- .invcalc(X = X, W = W, k = k)
                P <- W - W %*% X %*% stXWX %*% crossprod(X, W)
                RSS <- crossprod(Y, P) %*% Y
                tau2 <- tau2 * RSS/(k - p)
                change <- abs(tau2.old - tau2)
                if (iter > con$maxiter) {
                  conv <- 0
                  break
                }
            }
            if (conv == 0L) 
                stop("Algorithm did not converge.")
        }
        if (method == "PM") {
            if (!allvipos) 
                stop("PM estimator cannot be used with non-positive sampling variances.")
            if (.QE.func(con$tau2.min, Y = Y, vi = vi, X = X, 
                k = k, objective = k - p) < 0) {
                tau2 <- con$tau2.min
            }
            else {
                tau2 <- try(uniroot(.QE.func, interval = c(con$tau2.min, 
                  con$tau2.max), tol = con$threshold, Y = Y, 
                  vi = vi, X = X, k = k, objective = k - p)$root, 
                  silent = TRUE)
                if (!is.numeric(tau2)) {
                  stop("Error in iterative search for tau2. Try increasing tau2.max or switch to another 'method'.")
                }
            }
        }
        if (is.element(method, c("ML", "REML", "EB"))) {
            conv <- 1
            change <- con$threshold + 1
            if (is.null(con$tau2.init)) {
                tau2 <- max(0, tau2)
            }
            else {
                tau2 <- con$tau2.init
            }
            while (change > con$threshold) {
                if (con$verbose) 
                  cat("Iteration", iter, "\ttau^2 =", round(tau2, 
                    digits), "\n")
                iter <- iter + 1
                tau2.old <- tau2
                wi <- 1/(vi + tau2)
                W <- .diag(wi)
                stXWX <- .invcalc(X = X, W = W, k = k)
                P <- W - W %*% X %*% stXWX %*% crossprod(X, W)
                if (method == "ML") {
                  PP <- P %*% P
                  adj <- (crossprod(Y, PP) %*% Y - sum(wi))/sum(wi^2)
                }
                if (method == "REML") {
                  PP <- P %*% P
                  adj <- (crossprod(Y, PP) %*% Y - .tr(P))/.tr(PP)
                }
                if (method == "EB") {
                  adj <- (crossprod(Y, P) %*% Y * k/(k - p) - 
                    k)/sum(wi)
                }
                adj <- adj * con$stepadj
                while (tau2 + adj < con$tau2.min) {
                  adj <- adj/2
                }
                tau2 <- tau2 + adj
                change <- abs(tau2.old - tau2)
                if (iter > con$maxiter) {
                  conv <- 0
                  break
                }
            }
            if (conv == 0L) 
                stop("Fisher scoring algorithm did not converge. Try increasing the number of iterations,\n  adjusting the threshold, or use a different estimator for tau^2.")
            if (method == "ML") 
                se.tau2 <- sqrt(2/sum(wi^2))
            if (method == "REML") 
                se.tau2 <- sqrt(2/.tr(PP))
            if (method == "EB") 
                se.tau2 <- sqrt((k/(k - p))^2/sum(wi)^2 * (2 * 
                  .tr(P %*% .diag(vi) %*% P %*% .diag(vi)) + 
                  4 * max(tau2, 0) * .tr(P %*% .diag(vi) %*% 
                    P) + 2 * max(tau2, 0)^2 * .tr(P %*% P)))
        }
        tau2 <- max(con$tau2.min, c(tau2))
        if (con$verbose && is.element(method, c("ML", "REML", 
            "EB"))) 
            cat("Fisher scoring algorithm converged after", iter, 
                "iterations.", "\n")
    }
    else {
        if (method == "HS") {
            wi <- 1/vi
            W <- .diag(wi)
            stXWX <- .invcalc(X = X, W = W, k = k)
            P <- W - W %*% X %*% stXWX %*% crossprod(X, W)
            se.tau2 <- sqrt(1/sum(wi)^2 * (2 * (k - p) + 4 * 
                max(tau2, 0) * .tr(P) + 2 * max(tau2, 0)^2 * 
                .tr(P %*% P)))
        }
        if (method == "HE") {
            stXX <- .invcalc(X = X, W = diag(k), k = k)
            P <- diag(k) - X %*% tcrossprod(stXX, X)
            RSS <- crossprod(Y, P) %*% Y
            trPV <- .tr(P %*% .diag(vi))
            se.tau2 <- sqrt(1/(k - p)^2 * (2 * .tr(P %*% .diag(vi) %*% 
                P %*% .diag(vi)) + 4 * max(tau2, 0) * trPV + 
                2 * max(tau2, 0)^2 * (k - p)))
        }
        if (method == "DL") {
            wi <- 1/vi
            W <- .diag(wi)
            stXWX <- .invcalc(X = X, W = W, k = k)
            P <- W - W %*% X %*% stXWX %*% crossprod(X, W)
            trP <- .tr(P)
            se.tau2 <- sqrt(1/trP^2 * (2 * (k - p) + 4 * max(tau2, 
                0) * trP + 2 * max(tau2, 0)^2 * .tr(P %*% P)))
        }
        if (method == "SJ") {
            if (is.null(con$tau2.init)) {
                tau2.0 <- var(yi) * (k - 1)/k
            }
            else {
                tau2.0 <- con$tau2.init
            }
            wi <- 1/(vi + tau2.0)
            W <- .diag(wi)
            stXWX <- .invcalc(X = X, W = W, k = k)
            P <- W - W %*% X %*% stXWX %*% crossprod(X, W)
            se.tau2 <- sqrt(tau2.0^2/(k - p)^2 * (2 * .tr(P %*% 
                .diag(vi) %*% P %*% .diag(vi)) + 4 * max(tau2, 
                0) * .tr(P %*% .diag(vi) %*% P) + 2 * max(tau2, 
                0)^2 * .tr(P %*% P)))
        }
        if (method == "ML") {
            wi <- 1/(vi + tau2)
            se.tau2 <- sqrt(2/sum(wi^2))
        }
        if (method == "REML") {
            wi <- 1/(vi + tau2)
            W <- .diag(wi)
            stXWX <- .invcalc(X = X, W = W, k = k)
            P <- W - W %*% X %*% stXWX %*% crossprod(X, W)
            se.tau2 <- sqrt(2/.tr(P %*% P))
        }
        if (method == "EB") {
            wi <- 1/(vi + tau2)
            W <- .diag(wi)
            stXWX <- .invcalc(X = X, W = W, k = k)
            P <- W - W %*% X %*% stXWX %*% crossprod(X, W)
            se.tau2 <- sqrt((k/(k - p))^2/sum(wi)^2 * (2 * .tr(P %*% 
                .diag(vi) %*% P %*% .diag(vi)) + 4 * tau2 * .tr(P %*% 
                .diag(vi) %*% P) + 2 * tau2^2 * .tr(P %*% P)))
        }
    }
    if (method == "FE") {
        tau2 <- 0
        if (!allvipos && weighted) 
            stop("Weighted estimation cannot be used with a fixed-effects\n  model when there are non-positive sampling variances.")
    }
    if (allvipos) {
        wi <- 1/vi
        W <- .diag(wi)
        stXWX <- .invcalc(X = X, W = W, k = k)
        P <- W - W %*% X %*% stXWX %*% crossprod(X, W)
        QE <- max(0, c(crossprod(Y, P) %*% Y))
        QEp <- ifelse(k - p >= 1, pchisq(QE, df = k - p, lower.tail = FALSE), 
            1)
        vi.avg <- (k - p)/.tr(P)
        I2 <- 100 * tau2/(vi.avg + tau2)
        H2 <- tau2/vi.avg + 1
    }
    wi <- 1/(vi + tau2)
    W <- .diag(wi)
    if (weighted) {
        stXWX <- .invcalc(X = X, W = W, k = k)
        b <- stXWX %*% crossprod(X, W) %*% Y
        vb <- stXWX
        P <- W - W %*% X %*% stXWX %*% crossprod(X, W)
        RSS.f <- crossprod(Y, P) %*% Y
        if (robust) {
            resid <- c(Y - X %*% b)
            vb <- vb %*% t(X) %*% W %*% diag(resid^2) %*% W %*% 
                X %*% vb
            vb <- vb * k/(k - p)
        }
        if (knha) {
            if (RSS.f <= .Machine$double.eps) {
                s2w <- 1
            }
            else {
                s2w <- c(RSS.f)/(k - p)
            }
            vb <- s2w * vb
            if (method == "FE") 
                warning("Knapp & Hartung (2003) method is not meant to be used in the context of fixed-effects models.")
        }
        if (robust) {
            QM <- c(t(b)[btt] %*% chol2inv(chol(vb[btt, btt])) %*% 
                b[btt])
        }
        else {
            if (length(bntt) == 0L) {
                QM <- c(sum(wi * yi^2) - RSS.f)/s2w
            }
            else {
                Xr <- X[, bntt, drop = FALSE]
                stXWX <- .invcalc(X = Xr, W = W, k = k)
                P <- W - W %*% Xr %*% stXWX %*% crossprod(Xr, 
                  W)
                RSS.r <- crossprod(Y, P) %*% Y
                QM <- c(RSS.r - RSS.f)/s2w
            }
        }
    }
    else {
        stXX <- .invcalc(X = X, W = diag(k), k = k)
        b <- stXX %*% crossprod(X, Y)
        vb <- tcrossprod(stXX, X) %*% .diag(vi + tau2) %*% X %*% 
            stXX
        P <- W - W %*% X %*% tcrossprod(stXX, X) - X %*% stXX %*% 
            crossprod(X, W) + X %*% stXX %*% crossprod(X, W) %*% 
            X %*% tcrossprod(stXX, X)
        RSS.f <- crossprod(Y, P) %*% Y
        if (robust) {
            resid <- c(Y - X %*% b)
            vb <- stXX %*% t(X) %*% diag(resid^2) %*% X %*% stXX
            vb <- vb * k/(k - p)
        }
        if (knha) {
            stXWX <- .invcalc(X = X, W = W, k = k)
            P <- W - W %*% X %*% stXWX %*% crossprod(X, W)
            RSS.k <- c(crossprod(Y, P) %*% Y)
            if (RSS.k <= .Machine$double.eps) {
                s2w <- 1
            }
            else {
                s2w <- c(RSS.k)/(k - p)
            }
            vb <- s2w * vb
            if (method == "FE") 
                warning("Knapp & Hartung (2003) method is not meant to be used in the context of fixed-effects models.")
        }
        QM <- c(t(b)[btt] %*% chol2inv(chol(vb[btt, btt])) %*% 
            b[btt])
    }
    se <- sqrt(diag(vb))
    names(se) <- NULL
    zval <- c(b/se)
    if (knha || robust) {
        QM <- QM/m
        QMp <- pf(QM, df1 = m, df2 = k - p, lower.tail = FALSE)
        pval <- 2 * pt(abs(zval), df = k - p, lower.tail = FALSE)
        crit <- qt(alpha/2, df = k - p, lower.tail = FALSE)
    }
    else {
        QMp <- pchisq(QM, df = m, lower.tail = FALSE)
        pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)
        crit <- qnorm(alpha/2, lower.tail = FALSE)
    }
    ci.lb <- c(b - crit * se)
    ci.ub <- c(b + crit * se)
    parms <- p + ifelse(method == "FE", 0, 1)
    ll.ML <- -1/2 * (k) * log(2 * base::pi) - 1/2 * sum(log(vi + 
        tau2)) - 1/2 * RSS.f
    ll.REML <- -1/2 * (k - p) * log(2 * base::pi) + ifelse(con$REML2, 
        0, 1/2 * determinant(crossprod(X, X), logarithm = TRUE)$modulus) - 
        1/2 * sum(log(vi + tau2)) - 1/2 * determinant(crossprod(X, 
        W) %*% X, logarithm = TRUE)$modulus - 1/2 * RSS.f
    if (k - p > 0L) {
        dev.ML <- -2 * (ll.ML - sum(dnorm(yi, mean = yi, sd = sqrt(vi), 
            log = TRUE)))
    }
    else {
        dev.ML <- 0
    }
    AIC.ML <- -2 * ll.ML + 2 * parms
    BIC.ML <- -2 * ll.ML + parms * log(k)
    dev.REML <- -2 * (ll.REML - 0)
    AIC.REML <- -2 * ll.REML + 2 * parms
    BIC.REML <- -2 * ll.REML + parms * log(k - p)
    fit.stats <- matrix(c(ll.ML, dev.ML, AIC.ML, BIC.ML, ll.REML, 
        dev.REML, AIC.REML, BIC.REML), ncol = 2, byrow = FALSE)
    dimnames(fit.stats) <- list(c("ll", "dev", "AIC", "BIC"), 
        c("ML", "REML"))
    fit.stats <- data.frame(fit.stats)
    p.eff <- p
    k.eff <- k
    res <- list(b, se, zval, pval, ci.lb, ci.ub, vb, tau2, se.tau2, 
        k, k.f, k.eff, p, p.eff, parms, m, QE, QEp, QM, QMp, 
        I2, H2, int.only, int.incl, allvipos, yi, vi, X, yi.f, 
        vi.f, X.f, ai.f, bi.f, ci.f, di.f, x1i.f, x2i.f, t1i.f, 
        t2i.f, ni, ni.f, ids, not.na, slab, slab.null, measure, 
        method, weighted, knha, robust, s2w, btt, intercept, 
        digits, level, con, add, to, drop00, fit.stats)
    names(res) <- c("b", "se", "zval", "pval", "ci.lb", "ci.ub", 
        "vb", "tau2", "se.tau2", "k", "k.f", "k.eff", "p", "p.eff", 
        "parms", "m", "QE", "QEp", "QM", "QMp", "I2", "H2", "int.only", 
        "int.incl", "allvipos", "yi", "vi", "X", "yi.f", "vi.f", 
        "X.f", "ai.f", "bi.f", "ci.f", "di.f", "x1i.f", "x2i.f", 
        "t1i.f", "t2i.f", "ni", "ni.f", "ids", "not.na", "slab", 
        "slab.null", "measure", "method", "weighted", "knha", 
        "robust", "s2w", "btt", "intercept", "digits", "level", 
        "control", "add", "to", "drop00", "fit.stats")
    class(res) <- c("rma.uni", "rma")
    return(res)
}
rstandard.rma.mh <-
function (model, digits = model$digits, ...) 
{
    if (!is.element("rma.mh", class(model))) 
        stop("Argument 'model' must be an object of class \"rma.mh\".")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    x <- model
    ei <- c(x$yi.f - x$b)
    ei[abs(ei) < 100 * .Machine$double.eps] <- 0
    sei <- sqrt(x$vi.f)
    zi <- ei/sei
    if (na.act == "na.omit") {
        out <- list(resid = ei[x$not.na.yivi], se = sei[x$not.na.yivi], 
            z = zi[x$not.na.yivi])
        out$slab <- x$slab[x$not.na.yivi]
    }
    if (na.act == "na.exclude" || na.act == "na.pass") {
        out <- list(resid = ei, se = sei, z = zi)
        out$slab <- x$slab
    }
    if (na.act == "na.fail" && any(!x$not.na.yivi)) 
        stop("Missing values in results.")
    out$digits <- digits
    class(out) <- c("list.rma")
    return(out)
}
rstandard.rma.peto <-
function (model, digits = model$digits, ...) 
{
    if (!is.element("rma.peto", class(model))) 
        stop("Argument 'model' must be an object of class \"rma.peto\".")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    x <- model
    ei <- c(x$yi.f - x$b)
    ei[abs(ei) < 100 * .Machine$double.eps] <- 0
    sei <- sqrt(x$vi.f)
    zi <- ei/sei
    if (na.act == "na.omit") {
        out <- list(resid = ei[x$not.na.yivi], se = sei[x$not.na.yivi], 
            z = zi[x$not.na.yivi])
        out$slab <- x$slab[x$not.na.yivi]
    }
    if (na.act == "na.exclude" || na.act == "na.pass") {
        out <- list(resid = ei, se = sei, z = zi)
        out$slab <- x$slab
    }
    if (na.act == "na.fail" && any(!x$not.na.yivi)) 
        stop("Missing values in results.")
    out$digits <- digits
    class(out) <- c("list.rma")
    return(out)
}
rstandard.rma.uni <-
function (model, digits = model$digits, ...) 
{
    if (!is.element("rma.uni", class(model))) 
        stop("Argument 'model' must be an object of class \"rma.uni\".")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    x <- model
    V <- .diag(x$vi + x$tau2)
    if (x$weighted) {
        wi <- 1/(x$vi + x$tau2)
        W <- .diag(wi)
        stXWX <- .invcalc(X = x$X, W = W, k = x$k)
        H <- x$X %*% stXWX %*% crossprod(x$X, W)
    }
    else {
        stXX <- .invcalc(X = x$X, W = diag(x$k), k = x$k)
        H <- x$X %*% tcrossprod(stXX, x$X)
    }
    ImH <- diag(x$k) - H
    ei <- ImH %*% cbind(x$yi)
    ei[abs(ei) < 100 * .Machine$double.eps] <- 0
    ve <- ImH %*% tcrossprod(V, ImH)
    sei <- sqrt(diag(ve))
    resid <- rep(NA, x$k.f)
    seresid <- rep(NA, x$k.f)
    stanres <- rep(NA, x$k.f)
    resid[x$not.na] <- ei
    seresid[x$not.na] <- sei
    stanres[x$not.na] <- ei/sei
    if (na.act == "na.omit") {
        out <- list(resid = resid[x$not.na], se = seresid[x$not.na], 
            z = stanres[x$not.na])
        out$slab <- x$slab[x$not.na]
    }
    if (na.act == "na.exclude" || na.act == "na.pass") {
        out <- list(resid = resid, se = seresid, z = stanres)
        out$slab <- x$slab
    }
    if (na.act == "na.fail" && any(!x$not.na)) 
        stop("Missing values in results.")
    out$digits <- digits
    class(out) <- c("list.rma")
    return(out)
}
rstudent.rma.mh <-
function (model, digits = model$digits, ...) 
{
    if (!is.element("rma.mh", class(model))) 
        stop("Argument 'model' must be an object of class \"rma.mh\".")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    x <- model
    delpred <- rep(NA, x$k.f)
    vdelpred <- rep(NA, x$k.f)
    o.warn <- getOption("warn")
    on.exit(options(warn = o.warn))
    options(warn = -1)
    for (i in seq.int(x$k.f)[x$not.na]) {
        if (is.element(x$measure, c("RR", "OR", "RD"))) {
            res <- try(rma.mh(ai = x$ai.f[-i], bi = x$bi.f[-i], 
                ci = x$ci.f[-i], di = x$di.f[-i], measure = x$measure, 
                add = x$add, to = x$to, ...), silent = TRUE)
        }
        else {
            res <- try(rma.mh(x1i = x$x1i.f[-i], x2i = x$x2i.f[-i], 
                t1i = x$t1i.f[-i], t2i = x$t2i.f[-i], measure = x$measure, 
                add = x$add, to = x$to, ...), silent = TRUE)
        }
        if (inherits(res, "try-error")) 
            next
        delpred[i] <- res$b
        vdelpred[i] <- res$vb
    }
    delresid <- x$yi.f - delpred
    delresid[abs(delresid) < 100 * .Machine$double.eps] <- 0
    sedelresid <- sqrt(x$vi.f + vdelpred)
    standelres <- delresid/sedelresid
    if (na.act == "na.omit") {
        out <- list(resid = delresid[x$not.na.yivi], se = sedelresid[x$not.na.yivi], 
            z = standelres[x$not.na.yivi])
        out$slab <- x$slab[x$not.na.yivi]
    }
    if (na.act == "na.exclude" || na.act == "na.pass") {
        out <- list(resid = delresid, se = sedelresid, z = standelres)
        out$slab <- x$slab
    }
    if (na.act == "na.fail" && any(!x$not.na.yivi)) 
        stop("Missing values in results.")
    out$digits <- digits
    class(out) <- c("list.rma")
    return(out)
}
rstudent.rma.peto <-
function (model, digits = model$digits, ...) 
{
    if (!is.element("rma.peto", class(model))) 
        stop("Argument 'model' must be an object of class \"rma.peto\".")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    x <- model
    delpred <- rep(NA, x$k.f)
    vdelpred <- rep(NA, x$k.f)
    o.warn <- getOption("warn")
    on.exit(options(warn = o.warn))
    options(warn = -1)
    for (i in seq.int(x$k.f)[x$not.na]) {
        res <- try(rma.peto(ai = x$ai.f[-i], bi = x$bi.f[-i], 
            ci = x$ci.f[-i], di = x$di.f[-i], add = x$add, to = x$to, 
            ...), silent = TRUE)
        if (inherits(res, "try-error")) 
            next
        delpred[i] <- res$b
        vdelpred[i] <- res$vb
    }
    delresid <- x$yi.f - delpred
    delresid[abs(delresid) < 100 * .Machine$double.eps] <- 0
    sedelresid <- sqrt(x$vi.f + vdelpred)
    standelres <- delresid/sedelresid
    if (na.act == "na.omit") {
        out <- list(resid = delresid[x$not.na.yivi], se = sedelresid[x$not.na.yivi], 
            z = standelres[x$not.na.yivi])
        out$slab <- x$slab[x$not.na.yivi]
    }
    if (na.act == "na.exclude" || na.act == "na.pass") {
        out <- list(resid = delresid, se = sedelresid, z = standelres)
        out$slab <- x$slab
    }
    if (na.act == "na.fail" && any(!x$not.na.yivi)) 
        stop("Missing values in results.")
    out$digits <- digits
    class(out) <- c("list.rma")
    return(out)
}
rstudent.rma.uni <-
function (model, digits = model$digits, ...) 
{
    if (!is.element("rma.uni", class(model))) 
        stop("Argument 'model' must be an object of class \"rma.uni\".")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    x <- model
    tau2.del <- rep(NA, x$k.f)
    delpred <- rep(NA, x$k.f)
    vdelpred <- rep(NA, x$k.f)
    o.warn <- getOption("warn")
    on.exit(options(warn = o.warn))
    options(warn = -1)
    for (i in seq.int(x$k.f)[x$not.na]) {
        res <- try(rma(x$yi.f[-i], x$vi.f[-i], mods = cbind(x$X.f[-i, 
            ]), method = x$method, weighted = x$weighted, intercept = FALSE, 
            knha = x$knha, control = x$control, ...), silent = TRUE)
        if (inherits(res, "try-error")) 
            next
        tau2.del[i] <- res$tau2
        Xi <- matrix(x$X.f[i, ], nrow = 1)
        delpred[i] <- Xi %*% res$b
        vdelpred[i] <- Xi %*% tcrossprod(res$vb, Xi)
    }
    delresid <- x$yi.f - delpred
    delresid[abs(delresid) < 100 * .Machine$double.eps] <- 0
    sedelresid <- sqrt(x$vi.f + vdelpred + tau2.del)
    standelres <- delresid/sedelresid
    if (na.act == "na.omit") {
        out <- list(resid = delresid[x$not.na], se = sedelresid[x$not.na], 
            z = standelres[x$not.na])
        out$slab <- x$slab[x$not.na]
    }
    if (na.act == "na.exclude" || na.act == "na.pass") {
        out <- list(resid = delresid, se = sedelresid, z = standelres)
        out$slab <- x$slab
    }
    if (na.act == "na.fail" && any(!x$not.na)) 
        stop("Missing values in results.")
    out$digits <- digits
    class(out) <- c("list.rma")
    return(out)
}
summary.escalc <-
function (object, out.names = c("sei", "zi", "ci.lb", "ci.ub"), 
    var.names, append = TRUE, replace = TRUE, level = 95, digits, 
    transf = FALSE, ...) 
{
    if (!is.element("escalc", class(object))) 
        stop("Argument 'object' must be an object of class \"escalc\".")
    x <- object
    alpha <- (100 - level)/100
    crit <- qnorm(alpha/2, lower.tail = FALSE)
    if (length(out.names) != 4) 
        stop("Argument out.names must be of length 4.")
    if (any(out.names != make.names(out.names, unique = TRUE))) {
        out.names <- make.names(out.names, unique = TRUE)
        warning(paste("Argument 'out.names' does not contain syntactically valid variable names.\n  Variable names adjusted to: out.names = c('", 
            out.names[1], "', '", out.names[2], "', '", out.names[3], 
            "', '", out.names[4], "').", sep = ""))
    }
    if (missing(var.names)) {
        if (!is.null(attr(x, "var.names"))) {
            yi.name <- attr(x, "var.names")[1]
            vi.name <- attr(x, "var.names")[2]
        }
        else {
            if (!is.element("yi", names(x))) 
                stop("Cannot determine name of the 'yi' variable.")
            if (!is.element("vi", names(x))) 
                stop("Cannot determine name of the 'vi' variable.")
            yi.name <- "yi"
            vi.name <- "vi"
        }
    }
    else {
        if (length(var.names) != 2) 
            stop("Argument var.names must be of length 2.")
        if (any(var.names != make.names(var.names, unique = TRUE))) {
            var.names <- make.names(var.names, unique = TRUE)
            warning(paste("Argument 'var.names' does not contain syntactically valid variable names.\n  Variable names adjusted to: var.names = c('", 
                var.names[1], "', '", var.names[2], "').", sep = ""))
        }
        yi.name <- var.names[1]
        vi.name <- var.names[2]
    }
    yi <- x[[yi.name]]
    vi <- x[[vi.name]]
    if (is.null(yi) || is.null(vi)) 
        stop(paste("Cannot find variables '", yi.name, "' and/or '", 
            vi.name, "' in the data frame.", sep = ""))
    sei <- sqrt(vi)
    zi <- yi/sei
    if (is.function(transf)) {
        yi <- mapply(transf, yi, ...)
        ci.lb <- mapply(transf, yi - crit * sei, ...)
        ci.ub <- mapply(transf, yi + crit * sei, ...)
        attr(x, "transf") <- TRUE
    }
    else {
        ci.lb <- yi - crit * sei
        ci.ub <- yi + crit * sei
        attr(x, "transf") <- FALSE
    }
    x[[yi.name]] <- yi
    x[[vi.name]] <- vi
    if (append) {
        dat <- data.frame(x)
        if (replace) {
            dat[[out.names[1]]] <- sei
            dat[[out.names[2]]] <- zi
            dat[[out.names[3]]] <- ci.lb
            dat[[out.names[4]]] <- ci.ub
        }
        else {
            if (is.element(out.names[1], names(dat))) {
                is.na.sei <- is.na(dat[[out.names[1]]])
                dat[[out.names[1]]][is.na.sei] <- sei[is.na.sei]
            }
            else {
                dat[[out.names[1]]] <- sei
            }
            if (is.element(out.names[2], names(dat))) {
                is.na.zi <- is.na(dat[[out.names[2]]])
                dat[[out.names[2]]][is.na.zi] <- zi[is.na.zi]
            }
            else {
                dat[[out.names[2]]] <- zi
            }
            if (is.element(out.names[3], names(dat))) {
                is.na.ci.lb <- is.na(dat[[out.names[3]]])
                dat[[out.names[3]]][is.na.ci.lb] <- ci.lb[is.na.ci.lb]
            }
            else {
                dat[[out.names[3]]] <- ci.lb
            }
            if (is.element(out.names[4], names(dat))) {
                is.na.ci.ub <- is.na(dat[[out.names[4]]])
                dat[[out.names[4]]][is.na.ci.ub] <- ci.ub[is.na.ci.ub]
            }
            else {
                dat[[out.names[4]]] <- ci.ub
            }
        }
    }
    else {
        dat <- data.frame(yi, vi, sei, zi, ci.lb, ci.ub)
        names(dat) <- c(yi.name, vi.name, out.names)
    }
    if (!missing(var.names)) {
        attr(dat, "var.names") <- var.names
    }
    else {
        attr(dat, "var.names") <- c(yi.name, vi.name)
    }
    if (!missing(digits)) {
        attr(dat, "digits") <- digits
    }
    else {
        attr(dat, "digits") <- attr(x, "digits")
    }
    if (is.null(attr(dat, "digits"))) 
        attr(dat, "digits") <- 4
    attr(dat, "out.names") <- out.names
    class(dat) <- c("summary.escalc", "data.frame")
    return(dat)
}
summary.rma <-
function (object, digits = object$digits, showfit = TRUE, ...) 
{
    if (!is.element("rma", class(object))) 
        stop("Argument 'object' must be an object of class \"rma\".")
    class(object) <- c("summary.rma", class(object))
    return(object)
}
transf.abt <-
function (xi, ...) 
{
    zi <- -log(1 - xi)
    return(c(zi))
}
transf.ahw <-
function (xi, ...) 
{
    zi <- 1 - (1 - xi)^(1/3)
    return(c(zi))
}
transf.arcsin <-
function (xi, ...) 
{
    zi <- asin(sqrt(xi))
    return(c(zi))
}
transf.exp.int <-
function (xi, targs = NULL, ...) 
{
    if (is.null(targs$tau2)) 
        targs$tau2 <- 0
    if (is.null(targs$lower)) 
        targs$lower <- xi - 5 * sqrt(targs$tau2)
    if (is.null(targs$upper)) 
        targs$upper <- xi + 5 * sqrt(targs$tau2)
    toint <- function(zval, xi, tau2) {
        exp(zval) * dnorm(zval, mean = xi, sd = sqrt(tau2))
    }
    cfunc <- function(xi, tau2, lower, upper) {
        integrate(toint, lower = lower, upper = upper, xi = xi, 
            tau2 = tau2)$value
    }
    zi <- mapply(xi, FUN = cfunc, tau2 = targs$tau2, lower = targs$lower, 
        upper = targs$upper)
    return(c(zi))
}
transf.iabt <-
function (xi, ...) 
{
    zi <- 1 - exp(-xi)
    zi <- ifelse(is.nan(zi), NA, zi)
    zi[xi < 0] <- 0
    return(c(zi))
}
transf.iahw <-
function (xi, ...) 
{
    zi <- 1 - (1 - xi)^3
    zi <- ifelse(is.nan(zi), NA, zi)
    zi[xi > 1] <- 1
    zi[xi < 0] <- 0
    return(c(zi))
}
transf.iarcsin <-
function (xi, ...) 
{
    zi <- sin(xi)^2
    zi[xi < 0] <- 0
    zi[xi > asin(1)] <- 1
    return(c(zi))
}
transf.iirft <-
function (xi, ti, ...) 
{
    zi <- (1/ti - 8 * xi^2 + 16 * ti * xi^4)/(16 * xi^2 * ti)
    zi <- ifelse(is.nan(zi), NA, zi)
    zi[xi < transf.irft(0, ti)] <- 0
    zi[zi <= .Machine$double.eps] <- 0
    return(c(zi))
}
transf.ilogit <-
function (xi, ...) 
{
    zi <- exp(xi)/(1 + exp(xi))
    zi[xi == -Inf] <- 0
    zi[xi == Inf] <- 1
    zi[is.nan(zi) & (xi > 0.5)] <- 1
    return(c(zi))
}
transf.ilogit.int <-
function (xi, targs = NULL, ...) 
{
    if (is.null(targs$tau2)) 
        targs$tau2 <- 0
    if (is.null(targs$lower)) 
        targs$lower <- xi - 5 * sqrt(targs$tau2)
    if (is.null(targs$upper)) 
        targs$upper <- xi + 5 * sqrt(targs$tau2)
    toint <- function(zval, xi, tau2) {
        zi <- exp(zval)/(1 + exp(zval))
        zi[xi == -Inf] <- 0
        zi[xi == Inf] <- 1
        zi[is.nan(zi) & (xi > 0.5)] <- 1
        zi * dnorm(zval, mean = xi, sd = sqrt(tau2))
    }
    cfunc <- function(xi, tau2, lower, upper) {
        integrate(toint, lower = lower, upper = upper, xi = xi, 
            tau2 = tau2)$value
    }
    zi <- mapply(xi, FUN = cfunc, tau2 = targs$tau2, lower = targs$lower, 
        upper = targs$upper)
    return(c(zi))
}
transf.ipft <-
function (xi, ni, ...) 
{
    zi <- 1/2 * (1 - sign(cos(2 * xi)) * sqrt(1 - (sin(2 * xi) + 
        (sin(2 * xi) - 1/sin(2 * xi))/ni)^2))
    zi <- ifelse(is.nan(zi), NA, zi)
    zi[xi > transf.pft(1, ni)] <- 1
    zi[xi < transf.pft(0, ni)] <- 0
    return(c(zi))
}
transf.ipft.hm <-
function (xi, targs, ...) 
{
    ni <- 1/(mean(1/targs$ni, na.rm = TRUE))
    zi <- 1/2 * (1 - sign(cos(2 * xi)) * sqrt(1 - (sin(2 * xi) + 
        (sin(2 * xi) - 1/sin(2 * xi))/ni)^2))
    zi <- ifelse(is.nan(zi), NA, zi)
    zi[xi > transf.pft(1, ni)] <- 1
    zi[xi < transf.pft(0, ni)] <- 0
    return(c(zi))
}
transf.irft <-
function (xi, ti, ...) 
{
    zi <- 1/2 * (sqrt(xi) + sqrt(xi + 1/ti))
    return(c(zi))
}
transf.isqrt <-
function (xi, ...) 
{
    zi <- xi * xi
    zi[xi < 0] <- 0
    return(c(zi))
}
transf.logit <-
function (xi, ...) 
{
    zi <- log(xi/(1 - xi))
    return(c(zi))
}
transf.pft <-
function (xi, ni, ...) 
{
    xi <- xi * ni
    zi <- 1/2 * (asin(sqrt(xi/(ni + 1))) + asin(sqrt((xi + 1)/(ni + 
        1))))
    return(c(zi))
}
transf.rtoz <-
function (xi, ...) 
{
    zi <- 1/2 * log((1 + xi)/(1 - xi))
    return(c(zi))
}
transf.ztor <-
function (xi, ...) 
{
    zi <- (exp(2 * xi) - 1)/(exp(2 * xi) + 1)
    zi[xi == -Inf] <- -1
    zi[xi == Inf] <- 1
    zi[is.nan(zi) & (xi > 0)] <- 1
    return(c(zi))
}
transf.ztor.int <-
function (xi, targs = NULL, ...) 
{
    if (is.null(targs$tau2)) 
        targs$tau2 <- 0
    if (is.null(targs$lower)) 
        targs$lower <- xi - 5 * sqrt(targs$tau2)
    if (is.null(targs$upper)) 
        targs$upper <- xi + 5 * sqrt(targs$tau2)
    toint <- function(zval, xi, tau2) {
        zi <- (exp(2 * zval) - 1)/(exp(2 * zval) + 1)
        zi[xi == -Inf] <- -1
        zi[xi == Inf] <- 1
        zi[is.nan(zi) & (xi > 0)] <- 1
        zi * dnorm(zval, mean = xi, sd = sqrt(tau2))
    }
    cfunc <- function(xi, tau2, lower, upper) {
        integrate(toint, lower = lower, upper = upper, xi = xi, 
            tau2 = tau2)$value
    }
    zi <- mapply(xi, FUN = cfunc, tau2 = targs$tau2, lower = targs$lower, 
        upper = targs$upper)
    return(c(zi))
}
trimfill <-
function (x, ...) 
UseMethod("trimfill")
trimfill.rma.uni <-
function (x, side, estimator = "L0", maxiter = 50, verbose = FALSE, 
    ...) 
{
    if (!is.element("rma.uni", class(x))) 
        stop("Argument 'x' must be an object of class \"rma.uni\".")
    if (!x$int.only) 
        stop("Trim-and-fill method only applicable for models without moderators.")
    if (missing(side)) 
        side <- NULL
    estimator <- match.arg(estimator, c("L0", "R0"))
    if (x$k == 1) 
        stop("Stopped because k = 1.")
    yi <- x$yi
    vi <- x$vi
    ni <- x$ni
    if (is.null(side)) {
        res <- rma(yi, vi, mods = sqrt(vi), intercept = TRUE, 
            method = x$method, weighted = x$weighted, ...)
        if (res$b[2] < 0) {
            side <- "right"
        }
        else {
            side <- "left"
        }
    }
    else {
        side <- match.arg(side, c("left", "right"))
    }
    if (side == "right") {
        yi <- -1 * yi
    }
    idix <- sort(yi, index.return = TRUE)$ix
    yi <- yi[idix]
    vi <- vi[idix]
    ni <- ni[idix]
    k <- length(yi)
    k0.sav <- -1
    k0 <- 0
    iter <- 0
    while (abs(k0 - k0.sav) > 0) {
        k0.sav <- k0
        iter <- iter + 1
        if (iter > maxiter) 
            stop("Trim and fill algorithm did not converge.")
        yi.t <- yi[1:(k - k0)]
        vi.t <- vi[1:(k - k0)]
        res <- rma(yi.t, vi.t, intercept = TRUE, method = x$method, 
            weighted = x$weighted, ...)
        b <- c(res$b)
        yi.c <- yi - b
        yi.c.r <- rank(abs(yi.c), ties.method = "first")
        yi.c.r.s <- sign(yi.c) * yi.c.r
        if (estimator == "L0") {
            Sr <- sum(yi.c.r.s[yi.c.r.s > 0])
            k0 <- round((4 * Sr - k * (k + 1))/(2 * k - 1))
        }
        if (estimator == "R0") {
            k0 <- (k - max(-1 * yi.c.r.s[yi.c.r.s < 0])) - 1
        }
        k0 <- max(0, k0)
        if (verbose) 
            cat("Iteration:", iter, "\tmissing =", k0, "\t  b =", 
                ifelse(side == "right", -1 * b, b), "\n")
    }
    if (k0 > 0) {
        if (side == "right") {
            yi.c <- -1 * (yi.c - b)
        }
        else {
            yi.c <- yi.c - b
        }
        yi.fill <- c(x$yi.f, -1 * yi.c[(k - k0 + 1):k])
        vi.fill <- c(x$vi.f, vi[(k - k0 + 1):k])
        ni.fill <- c(x$ni.f, ni[(k - k0 + 1):k])
        attr(yi.fill, "measure") <- x$measure
        res <- rma(yi.fill, vi.fill, ni = ni.fill, intercept = TRUE, 
            method = x$method, weighted = x$weighted, ...)
        res$fill <- c(rep(0, k), rep(1, k0))
        res$k0 <- k0
        res$side <- side
        class(res) <- c("rma.uni.trimfill", class(res))
        res$ids <- c(x$ids, (x$k.f + 1):(x$k.f + k0))
        if (x$slab.null) {
            res$slab <- c(paste("Study", x$ids), paste("Filled", 
                seq.int(k0)))
            res$slab.null <- FALSE
        }
        else {
            res$slab <- c(x$slab, paste("Filled", seq.int(k0)))
            res$slab.null <- FALSE
        }
        return(res)
    }
    else {
        res <- x
        res$fill <- c(rep(0, k), rep(1, k0))
        res$k0 <- k0
        res$side <- side
        class(res) <- c("rma.uni.trimfill", class(res))
        return(res)
    }
}
vcov.rma <-
function (object, ...) 
{
    if (!is.element("rma", class(object))) 
        stop("Argument 'object' must be an object of class \"rma\".")
    return(object$vb)
}
weights.rma.mh <-
function (object, ...) 
{
    if (!is.element("rma.mh", class(object))) 
        stop("Argument 'object' must be an object of class \"rma.mh\".")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    x <- object
    weight <- rep(NA, x$k.f)
    if (is.element(x$measure, c("RR", "OR", "RD"))) {
        Ni <- x$ai + x$bi + x$ci + x$di
    }
    else {
        Ti <- x$t1i + x$t2i
    }
    if (x$measure == "OR") 
        wi <- x$bi * x$ci/Ni
    if (x$measure == "RR") 
        wi <- x$ci * (x$ai + x$bi)/Ni
    if (x$measure == "RD") 
        wi <- (x$ai + x$bi) * (x$ci + x$di)/Ni
    if (x$measure == "IRR") 
        wi <- x$x2i * (x$t1i)/Ti
    weight[x$not.na] <- wi/sum(wi) * 100
    names(weight) <- x$slab
    if (na.act == "na.omit") 
        weight <- weight[x$not.na]
    if (na.act == "na.fail" && any(!x$not.na)) 
        stop("Missing values in weights.")
    return(weight)
}
weights.rma.peto <-
function (object, ...) 
{
    if (!is.element("rma.peto", class(object))) 
        stop("Argument 'object' must be an object of class \"rma.peto\".")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    x <- object
    weight <- rep(NA, x$k.f)
    n1i <- x$ai + x$bi
    n2i <- x$ci + x$di
    Ni <- x$ai + x$bi + x$ci + x$di
    xt <- x$ai + x$ci
    yt <- x$bi + x$di
    Vi <- xt * yt * (n1i/Ni) * (n2i/Ni)/(Ni - 1)
    weight[x$not.na] <- Vi/sum(Vi) * 100
    names(weight) <- x$slab
    if (na.act == "na.omit") 
        weight <- weight[x$not.na]
    if (na.act == "na.fail" && any(!x$not.na)) 
        stop("Missing values in weights.")
    return(weight)
}
weights.rma.uni <-
function (object, ...) 
{
    if (!is.element("rma.uni", class(object))) 
        stop("Argument 'object' must be an object of class \"rma.uni\".")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    x <- object
    weight <- rep(NA, x$k.f)
    if (x$weighted) {
        wi <- 1/(x$vi + x$tau2)
        weight[x$not.na] <- wi/sum(wi) * 100
    }
    else {
        weight[x$not.na] <- 1/x$k * 100
    }
    names(weight) <- x$slab
    if (na.act == "na.omit") 
        weight <- weight[x$not.na]
    if (na.act == "na.fail" && any(!x$not.na)) 
        stop("Missing values in weights.")
    return(weight)
}
