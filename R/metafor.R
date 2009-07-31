.onLoad <-
function (lib, pkg) 
{
    loadmsg <- "Loading 'metafor' package. Type: help(metafor)\nfor an overview and introduction to the package."
    packageStartupMessage(loadmsg, domain = NULL, appendLF = TRUE)
}
addpoly <-
function (x, ...) 
UseMethod("addpoly")
addpoly.default <-
function (x, vi, sei, row = -1, level = 95, digits = 2, annotate = TRUE, 
    mlab = NULL, transf = FALSE, atransf = FALSE, targs = NULL, 
    col = "black", efac = 1, cex = NULL, ...) 
{
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail"))) 
        stop("Unknwn 'na.action' specified under options().")
    yi <- x
    if (missing(vi)) 
        vi <- sei^2
    yivi.na <- is.na(cbind(yi, vi))
    if (sum(yivi.na) > 0) {
        not.na <- apply(yivi.na, MARGIN = 1, sum) == 0
        if (na.act == "na.omit") {
            yi <- yi[not.na]
            vi <- vi[not.na]
            mlab <- mlab[not.na]
        }
        if (na.act == "na.fail") 
            stop("Missing values in results.")
    }
    k <- length(yi)
    alpha <- (100 - level)/100
    ci.lb <- yi - qnorm(1 - alpha/2) * sqrt(vi)
    ci.ub <- yi + qnorm(1 - alpha/2) * sqrt(vi)
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
    par.usr <- par("usr")
    height <- par.usr[4] - par.usr[3]
    cex.adj <- min(1, 20/height)
    xlim <- par.usr[1:2]
    if (is.null(cex)) 
        cex <- par("cex") * cex.adj
    if (length(row) == 1) 
        row <- row:(row - k + 1)
    for (i in 1:k) {
        polygon(x = c(ci.lb[i], yi[i], ci.ub[i], yi[i]), y = c(row[i], 
            row[i] + (height/100) * cex * efac, row[i], row[i] - 
                (height/100) * cex * efac), col = col, ...)
        if (annotate) {
            if (is.function(atransf)) {
                if (is.null(targs)) {
                  text(x = xlim[2], row[i], labels = paste(formatC(sapply(yi.ut[i], 
                    atransf), digits = digits, format = "f", 
                    flag = " "), "[", formatC(sapply(ci.lb.ut[i], 
                    atransf), digits = digits, format = "f", 
                    flag = " "), ",", formatC(sapply(ci.ub.ut[i], 
                    atransf), digits = digits, format = "f", 
                    flag = " "), "]"), pos = 2, cex = cex, ...)
                }
                else {
                  text(x = xlim[2], row[i], labels = paste(formatC(sapply(yi.ut[i], 
                    atransf, targs), digits = digits, format = "f", 
                    flag = " "), "[", formatC(sapply(ci.lb.ut[i], 
                    atransf, targs), digits = digits, format = "f", 
                    flag = " "), ",", formatC(sapply(ci.ub.ut[i], 
                    atransf, targs), digits = digits, format = "f", 
                    flag = " "), "]"), pos = 2, cex = cex, ...)
                }
            }
            else {
                text(x = xlim[2], row[i], labels = paste(formatC(yi[i], 
                  digits = digits, format = "f", flag = " "), 
                  "[", formatC(ci.lb[i], digits = digits, format = "f", 
                    flag = " "), ",", formatC(ci.ub[i], digits = digits, 
                    format = "f", flag = " "), "]"), pos = 2, 
                  cex = cex)
            }
        }
        if (!is.null(mlab)) {
            text(xlim[1], row[i], mlab[i], pos = 4, cex = cex, 
                ...)
        }
    }
}
addpoly.rma <-
function (x, row = -2, level = x$level, digits = 2, annotate = TRUE, 
    mlab = NULL, transf = FALSE, atransf = FALSE, targs = NULL, 
    col = "black", efac = 1, cex = NULL, ...) 
{
    if (!is.element("rma", class(x))) 
        stop("Argument 'x' must be an object of class \"rma\".")
    if (!x$int.only) 
        stop("The model should not contain moderators.")
    if (is.null(mlab)) 
        mlab <- ifelse((x$method == "FE"), "FE Model", "RE Model")
    addpoly(x$b, vi = x$vb, row = row, level = level, digits = digits, 
        annotate = annotate, mlab = mlab, transf = transf, atransf = atransf, 
        col = col, targs = targs, efac = efac, cex = cex, ...)
}
blup <-
function (x, ...) 
UseMethod("blup")
blup.rma.uni <-
function (x, level = x$level, digits = x$digits, transf = FALSE, 
    targs = NULL, ...) 
{
    if (!is.element("rma.uni", class(x))) 
        stop("Argument 'x' must be an object of class \"rma.uni\".")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail"))) 
        stop("Unknwn 'na.action' specified under options().")
    alpha <- (100 - level)/100
    if (!x$knha) {
        crit <- qnorm(1 - alpha/2)
    }
    else {
        crit <- qt(1 - alpha/2, df = x$k - x$p)
    }
    pred <- rep(NA, x$k.f)
    vpred <- rep(NA, x$k.f)
    li <- x$tau2/(x$tau2 + x$vi.f)
    for (i in (1:x$k.f)[x$not.na]) {
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
    if (na.act == "na.omit") {
        out <- list(pred = pred[x$not.na], se = se[x$not.na], 
            pi.lb = pi.lb[x$not.na], pi.ub = pi.ub[x$not.na])
        out$slab <- x$slab[x$not.na]
    }
    if (na.act == "na.exclude") {
        out <- list(pred = pred, se = se, pi.lb = pi.lb, pi.ub = pi.ub)
        out$slab <- x$slab
    }
    if (na.act == "na.fail") 
        stop("Missing values in results.")
    out$digits <- digits
    class(out) <- c("rma.list")
    return(out)
}
cint <-
function (object, ...) 
UseMethod("cint")
cint.rma.uni <-
function (object, fixed = FALSE, random = TRUE, level = object$level, 
    digits = object$digits, control = list(), ...) 
{
    if (!is.element("rma.uni", class(object))) 
        stop("Argument 'object' must be an object of class \"rma.uni\".")
    x <- object
    if (random) {
        .invcalc <- function(X, W, k) {
            wX <- sqrt(W) %*% X
            res.qrs <- qr.solve(wX, diag(k))
            res.qrs %*% t(res.qrs)
        }
        QE.func <- function(tau2val, Y, vi, X, objective, verbose) {
            wi <- 1/(vi + tau2val)
            W <- diag(wi)
            stXWX <- .invcalc(X = X, W = W, k = x$k)
            P <- W - W %*% X %*% stXWX %*% crossprod(X, W)
            RSS <- crossprod(Y, P) %*% Y
            if (verbose) 
                print(c(RSS - objective))
            RSS - objective
        }
        alpha <- (100 - level)/100
        crit.u <- qchisq(1 - alpha/2, x$k - x$p)
        crit.l <- qchisq(alpha/2, x$k - x$p)
        con <- list(tol = .Machine$double.eps^0.25, maxiter = 1000, 
            tau2.min = x$control$tau2.min, tau2.max = 50, verbose = FALSE)
        con[pmatch(names(control), names(con))] <- control
        status.lb <- 1
        status.ub <- 1
        conv <- 1
        if (QE.func(con$tau2.min, Y = cbind(x$yi), vi = x$vi, 
            X = x$X, objective = 0, verbose = FALSE) < crit.l) {
            tau2.lb <- NA
            tau2.ub <- NA
        }
        else {
            if (QE.func(con$tau2.min, Y = cbind(x$yi), vi = x$vi, 
                X = x$X, objective = 0, verbose = FALSE) > crit.u) {
                tau2.lb <- try(uniroot(QE.func, interval = c(con$tau2.min, 
                  con$tau2.max), tol = con$tol, maxiter = con$maxiter, 
                  Y = cbind(x$yi), vi = x$vi, X = x$X, objective = crit.u, 
                  verbose = con$verbose)$root, silent = TRUE)
                if (!is.numeric(tau2.lb)) {
                  tau2.lb <- NA
                  status.lb <- 0
                  conv <- 0
                }
            }
            else {
                tau2.lb <- con$tau2.min
            }
            tau2.ub <- try(uniroot(QE.func, interval = c(tau2.lb, 
                con$tau2.max), tol = con$tol, maxiter = con$maxiter, 
                Y = cbind(x$yi), vi = x$vi, X = x$X, objective = crit.l, 
                verbose = con$verbose)$root, silent = TRUE)
            if (is.numeric(tau2.ub) == FALSE) {
                tau2.ub <- NA
                status.ub <- 0
                conv <- 0
            }
        }
        if (status.lb == 0) {
            warning("Error in iterative search for the lower bound.")
        }
        if (status.ub == 0) {
            warning("Error in iterative search for the upper bound.")
        }
        if (conv == 0) {
            stop("Try increasing tau2.max (via the 'control' argument).")
        }
        if (x$int.only) {
            wi <- 1/x$vi
            s2 <- (x$k - 1) * sum(wi)/(sum(wi)^2 - sum(wi^2))
            I2.lb <- tau2.lb/(tau2.lb + s2) * 100
            I2.ub <- tau2.ub/(tau2.ub + s2) * 100
        }
        else {
            I2.lb <- NA
            I2.ub <- NA
        }
        if (is.na(tau2.lb) && is.na(tau2.lb)) {
            cat("The upper and lower bound both fall below ", 
                con$tau2.min, ".\nThe CI is therefore equal to the null set.\n\n", 
                sep = "")
        }
        tau2 <- round(c(x$tau2, tau2.lb, tau2.ub), digits)
        tau <- round(sqrt(c(ifelse(x$tau2 >= 0, x$tau2, NA), 
            ifelse(tau2.lb >= 0, tau2.lb, NA), ifelse(tau2.ub >= 
                0, tau2.ub, NA))), digits)
        I2 <- round(c(x$I2, I2.lb, I2.ub), digits)
        if (x$int.only) {
            res.random <- rbind(tau2, tau, I2)
            dimnames(res.random)[[1]] <- c("tau^2", "tau", "I^2(%)")
        }
        else {
            res.random <- rbind(tau2, tau)
            dimnames(res.random)[[1]] <- c("tau^2", "tau")
        }
        dimnames(res.random)[[2]] <- c("estimate", "ci.lb", "ci.ub")
    }
    if (fixed) {
        alpha <- (100 - level)/100
        if (x$knha) {
            crit <- qt(1 - alpha/2, df = x$k - x$p)
        }
        else {
            crit <- qnorm(1 - alpha/2)
        }
        ci.lb <- c(x$b - crit * x$se)
        ci.ub <- c(x$b + crit * x$se)
        res.fixed <- round(cbind(x$b, ci.lb, ci.ub), digits)
        dimnames(res.fixed)[[2]] <- c("estimate", "ci.lb", "ci.ub")
    }
    if (fixed && random) {
        res <- list(fixed = data.frame(res.fixed), random = data.frame(res.random))
        return(res)
    }
    if (fixed) 
        return(data.frame(res.fixed))
    if (random) 
        return(data.frame(res.random))
}
coef.rma <-
function (object, ...) 
{
    if (!is.element("rma", class(object))) 
        stop("Argument 'object' must be an object of class \"rma\".")
    x <- object
    res.table <- cbind(x$b, x$se, x$zval, x$pval, x$ci.lb, x$ci.ub)
    dimnames(res.table)[[2]] <- c("estimate", "se", "zval", "pval", 
        "ci.lb", "ci.ub")
    if (is.element("rma.uni", class(x)) && x$knha) {
        dimnames(res.table)[[2]][3] <- c("tval")
    }
    res.table <- data.frame(res.table)
    return(res.table)
}
escalc <-
function (measure, ai, bi, ci, di, n1i, n2i, m1i, m2i, sd1i, 
    sd2i, xi, mi, ri, ni, data = NULL, add = 1/2, to = "only0", 
    vtype = "LS") 
{
    if (!is.element(measure, c("MD", "SMD", "RR", "OR", "PETO", 
        "RD", "AS", "PR", "PLN", "PLO", "PAS", "PFT", "COR", 
        "ZCOR"))) 
        stop("Unknown 'measure' specified.")
    if (!is.element(to, c("all", "only0", "if0all", "none"))) 
        stop("Unknown 'to' argument specified.")
    if (!is.element(vtype, c("UB", "LS", "HS"))) 
        stop("Unknown 'vtype' argument specified.")
    if (is.null(data)) {
        data <- sys.frame(sys.parent())
    }
    else {
        if (!is.data.frame(data)) {
            data <- data.frame(data)
        }
    }
    mf <- match.call()
    if (is.element(measure, c("MD", "SMD"))) {
        mf.m1i <- mf[[match("m1i", names(mf))]]
        mf.m2i <- mf[[match("m2i", names(mf))]]
        mf.sd1i <- mf[[match("sd1i", names(mf))]]
        mf.sd2i <- mf[[match("sd2i", names(mf))]]
        mf.n1i <- mf[[match("n1i", names(mf))]]
        mf.n2i <- mf[[match("n2i", names(mf))]]
        m1i <- eval(mf.m1i, data)
        m2i <- eval(mf.m2i, data)
        sd1i <- eval(mf.sd1i, data)
        sd2i <- eval(mf.sd2i, data)
        n1i <- eval(mf.n1i, data)
        n2i <- eval(mf.n2i, data)
        if (measure == "MD") {
            yi <- m1i - m2i
            vi <- sd1i^2/n1i + sd2i^2/n2i
        }
        if (measure == "SMD") {
            cNm2ifunc <- function(Nm2i) {
                cNm2i <- gamma(Nm2i/2)/(sqrt(Nm2i/2) * gamma((Nm2i - 
                  1)/2))
                isna <- is.na(cNm2i)
                cNm2i[isna] <- 1 - 3/(4 * Nm2i[isna] - 1)
                cNm2i
            }
            Nm2i <- n1i + n2i - 2
            warn.before <- getOption("warn")
            options(warn = -1)
            cNm2i <- cNm2ifunc(Nm2i)
            options(warn = warn.before)
            nti <- (n1i * n2i)/(n1i + n2i)
            yi <- cNm2i * (m1i - m2i)/sqrt(((n1i - 1) * sd1i^2 + 
                (n2i - 1) * sd2i^2)/(n1i + n2i - 2))
            if (vtype == "UB") {
                vi <- 1/nti + (1 - (Nm2i - 2)/(Nm2i * cNm2i^2)) * 
                  yi^2
            }
            if (vtype == "LS") {
                vi <- 1/nti + yi^2/(2 * (n1i + n2i))
            }
            if (vtype == "HS") {
                myi <- sum((n1i + n2i) * yi)/sum(n1i + n2i)
                vi <- 1/nti + myi^2/(2 * (n1i + n2i))
            }
        }
    }
    if (is.element(measure, c("RR", "OR", "RD", "AS", "PETO"))) {
        mf.ai <- mf[[match("ai", names(mf))]]
        mf.bi <- mf[[match("bi", names(mf))]]
        mf.ci <- mf[[match("ci", names(mf))]]
        mf.di <- mf[[match("di", names(mf))]]
        mf.n1i <- mf[[match("n1i", names(mf))]]
        mf.n2i <- mf[[match("n2i", names(mf))]]
        ai <- eval(mf.ai, data)
        bi <- eval(mf.bi, data)
        ci <- eval(mf.ci, data)
        di <- eval(mf.di, data)
        n1i <- eval(mf.n1i, data)
        n2i <- eval(mf.n2i, data)
        if (is.null(bi)) {
            bi <- n1i - ai
        }
        if (is.null(di)) {
            di <- n2i - ci
        }
        if (to == "all") {
            ai <- ai + add
            ci <- ci + add
            bi <- bi + add
            di <- di + add
        }
        if (to == "only0") {
            id0 <- c(ai == 0 | ci == 0 | bi == 0 | di == 0)
            id0[is.na(id0)] <- FALSE
            ai[id0] <- ai[id0] + add
            ci[id0] <- ci[id0] + add
            bi[id0] <- bi[id0] + add
            di[id0] <- di[id0] + add
        }
        if (to == "if0all") {
            id0 <- c(ai == 0 | ci == 0 | bi == 0 | di == 0)
            id0[is.na(id0)] <- FALSE
            if (sum(id0) > 0) {
                ai <- ai + add
                ci <- ci + add
                bi <- bi + add
                di <- di + add
            }
        }
        n1i <- ai + bi
        n2i <- ci + di
        p1 <- ai/n1i
        p2 <- ci/n2i
        if (measure == "RR") {
            yi <- log(p1) - log(p2)
            vi <- 1/ai - 1/n1i + 1/ci - 1/n2i
        }
        if (measure == "OR") {
            yi <- log(p1/(1 - p1)) - log(p2/(1 - p2))
            vi <- 1/ai + 1/bi + 1/ci + 1/di
        }
        if (measure == "PETO") {
            xt <- ai + ci
            yt <- bi + di
            Ni <- ai + ci + bi + di
            Oi <- ai
            Ei <- xt * n1i/Ni
            Vi <- xt * yt * (n1i/Ni) * (n2i/Ni)/(Ni - 1)
            yi <- (ai - Ei)/Vi
            vi <- 1/Vi
        }
        if (measure == "RD") {
            yi <- p1 - p2
            vi <- p1 * (1 - p1)/n1i + p2 * (1 - p2)/n2i
        }
        if (measure == "AS") {
            yi <- asin(sqrt(p1)) - asin(sqrt(p2))
            vi <- 1/(4 * n1i) + 1/(4 * n2i)
        }
    }
    if (is.element(measure, c("PR", "PLN", "PLO", "PAS", "PFT"))) {
        mf.xi <- mf[[match("xi", names(mf))]]
        mf.mi <- mf[[match("mi", names(mf))]]
        mf.ni <- mf[[match("ni", names(mf))]]
        xi <- eval(mf.xi, data)
        mi <- eval(mf.mi, data)
        ni <- eval(mf.ni, data)
        if (is.null(mi)) {
            mi <- ni - xi
        }
        if (to == "all") {
            xi <- xi + add
            mi <- mi + add
        }
        if (to == "only0") {
            id0 <- c(xi == 0 | mi == 0)
            id0[is.na(id0)] <- FALSE
            xi[id0] <- xi[id0] + add
            mi[id0] <- mi[id0] + add
        }
        if (to == "if0all") {
            id0 <- c(xi == 0 | mi == 0)
            id0[is.na(id0)] <- FALSE
            if (sum(id0) > 0) {
                xi <- xi + add
                mi <- mi + add
            }
        }
        ni <- xi + mi
        pri <- xi/ni
        if (measure == "PR") {
            yi <- pri
            vi <- pri * (1 - pri)/ni
        }
        if (measure == "PLN") {
            yi <- log(pri)
            vi <- 1/xi - 1/ni
        }
        if (measure == "PLO") {
            yi <- log(pri/(1 - pri))
            vi <- 1/xi + 1/mi
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
    if (is.element(measure, c("COR", "ZCOR"))) {
        mf.ri <- mf[[match("ri", names(mf))]]
        mf.ni <- mf[[match("ni", names(mf))]]
        ri <- eval(mf.ri, data)
        ni <- eval(mf.ni, data)
        if (measure == "COR") {
            yi <- ri + ri * (1 - ri^2)/(2 * (ni - 4))
            if (vtype == "UB") {
                vi <- yi^2 - 1 + (ni - 3)/(ni - 2) * ((1 - ri^2) + 
                  2 * (1 - ri^2)^2/ni + 8 * (1 - ri^2)^3/(ni * 
                  (ni + 2)) + 48 * (1 - ri^2)^4/(ni * (ni + 2) * 
                  (ni + 4)))
            }
            if (vtype == "LS") {
                vi <- (1 - ri^2)^2/(ni - 1)
            }
            if (vtype == "HS") {
                mr <- sum(ni * ri)/sum(ni)
                vi <- (1 - mr^2)^2/(ni - 1)
            }
        }
        if (measure == "ZCOR") {
            yi <- 1/2 * log((1 + ri)/(1 - ri))
            vi <- 1/(ni - 3)
        }
    }
    if (sum(is.infinite(c(yi, vi)) == TRUE) > 0) {
        warning("Some yi and/or vi equal to +-Inf. Recoded to NAs.")
        k <- length(yi)
        inf.ids <- (1:k)[is.infinite(yi) == TRUE | is.infinite(vi) == 
            TRUE]
        yi[inf.ids] <- NA
        vi[inf.ids] <- NA
    }
    dat <- data.frame(cbind(yi, vi))
    return(dat)
}
fitstats <-
function (x, ...) 
UseMethod("fitstats")
fitstats.rma <-
function (x, REML = NULL, ...) 
{
    if (!is.element("rma", class(x))) 
        stop("Argument 'x' must be an object of class \"rma\".")
    if (is.null(REML)) {
        if (x$method == "REML") {
            REML <- TRUE
        }
        else {
            REML <- FALSE
        }
    }
    if (REML) {
        out <- cbind(x$fit.stats$REML)
        dimnames(out)[[1]] <- c("Log-Likelihood: ", "Deviance (-2RLL): ", 
            "AIC: ", "BIC: ")
        dimnames(out)[[2]] <- c("REML")
    }
    else {
        out <- cbind(x$fit.stats$ML)
        dimnames(out)[[1]] <- c("Log-Likelihood: ", "Deviance (-2LL): ", 
            "AIC: ", "BIC: ")
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
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail"))) 
        stop("Unknwn 'na.action' specified under options().")
    out <- c(object$X.f %*% object$b)
    names(out) <- object$slab
    not.na <- !is.na(out)
    if (na.act == "na.omit") {
        out <- out[not.na]
    }
    if (na.act == "na.fail") 
        stop("Missing values in results.")
    return(out)
}
forest <-
function (x, ...) 
UseMethod("forest")
forest.default <-
function (x, vi, sei, annotate = TRUE, xlim = NULL, alim = NULL, 
    ylim = NULL, at = NULL, steps = 5, level = 95, digits = 2, 
    refline = 0, xlab = NULL, slab = NULL, ilab = NULL, ilab.xpos = NULL, 
    ilab.pos = NULL, subset = NULL, transf = FALSE, atransf = FALSE, 
    targs = NULL, addrows = 0, efac = 1, pch = 15, psize = NULL, 
    cex = NULL, cex.lab = NULL, cex.axis = NULL, ...) 
{
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail"))) 
        stop("Unknwn 'na.action' specified under options().")
    yi <- x
    if (missing(vi)) 
        vi <- sei^2
    if (length(yi) != length(vi)) 
        stop("Length of yi and vi (or sei) is not the same.")
    k <- length(yi)
    if (is.null(slab)) 
        slab <- paste("Study ", 1:k)
    if (is.vector(ilab)) 
        ilab <- cbind(ilab)
    if (is.null(psize)) {
        wi <- 1/vi
        wi[is.infinite(wi)] <- 2 * max(wi, na.rm = TRUE)
        psize <- wi/sum(wi, na.rm = TRUE)
        psize <- (psize - min(psize, na.rm = TRUE))/(max(psize, 
            na.rm = TRUE) - min(psize, na.rm = TRUE))
        psize <- (psize * 0.9) + 0.6
    }
    else {
        if (length(psize) == 1) 
            psize <- rep(psize, k)
    }
    if (!is.null(subset)) {
        yi <- yi[subset]
        vi <- vi[subset]
        slab <- slab[subset]
        psize <- psize[subset]
        ilab <- ilab[subset, , drop = FALSE]
    }
    yivi.na <- is.na(cbind(yi, vi))
    if (sum(yivi.na) > 0) {
        not.na <- apply(yivi.na, MARGIN = 1, sum) == 0
        if (na.act == "na.omit") {
            yi <- yi[not.na]
            vi <- vi[not.na]
            slab <- slab[not.na]
            psize <- psize[not.na]
            ilab <- ilab[not.na, , drop = FALSE]
        }
        if (na.act == "na.fail") 
            stop("Missing values in results.")
    }
    k <- length(yi)
    alpha <- (100 - level)/100
    ci.lb <- yi - qnorm(1 - alpha/2) * sqrt(vi)
    ci.ub <- yi + qnorm(1 - alpha/2) * sqrt(vi)
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
    rng <- max(ci.ub, na.rm = TRUE) - min(ci.lb, na.rm = TRUE)
    if (is.null(xlim)) {
        xlim <- c(min(ci.lb, na.rm = TRUE) - rng * 1.2, max(ci.ub, 
            na.rm = TRUE) + rng * 1.2)
        xlim <- round(xlim, digits)
    }
    if (is.null(alim)) {
        if (is.null(at)) {
            alim <- c(min(ci.lb, na.rm = TRUE) - rng * 0.2, max(ci.ub, 
                na.rm = TRUE) + rng * 0.2)
            alim <- round(alim, digits)
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
    addrows <- round(addrows)
    addrows[addrows < 0] <- 0
    if (addrows > 0) 
        addrows <- addrows + 1
    if (is.null(ylim)) {
        ylim <- c(0.5 - addrows, k + 1)
    }
    else {
        ylim <- sort(ylim)
    }
    if (is.null(at)) {
        at <- seq(alim[1], alim[2], length = steps)
    }
    else {
        at[at < alim[1]] <- alim[1]
        at[at > alim[2]] <- alim[2]
        at <- unique(at)
    }
    at.lab <- at
    if (is.function(atransf)) {
        if (is.null(targs)) {
            at.lab <- formatC(sapply(at.lab, atransf), digits = digits, 
                format = "f")
        }
        else {
            at.lab <- formatC(sapply(at.lab, atransf, targs), 
                digits = digits, format = "f")
        }
    }
    else {
        at.lab <- round(at, digits)
    }
    par.mar <- par("mar")
    par.mar.adj <- par.mar - c(0, 3, 0, 1)
    par.mar.adj[par.mar.adj < 1] <- 1
    par(mar = par.mar.adj)
    on.exit(par(mar = par.mar))
    plot(NA, NA, xlim = xlim, ylim = ylim, xlab = "", ylab = "", 
        yaxt = "n", xaxt = "n", xaxs = "i", bty = "n", ...)
    abline(h = k + 1, ...)
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
        cex.lab <- par("cex.lab") * cex.adj
    if (is.null(cex.axis)) 
        cex.axis <- par("cex.axis") * cex.adj
    axis(side = 1, at = at, labels = at.lab, cex.axis = cex.axis, 
        ...)
    if (!is.null(xlab)) 
        mtext(xlab, side = 1, at = min(at) + (max(at) - min(at))/2, 
            line = 2.75, cex = cex.lab, ...)
    for (i in 1:k) {
        if (is.na(yi[i]) || is.na(vi)[i]) 
            next
        segments(max(ci.lb[i], alim[1]), i, min(ci.ub[i], alim[2]), 
            i, ...)
        if (ci.lb[i] >= alim[1]) {
            segments(ci.lb[i], i - (k/100) * cex * efac, ci.lb[i], 
                i + (k/100) * cex * efac, ...)
        }
        else {
            polygon(x = c(alim[1], alim[1] + (1.4/100) * cex * 
                (xlim[2] - xlim[1]), alim[1] + (1.4/100) * cex * 
                (xlim[2] - xlim[1]), alim[1]), y = c(i, i + (k/100) * 
                cex * efac, i - (k/100) * cex * efac, i), col = "black", 
                ...)
        }
        if (ci.ub[i] <= alim[2]) {
            segments(ci.ub[i], i - (k/100) * cex * efac, ci.ub[i], 
                i + (k/100) * cex * efac, ...)
        }
        else {
            polygon(x = c(alim[2], alim[2] - (1.4/100) * cex * 
                (xlim[2] - xlim[1]), alim[2] - (1.4/100) * cex * 
                (xlim[2] - xlim[1]), alim[2]), y = c(i, i + (k/100) * 
                cex * efac, i - (k/100) * cex * efac, i), col = "black", 
                ...)
        }
    }
    text(xlim[1], 1:k, slab, pos = 4, cex = cex, ...)
    if (!is.null(ilab)) {
        for (l in 1:dim(ilab)[2]) {
            text(ilab.xpos[l], 1:k, ilab[, l], offset = 0, pos = ilab.pos[l], 
                cex = cex, ...)
        }
    }
    if (is.numeric(refline)) 
        segments(refline, ylim[1] - 1, refline, (k + 1), lty = "dotted", 
            ...)
    if (annotate) {
        if (is.function(atransf)) {
            if (is.null(targs)) {
                text(x = xlim[2], 1:k, labels = paste(formatC(sapply(yi, 
                  atransf), digits = digits, format = "f", flag = " "), 
                  "[", formatC(sapply(ci.lb, atransf), digits = digits, 
                    format = "f", flag = " "), ",", formatC(sapply(ci.ub, 
                    atransf), digits = digits, format = "f", 
                    flag = " "), "]"), pos = 2, cex = cex, ...)
            }
            else {
                text(x = xlim[2], 1:k, labels = paste(formatC(sapply(yi, 
                  atransf, targs), digits = digits, format = "f", 
                  flag = " "), "[", formatC(sapply(ci.lb, atransf, 
                  targs), digits = digits, format = "f", flag = " "), 
                  ",", formatC(sapply(ci.ub, atransf, targs), 
                    digits = digits, format = "f", flag = " "), 
                  "]"), pos = 2, cex = cex, ...)
            }
        }
        else {
            text(x = xlim[2], 1:k, labels = paste(formatC(yi, 
                digits = digits, format = "f", flag = " "), "[", 
                formatC(ci.lb, digits = digits, format = "f", 
                  flag = " "), ",", formatC(ci.ub, digits = digits, 
                  format = "f", flag = " "), "]"), pos = 2, cex = cex, 
                ...)
        }
    }
    points(yi, 1:k, pch = pch, cex = cex * psize, ...)
    if (addrows > 0) 
        abline(h = 0, ...)
    invisible()
}
forest.rma <-
function (x, annotate = TRUE, addfit = TRUE, xlim = NULL, alim = NULL, 
    ylim = NULL, at = NULL, steps = 5, level = x$level, digits = 2, 
    refline = 0, xlab = NULL, slab = NULL, mlab = NULL, ilab = NULL, 
    ilab.xpos = NULL, ilab.pos = NULL, order = NULL, transf = FALSE, 
    atransf = FALSE, targs = NULL, addrows = 0, efac = 1, pch = 15, 
    psize = NULL, col = "darkgray", border = "darkgray", cex = NULL, 
    cex.lab = NULL, cex.axis = NULL, ...) 
{
    if (!is.element("rma", class(x))) 
        stop("Argument 'x' must be an object of class \"rma\".")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail"))) 
        stop("Unknwn 'na.action' specified under options().")
    yi <- x$yi.f
    vi <- x$vi.f
    X <- x$X.f
    ids <- x$ids
    if (is.null(slab)) {
        if (x$slab.null) {
            slab <- paste("Study ", ids)
        }
        else {
            slab <- x$slab
        }
    }
    if (is.vector(ilab)) 
        ilab <- cbind(ilab)
    k <- length(yi)
    options(na.action = "na.exclude")
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
    options(na.action = na.act)
    if (is.null(psize)) {
        wi <- 1/vi
        wi[is.infinite(wi)] <- 2 * max(wi, na.rm = TRUE)
        psize <- wi/sum(wi, na.rm = TRUE)
        psize <- (psize - min(psize, na.rm = TRUE))/(max(psize, 
            na.rm = TRUE) - min(psize, na.rm = TRUE))
        psize <- (psize * 0.9) + 0.6
    }
    else {
        if (length(psize) == 1) 
            psize <- rep(psize, k)
    }
    if (!is.null(order)) {
        sort.ids <- 1:k
        if (length(order) == k) {
            sort.ids <- order
        }
        else {
            if (order == "obs") 
                sort.ids <- order(yi)
            if (order == "fit") 
                sort.ids <- order(pred)
            if (order == "prec") 
                sort.ids <- order(vi, yi)
            if (order == "resid") 
                sort.ids <- order(abs(yi - pred), yi)
        }
        yi <- yi[sort.ids]
        vi <- vi[sort.ids]
        X <- X[sort.ids, , drop = FALSE]
        ids <- ids[sort.ids]
        slab <- slab[sort.ids]
        pred <- pred[sort.ids]
        pred.ci.lb <- pred.ci.lb[sort.ids]
        pred.ci.ub <- pred.ci.ub[sort.ids]
        psize <- psize[sort.ids]
        ilab <- ilab[sort.ids, , drop = FALSE]
    }
    yiviX.na <- is.na(cbind(yi, vi, X))
    if (sum(yiviX.na) > 0) {
        not.na <- apply(yiviX.na, MARGIN = 1, sum) == 0
        if (na.act == "na.omit") {
            yi <- yi[not.na]
            vi <- vi[not.na]
            X <- X[not.na, , drop = FALSE]
            ids <- ids[not.na]
            slab <- slab[not.na]
            pred <- pred[not.na]
            pred.ci.lb <- pred.ci.lb[not.na]
            pred.ci.ub <- pred.ci.ub[not.na]
            psize <- psize[not.na]
            ilab <- ilab[not.na, , drop = FALSE]
        }
        if (na.act == "na.fail") 
            stop("Missing values in results.")
    }
    k <- length(yi)
    alpha <- (100 - level)/100
    ci.lb <- yi - qnorm(1 - alpha/2) * sqrt(vi)
    ci.ub <- yi + qnorm(1 - alpha/2) * sqrt(vi)
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
    rng <- max(ci.ub, na.rm = TRUE) - min(ci.lb, na.rm = TRUE)
    if (is.null(xlim)) {
        xlim <- c(min(ci.lb, na.rm = TRUE) - rng * 1.2, max(ci.ub, 
            na.rm = TRUE) + rng * 1.2)
        xlim <- round(xlim, digits)
    }
    if (is.null(alim)) {
        if (is.null(at)) {
            alim <- c(min(ci.lb, na.rm = TRUE) - rng * 0.2, max(ci.ub, 
                na.rm = TRUE) + rng * 0.2)
            alim <- round(alim, digits)
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
    addrows <- round(addrows)
    addrows[addrows < 0] <- 0
    if (x$int.only && addfit) {
        addrows <- addrows + 2
    }
    else {
        if (addrows > 0) 
            addrows <- addrows + 1
    }
    if (is.null(ylim)) {
        ylim <- c(0.5 - addrows, k + 1)
    }
    else {
        ylim <- sort(ylim)
    }
    if (is.null(at)) {
        at <- seq(alim[1], alim[2], length = steps)
    }
    else {
        at[at < alim[1]] <- alim[1]
        at[at > alim[2]] <- alim[2]
        at <- unique(at)
    }
    at.lab <- at
    if (is.function(atransf)) {
        if (is.null(targs)) {
            at.lab <- formatC(sapply(at.lab, atransf), digits = digits, 
                format = "f")
        }
        else {
            at.lab <- formatC(sapply(at.lab, atransf, targs), 
                digits = digits, format = "f")
        }
    }
    else {
        at.lab <- round(at, digits)
    }
    par.mar <- par("mar")
    par.mar.adj <- par.mar - c(0, 3, 0, 1)
    par.mar.adj[par.mar.adj < 1] <- 1
    par(mar = par.mar.adj)
    on.exit(par(mar = par.mar))
    plot(NA, NA, xlim = xlim, ylim = ylim, xlab = "", ylab = "", 
        yaxt = "n", xaxt = "n", xaxs = "i", bty = "n", ...)
    abline(h = k + 1, ...)
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
        cex.lab <- par("cex.lab") * cex.adj
    if (is.null(cex.axis)) 
        cex.axis <- par("cex.axis") * cex.adj
    if (addfit && !x$int.only) {
        for (i in 1:k) {
            if (is.na(pred[i])) 
                next
            if ((pred.ci.lb[i] > alim[1]) & (pred.ci.ub[i] < 
                alim[2])) 
                polygon(x = c(pred.ci.lb[i], pred[i], pred.ci.ub[i], 
                  pred[i]), y = c(i, i + (height/100) * cex * 
                  efac, i, i - (height/100) * cex * efac), col = col, 
                  border = border, ...)
        }
    }
    if (addfit && x$int.only) {
        b <- x$b
        b.ci.lb <- x$ci.lb
        b.ci.ub <- x$ci.ub
        if (is.function(transf)) {
            if (is.null(targs)) {
                b <- sapply(b, transf)
                b.ci.lb <- sapply(b.ci.lb, transf)
                b.ci.ub <- sapply(b.ci.ub, transf)
            }
            else {
                b <- sapply(b, transf, targs)
                b.ci.lb <- sapply(b.ci.lb, transf, targs)
                b.ci.ub <- sapply(b.ci.ub, transf, targs)
            }
        }
        polygon(x = c(b.ci.lb, b, b.ci.ub, b), y = c(-1, -1 + 
            (height/100) * cex * efac, -1, -1 - (height/100) * 
            cex * efac), col = "black", ...)
        if (annotate) {
            if (is.function(atransf)) {
                if (is.null(targs)) {
                  text(x = xlim[2], -1, labels = paste(formatC(sapply(b, 
                    atransf), digits = digits, format = "f", 
                    flag = " "), "[", formatC(sapply(b.ci.lb, 
                    atransf), digits = digits, format = "f", 
                    flag = " "), ",", formatC(sapply(b.ci.ub, 
                    atransf), digits = digits, format = "f", 
                    flag = " "), "]"), pos = 2, cex = cex, ...)
                }
                else {
                  text(x = xlim[2], -1, labels = paste(formatC(sapply(b, 
                    atransf, targs), digits = digits, format = "f", 
                    flag = " "), "[", formatC(sapply(b.ci.lb, 
                    atransf, targs), digits = digits, format = "f", 
                    flag = " "), ",", formatC(sapply(b.ci.ub, 
                    atransf, targs), digits = digits, format = "f", 
                    flag = " "), "]"), pos = 2, cex = cex, ...)
                }
            }
            else {
                text(x = xlim[2], -1, labels = paste(formatC(b, 
                  digits = digits, format = "f", flag = " "), 
                  "[", formatC(b.ci.lb, digits = digits, format = "f", 
                    flag = " "), ",", formatC(b.ci.ub, digits = digits, 
                    format = "f", flag = " "), "]"), pos = 2, 
                  cex = cex, ...)
            }
        }
        if (is.null(mlab)) 
            mlab <- ifelse((x$method == "FE"), "FE Model", "RE Model")
        text(xlim[1], -1, mlab, pos = 4, cex = cex, ...)
    }
    axis(side = 1, at = at, labels = at.lab, cex.axis = cex.axis, 
        ...)
    if (!is.null(xlab)) 
        mtext(xlab, side = 1, at = min(at) + (max(at) - min(at))/2, 
            line = 2.75, cex = cex.lab, ...)
    for (i in 1:k) {
        if (is.na(yi[i]) || is.na(vi)[i]) 
            next
        segments(max(ci.lb[i], alim[1]), i, min(ci.ub[i], alim[2]), 
            i, ...)
        if (ci.lb[i] >= alim[1]) {
            segments(ci.lb[i], i - (k/100) * cex * efac, ci.lb[i], 
                i + (k/100) * cex * efac, ...)
        }
        else {
            polygon(x = c(alim[1], alim[1] + (1.4/100) * cex * 
                (xlim[2] - xlim[1]), alim[1] + (1.4/100) * cex * 
                (xlim[2] - xlim[1]), alim[1]), y = c(i, i + (k/100) * 
                cex * efac, i - (k/100) * cex * efac, i), col = "black", 
                ...)
        }
        if (ci.ub[i] <= alim[2]) {
            segments(ci.ub[i], i - (k/100) * cex * efac, ci.ub[i], 
                i + (k/100) * cex * efac, ...)
        }
        else {
            polygon(x = c(alim[2], alim[2] - (1.4/100) * cex * 
                (xlim[2] - xlim[1]), alim[2] - (1.4/100) * cex * 
                (xlim[2] - xlim[1]), alim[2]), y = c(i, i + (k/100) * 
                cex * efac, i - (k/100) * cex * efac, i), col = "black", 
                ...)
        }
    }
    text(xlim[1], 1:k, slab, pos = 4, cex = cex, ...)
    if (!is.null(ilab)) {
        for (l in 1:dim(ilab)[2]) {
            text(ilab.xpos[l], 1:k, ilab[, l], offset = 0, pos = ilab.pos[l], 
                cex = cex, ...)
        }
    }
    if (is.numeric(refline)) 
        segments(refline, ylim[1] - 1, refline, (k + 1), lty = "dotted", 
            ...)
    if (annotate) {
        if (is.function(atransf)) {
            if (is.null(targs)) {
                text(x = xlim[2], 1:k, labels = paste(formatC(sapply(yi, 
                  atransf), digits = digits, format = "f", flag = " "), 
                  "[", formatC(sapply(ci.lb, atransf), digits = digits, 
                    format = "f", flag = " "), ",", formatC(sapply(ci.ub, 
                    atransf), digits = digits, format = "f", 
                    flag = " "), "]"), pos = 2, cex = cex, ...)
            }
            else {
                text(x = xlim[2], 1:k, labels = paste(formatC(sapply(yi, 
                  atransf, targs), digits = digits, format = "f", 
                  flag = " "), "[", formatC(sapply(ci.lb, atransf, 
                  targs), digits = digits, format = "f", flag = " "), 
                  ",", formatC(sapply(ci.ub, atransf, targs), 
                    digits = digits, format = "f", flag = " "), 
                  "]"), pos = 2, cex = cex, ...)
            }
        }
        else {
            text(x = xlim[2], 1:k, labels = paste(formatC(yi, 
                digits = digits, format = "f", flag = " "), "[", 
                formatC(ci.lb, digits = digits, format = "f", 
                  flag = " "), ",", formatC(ci.ub, digits = digits, 
                  format = "f", flag = " "), "]"), pos = 2, cex = cex, 
                ...)
        }
    }
    points(yi, 1:k, pch = pch, cex = cex * psize, ...)
    if (addrows > 0) 
        abline(h = 0, ...)
    invisible()
}
fsn <-
function (yi, vi, sei, data = NULL, subset = NULL, alpha = 0.05) 
{
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail"))) 
        stop("Unknwn 'na.action' specified under options().")
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
    yi <- eval(mf.yi, data)
    vi <- eval(mf.vi, data)
    sei <- eval(mf.sei, data)
    subset <- eval(mf.subset, data)
    if (missing(vi)) 
        vi <- sei^2
    if (length(yi) != length(vi)) 
        stop("Length of yi and vi (or sei) is not the same.")
    if (!is.null(subset)) {
        yi <- yi[subset]
        vi <- vi[subset]
    }
    yivi.na <- is.na(cbind(yi, vi))
    if (sum(yivi.na) > 0) {
        not.na <- apply(yivi.na, MARGIN = 1, sum) == 0
        if (na.act == "na.omit" || na.act == "na.exclude") {
            yi <- yi[not.na]
            vi <- vi[not.na]
        }
        if (na.act == "na.fail") 
            stop("Missing values in results.")
    }
    k <- length(yi)
    zi <- yi/sqrt(vi)
    z.avg <- sum(zi)/sqrt(k)
    fsnum <- max(0, k * (z.avg/qnorm(1 - alpha))^2 - k)
    fsnum <- data.frame(ceiling(fsnum), row.names = c("Fail-safe N:"))
    dimnames(fsnum)[[2]] <- " "
    return(fsnum)
}
funnel <-
function (x, ...) 
UseMethod("funnel")
funnel.rma <-
function (x, xlim = NULL, ylim = NULL, xlab = NULL, ylab = "Standard Error", 
    steps = 5, level = x$level, digits = 3, addtau2 = FALSE, 
    type = "rstandard", back = "lightgray", shade = "white", 
    hlines = "white", refline = NULL, pch = 19, pch.fill = 21, 
    ...) 
{
    if (!is.element("rma", class(x))) 
        stop("Argument 'x' must be an object of class \"rma\".")
    if (!is.element(type, c("rstandard", "rstudent"))) 
        stop("Argument 'type' must be an either 'rstandard' or 'rstudent'.")
    if (x$int.only) {
        if (is.null(refline)) 
            refline <- x$b
        if (addtau2) {
            tau2 <- x$tau2
        }
        else {
            tau2 <- 0
        }
        yi <- x$yi
        vi <- x$vi
        sei <- sqrt(vi)
        if (is.null(xlab)) 
            xlab <- "Observed Outcome"
    }
    else {
        if (is.null(refline)) 
            refline <- 0
        tau2 <- 0
        if (type == "rstandard") {
            res <- rstandard(x)
        }
        else {
            res <- rstudent(x)
        }
        not.na <- !is.na(res$resid)
        yi <- res$resid[not.na]
        sei <- res$se[not.na]
        if (is.null(xlab)) 
            xlab <- "Residual Value"
    }
    if (is.null(ylim)) {
        ylim <- c(0, max(sei) * 1)
    }
    else {
        ylim <- sort(ylim)
        if (ylim[1] < 0 || ylim[2] < 0) 
            stop("Both limits for the y-axis must be >= 0.")
    }
    alpha <- (100 - level)/100
    alpha.min <- min(alpha)
    x.lb.bot <- refline - qnorm(1 - alpha.min/2) * sqrt(ylim[2]^2 + 
        tau2)
    x.ub.bot <- refline + qnorm(1 - alpha.min/2) * sqrt(ylim[2]^2 + 
        tau2)
    x.lb.top <- refline - qnorm(1 - alpha.min/2) * sqrt(ylim[1]^2 + 
        tau2)
    x.ub.top <- refline + qnorm(1 - alpha.min/2) * sqrt(ylim[1]^2 + 
        tau2)
    if (is.null(xlim)) {
        xlim <- c(min(x.lb.bot, min(yi)), max(x.ub.bot, max(yi)))
        rxlim <- xlim[2] - xlim[1]
        xlim[1] <- xlim[1] - (rxlim * 0.1)
        xlim[2] <- xlim[2] + (rxlim * 0.1)
    }
    else {
        xlim <- sort(xlim)
    }
    plot(NA, NA, xlim = xlim, ylim = max(sei) - c(ylim[2], ylim[1]), 
        xlab = xlab, ylab = ylab, xaxt = "n", yaxt = "n", bty = "n", 
        ...)
    par.usr <- par("usr")
    rect(par.usr[1], par.usr[3], par.usr[2], par.usr[4], col = back, 
        border = NA, ...)
    axis(side = 2, at = max(sei) - seq(ylim[2], ylim[1], length = steps), 
        labels = formatC(seq(ylim[2], ylim[1], length = steps), 
            digits = digits, format = "f"), ...)
    abline(h = max(sei) - seq(ylim[2], ylim[1], length = steps), 
        col = hlines, ...)
    avals <- length(alpha)
    rylim <- ylim[2] - ylim[1]
    ylim[1] <- max(0, ylim[1] - (rylim * 0.1))
    ylim[2] <- ylim[2] + (rylim * 0.1)
    if (x$method == "FE") {
        for (m in avals:1) {
            x.lb.bot <- refline - qnorm(1 - alpha[m]/2) * sqrt(ylim[2]^2 + 
                tau2)
            x.ub.bot <- refline + qnorm(1 - alpha[m]/2) * sqrt(ylim[2]^2 + 
                tau2)
            x.lb.top <- refline - qnorm(1 - alpha[m]/2) * sqrt(ylim[1]^2 + 
                tau2)
            x.ub.top <- refline + qnorm(1 - alpha[m]/2) * sqrt(ylim[1]^2 + 
                tau2)
            polygon(c(x.lb.bot, x.lb.top, x.ub.top, x.ub.bot), 
                c(max(sei) - ylim[2], max(sei) - ylim[1], max(sei) - 
                  ylim[1], max(sei) - ylim[2]), border = NA, 
                col = shade[m], ...)
            segments(refline, max(sei) - ylim[1], refline, max(sei) - 
                ylim[2], ...)
            segments(x.lb.bot, max(sei) - ylim[2], x.lb.top, 
                max(sei) - ylim[1], lty = "dotted", ...)
            segments(x.ub.bot, max(sei) - ylim[2], x.ub.top, 
                max(sei) - ylim[1], lty = "dotted", ...)
        }
    }
    else {
        for (m in avals:1) {
            vi.vals <- seq(ylim[1]^2, ylim[2]^2, length = 100)
            ci.left <- refline - qnorm(1 - alpha[m]/2) * sqrt(vi.vals + 
                tau2)
            ci.right <- refline + qnorm(1 - alpha[m]/2) * sqrt(vi.vals + 
                tau2)
            lvi <- length(vi.vals)
            polygon(c(ci.left[lvi:1], ci.right), c(max(sei) - 
                sqrt(vi.vals)[lvi:1], max(sei) - sqrt(vi.vals)), 
                border = NA, col = shade[m], ...)
            segments(refline, max(sei), refline, max(sei) - ylim[2], 
                ...)
            lines(ci.left, max(sei) - sqrt(vi.vals), lty = "dotted", 
                ...)
            lines(ci.right, max(sei) - sqrt(vi.vals), lty = "dotted", 
                ...)
        }
    }
    box(bty = "l", ...)
    axis(side = 1, ...)
    points(yi, max(sei) - sei, pch = pch, ...)
    if (is.element("rma.uni.trimfill", class(x))) 
        points(yi[x$fill == 1], max(sei) - sei[x$fill == 1], 
            pch = pch.fill, bg = "white", ...)
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
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail"))) 
        stop("Unknwn 'na.action' specified under options().")
    .invcalc <- function(X, W, k) {
        wX <- sqrt(W) %*% X
        res.qrs <- qr.solve(wX, diag(k))
        res.qrs %*% t(res.qrs)
    }
    x <- model
    if (x$weighted) {
        wi <- 1/(x$vi + x$tau2)
        W <- diag(wi)
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
    if (na.act == "na.omit") {
        hii <- hii[x$not.na]
    }
    if (na.act == "na.fail") 
        stop("Missing values in results.")
    return(hii)
}
influence.rma.uni <-
function (model, digits = model$digits, ...) 
{
    if (!is.element("rma.uni", class(model))) 
        stop("Argument 'model' must be an object of class \"rma.uni\".")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail"))) 
        stop("Unknwn 'na.action' specified under options().")
    x <- model
    tau2.del <- rep(NA, x$k.f)
    delpred <- rep(NA, x$k.f)
    vdelpred <- rep(NA, x$k.f)
    QE.del <- rep(NA, x$k.f)
    dfbetas <- matrix(NA, nrow = x$k.f, ncol = length(x$b))
    cooks.d <- rep(NA, x$k.f)
    covratio <- rep(NA, x$k.f)
    detx <- det(x$vb)
    pred <- x$X.f %*% x$b
    if (x$weighted) {
        wi <- 1/(x$vi + x$tau2)
        W <- diag(wi)
        svb <- crossprod(x$X, W) %*% x$X/x$s2w
    }
    else {
        svb <- solve(x$vb)
    }
    o.warn <- getOption("warn")
    on.exit(options(warn = o.warn))
    options(warn = -1)
    for (i in (1:x$k.f)[x$not.na]) {
        res <- try(rma(x$yi.f[-i], x$vi.f[-i], mods = cbind(x$X.f[-i, 
            ]), method = x$method, weighted = x$weighted, intercept = FALSE, 
            knha = x$knha, control = x$control, ...), silent = TRUE)
        if (class(res) == "try-error") 
            next
        tau2.del[i] <- res$tau2
        Xi <- matrix(x$X.f[i, ], nrow = 1)
        delpred[i] <- Xi %*% res$b
        vdelpred[i] <- Xi %*% tcrossprod(res$vb, Xi)
        QE.del[i] <- res$QE
        dfbeta <- x$b - res$b
        dfbetas[i, ] <- dfbeta/sqrt(diag(res$vb))
        cooks.d[i] <- (crossprod(dfbeta, svb) %*% dfbeta)
        covratio[i] <- det(res$vb)/detx
    }
    delresid <- x$yi.f - delpred
    sedelresid <- sqrt(x$vi.f + vdelpred + tau2.del)
    standelres <- delresid/sedelresid
    dffits <- (pred - delpred)/sqrt(vdelpred)
    options(na.action = "na.exclude")
    hii <- hatvalues(x)
    options(na.action = na.act)
    weight <- rep(NA, x$k.f)
    if (x$weighted) {
        weight[x$not.na] <- wi/sum(wi) * 100
    }
    else {
        weight[x$not.na] <- 1/x$k * 100
    }
    if (na.act == "na.omit") {
        inf <- cbind(standelres[x$not.na], dffits[x$not.na], 
            cooks.d[x$not.na], covratio[x$not.na], tau2.del[x$not.na], 
            QE.del[x$not.na], hii[x$not.na], weight[x$not.na])
        dfb <- cbind(dfbetas[x$not.na, ])
        out <- list(inf = inf, dfb = dfb, tau2 = x$tau2, QE = x$QE, 
            ids = x$ids[x$not.na], not.na = x$not.na[x$not.na], 
            k = x$k, p = x$p, digits = digits)
        dimnames(out$inf)[[1]] <- x$slab[x$not.na]
        dimnames(out$dfb)[[1]] <- x$slab[x$not.na]
    }
    if (na.act == "na.exclude") {
        inf <- cbind(standelres, dffits, cooks.d, covratio, tau2.del, 
            QE.del, hii, weight)
        dfb <- cbind(dfbetas)
        out <- list(inf = inf, dfb = dfb, tau2 = x$tau2, QE = x$QE, 
            ids = x$ids, not.na = x$not.na, k = x$k, p = x$p, 
            digits = digits)
        dimnames(out$inf)[[1]] <- x$slab
        dimnames(out$dfb)[[1]] <- x$slab
    }
    if (na.act == "na.fail") 
        stop("Missing values in results.")
    dimnames(out$dfb)[[2]] <- dimnames(x$b)[[1]]
    dimnames(out$inf)[[2]] <- c("rstudent", "dffits", "cook.d", 
        "cov.r", "tau2.del", "QE.del", "hat", "weight")
    out$inf <- data.frame(out$inf)
    out$dfb <- data.frame(out$dfb)
    class(out) <- "rma.uni.infl"
    return(out)
}
logLik.rma <-
function (object, REML = NULL, ...) 
{
    if (!is.element("rma", class(object))) 
        stop("Argument 'object' must be an object of class \"rma\".")
    if (is.null(REML)) {
        if (object$method == "REML") {
            REML <- TRUE
        }
        else {
            REML <- FALSE
        }
    }
    if (REML) {
        out <- object$fit.stats$REML[1]
        names(out) <- c("ll (REML)")
    }
    else {
        out <- object$fit.stats$ML[1]
        names(out) <- c("ll (ML)")
    }
    return(out)
}
plot.rma.mh <-
function (x, qqplot = FALSE, ...) 
{
    if (!is.element("rma.mh", class(x))) 
        stop("Argument 'x' must be an object of class \"rma.mh\".")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail"))) 
        stop("Unknwn 'na.action' specified under options().")
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
        options(na.action = "na.exclude")
        z <- rstandard(x)$z
        options(na.action = na.act)
        not.na <- !is.na(z)
        if (na.act == "na.omit") {
            z <- z[not.na]
            ids <- x$ids[not.na]
            not.na <- not.na[not.na]
        }
        if (na.act == "na.exclude") 
            ids <- x$ids
        k <- length(z)
        plot(NA, NA, xlim = c(1, k), ylim = c(min(z, -2, na.rm = TRUE), 
            max(z, 2, na.rm = TRUE)), xaxt = "n", xlab = "Study", 
            ylab = "", bty = "l", ...)
        lines((1:k)[not.na], z[not.na], col = "lightgray", ...)
        lines(1:k, z, ...)
        points(1:k, z, pch = 21, bg = "black", ...)
        axis(side = 1, at = 1:k, label = ids, ...)
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
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail"))) 
        stop("Unknwn 'na.action' specified under options().")
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
        options(na.action = "na.exclude")
        z <- rstandard(x)$z
        options(na.action = na.act)
        not.na <- !is.na(z)
        if (na.act == "na.omit") {
            z <- z[not.na]
            ids <- x$ids[not.na]
            not.na <- not.na[not.na]
        }
        if (na.act == "na.exclude") 
            ids <- x$ids
        k <- length(z)
        plot(NA, NA, xlim = c(1, k), ylim = c(min(z, -2, na.rm = TRUE), 
            max(z, 2, na.rm = TRUE)), xaxt = "n", xlab = "Study", 
            ylab = "", bty = "l", ...)
        lines((1:k)[not.na], z[not.na], col = "lightgray", ...)
        lines(1:k, z, ...)
        points(1:k, z, pch = 21, bg = "black", ...)
        axis(side = 1, at = 1:k, label = ids, ...)
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
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail"))) 
        stop("Unknwn 'na.action' specified under options().")
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
            options(na.action = "na.exclude")
            z <- rstandard(x)$z
            options(na.action = na.act)
            not.na <- !is.na(z)
            if (na.act == "na.omit") {
                z <- z[not.na]
                ids <- x$ids[not.na]
                not.na <- not.na[not.na]
            }
            if (na.act == "na.exclude") 
                ids <- x$ids
            k <- length(z)
            plot(NA, NA, xlim = c(1, k), ylim = c(min(z, -2, 
                na.rm = TRUE), max(z, 2, na.rm = TRUE)), xaxt = "n", 
                xlab = "Study", ylab = "", bty = "l", ...)
            lines((1:k)[not.na], z[not.na], col = "lightgray", 
                ...)
            lines(1:k, z, ...)
            points(1:k, z, pch = 21, bg = "black", ...)
            axis(side = 1, at = 1:k, label = ids, ...)
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
        options(na.action = "na.exclude")
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
            options(na.action = "na.exclude")
            z <- rstandard(x)$z
            options(na.action = na.act)
            not.na <- !is.na(z)
            if (na.act == "na.omit") {
                z <- z[not.na]
                ids <- x$ids[not.na]
                not.na <- not.na[not.na]
            }
            if (na.act == "na.exclude") {
                z <- z
                ids <- x$ids
            }
            k <- length(z)
            plot(NA, NA, xlim = c(1, k), ylim = c(min(z, -2, 
                na.rm = TRUE), max(z, 2, na.rm = TRUE)), xaxt = "n", 
                xlab = "Study", ylab = "", bty = "l", ...)
            lines((1:k)[not.na], z[not.na], col = "lightgray", 
                ...)
            lines(1:k, z, ...)
            points(1:k, z, pch = 21, bg = "black", ...)
            axis(side = 1, at = 1:k, label = ids, ...)
            abline(h = 0, lty = "dashed", ...)
            abline(h = c(qnorm(0.025), qnorm(0.975)), lty = "dotted", 
                ...)
            title("Standardized Residuals", ...)
        }
    }
    invisible()
}
plot.rma.uni.infl <-
function (x, plotdfb = FALSE, dfbnew = FALSE, pch = 21, bg = "black", 
    bg.infl = "red", col.na = "lightgray", ...) 
{
    if (class(x) != "rma.uni.infl") 
        stop("Argument 'x' must be an object of class \"rma.uni.infl\".")
    ids <- x$ids
    lids <- length(ids)
    not.na <- x$not.na
    ids.infl <- abs(x$inf$dffits) > 3 * sqrt(x$p/(x$k - x$p)) | 
        pchisq(x$inf$cook.d, df = x$p) > 0.5 | x$inf$hat > 3 * 
        x$p/x$k | apply(abs(x$dfb) > 1, 1, any)
    par.mfrow <- par("mfrow")
    par(mfrow = c(4, 2))
    on.exit(par(mfrow = par.mfrow))
    par.mar <- par("mar")
    par.mar.adj <- par.mar - c(2, 2, 2, 1)
    par.mar.adj[par.mar.adj < 1] <- 1
    par(mar = par.mar.adj)
    on.exit(par(mar = par.mar), add = TRUE)
    plot(NA, NA, xlim = c(1, lids), ylim = c(min(x$inf$rstudent, 
        -2, na.rm = TRUE), max(x$inf$rstudent, 2, na.rm = TRUE)), 
        xaxt = "n", main = "rstudent", xlab = "", ylab = "", 
        ...)
    axis(side = 1, at = 1:lids, label = ids, xlab = "", ...)
    abline(h = 0, lty = "dashed", ...)
    abline(h = c(qnorm(0.025), qnorm(0.975)), lty = "dotted", 
        ...)
    lines((1:lids)[x$not.na], x$inf$rstudent[x$not.na], col = col.na, 
        ...)
    lines(1:lids, x$inf$rstudent, ...)
    points(1:lids, x$inf$rstudent, pch = pch, bg = bg, ...)
    points((1:lids)[ids.infl], x$inf$rstudent[ids.infl], bg = bg.infl, 
        pch = pch, ...)
    plot(NA, NA, xlim = c(1, lids), ylim = range(x$inf$dffits, 
        na.rm = TRUE), xaxt = "n", main = "dffits", xlab = "", 
        ylab = "", ...)
    axis(side = 1, at = 1:lids, label = ids, xlab = "", ...)
    abline(h = 0, lty = "dashed", ...)
    abline(h = 3 * sqrt(x$p/(x$k - x$p)), lty = "dotted", ...)
    abline(h = -3 * sqrt(x$p/(x$k - x$p)), lty = "dotted", ...)
    lines((1:lids)[x$not.na], x$inf$dffits[x$not.na], col = col.na, 
        ...)
    lines(1:lids, x$inf$dffits, ...)
    points(1:lids, x$inf$dffits, pch = pch, bg = bg, ...)
    points((1:lids)[ids.infl], x$inf$dffits[ids.infl], bg = bg.infl, 
        pch = pch, ...)
    plot(NA, NA, xlim = c(1, lids), ylim = range(x$inf$cook.d, 
        na.rm = TRUE), xaxt = "n", main = "cook.d", xlab = "", 
        ylab = "", ...)
    axis(side = 1, at = 1:lids, label = ids, xlab = "", ...)
    abline(h = qchisq(0.5, df = x$p), lty = "dotted", ...)
    lines((1:lids)[x$not.na], x$inf$cook.d[x$not.na], col = col.na, 
        ...)
    lines(1:lids, x$inf$cook.d, ...)
    points(1:lids, x$inf$cook.d, pch = pch, bg = bg, ...)
    points((1:lids)[ids.infl], x$inf$cook.d[ids.infl], bg = bg.infl, 
        pch = pch, ...)
    plot(NA, NA, xlim = c(1, lids), ylim = range(x$inf$cov.r, 
        na.rm = TRUE), xaxt = "n", main = "cov.r", xlab = "", 
        ylab = "", ...)
    axis(side = 1, at = 1:lids, label = ids, xlab = "", ...)
    abline(h = 1, lty = "dashed", ...)
    lines((1:lids)[x$not.na], x$inf$cov.r[x$not.na], col = col.na, 
        ...)
    lines(1:lids, x$inf$cov.r, ...)
    points(1:lids, x$inf$cov.r, pch = pch, bg = bg, ...)
    points((1:lids)[ids.infl], x$inf$cov.r[ids.infl], bg = bg.infl, 
        pch = pch, ...)
    plot(NA, NA, xlim = c(1, lids), ylim = range(x$inf$tau2.del, 
        na.rm = TRUE), xaxt = "n", main = "tau2.del", xlab = "", 
        ylab = "", ...)
    axis(side = 1, at = 1:lids, label = ids, xlab = "", ...)
    abline(h = x$tau2, lty = "dashed", ...)
    lines((1:lids)[x$not.na], x$inf$tau2.del[x$not.na], col = col.na, 
        ...)
    lines(1:lids, x$inf$tau2.del, ...)
    points(1:lids, x$inf$tau2.del, pch = pch, bg = bg, ...)
    points((1:lids)[ids.infl], x$inf$tau2.del[ids.infl], bg = bg.infl, 
        pch = pch, ...)
    plot(NA, NA, xlim = c(1, lids), ylim = range(x$inf$QE.del, 
        na.rm = TRUE), xaxt = "n", main = "QE.del", xlab = "", 
        ylab = "", ...)
    axis(side = 1, at = 1:lids, label = ids, xlab = "", ...)
    abline(h = x$QE, lty = "dashed", ...)
    lines((1:lids)[x$not.na], x$inf$QE.del[x$not.na], col = col.na, 
        ...)
    lines(1:lids, x$inf$QE.del, ...)
    points(1:lids, x$inf$QE.del, pch = pch, bg = bg, ...)
    points((1:lids)[ids.infl], x$inf$QE.del[ids.infl], bg = bg.infl, 
        pch = pch, ...)
    plot(NA, NA, xlim = c(1, lids), ylim = c(0, max(x$inf$hat, 
        na.rm = TRUE)), xaxt = "n", main = "hat", xlab = "", 
        ylab = "", ...)
    axis(side = 1, at = 1:lids, label = ids, xlab = "", ...)
    abline(h = x$p/x$k, lty = "dashed", ...)
    abline(h = 3 * x$p/x$k, lty = "dotted", ...)
    lines((1:lids)[x$not.na], x$inf$hat[x$not.na], col = col.na, 
        ...)
    lines(1:lids, x$inf$hat, ...)
    points(1:lids, x$inf$hat, pch = pch, bg = bg, ...)
    points((1:lids)[ids.infl], x$inf$hat[ids.infl], bg = bg.infl, 
        pch = pch, ...)
    plot(NA, NA, xlim = c(1, lids), ylim = c(0, max(x$inf$weight, 
        na.rm = TRUE)), xaxt = "n", main = "weight", xlab = "", 
        ylab = "", ...)
    axis(side = 1, at = 1:lids, label = ids, xlab = "", ...)
    abline(h = 100/x$k, lty = "dashed", ...)
    lines((1:lids)[x$not.na], x$inf$weight[x$not.na], col = col.na, 
        ...)
    lines(1:lids, x$inf$weight, ...)
    points(1:lids, x$inf$weight, pch = pch, bg = bg, ...)
    points((1:lids)[ids.infl], x$inf$weight[ids.infl], bg = bg.infl, 
        pch = pch, ...)
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
            par.ask <- par("ask")
            par(ask = TRUE)
        }
        par(mfrow = c(x$p, 1))
        for (i in 1:x$p) {
            plot(NA, NA, xlim = c(1, lids), ylim = range(x$dfb[, 
                i], na.rm = TRUE), xaxt = "n", main = paste("dfb: ", 
                dimnames(x$dfb)[[2]][i]), xlab = "", ylab = "", 
                ...)
            axis(side = 1, at = 1:lids, label = ids, xlab = "", 
                ...)
            abline(h = 0, lty = "dashed", ...)
            abline(h = 1, lty = "dotted", ...)
            abline(h = -1, lty = "dotted", ...)
            lines((1:lids)[x$not.na], x$dfb[x$not.na, i], col = col.na, 
                ...)
            lines(1:lids, x$dfb[, i], ...)
            points(1:lids, x$dfb[, i], pch = pch, bg = bg, ...)
            points((1:lids)[ids.infl], x$dfb[ids.infl, i], bg = bg.infl, 
                pch = pch, ...)
        }
        if (!dfbnew) {
            par(ask = par.ask)
        }
    }
    invisible()
}
predict.rma.uni <-
function (object, newmods = NULL, level = object$level, digits = object$digits, 
    transf = FALSE, targs = NULL, ...) 
{
    if (!is.element("rma.uni", class(object))) 
        stop("Argument 'object' must be an object of class \"rma.uni\".")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail"))) 
        stop("Unknwn 'na.action' specified under options().")
    x <- object
    alpha <- (100 - level)/100
    if (x$knha) {
        crit <- qt(1 - alpha/2, df = x$k - x$p)
    }
    else {
        crit <- qnorm(1 - alpha/2)
    }
    if (x$int.only && !is.null(newmods)) 
        stop("Cannot specify new moderator values for models without moderators.")
    if (is.null(newmods)) {
        if (x$int.only) {
            knew <- 1
            Xnew <- cbind(1)
        }
        else {
            knew <- x$k.f
            Xnew <- x$X.f
        }
    }
    else {
        if (x$intercept && x$p == 2) {
            knew <- length(newmods)
            Xnew <- cbind(c(newmods))
        }
        else {
            if (is.vector(newmods) || nrow(newmods) == 1) {
                knew <- 1
                Xnew <- rbind(newmods)
            }
            else {
                knew <- dim(newmods)[1]
                Xnew <- cbind(newmods)
            }
        }
        if (x$intercept) {
            Xnew <- cbind(rep(1, knew), Xnew)
        }
    }
    pred <- rep(NA, knew)
    vpred <- rep(NA, knew)
    for (i in 1:knew) {
        Xinew <- matrix(Xnew[i, ], nrow = 1)
        pred[i] <- Xinew %*% x$b
        vpred[i] <- Xinew %*% tcrossprod(x$vb, Xinew)
    }
    se <- sqrt(vpred)
    ci.lb <- pred - crit * se
    ci.ub <- pred + crit * se
    cr.lb <- pred - crit * sqrt(vpred + x$tau2)
    cr.ub <- pred + crit * sqrt(vpred + x$tau2)
    if (is.function(transf)) {
        if (is.null(targs)) {
            pred <- sapply(pred, transf)
            se <- rep(NA, knew)
            ci.lb <- sapply(ci.lb, transf)
            ci.ub <- sapply(ci.ub, transf)
            cr.lb <- sapply(cr.lb, transf)
            cr.ub <- sapply(cr.ub, transf)
        }
        else {
            pred <- sapply(pred, transf, targs)
            se <- rep(NA, knew)
            ci.lb <- sapply(ci.lb, transf, targs)
            ci.ub <- sapply(ci.ub, transf, targs)
            cr.lb <- sapply(cr.lb, transf, targs)
            cr.ub <- sapply(cr.ub, transf, targs)
        }
    }
    if (is.null(newmods) && !x$int.only) {
        slab <- x$slab
    }
    else {
        slab <- 1:knew
    }
    if (x$int.only) 
        slab <- ""
    if (na.act == "na.omit") {
        not.na <- !is.na(pred)
        out <- list(pred = pred[not.na], se = se[not.na], ci.lb = ci.lb[not.na], 
            ci.ub = ci.ub[not.na], cr.lb = cr.lb[not.na], cr.ub = cr.ub[not.na])
        out$slab <- slab[not.na]
    }
    if (na.act == "na.exclude") {
        out <- list(pred = pred, se = se, ci.lb = ci.lb, ci.ub = ci.ub, 
            cr.lb = cr.lb, cr.ub = cr.ub)
        out$slab <- slab
    }
    if (na.act == "na.fail") 
        stop("Missing values in results.")
    if (x$method == "FE") {
        out$cr.lb <- NULL
        out$cr.ub <- NULL
    }
    out$digits <- digits
    class(out) <- c("rma.list")
    return(out)
}
print.rma.list <-
function (x, digits = x$digits, ...) 
{
    if (class(x) != "rma.list") 
        stop("Argument 'x' must be an object of class \"rma.list\".")
    force(digits)
    attr(x, "class") <- NULL
    out <- x[1:(which(names(x) == "slab") - 1)]
    out <- data.frame(out, row.names = x$slab)
    out <- apply(out, 2, round, digits)
    print(out)
}
print.rma.mh <-
function (x, digits = x$digits, showfit = FALSE, ...) 
{
    if (!is.element("rma.mh", class(x))) 
        stop("Argument 'x' must be an object of class \"rma.mh\".")
    cat("\n")
    cat("Fixed-Effects Model (k = ", x$k, ")", sep = "")
    if (showfit) {
        cat("\n")
        fs <- cbind(round(x$fit.stats$ML, digits = digits))
        dimnames(fs)[[1]] <- c("Log-Likelihood: ", "Deviance (-2LL): ", 
            "AIC: ", "BIC: ")
        dimnames(fs)[[2]] <- "ML"
        cat("\n")
        print(fs)
        cat("\n")
    }
    else {
        cat("\n\n")
    }
    cat("Test for Heterogeneity: \n\n")
    cat("Q(df = ", x$k.yi - 1, ") = ", round(x$QE, digits), ", p-val = ", 
        round(x$QEp, digits), sep = "")
    if (x$measure == "OR" || x$measure == "RR") {
        res.table <- c(x$b, x$se, x$zval, x$pval, x$ci.lb, x$ci.ub)
        res.table.exp <- c(exp(x$b), exp(x$ci.lb), exp(x$ci.ub))
        if (!is.na(x$b)) {
            res.table <- formatC(res.table, digits = digits, 
                format = "f")
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
                cat("Cochran-Mantel-Haenszel Test:     CMH Test not defined for these data \n", 
                  sep = "")
            }
            else {
                cat("Cochran-Mantel-Haenszel Test:     CMH = ", 
                  formatC(x$CMH, digits, format = "f"), ", df = 1, p-val = ", 
                  formatC(x$CMHp, digits, format = "f"), "\n", 
                  sep = "")
            }
            if (is.na(x$TAp)) {
                cat("Tarone's Test for Heterogeneity:  Tarone's Test not defined for these data \n\n", 
                  sep = "")
            }
            else {
                cat("Tarone's Test for Heterogeneity:  X^2 = ", 
                  formatC(x$TA, digits, format = "f"), ", df = ", 
                  x$k.pos - 1, ", p-val = ", formatC(x$TAp, digits, 
                    format = "f"), "\n\n", sep = "")
            }
        }
    }
    else {
        res.table <- c(x$b, x$se, x$zval, x$pval, x$ci.lb, x$ci.ub)
        if (!is.na(x$b)) {
            res.table <- formatC(res.table, digits = digits, 
                format = "f")
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
    cat("\n")
    cat("Fixed-Effects Model (k = ", x$k, ")", sep = "")
    if (showfit) {
        cat("\n")
        fs <- cbind(round(x$fit.stats$ML, digits = digits))
        dimnames(fs)[[1]] <- c("Log-Likelihood: ", "Deviance (-2LL): ", 
            "AIC: ", "BIC: ")
        dimnames(fs)[[2]] <- "ML"
        cat("\n")
        print(fs)
        cat("\n")
    }
    else {
        cat("\n\n")
    }
    cat("Test for Heterogeneity: \n\n")
    cat("Q(df = ", x$k.yi - 1, ") = ", round(x$QE, digits), ", p-val = ", 
        round(x$QEp, digits), sep = "")
    res.table <- c(x$b, x$se, x$zval, x$pval, x$ci.lb, x$ci.ub)
    res.table.exp <- c(exp(x$b), exp(x$ci.lb), exp(x$ci.ub))
    if (!is.na(x$b)) {
        res.table <- formatC(res.table, digits = digits, format = "f")
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
function (x, digits = x$digits, showfit = FALSE, signif.legend = TRUE, 
    ...) 
{
    if (!is.element("rma.uni", class(x))) 
        stop("Argument 'x' must be an object of class \"rma.uni\".")
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
            fs <- cbind(round(x$fit.stats$REML, digits = digits))
            dimnames(fs)[[1]] <- c("Log-Likelihood: ", "Deviance (-2RLL): ", 
                "AIC: ", "BIC: ")
            dimnames(fs)[[2]] <- c("REML")
        }
        else {
            fs <- cbind(round(x$fit.stats$ML, digits = digits))
            dimnames(fs)[[1]] <- c("Log-Likelihood: ", "Deviance (-2LL): ", 
                "AIC: ", "BIC: ")
            dimnames(fs)[[2]] <- c("ML")
        }
        cat("\n")
        print(fs)
        cat("\n")
    }
    else {
        cat("\n\n")
    }
    if (x$method != "FE") {
        if (x$int.only) {
            if (x$method == "ML" || x$method == "REML") {
                cat("tau^2 (estimate of total amount of heterogeneity): ", 
                  round(x$tau2, digits), " (SE = ", round(x$se.tau2, 
                    digits), ")", "\n", sep = "")
            }
            else {
                cat("tau^2 (estimate of total amount of heterogeneity): ", 
                  round(x$tau2, digits), "\n", sep = "")
            }
            cat("tau (sqrt of the estimate of total heterogeneity): ", 
                round(sqrt(ifelse(x$tau2 >= 0, x$tau2, NA)), 
                  digits), "\n", sep = "")
            cat("I^2 (% of total variability due to heterogeneity): ", 
                round(x$I2, 1), "%", sep = "")
        }
        else {
            if (x$method == "ML" || x$method == "REML") {
                cat("tau^2 (estimate of residual amount of heterogeneity): ", 
                  round(x$tau2, digits), " (SE = ", round(x$se.tau2, 
                    digits), ")", "\n", sep = "")
            }
            else {
                cat("tau^2 (estimate of residual amount of heterogeneity): ", 
                  round(x$tau2, digits), "\n", sep = "")
            }
            cat("tau (sqrt of the estimate of residual heterogeneity): ", 
                round(sqrt(ifelse(x$tau2 >= 0, x$tau2, NA)), 
                  digits), sep = "")
        }
        cat("\n\n")
    }
    if (x$int.only) {
        cat("Test for Heterogeneity: \n\n")
        cat("Q(df = ", x$k - x$p, ") = ", round(x$QE, digits), 
            ", p-val = ", round(x$QEp, digits), "\n\n", sep = "")
    }
    else {
        cat("Test for Residual Heterogeneity: \n\n")
        cat("QE(df = ", x$k - x$p, ") = ", round(x$QE, digits), 
            ", p-val = ", round(x$QEp, digits), "\n\n", sep = "")
    }
    if (x$p > 1) {
        cat("Test of Moderators (coefficient(s) ", paste(x$btt, 
            collapse = ","), "): \n\n", sep = "")
        if (!x$knha) {
            cat("QM(df = ", x$m, ") = ", round(x$QM, digits), 
                ", p-val = ", round(x$QMp, digits), "\n\n", sep = "")
        }
        else {
            cat("F(df1 = ", x$m, ", df2 = ", x$k - x$p, ") = ", 
                round(x$QM/x$m, digits), ", p-val = ", round(x$QMp, 
                  digits), "\n\n", sep = "")
        }
    }
    if (x$int.only) {
        res.table <- c(x$b, x$se, x$zval, x$pval, x$ci.lb, x$ci.ub)
        names(res.table) <- c("estimate", "se", "zval", "pval", 
            "ci.lb", "ci.ub")
        if (x$knha) {
            names(res.table)[3] <- c("tval")
        }
        res.table <- formatC(res.table, digits = digits, format = "f")
        signif <- symnum(x$pval, corr = FALSE, na = FALSE, cutpoints = c(0, 
            0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", 
            "*", ".", " "))
        res.table <- c(formatC(res.table, digits = digits, format = "f"), 
            signif)
        names(res.table)[7] <- ""
    }
    else {
        res.table <- cbind(x$b, x$se, x$zval, x$pval, x$ci.lb, 
            x$ci.ub)
        dimnames(res.table)[[2]] <- c("estimate", "se", "zval", 
            "pval", "ci.lb", "ci.ub")
        if (x$knha) {
            dimnames(res.table)[[2]][3] <- c("tval")
        }
        signif <- symnum(x$pval, corr = FALSE, na = FALSE, cutpoints = c(0, 
            0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", 
            "*", ".", " "))
        res.table <- cbind(formatC(res.table, digits = digits, 
            format = "f"), signif)
        dimnames(res.table)[[2]][7] <- ""
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
    if (signif.legend == TRUE) {
        cat("---\nSignif. codes: ", attr(signif, "legend"), "\n\n")
    }
    invisible()
}
print.rma.uni.infl <-
function (x, digits = x$digits, ...) 
{
    if (class(x) != "rma.uni.infl") 
        stop("Argument 'x' must be an object of class \"rma.uni.infl\".")
    x <- list(inf = round(x$inf, digits), dfb = round(x$dfb, 
        digits))
    print(x)
}
qqnorm.rma.mh <-
function (y, type = "rstandard", pch = 19, ...) 
{
    if (!is.element("rma.mh", class(y))) 
        stop("Argument 'y' must be an object of class \"rma.mh\".")
    if (!is.element(type, c("rstandard", "rstudent"))) 
        stop("Argument 'type' must be an either 'rstandard' or 'rstudent'.")
    if (type == "rstandard") {
        res <- rstandard(y)
        not.na <- !is.na(res$z)
        zi <- res$z[not.na]
    }
    else {
        res <- rstudent(y)
        not.na <- !is.na(res$z)
        zi <- res$z[not.na]
    }
    qqnorm(zi, pch = pch, bty = "l", ...)
    abline(a = 0, b = 1, lty = "solid", ...)
    invisible()
}
qqnorm.rma.peto <-
function (y, type = "rstandard", pch = 19, ...) 
{
    if (!is.element("rma.peto", class(y))) 
        stop("Argument 'y' must be an y of class \"rma.peto\".")
    if (!is.element(type, c("rstandard", "rstudent"))) 
        stop("Argument 'type' must be an either 'rstandard' or 'rstudent'.")
    if (type == "rstandard") {
        res <- rstandard(y)
        not.na <- !is.na(res$z)
        zi <- res$z[not.na]
    }
    else {
        res <- rstudent(y)
        not.na <- !is.na(res$z)
        zi <- res$z[not.na]
    }
    qqnorm(zi, pch = pch, bty = "l", ...)
    abline(a = 0, b = 1, lty = "solid", ...)
    invisible()
}
qqnorm.rma.uni <-
function (y, type = "rstandard", pch = 19, envelope = TRUE, level = y$level, 
    reps = 1000, smooth = TRUE, bass = 0, ...) 
{
    if (!is.element("rma.uni", class(y))) 
        stop("Argument 'y' must be an y of class \"rma.uni\".")
    if (!is.element(type, c("rstandard", "rstudent"))) 
        stop("Argument 'type' must be an either 'rstandard' or 'rstudent'.")
    .invcalc <- function(X, W, k) {
        wX <- sqrt(W) %*% X
        res.qrs <- qr.solve(wX, diag(k))
        res.qrs %*% t(res.qrs)
    }
    if (type == "rstandard") {
        res <- rstandard(y)
        not.na <- !is.na(res$z)
        zi <- res$z[not.na]
    }
    else {
        res <- rstudent(y)
        not.na <- !is.na(res$z)
        zi <- res$z[not.na]
    }
    qqnorm(zi, pch = pch, bty = "l", ...)
    abline(a = 0, b = 1, lty = "solid", ...)
    if (envelope) {
        alpha <- (100 - level)/100
        x <- y
        dat <- matrix(rnorm(x$k * reps), nrow = x$k, ncol = reps)
        if (x$weighted) {
            wi <- 1/(x$vi + x$tau2)
            W <- diag(wi)
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
        lb <- apply(ei, 1, quantile, (alpha/2)/x$k)
        ub <- apply(ei, 1, quantile, 1 - (alpha/2)/x$k)
        temp <- qqnorm(lb, plot.it = FALSE)
        if (smooth) 
            temp <- supsmu(temp$x, temp$y, bass = bass)
        lines(temp$x, temp$y, lty = "dotted", ...)
        temp <- qqnorm(ub, plot.it = FALSE)
        if (smooth) 
            temp <- supsmu(temp$x, temp$y, bass = bass)
        lines(temp$x, temp$y, lty = "dotted", ...)
    }
    invisible()
}
radial <-
function (x, ...) 
UseMethod("radial")
radial.rma <-
function (x, center = FALSE, xlim = NULL, zlim = NULL, xlab = NULL, 
    zlab = NULL, atz = NULL, aty = NULL, steps = 7, level = x$level, 
    digits = 2, back = "lightgray", transf = FALSE, targs = NULL, 
    pch = 19, arc.res = 100, cex = NULL, ...) 
{
    if (!is.element("rma", class(x))) 
        stop("Argument 'x' must be an object of class \"rma\".")
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
    if (is.null(xlim)) {
        xlims <- c(0, (1.3 * max(xi)))
    }
    else {
        xlims <- sort(xlim)
    }
    ci.xpos <- xlims[2] + 0.12 * (xlims[2] - xlims[1])
    ya.xpos <- xlims[2] + 0.14 * (xlims[2] - xlims[1])
    xaxismax <- xlims[2]
    if (is.null(zlim)) {
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
    if (is.null(xlab)) {
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
    if (is.null(cex)) 
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
    if (is.null(zlab)) {
        if (center) {
            if (x$method == "FE") {
                mtext(expression(z[i] == frac(y[i] - hat(theta), 
                  v[i])), side = 2, line = par.mar.adj[2] - 1, 
                  at = 0, adj = 0, las = 1, cex = cex, ...)
            }
            else {
                mtext(expression(z[i] == frac(y[i] - hat(mu), 
                  sqrt(v[i] + tau^2))), side = 2, line = par.mar.adj[2] - 
                  1, adj = 0, at = 0, las = 1, cex = cex, ...)
            }
        }
        else {
            if (x$method == "FE") {
                mtext(expression(z[i] == frac(y[i], v[i])), side = 2, 
                  line = par.mar.adj[2] - 2, at = 0, adj = 0, 
                  las = 1, cex = cex, ...)
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
    if (length(arc.res) == 1) 
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
    for (i in 1:length(atyis)) {
        xis[i] <- sqrt(len^2/(1 + (atyis[i]/asp.rat)^2))
        zis[i] <- xis[i] * atyis[i]
    }
    ids <- zis > zlims[1] & zis < zlims[2]
    lines(xis[ids], zis[ids], ...)
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
    for (i in 1:length(atyis)) {
        xis.l[i] <- sqrt(len.l^2/(1 + (atyis[i]/asp.rat)^2))
        zis.l[i] <- xis.l[i] * atyis[i]
        xis.u[i] <- sqrt(len.u^2/(1 + (atyis[i]/asp.rat)^2))
        zis.u[i] <- xis.u[i] * atyis[i]
    }
    ids <- zis.l > zlims[1] & zis.u > zlims[1] & zis.l < zlims[2] & 
        zis.u < zlims[2]
    if (sum(ids) > 0) 
        segments(xis.l[ids], zis.l[ids], xis.u[ids], (xis.u * 
            atyis)[ids], ...)
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
    for (i in 1:length(atyis)) {
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
    ids <- zis > zlims[1] & zis < zlims[2]
    if (sum(ids) > 0) 
        text(xis[ids], zis[ids], formatC(atyis.lab[ids], digits = digits, 
            format = "f"), pos = 4, cex = cex, ...)
    atyis <- seq(ci.lb, ci.ub, length = arc.res[2])
    len <- ci.xpos
    xis <- rep(NA, length(atyis))
    zis <- rep(NA, length(atyis))
    for (i in 1:length(atyis)) {
        xis[i] <- sqrt(len^2/(1 + (atyis[i]/asp.rat)^2))
        zis[i] <- xis[i] * atyis[i]
    }
    ids <- zis > zlims[1] & zis < zlims[2]
    if (sum(ids) > 0) 
        lines(xis[ids], zis[ids], ...)
    atyis <- c(ci.lb, b, ci.ub)
    len.l <- ci.xpos - 0.007 * (xlims[2] - xlims[1])
    len.u <- ci.xpos + 0.007 * (xlims[2] - xlims[1])
    xis.l <- rep(NA, 3)
    zis.l <- rep(NA, 3)
    xis.u <- rep(NA, 3)
    zis.u <- rep(NA, 3)
    for (i in 1:length(atyis)) {
        xis.l[i] <- sqrt(len.l^2/(1 + (atyis[i]/asp.rat)^2))
        zis.l[i] <- xis.l[i] * atyis[i]
        xis.u[i] <- sqrt(len.u^2/(1 + (atyis[i]/asp.rat)^2))
        zis.u[i] <- xis.u[i] * atyis[i]
    }
    ids <- zis.l > zlims[1] & zis.u > zlims[1] & zis.l < zlims[2] & 
        zis.u < zlims[2]
    if (sum(ids) > 0) 
        segments(xis.l[ids], zis.l[ids], xis.u[ids], (xis.u * 
            atyis)[ids], ...)
    par(xpd = par.xpd)
    points(xi, zi, pch = pch, cex = cex, ...)
    invisible()
}
ranktest <-
function (yi, vi, sei, data = NULL, subset = NULL) 
{
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail"))) 
        stop("Unknwn 'na.action' specified under options().")
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
    yi <- eval(mf.yi, data)
    vi <- eval(mf.vi, data)
    sei <- eval(mf.sei, data)
    subset <- eval(mf.subset, data)
    if (missing(vi)) 
        vi <- sei^2
    if (length(yi) != length(vi)) 
        stop("Length of yi and vi (or sei) is not the same.")
    if (!is.null(subset)) {
        yi <- yi[subset]
        vi <- vi[subset]
    }
    yivi.na <- is.na(cbind(yi, vi))
    if (sum(yivi.na) > 0) {
        not.na <- apply(yivi.na, MARGIN = 1, sum) == 0
        if (na.act == "na.omit" || na.act == "na.exclude") {
            yi <- yi[not.na]
            vi <- vi[not.na]
        }
        if (na.act == "na.fail") 
            stop("Missing values in results.")
    }
    res <- rma(yi, vi, method = "FE")
    b <- res$b
    vb <- res$vb
    vi.star <- vi - c(vb)
    yi.star <- (yi - c(b))/sqrt(vi.star)
    cor.test(yi.star, vi, method = "kendall", exact = TRUE)
}
residuals.rma <-
function (object, ...) 
{
    if (!is.element("rma", class(object))) 
        stop("Argument 'object' must be an object of class \"rma\".")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail"))) 
        stop("Unknwn 'na.action' specified under options().")
    x <- object
    out <- c(x$yi.f - x$X.f %*% x$b)
    names(out) <- x$slab
    if (na.act == "na.omit") {
        out <- na.omit(out)
        attr(out, "na.action") <- NULL
    }
    if (na.act == "na.fail") 
        stop("Missing values in results.")
    return(out)
}
rma <-
function (yi, vi, sei, ai, bi, ci, di, n1i, n2i, m1i, m2i, sd1i, 
    sd2i, xi, mi, ri, ni, mods = NULL, data = NULL, intercept = TRUE, 
    slab = NULL, subset = NULL, measure = "GEN", add = 1/2, to = "only0", 
    vtype = "LS", method = "REML", weighted = TRUE, level = 95, 
    digits = 4, btt = NULL, tau2 = NULL, knha = FALSE, control = list()) 
{
    if (!is.element(measure, c("GEN", "MD", "SMD", "RR", "OR", 
        "PETO", "RD", "AS", "PR", "PLN", "PLO", "PAS", "PFT", 
        "COR", "ZCOR"))) 
        stop("Unknown 'measure' specified.")
    if (!is.element(method, c("FE", "HS", "HE", "DL",
        "SJ", "ML", "REML", "EB"))) 
        stop("Unknown 'method' specified.")
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
    slab <- eval(mf.slab, data)
    subset <- eval(mf.subset, data)
    mods <- eval(mf.mods, data)
    if (measure == "GEN") {
        mf.yi <- mf[[match("yi", names(mf))]]
        mf.vi <- mf[[match("vi", names(mf))]]
        mf.sei <- mf[[match("sei", names(mf))]]
        yi <- eval(mf.yi, data)
        vi <- eval(mf.vi, data)
        sei <- eval(mf.sei, data)
        if (is.null(vi)) {
            vi <- sei^2
        }
        if (length(yi) != length(vi)) 
            stop("Length of yi and vi (or sei) is not the same.")
    }
    else {
        if (is.element(measure, c("RR", "OR", "PETO", "RD", "AS"))) {
            mf.ai <- mf[[match("ai", names(mf))]]
            mf.bi <- mf[[match("bi", names(mf))]]
            mf.ci <- mf[[match("ci", names(mf))]]
            mf.di <- mf[[match("di", names(mf))]]
            mf.n1i <- mf[[match("n1i", names(mf))]]
            mf.n2i <- mf[[match("n2i", names(mf))]]
            ai <- eval(mf.ai, data)
            bi <- eval(mf.bi, data)
            ci <- eval(mf.ci, data)
            di <- eval(mf.di, data)
            n1i <- eval(mf.n1i, data)
            n2i <- eval(mf.n2i, data)
            if (is.null(bi)) {
                bi <- n1i - ai
            }
            if (is.null(di)) {
                di <- n2i - ci
            }
            dat <- escalc(measure, ai = ai, bi = bi, ci = ci, 
                di = di, add = add, to = to)
        }
        if (is.element(measure, c("MD", "SMD"))) {
            mf.m1i <- mf[[match("m1i", names(mf))]]
            mf.m2i <- mf[[match("m2i", names(mf))]]
            mf.sd1i <- mf[[match("sd1i", names(mf))]]
            mf.sd2i <- mf[[match("sd2i", names(mf))]]
            mf.n1i <- mf[[match("n1i", names(mf))]]
            mf.n2i <- mf[[match("n2i", names(mf))]]
            m1i <- eval(mf.m1i, data)
            m2i <- eval(mf.m2i, data)
            sd1i <- eval(mf.sd1i, data)
            sd2i <- eval(mf.sd2i, data)
            n1i <- eval(mf.n1i, data)
            n2i <- eval(mf.n2i, data)
            dat <- escalc(measure, m1i = m1i, m2i = m2i, sd1i = sd1i, 
                sd2i = sd2i, n1i = n1i, n2i = n2i, vtype = vtype)
        }
        if (is.element(measure, c("PR", "PLN", "PLO", "PAS", 
            "PFT"))) {
            mf.xi <- mf[[match("xi", names(mf))]]
            mf.mi <- mf[[match("mi", names(mf))]]
            mf.ni <- mf[[match("ni", names(mf))]]
            xi <- eval(mf.xi, data)
            mi <- eval(mf.mi, data)
            ni <- eval(mf.ni, data)
            if (is.null(mi)) {
                mi <- ni - xi
            }
            dat <- escalc(measure, xi = xi, mi = mi, add = add, 
                to = to)
        }
        if (is.element(measure, c("COR", "ZCOR"))) {
            mf.ri <- mf[[match("ri", names(mf))]]
            mf.ni <- mf[[match("ni", names(mf))]]
            ri <- eval(mf.ri, data)
            ni <- eval(mf.ni, data)
            dat <- escalc(measure, ri = ri, ni = ni, vtype = vtype)
        }
        yi <- dat$yi
        vi <- dat$vi
    }
    if (is.vector(mods)) 
        mods <- cbind(mods)
    if (is.data.frame(mods)) 
        mods <- as.matrix(mods)
    k <- length(yi)
    ids <- 1:k
    if (is.null(slab)) {
        slab.null <- TRUE
        slab <- 1:k
    }
    else {
        if (length(slab) != unique(length(slab))) 
            stop("Study labels must be unique.")
        if (length(slab) != length(yi)) 
            stop("Study labels not of same length as data.")
        slab.null <- FALSE
    }
    if (!is.null(subset)) {
        yi <- yi[subset]
        vi <- vi[subset]
        mods <- mods[subset, , drop = FALSE]
        slab <- slab[subset]
        ids <- ids[subset]
        k <- length(yi)
    }
    yi.f <- yi
    vi.f <- vi
    mods.f <- mods
    k.f <- k
    YVM.na <- is.na(cbind(yi, vi, mods))
    if (sum(YVM.na) > 0) {
        na.act <- getOption("na.action")
        if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail"))) 
            stop("Unknwn 'na.action' specified under options().")
        not.na <- apply(YVM.na, MARGIN = 1, sum) == 0
        if (na.act == "na.omit" || na.act == "na.exclude") {
            yi <- yi[not.na]
            vi <- vi[not.na]
            mods <- mods[not.na, , drop = FALSE]
            k <- length(yi)
            warning("Cases with NAs omitted from model fitting.")
        }
        if (na.act == "na.fail") 
            stop("Missing values in data.")
    }
    else {
        not.na <- rep(TRUE, k)
    }
    if (k <= 1) 
        stop("Processing terminated since k <= 1.")
    if (sum(vi <= 0) > 0) {
        allvipos <- FALSE
        vi[vi <= 0] <- 0
        warning("There are outcomes with non-positive sampling variances.")
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
    p <- dim(X)[2]
    if (method == "FE") {
        if (p > k) {
            stop("The number of parameters to be estimated is larger than the number of observations.")
        }
    }
    else {
        if (!is.numeric(tau2)) {
            if (p + 1 > k) {
                stop("The number of parameters to be estimated is larger than the number of observations.")
            }
        }
        else {
            if (p > k) {
                stop("The number of parameters to be estimated is larger than the number of observations.")
            }
        }
    }
    if ((p == 1) && (sum(sapply(X, identical, 1)) == k)) {
        int.only <- TRUE
    }
    else {
        int.only <- FALSE
    }
    if (is.null(btt)) {
        if (p > 1) {
            if (intercept) {
                btt <- 2:p
            }
            else {
                btt <- 1:p
            }
        }
        else {
            btt <- 1
        }
    }
    else {
        btt <- btt[(btt >= 1) & (btt <= p)]
        btt <- unique(round(btt))
        if (length(intersect(btt, 1:p)) == 0) {
            stop("Non-existent coefficients specified with 'btt'.")
        }
    }
    bntt <- setdiff(1:p, btt)
    m <- length(btt)
    tr <- function(X) sum(diag(X))
    .invcalc <- function(X, W, k) {
        wX <- sqrt(W) %*% X
        res.qrs <- qr.solve(wX, diag(k))
        res.qrs %*% t(res.qrs)
    }
    con <- list(tau2.init = NULL, tau2.min = 0, tau2.max = 50, 
        threshold = 10^-5, maxit = 50, verbose = FALSE)
    con[pmatch(names(control), names(con))] <- control
    se.tau2 <- I2 <- QE <- QEp <- NA
    s2w <- 1
    Y <- as.matrix(yi)
    alpha <- (100 - level)/100
    if (!is.numeric(tau2)) {
        if (method == "HS") {
            if (!allvipos) 
                stop("HS estimator cannot be used with non-positive sampling variances.")
            wi <- 1/vi
            W <- diag(wi)
            stXWX <- .invcalc(X = X, W = W, k = k)
            P <- W - W %*% X %*% stXWX %*% crossprod(X, W)
            RSS <- crossprod(Y, P) %*% Y
            tau2 <- RSS/sum(wi) - k/sum(wi)
        }
        if (is.element(method, c("HE", "ML", "REML", "EB"))) {
            stXX <- .invcalc(X = X, W = diag(k), k = k)
            P <- diag(k) - X %*% tcrossprod(stXX, X)
            RSS <- crossprod(Y, P) %*% Y
            tau2 <- (RSS - tr(P %*% diag(vi)))/(k - p)
        }
        if (method == "DL") {
            if (!allvipos) 
                stop("DL estimator cannot be used with non-positive sampling variances.")
            wi <- 1/vi
            W <- diag(wi)
            stXWX <- .invcalc(X = X, W = W, k = k)
            P <- W - W %*% X %*% stXWX %*% crossprod(X, W)
            RSS <- crossprod(Y, P) %*% Y
            tau2 <- (RSS - (k - p))/tr(P)
        }
        if (method == "SJ") {
            if (is.null(con$tau2.init)) {
                tau2 <- var(yi) * (k - 1)/k
            }
            else {
                tau2 <- con$tau2.init
            }
            wi <- 1/(vi + tau2)
            W <- diag(wi)
            stXWX <- .invcalc(X = X, W = W, k = k)
            P <- W - W %*% X %*% stXWX %*% crossprod(X, W)
            RSS <- crossprod(Y, P) %*% Y
            tau2 <- tau2 * RSS/(k - p)
        }
        if (is.element(method, c("ML", "REML", "EB"))) {
            conv <- 1
            change <- con$threshold + 1
            iter <- 0
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
                W <- diag(wi)
                stXWX <- .invcalc(X = X, W = W, k = k)
                P <- W - W %*% X %*% stXWX %*% crossprod(X, W)
                if (method == "REML") {
                  PP <- P %*% P
                  adj <- 1/tr(PP) * (crossprod(Y, PP) %*% Y - 
                    tr(P))
                }
                if (method == "ML") {
                  PP <- P %*% P
                  adj <- 1/sum(wi^2) * (crossprod(Y, PP) %*% 
                    Y - sum(wi))
                }
                if (method == "EB") {
                  adj <- 1/sum(wi) * (crossprod(Y, P) %*% Y * 
                    k/(k - p) - k)
                }
                while (tau2 + adj < con$tau2.min) {
                  adj <- adj/2
                }
                tau2 <- tau2 + adj
                change <- abs(tau2.old - tau2)
                if (iter > con$maxit) {
                  conv <- 0
                  break
                }
            }
            if (conv == 0) 
                stop("Fisher scoring algorithm did not converge. Try increasing the number of iterations (maxit), adjust the threshold (threshold), or use a different estimator for tau^2.")
            if (method == "ML") {
                se.tau2 <- sqrt(2/sum(wi^2))
            }
            if (method == "REML") {
                se.tau2 <- sqrt(2/tr(PP))
            }
        }
        tau2 <- max(con$tau2.min, c(tau2))
    }
    if (method == "FE") {
        tau2 <- 0
        if (!allvipos && weighted) 
            stop("Weighted estimation cannot be used with a fixed-effects\n  model when there are non-positive sampling variances.")
    }
    if (allvipos) {
        wi <- 1/vi
        W <- diag(wi)
        stXWX <- .invcalc(X = X, W = W, k = k)
        P <- W - W %*% X %*% stXWX %*% crossprod(X, W)
        QE <- c(crossprod(Y, P) %*% Y)
        QEp <- 1 - pchisq(QE, df = k - p)
        if (int.only) {
            sumwi <- sum(wi)
            vi.avg <- (k - 1)/(sumwi - sum(wi^2)/sumwi)
            I2 <- 100 * tau2/(vi.avg + tau2)
        }
    }
    wi <- 1/(vi + tau2)
    W <- diag(wi)
    if (weighted) {
        stXWX <- .invcalc(X = X, W = W, k = k)
        b <- stXWX %*% crossprod(X, W) %*% Y
        vb <- stXWX
        P <- W - W %*% X %*% stXWX %*% crossprod(X, W)
        RSS.f <- crossprod(Y, P) %*% Y
        if (knha) {
            s2w <- c(RSS.f)/(k - p)
            vb <- s2w * vb
            if (method == "FE") 
                warning("The Knapp & Hartung (2003) method is not meant to be used in the context of fixed-effects models.")
        }
        if (length(bntt) == 0) {
            QM <- c(sum(wi * yi^2) - RSS.f)/s2w
        }
        else {
            Xr <- X[, bntt, drop = FALSE]
            stXWX <- .invcalc(X = Xr, W = W, k = k)
            P <- W - W %*% Xr %*% stXWX %*% crossprod(Xr, W)
            RSS.r <- crossprod(Y, P) %*% Y
            QM <- c(RSS.r - RSS.f)/s2w
        }
    }
    else {
        stXX <- .invcalc(X = X, W = diag(k), k = k)
        b <- stXX %*% crossprod(X, Y)
        vb <- tcrossprod(stXX, X) %*% diag(vi + tau2) %*% X %*% 
            stXX
        P <- W - W %*% X %*% tcrossprod(stXX, X) - X %*% stXX %*% 
            crossprod(X, W) + X %*% stXX %*% crossprod(X, W) %*% 
            X %*% tcrossprod(stXX, X)
        RSS.f <- crossprod(Y, P) %*% Y
        if (knha) {
            stXWX <- .invcalc(X = X, W = W, k = k)
            P <- W - W %*% X %*% stXWX %*% crossprod(X, W)
            s2w <- c(crossprod(Y, P) %*% Y)/(k - p)
            vb <- s2w * vb
            if (method == "FE") 
                warning("The Knapp & Hartung (2003) method is not meant to be used in the context of fixed-effects models.")
        }
        QM <- t(b)[btt] %*% solve(vb[btt, btt]) %*% b[btt]
    }
    se <- sqrt(diag(vb))
    names(se) <- NULL
    zval <- c(b/se)
    if (knha) {
        QMp <- 1 - pf(QM/m, df1 = m, df2 = k - p)
        pval <- 2 * (1 - pt(abs(zval), df = k - p))
        crit <- qt(1 - alpha/2, df = k - p)
    }
    else {
        QMp <- 1 - pchisq(QM, df = m)
        pval <- 2 * (1 - pnorm(abs(zval)))
        crit <- qnorm(1 - alpha/2)
    }
    ci.lb <- c(b - crit * se)
    ci.ub <- c(b + crit * se)
    ll.ML <- -1/2 * (k) * log(2 * get("pi", pos = "package:base")) - 
        1/2 * sum(log(vi + tau2)) - 1/2 * RSS.f
    ll.REML <- -1/2 * (k - p) * log(2 * get("pi", pos = "package:base")) - 
        1/2 * sum(log(vi + tau2)) - 1/2 * log(det(crossprod(X, 
        W) %*% X)) - 1/2 * RSS.f
    dev.ML <- -2 * ll.ML
    dev.REML <- -2 * ll.REML
    AIC.ML <- -2 * ll.ML + 2 * (p + (if (method == "FE") 
        0
    else 1))
    BIC.ML <- -2 * ll.ML + (p + (if (method == "FE") 
        0
    else 1)) * log(k)
    AIC.REML <- -2 * ll.REML + 2 * (p + (if (method == "FE") 
        0
    else 1))
    BIC.REML <- -2 * ll.REML + (p + (if (method == "FE") 
        0
    else 1)) * log(k - p)
    fit.stats <- matrix(c(ll.ML, dev.ML, AIC.ML, BIC.ML, ll.REML, 
        dev.REML, AIC.REML, BIC.REML), ncol = 2, byrow = FALSE)
    dimnames(fit.stats) <- list(c("ll", "dev", "AIC", "BIC"), 
        c("ML", "REML"))
    fit.stats <- data.frame(fit.stats)
    res <- list(b, se, zval, pval, ci.lb, ci.ub, vb, tau2, se.tau2, 
        k, k.f, p, m, QE, QEp, QM, QMp, I2, int.only, yi, vi, 
        X, yi.f, vi.f, X.f, ids, not.na, slab, slab.null, measure, 
        method, weighted, knha, s2w, btt, intercept, digits, 
        level, con, fit.stats)
    names(res) <- c("b", "se", "zval", "pval", "ci.lb", "ci.ub", 
        "vb", "tau2", "se.tau2", "k", "k.f", "p", "m", "QE", 
        "QEp", "QM", "QMp", "I2", "int.only", "yi", "vi", "X", 
        "yi.f", "vi.f", "X.f", "ids", "not.na", "slab", "slab.null", 
        "measure", "method", "weighted", "knha", "s2w", "btt", 
        "intercept", "digits", "level", "control", "fit.stats")
    class(res) <- c("rma.uni", "rma")
    return(res)
}
rma.mh <-
function (ai, bi, ci, di, n1i, n2i, data = NULL, slab = NULL, 
    subset = NULL, measure = "OR", add = c(1/2, 0), to = c("only0", 
        "none"), level = 95, digits = 4) 
{
    if (!is.element(measure, c("OR", "RR", "RD"))) 
        stop("Mantel-Haenszel method can only be used with measures OR, RR, and RD.")
    if (length(add) != 2) 
        stop("Argument 'add' should specify two values (see 'help(rma.mh)').")
    if (length(to) != 2) 
        stop("Argument 'to' should specify two values (see 'help(rma.mh)').")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail"))) 
        stop("Unknwn 'na.action' specified under options().")
    if (!is.element(to[1], c("all", "only0", "if0all", "none"))) 
        stop("Unknown 'to' argument specified.")
    if (!is.element(to[2], c("all", "only0", "if0all", "none"))) 
        stop("Unknown 'to' argument specified.")
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
    slab <- eval(mf.slab, data)
    subset <- eval(mf.subset, data)
    mf.ai <- mf[[match("ai", names(mf))]]
    mf.bi <- mf[[match("bi", names(mf))]]
    mf.ci <- mf[[match("ci", names(mf))]]
    mf.di <- mf[[match("di", names(mf))]]
    mf.n1i <- mf[[match("n1i", names(mf))]]
    mf.n2i <- mf[[match("n2i", names(mf))]]
    ai <- eval(mf.ai, data)
    bi <- eval(mf.bi, data)
    ci <- eval(mf.ci, data)
    di <- eval(mf.di, data)
    n1i <- eval(mf.n1i, data)
    n2i <- eval(mf.n2i, data)
    if (is.null(bi)) {
        bi <- n1i - ai
    }
    if (is.null(di)) {
        di <- n2i - ci
    }
    k <- length(ai)
    ids <- 1:k
    if (is.null(slab)) {
        slab.null <- TRUE
        slab <- 1:k
    }
    else {
        if (length(slab) != unique(length(slab))) 
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
        slab <- slab[subset]
        ids <- ids[subset]
        k <- length(ai)
    }
    k.f <- k
    dat <- escalc(measure, ai = ai, bi = bi, ci = ci, di = di, 
        add = add[1], to = to[1])
    yi <- dat$yi
    vi <- dat$vi
    ai.f <- ai
    bi.f <- bi
    ci.f <- ci
    di.f <- di
    yi.f <- yi
    vi.f <- vi
    aibicidi.na <- is.na(cbind(ai, bi, ci, di))
    if (sum(aibicidi.na) > 0) {
        not.na <- apply(aibicidi.na, MARGIN = 1, sum) == 0
        if (na.act == "na.omit" || na.act == "na.exclude") {
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
    if (k <= 1) 
        stop("Processing terminated since k <= 1.")
    yivi.na <- is.na(cbind(yi, vi))
    if (sum(yivi.na) > 0) {
        not.na.yivi <- apply(yivi.na, MARGIN = 1, sum) == 0
        if (na.act == "na.omit" || na.act == "na.exclude") {
            yi <- yi[not.na.yivi]
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
    if (to[2] == "all") {
        ai <- ai + add[2]
        bi <- bi + add[2]
        ci <- ci + add[2]
        di <- di + add[2]
    }
    if (to[2] == "only0") {
        id0 <- c(ai == 0 | bi == 0 | ci == 0 | di == 0)
        ai[id0] <- ai[id0] + add[2]
        bi[id0] <- bi[id0] + add[2]
        ci[id0] <- ci[id0] + add[2]
        di[id0] <- di[id0] + add[2]
    }
    if (to[2] == "if0all") {
        id0 <- c(ai == 0 | bi == 0 | ci == 0 | di == 0)
        if (sum(id0) > 0) {
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
            pval <- 2 * (1 - pnorm(abs(zval)))
            ci.lb <- b - qnorm(1 - alpha/2) * se
            ci.ub <- b + qnorm(1 - alpha/2) * se
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
            b <- NA
            se <- NA
            zval <- NA
            pval <- NA
            ci.lb <- NA
            ci.ub <- NA
        }
        else {
            b <- log(sum(ai * (n2i/Ni))/sum(ci * (n1i/Ni)))
            se <- sqrt(sum(((n1i/Ni) * (n2i/Ni) * (ai + ci) - 
                (ai/Ni) * ci))/(R * S))
            zval <- b/se
            pval <- 2 * (1 - pnorm(abs(zval)))
            ci.lb <- b - qnorm(1 - alpha/2) * se
            ci.ub <- b + qnorm(1 - alpha/2) * se
        }
        names(b) <- "intrcpt"
        vb <- matrix(se^2, dimnames = list("intrcpt", "intrcpt"))
        CO <- COp <- CMH <- CMHp <- BD <- BDp <- TA <- TAp <- k.pos <- NA
    }
    if (measure == "RD") {
        b <- sum(ai * (n2i/Ni) - ci * (n1i/Ni))/sum(n1i * (n2i/Ni))
        se <- sqrt(sum(((ai/Ni^2) * bi * (n2i^2/n1i) + (ci/Ni^2) * 
            di * (n1i^2/n2i)))/sum(n1i * (n2i/Ni))^2)
        zval <- b/se
        pval <- 2 * (1 - pnorm(abs(zval)))
        ci.lb <- b - qnorm(1 - alpha/2) * se
        ci.ub <- b + qnorm(1 - alpha/2) * se
        names(b) <- "intrcpt"
        vb <- matrix(se^2, dimnames = list("intrcpt", "intrcpt"))
        CO <- COp <- CMH <- CMHp <- BD <- BDp <- TA <- TAp <- k.pos <- NA
    }
    wi <- 1/vi
    QE <- sum(wi * (yi - b)^2)
    QEp <- 1 - pchisq(QE, df = k.yi - 1)
    ll.ML <- -1/2 * (k.yi) * log(2 * get("pi", pos = "package:base")) - 
        1/2 * sum(log(vi)) - 1/2 * QE
    ll.REML <- -1/2 * (k.yi - 1) * log(2 * get("pi", pos = "package:base")) - 
        1/2 * sum(log(vi)) - 1/2 * log(sum(wi)) - 1/2 * QE
    dev.ML <- -2 * ll.ML
    dev.REML <- -2 * ll.REML
    AIC.ML <- -2 * ll.ML + 2
    BIC.ML <- -2 * ll.ML + log(k.yi)
    AIC.REML <- -2 * ll.REML + 2
    BIC.REML <- -2 * ll.REML + log(k.yi - 1)
    fit.stats <- matrix(c(ll.ML, dev.ML, AIC.ML, BIC.ML, ll.REML, 
        dev.REML, AIC.REML, BIC.REML), ncol = 2, byrow = FALSE)
    dimnames(fit.stats) <- list(c("ll", "dev", "AIC", "BIC"), 
        c("ML", "REML"))
    fit.stats <- data.frame(fit.stats)
    tau2 <- 0
    X.f <- cbind(rep(1, k.f))
    int.only <- TRUE
    method <- "FE"
    res <- list(b, se, zval, pval, ci.lb, ci.ub, vb, tau2, k, 
        k.f, k.yi, k.pos, QE, QEp, CO, COp, CMH, CMHp, BD, BDp, 
        TA, TAp, int.only, yi, vi, yi.f, vi.f, X.f, ai, bi, ci, 
        di, ai.f, bi.f, ci.f, di.f, ids, not.na, not.na.yivi, 
        slab, slab.null, measure, method, digits, level, add, 
        to, fit.stats)
    names(res) <- c("b", "se", "zval", "pval", "ci.lb", "ci.ub", 
        "vb", "tau2", "k", "k.f", "k.yi", "k.pos", "QE", "QEp", 
        "CO", "COp", "CMH", "CMHp", "BD", "BDp", "TA", "TAp", 
        "int.only", "yi", "vi", "yi.f", "vi.f", "X.f", "ai", 
        "bi", "ci", "di", "ai.f", "bi.f", "ci.f", "di.f", "ids", 
        "not.na", "not.na.yivi", "slab", "slab.null", "measure", 
        "method", "digits", "level", "add", "to", "fit.stats")
    class(res) <- c("rma.mh", "rma")
    return(res)
}
rma.peto <-
function (ai, bi, ci, di, n1i, n2i, data = NULL, slab = NULL, 
    subset = NULL, add = c(1/2, 0), to = c("only0", "none"), 
    level = 95, digits = 4) 
{
    if (length(add) != 2) 
        stop("Argument 'add' should specify two values (see 'help(rma.peto)').")
    if (length(to) != 2) 
        stop("Argument 'to' should specify two values (see 'help(rma.peto)').")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail"))) 
        stop("Unknwn 'na.action' specified under options().")
    if (!is.element(to[1], c("all", "only0", "if0all", "none"))) 
        stop("Unknown 'to' argument specified.")
    if (!is.element(to[2], c("all", "only0", "if0all", "none"))) 
        stop("Unknown 'to' argument specified.")
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
    slab <- eval(mf.slab, data)
    subset <- eval(mf.subset, data)
    mf.ai <- mf[[match("ai", names(mf))]]
    mf.bi <- mf[[match("bi", names(mf))]]
    mf.ci <- mf[[match("ci", names(mf))]]
    mf.di <- mf[[match("di", names(mf))]]
    mf.n1i <- mf[[match("n1i", names(mf))]]
    mf.n2i <- mf[[match("n2i", names(mf))]]
    ai <- eval(mf.ai, data)
    bi <- eval(mf.bi, data)
    ci <- eval(mf.ci, data)
    di <- eval(mf.di, data)
    n1i <- eval(mf.n1i, data)
    n2i <- eval(mf.n2i, data)
    if (is.null(bi)) {
        bi <- n1i - ai
    }
    if (is.null(di)) {
        di <- n2i - ci
    }
    k <- length(ai)
    ids <- 1:k
    if (is.null(slab)) {
        slab.null <- TRUE
        slab <- 1:k
    }
    else {
        if (length(slab) != unique(length(slab))) 
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
        slab <- slab[subset]
        ids <- ids[subset]
        k <- length(ai)
    }
    k.f <- k
    dat <- escalc(measure = "PETO", ai = ai, bi = bi, ci = ci, 
        di = di, add = add[1], to = to[1])
    yi <- dat$yi
    vi <- dat$vi
    ai.f <- ai
    bi.f <- bi
    ci.f <- ci
    di.f <- di
    yi.f <- yi
    vi.f <- vi
    aibicidi.na <- is.na(cbind(ai, bi, ci, di))
    if (sum(aibicidi.na) > 0) {
        not.na <- apply(aibicidi.na, MARGIN = 1, sum) == 0
        if (na.act == "na.omit" || na.act == "na.exclude") {
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
    if (k <= 1) 
        stop("Processing terminated since k <= 1.")
    yivi.na <- is.na(cbind(yi, vi))
    if (sum(yivi.na) > 0) {
        not.na.yivi <- apply(yivi.na, MARGIN = 1, sum) == 0
        if (na.act == "na.omit" || na.act == "na.exclude") {
            yi <- yi[not.na.yivi]
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
    if (to[2] == "all") {
        ai <- ai + add[2]
        bi <- bi + add[2]
        ci <- ci + add[2]
        di <- di + add[2]
    }
    if (to[2] == "only0") {
        id0 <- c(ai == 0 | bi == 0 | ci == 0 | di == 0)
        ai[id0] <- ai[id0] + add[2]
        bi[id0] <- bi[id0] + add[2]
        ci[id0] <- ci[id0] + add[2]
        di[id0] <- di[id0] + add[2]
    }
    if (to[2] == "if0all") {
        id0 <- c(ai == 0 | bi == 0 | ci == 0 | di == 0)
        if (sum(id0) > 0) {
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
    if (sumVi == 0) 
        stop("All tables have either only events or no events at all. Peto's method cannot be used.")
    b <- sum(ai - Ei)/sumVi
    se <- sqrt(1/sumVi)
    zval <- b/se
    pval <- 2 * (1 - pnorm(abs(zval)))
    ci.lb <- b - qnorm(1 - alpha/2) * se
    ci.ub <- b + qnorm(1 - alpha/2) * se
    names(b) <- "intrcpt"
    vb <- matrix(se^2, dimnames = list("intrcpt", "intrcpt"))
    wi <- 1/vi
    QE <- sum(wi * (yi - b)^2)
    QEp <- 1 - pchisq(QE, df = k.yi - 1)
    ll.ML <- -1/2 * (k.yi) * log(2 * get("pi", pos = "package:base")) - 
        1/2 * sum(log(vi)) - 1/2 * QE
    ll.REML <- -1/2 * (k.yi - 1) * log(2 * get("pi", pos = "package:base")) - 
        1/2 * sum(log(vi)) - 1/2 * log(sum(wi)) - 1/2 * QE
    dev.ML <- -2 * ll.ML
    dev.REML <- -2 * ll.REML
    AIC.ML <- -2 * ll.ML + 2
    BIC.ML <- -2 * ll.ML + log(k.yi)
    AIC.REML <- -2 * ll.REML + 2
    BIC.REML <- -2 * ll.REML + log(k.yi - 1)
    fit.stats <- matrix(c(ll.ML, dev.ML, AIC.ML, BIC.ML, ll.REML, 
        dev.REML, AIC.REML, BIC.REML), ncol = 2, byrow = FALSE)
    dimnames(fit.stats) <- list(c("ll", "dev", "AIC", "BIC"), 
        c("ML", "REML"))
    fit.stats <- data.frame(fit.stats)
    tau2 <- 0
    X.f <- cbind(rep(1, k.f))
    int.only <- TRUE
    measure <- "PETO"
    method <- "FE"
    res <- list(b, se, zval, pval, ci.lb, ci.ub, vb, tau2, k, 
        k.f, k.yi, QE, QEp, int.only, yi, vi, yi.f, vi.f, X.f, 
        ai, bi, ci, di, ai.f, bi.f, ci.f, di.f, ids, not.na, 
        not.na.yivi, slab, slab.null, measure, method, digits, 
        level, add, to, fit.stats)
    names(res) <- c("b", "se", "zval", "pval", "ci.lb", "ci.ub", 
        "vb", "tau2", "k", "k.f", "k.yi", "QE", "QEp", "int.only", 
        "yi", "vi", "yi.f", "vi.f", "X.f", "ai", "bi", "ci", 
        "di", "ai.f", "bi.f", "ci.f", "di.f", "ids", "not.na", 
        "not.na.yivi", "slab", "slab.null", "measure", "method", 
        "digits", "level", "add", "to", "fit.stats")
    class(res) <- c("rma.peto", "rma")
    return(res)
}
rma.uni <-
function (yi, vi, sei, ai, bi, ci, di, n1i, n2i, m1i, m2i, sd1i, 
    sd2i, xi, mi, ri, ni, mods = NULL, data = NULL, intercept = TRUE, 
    slab = NULL, subset = NULL, measure = "GEN", add = 1/2, to = "only0", 
    vtype = "LS", method = "REML", weighted = TRUE, level = 95, 
    digits = 4, btt = NULL, tau2 = NULL, knha = FALSE, control = list()) 
{
    if (!is.element(measure, c("GEN", "MD", "SMD", "RR", "OR", 
        "PETO", "RD", "AS", "PR", "PLN", "PLO", "PAS", "PFT", 
        "COR", "ZCOR"))) 
        stop("Unknown 'measure' specified.")
    if (!is.element(method, c("FE", "HS", "HE", "DL",
        "SJ", "ML", "REML", "EB"))) 
        stop("Unknown 'method' specified.")
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
    slab <- eval(mf.slab, data)
    subset <- eval(mf.subset, data)
    mods <- eval(mf.mods, data)
    if (measure == "GEN") {
        mf.yi <- mf[[match("yi", names(mf))]]
        mf.vi <- mf[[match("vi", names(mf))]]
        mf.sei <- mf[[match("sei", names(mf))]]
        yi <- eval(mf.yi, data)
        vi <- eval(mf.vi, data)
        sei <- eval(mf.sei, data)
        if (is.null(vi)) {
            vi <- sei^2
        }
        if (length(yi) != length(vi)) 
            stop("Length of yi and vi (or sei) is not the same.")
    }
    else {
        if (is.element(measure, c("RR", "OR", "PETO", "RD", "AS"))) {
            mf.ai <- mf[[match("ai", names(mf))]]
            mf.bi <- mf[[match("bi", names(mf))]]
            mf.ci <- mf[[match("ci", names(mf))]]
            mf.di <- mf[[match("di", names(mf))]]
            mf.n1i <- mf[[match("n1i", names(mf))]]
            mf.n2i <- mf[[match("n2i", names(mf))]]
            ai <- eval(mf.ai, data)
            bi <- eval(mf.bi, data)
            ci <- eval(mf.ci, data)
            di <- eval(mf.di, data)
            n1i <- eval(mf.n1i, data)
            n2i <- eval(mf.n2i, data)
            if (is.null(bi)) {
                bi <- n1i - ai
            }
            if (is.null(di)) {
                di <- n2i - ci
            }
            dat <- escalc(measure, ai = ai, bi = bi, ci = ci, 
                di = di, add = add, to = to)
        }
        if (is.element(measure, c("MD", "SMD"))) {
            mf.m1i <- mf[[match("m1i", names(mf))]]
            mf.m2i <- mf[[match("m2i", names(mf))]]
            mf.sd1i <- mf[[match("sd1i", names(mf))]]
            mf.sd2i <- mf[[match("sd2i", names(mf))]]
            mf.n1i <- mf[[match("n1i", names(mf))]]
            mf.n2i <- mf[[match("n2i", names(mf))]]
            m1i <- eval(mf.m1i, data)
            m2i <- eval(mf.m2i, data)
            sd1i <- eval(mf.sd1i, data)
            sd2i <- eval(mf.sd2i, data)
            n1i <- eval(mf.n1i, data)
            n2i <- eval(mf.n2i, data)
            dat <- escalc(measure, m1i = m1i, m2i = m2i, sd1i = sd1i, 
                sd2i = sd2i, n1i = n1i, n2i = n2i, vtype = vtype)
        }
        if (is.element(measure, c("PR", "PLN", "PLO", "PAS", 
            "PFT"))) {
            mf.xi <- mf[[match("xi", names(mf))]]
            mf.mi <- mf[[match("mi", names(mf))]]
            mf.ni <- mf[[match("ni", names(mf))]]
            xi <- eval(mf.xi, data)
            mi <- eval(mf.mi, data)
            ni <- eval(mf.ni, data)
            if (is.null(mi)) {
                mi <- ni - xi
            }
            dat <- escalc(measure, xi = xi, mi = mi, add = add, 
                to = to)
        }
        if (is.element(measure, c("COR", "ZCOR"))) {
            mf.ri <- mf[[match("ri", names(mf))]]
            mf.ni <- mf[[match("ni", names(mf))]]
            ri <- eval(mf.ri, data)
            ni <- eval(mf.ni, data)
            dat <- escalc(measure, ri = ri, ni = ni, vtype = vtype)
        }
        yi <- dat$yi
        vi <- dat$vi
    }
    if (is.vector(mods)) 
        mods <- cbind(mods)
    if (is.data.frame(mods)) 
        mods <- as.matrix(mods)
    k <- length(yi)
    ids <- 1:k
    if (is.null(slab)) {
        slab.null <- TRUE
        slab <- 1:k
    }
    else {
        if (length(slab) != unique(length(slab))) 
            stop("Study labels must be unique.")
        if (length(slab) != length(yi)) 
            stop("Study labels not of same length as data.")
        slab.null <- FALSE
    }
    if (!is.null(subset)) {
        yi <- yi[subset]
        vi <- vi[subset]
        mods <- mods[subset, , drop = FALSE]
        slab <- slab[subset]
        ids <- ids[subset]
        k <- length(yi)
    }
    yi.f <- yi
    vi.f <- vi
    mods.f <- mods
    k.f <- k
    YVM.na <- is.na(cbind(yi, vi, mods))
    if (sum(YVM.na) > 0) {
        na.act <- getOption("na.action")
        if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail"))) 
            stop("Unknwn 'na.action' specified under options().")
        not.na <- apply(YVM.na, MARGIN = 1, sum) == 0
        if (na.act == "na.omit" || na.act == "na.exclude") {
            yi <- yi[not.na]
            vi <- vi[not.na]
            mods <- mods[not.na, , drop = FALSE]
            k <- length(yi)
            warning("Cases with NAs omitted from model fitting.")
        }
        if (na.act == "na.fail") 
            stop("Missing values in data.")
    }
    else {
        not.na <- rep(TRUE, k)
    }
    if (k <= 1) 
        stop("Processing terminated since k <= 1.")
    if (sum(vi <= 0) > 0) {
        allvipos <- FALSE
        vi[vi <= 0] <- 0
        warning("There are outcomes with non-positive sampling variances.")
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
    p <- dim(X)[2]
    if (method == "FE") {
        if (p > k) {
            stop("The number of parameters to be estimated is larger than the number of observations.")
        }
    }
    else {
        if (!is.numeric(tau2)) {
            if (p + 1 > k) {
                stop("The number of parameters to be estimated is larger than the number of observations.")
            }
        }
        else {
            if (p > k) {
                stop("The number of parameters to be estimated is larger than the number of observations.")
            }
        }
    }
    if ((p == 1) && (sum(sapply(X, identical, 1)) == k)) {
        int.only <- TRUE
    }
    else {
        int.only <- FALSE
    }
    if (is.null(btt)) {
        if (p > 1) {
            if (intercept) {
                btt <- 2:p
            }
            else {
                btt <- 1:p
            }
        }
        else {
            btt <- 1
        }
    }
    else {
        btt <- btt[(btt >= 1) & (btt <= p)]
        btt <- unique(round(btt))
        if (length(intersect(btt, 1:p)) == 0) {
            stop("Non-existent coefficients specified with 'btt'.")
        }
    }
    bntt <- setdiff(1:p, btt)
    m <- length(btt)
    tr <- function(X) sum(diag(X))
    .invcalc <- function(X, W, k) {
        wX <- sqrt(W) %*% X
        res.qrs <- qr.solve(wX, diag(k))
        res.qrs %*% t(res.qrs)
    }
    con <- list(tau2.init = NULL, tau2.min = 0, tau2.max = 50, 
        threshold = 10^-5, maxit = 50, verbose = FALSE)
    con[pmatch(names(control), names(con))] <- control
    se.tau2 <- I2 <- QE <- QEp <- NA
    s2w <- 1
    Y <- as.matrix(yi)
    alpha <- (100 - level)/100
    if (!is.numeric(tau2)) {
        if (method == "HS") {
            if (!allvipos) 
                stop("HS estimator cannot be used with non-positive sampling variances.")
            wi <- 1/vi
            W <- diag(wi)
            stXWX <- .invcalc(X = X, W = W, k = k)
            P <- W - W %*% X %*% stXWX %*% crossprod(X, W)
            RSS <- crossprod(Y, P) %*% Y
            tau2 <- RSS/sum(wi) - k/sum(wi)
        }
        if (is.element(method, c("HE", "ML", "REML", "EB"))) {
            stXX <- .invcalc(X = X, W = diag(k), k = k)
            P <- diag(k) - X %*% tcrossprod(stXX, X)
            RSS <- crossprod(Y, P) %*% Y
            tau2 <- (RSS - tr(P %*% diag(vi)))/(k - p)
        }
        if (method == "DL") {
            if (!allvipos) 
                stop("DL estimator cannot be used with non-positive sampling variances.")
            wi <- 1/vi
            W <- diag(wi)
            stXWX <- .invcalc(X = X, W = W, k = k)
            P <- W - W %*% X %*% stXWX %*% crossprod(X, W)
            RSS <- crossprod(Y, P) %*% Y
            tau2 <- (RSS - (k - p))/tr(P)
        }
        if (method == "SJ") {
            if (is.null(con$tau2.init)) {
                tau2 <- var(yi) * (k - 1)/k
            }
            else {
                tau2 <- con$tau2.init
            }
            wi <- 1/(vi + tau2)
            W <- diag(wi)
            stXWX <- .invcalc(X = X, W = W, k = k)
            P <- W - W %*% X %*% stXWX %*% crossprod(X, W)
            RSS <- crossprod(Y, P) %*% Y
            tau2 <- tau2 * RSS/(k - p)
        }
        if (is.element(method, c("ML", "REML", "EB"))) {
            conv <- 1
            change <- con$threshold + 1
            iter <- 0
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
                W <- diag(wi)
                stXWX <- .invcalc(X = X, W = W, k = k)
                P <- W - W %*% X %*% stXWX %*% crossprod(X, W)
                if (method == "REML") {
                  PP <- P %*% P
                  adj <- 1/tr(PP) * (crossprod(Y, PP) %*% Y - 
                    tr(P))
                }
                if (method == "ML") {
                  PP <- P %*% P
                  adj <- 1/sum(wi^2) * (crossprod(Y, PP) %*% 
                    Y - sum(wi))
                }
                if (method == "EB") {
                  adj <- 1/sum(wi) * (crossprod(Y, P) %*% Y * 
                    k/(k - p) - k)
                }
                while (tau2 + adj < con$tau2.min) {
                  adj <- adj/2
                }
                tau2 <- tau2 + adj
                change <- abs(tau2.old - tau2)
                if (iter > con$maxit) {
                  conv <- 0
                  break
                }
            }
            if (conv == 0) 
                stop("Fisher scoring algorithm did not converge. Try increasing the number of iterations (maxit), adjust the threshold (threshold), or use a different estimator for tau^2.")
            if (method == "ML") {
                se.tau2 <- sqrt(2/sum(wi^2))
            }
            if (method == "REML") {
                se.tau2 <- sqrt(2/tr(PP))
            }
        }
        tau2 <- max(con$tau2.min, c(tau2))
    }
    if (method == "FE") {
        tau2 <- 0
        if (!allvipos && weighted) 
            stop("Weighted estimation cannot be used with a fixed-effects\n  model when there are non-positive sampling variances.")
    }
    if (allvipos) {
        wi <- 1/vi
        W <- diag(wi)
        stXWX <- .invcalc(X = X, W = W, k = k)
        P <- W - W %*% X %*% stXWX %*% crossprod(X, W)
        QE <- c(crossprod(Y, P) %*% Y)
        QEp <- 1 - pchisq(QE, df = k - p)
        if (int.only) {
            sumwi <- sum(wi)
            vi.avg <- (k - 1)/(sumwi - sum(wi^2)/sumwi)
            I2 <- 100 * tau2/(vi.avg + tau2)
        }
    }
    wi <- 1/(vi + tau2)
    W <- diag(wi)
    if (weighted) {
        stXWX <- .invcalc(X = X, W = W, k = k)
        b <- stXWX %*% crossprod(X, W) %*% Y
        vb <- stXWX
        P <- W - W %*% X %*% stXWX %*% crossprod(X, W)
        RSS.f <- crossprod(Y, P) %*% Y
        if (knha) {
            s2w <- c(RSS.f)/(k - p)
            vb <- s2w * vb
            if (method == "FE") 
                warning("The Knapp & Hartung (2003) method is not meant to be used in the context of fixed-effects models.")
        }
        if (length(bntt) == 0) {
            QM <- c(sum(wi * yi^2) - RSS.f)/s2w
        }
        else {
            Xr <- X[, bntt, drop = FALSE]
            stXWX <- .invcalc(X = Xr, W = W, k = k)
            P <- W - W %*% Xr %*% stXWX %*% crossprod(Xr, W)
            RSS.r <- crossprod(Y, P) %*% Y
            QM <- c(RSS.r - RSS.f)/s2w
        }
    }
    else {
        stXX <- .invcalc(X = X, W = diag(k), k = k)
        b <- stXX %*% crossprod(X, Y)
        vb <- tcrossprod(stXX, X) %*% diag(vi + tau2) %*% X %*% 
            stXX
        P <- W - W %*% X %*% tcrossprod(stXX, X) - X %*% stXX %*% 
            crossprod(X, W) + X %*% stXX %*% crossprod(X, W) %*% 
            X %*% tcrossprod(stXX, X)
        RSS.f <- crossprod(Y, P) %*% Y
        if (knha) {
            stXWX <- .invcalc(X = X, W = W, k = k)
            P <- W - W %*% X %*% stXWX %*% crossprod(X, W)
            s2w <- c(crossprod(Y, P) %*% Y)/(k - p)
            vb <- s2w * vb
            if (method == "FE") 
                warning("The Knapp & Hartung (2003) method is not meant to be used in the context of fixed-effects models.")
        }
        QM <- t(b)[btt] %*% solve(vb[btt, btt]) %*% b[btt]
    }
    se <- sqrt(diag(vb))
    names(se) <- NULL
    zval <- c(b/se)
    if (knha) {
        QMp <- 1 - pf(QM/m, df1 = m, df2 = k - p)
        pval <- 2 * (1 - pt(abs(zval), df = k - p))
        crit <- qt(1 - alpha/2, df = k - p)
    }
    else {
        QMp <- 1 - pchisq(QM, df = m)
        pval <- 2 * (1 - pnorm(abs(zval)))
        crit <- qnorm(1 - alpha/2)
    }
    ci.lb <- c(b - crit * se)
    ci.ub <- c(b + crit * se)
    ll.ML <- -1/2 * (k) * log(2 * get("pi", pos = "package:base")) - 
        1/2 * sum(log(vi + tau2)) - 1/2 * RSS.f
    ll.REML <- -1/2 * (k - p) * log(2 * get("pi", pos = "package:base")) - 
        1/2 * sum(log(vi + tau2)) - 1/2 * log(det(crossprod(X, 
        W) %*% X)) - 1/2 * RSS.f
    dev.ML <- -2 * ll.ML
    dev.REML <- -2 * ll.REML
    AIC.ML <- -2 * ll.ML + 2 * (p + (if (method == "FE") 
        0
    else 1))
    BIC.ML <- -2 * ll.ML + (p + (if (method == "FE") 
        0
    else 1)) * log(k)
    AIC.REML <- -2 * ll.REML + 2 * (p + (if (method == "FE") 
        0
    else 1))
    BIC.REML <- -2 * ll.REML + (p + (if (method == "FE") 
        0
    else 1)) * log(k - p)
    fit.stats <- matrix(c(ll.ML, dev.ML, AIC.ML, BIC.ML, ll.REML, 
        dev.REML, AIC.REML, BIC.REML), ncol = 2, byrow = FALSE)
    dimnames(fit.stats) <- list(c("ll", "dev", "AIC", "BIC"), 
        c("ML", "REML"))
    fit.stats <- data.frame(fit.stats)
    res <- list(b, se, zval, pval, ci.lb, ci.ub, vb, tau2, se.tau2, 
        k, k.f, p, m, QE, QEp, QM, QMp, I2, int.only, yi, vi, 
        X, yi.f, vi.f, X.f, ids, not.na, slab, slab.null, measure, 
        method, weighted, knha, s2w, btt, intercept, digits, 
        level, con, fit.stats)
    names(res) <- c("b", "se", "zval", "pval", "ci.lb", "ci.ub", 
        "vb", "tau2", "se.tau2", "k", "k.f", "p", "m", "QE", 
        "QEp", "QM", "QMp", "I2", "int.only", "yi", "vi", "X", 
        "yi.f", "vi.f", "X.f", "ids", "not.na", "slab", "slab.null", 
        "measure", "method", "weighted", "knha", "s2w", "btt", 
        "intercept", "digits", "level", "control", "fit.stats")
    class(res) <- c("rma.uni", "rma")
    return(res)
}
rstandard.rma.mh <-
function (model, digits = model$digits, ...) 
{
    if (!is.element("rma.mh", class(model))) 
        stop("Argument 'model' must be an object of class \"rma.mh\".")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail"))) 
        stop("Unknwn 'na.action' specified under options().")
    x <- model
    e <- c(x$yi.f - x$b)
    e[abs(e) < 100 * .Machine$double.eps * median(abs(e), na.rm = TRUE)] <- 0
    se <- sqrt(x$vi.f)
    z <- e/se
    if (na.act == "na.omit") {
        out <- list(resid = e[x$not.na.yivi], se = se[x$not.na.yivi], 
            z = z[x$not.na.yivi])
        out$slab <- x$slab[x$not.na.yivi]
    }
    if (na.act == "na.exclude") {
        out <- list(resid = e, se = se, z = z)
        out$slab <- x$slab
    }
    if (na.act == "na.fail") 
        stop("Missing values in results.")
    out$digits <- digits
    class(out) <- c("rma.list")
    return(out)
}
rstandard.rma.peto <-
function (model, digits = model$digits, ...) 
{
    if (!is.element("rma.peto", class(model))) 
        stop("Argument 'model' must be an object of class \"rma.peto\".")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail"))) 
        stop("Unknwn 'na.action' specified under options().")
    x <- model
    e <- c(x$yi.f - x$b)
    e[abs(e) < 100 * .Machine$double.eps * median(abs(e), na.rm = TRUE)] <- 0
    se <- sqrt(x$vi.f)
    z <- e/se
    if (na.act == "na.omit") {
        out <- list(resid = e[x$not.na.yivi], se = se[x$not.na.yivi], 
            z = z[x$not.na.yivi])
        out$slab <- x$slab[x$not.na.yivi]
    }
    if (na.act == "na.exclude") {
        out <- list(resid = e, se = se, z = z)
        out$slab <- x$slab
    }
    if (na.act == "na.fail") 
        stop("Missing values in results.")
    out$digits <- digits
    class(out) <- c("rma.list")
    return(out)
}
rstandard.rma.uni <-
function (model, digits = model$digits, ...) 
{
    if (!is.element("rma.uni", class(model))) 
        stop("Argument 'model' must be an object of class \"rma.uni\".")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail"))) 
        stop("Unknwn 'na.action' specified under options().")
    .invcalc <- function(X, W, k) {
        wX <- sqrt(W) %*% X
        res.qrs <- qr.solve(wX, diag(k))
        res.qrs %*% t(res.qrs)
    }
    x <- model
    V <- diag(x$vi + x$tau2)
    if (x$weighted) {
        wi <- 1/(x$vi + x$tau2)
        W <- diag(wi)
        stXWX <- .invcalc(X = x$X, W = W, k = x$k)
        H <- x$X %*% stXWX %*% crossprod(x$X, W)
    }
    else {
        stXX <- .invcalc(X = x$X, W = diag(x$k), k = x$k)
        H <- x$X %*% tcrossprod(stXX, x$X)
    }
    ImH <- diag(x$k) - H
    e <- ImH %*% cbind(x$yi)
    e[abs(e) < 100 * .Machine$double.eps * median(abs(e), na.rm = TRUE)] <- 0
    ve <- ImH %*% tcrossprod(V, ImH)
    se <- sqrt(diag(ve))
    resid <- rep(NA, x$k.f)
    seresid <- rep(NA, x$k.f)
    stanres <- rep(NA, x$k.f)
    resid[x$not.na] <- e
    seresid[x$not.na] <- se
    stanres[x$not.na] <- e/se
    if (na.act == "na.omit") {
        out <- list(resid = resid[x$not.na], se = seresid[x$not.na], 
            z = stanres[x$not.na])
        out$slab <- x$slab[x$not.na]
    }
    if (na.act == "na.exclude") {
        out <- list(resid = resid, se = seresid, z = stanres)
        out$slab <- x$slab
    }
    if (na.act == "na.fail") 
        stop("Missing values in results.")
    out$digits <- digits
    class(out) <- c("rma.list")
    return(out)
}
rstudent.rma.mh <-
function (model, digits = model$digits, ...) 
{
    if (!is.element("rma.mh", class(model))) 
        stop("Argument 'model' must be an object of class \"rma.mh\".")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail"))) 
        stop("Unknwn 'na.action' specified under options().")
    x <- model
    delpred <- rep(NA, x$k.f)
    vdelpred <- rep(NA, x$k.f)
    o.warn <- getOption("warn")
    on.exit(options(warn = o.warn))
    options(warn = -1)
    for (i in (1:x$k.f)[x$not.na]) {
        res <- try(rma.mh(ai = x$ai.f[-i], bi = x$bi.f[-i], ci = x$ci.f[-i], 
            di = x$di.f[-i], measure = x$measure, add = x$add, 
            to = x$to, ...), silent = TRUE)
        if (class(res) == "try-error") 
            next
        delpred[i] <- res$b
        vdelpred[i] <- res$vb
    }
    delresid <- x$yi.f - delpred
    delresid[abs(delresid) < 100 * .Machine$double.eps * median(abs(delresid), 
        na.rm = TRUE)] <- 0
    sedelresid <- sqrt(x$vi.f + vdelpred)
    standelres <- delresid/sedelresid
    if (na.act == "na.omit") {
        out <- list(resid = delresid[x$not.na.yivi], se = sedelresid[x$not.na.yivi], 
            z = standelres[x$not.na.yivi])
        out$slab <- x$slab[x$not.na.yivi]
    }
    if (na.act == "na.exclude") {
        out <- list(resid = delresid, se = sedelresid, z = standelres)
        out$slab <- x$slab
    }
    if (na.act == "na.fail") 
        stop("Missing values in results.")
    out$digits <- digits
    class(out) <- c("rma.list")
    return(out)
}
rstudent.rma.peto <-
function (model, digits = model$digits, ...) 
{
    if (!is.element("rma.peto", class(model))) 
        stop("Argument 'model' must be an object of class \"rma.peto\".")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail"))) 
        stop("Unknwn 'na.action' specified under options().")
    x <- model
    delpred <- rep(NA, x$k.f)
    vdelpred <- rep(NA, x$k.f)
    o.warn <- getOption("warn")
    on.exit(options(warn = o.warn))
    options(warn = -1)
    for (i in (1:x$k.f)[x$not.na]) {
        res <- try(rma.peto(ai = x$ai.f[-i], bi = x$bi.f[-i], 
            ci = x$ci.f[-i], di = x$di.f[-i], add = x$add, to = x$to, 
            ...), silent = TRUE)
        if (class(res) == "try-error") 
            next
        delpred[i] <- res$b
        vdelpred[i] <- res$vb
    }
    delresid <- x$yi.f - delpred
    delresid[abs(delresid) < 100 * .Machine$double.eps * median(abs(delresid), 
        na.rm = TRUE)] <- 0
    sedelresid <- sqrt(x$vi.f + vdelpred)
    standelres <- delresid/sedelresid
    if (na.act == "na.omit") {
        out <- list(resid = delresid[x$not.na.yivi], se = sedelresid[x$not.na.yivi], 
            z = standelres[x$not.na.yivi])
        out$slab <- x$slab[x$not.na.yivi]
    }
    if (na.act == "na.exclude") {
        out <- list(resid = delresid, se = sedelresid, z = standelres)
        out$slab <- x$slab
    }
    if (na.act == "na.fail") 
        stop("Missing values in results.")
    out$digits <- digits
    class(out) <- c("rma.list")
    return(out)
}
rstudent.rma.uni <-
function (model, digits = model$digits, ...) 
{
    if (!is.element("rma.uni", class(model))) 
        stop("Argument 'model' must be an object of class \"rma.uni\".")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail"))) 
        stop("Unknwn 'na.action' specified under options().")
    x <- model
    tau2.del <- rep(NA, x$k.f)
    delpred <- rep(NA, x$k.f)
    vdelpred <- rep(NA, x$k.f)
    o.warn <- getOption("warn")
    on.exit(options(warn = o.warn))
    options(warn = -1)
    for (i in (1:x$k.f)[x$not.na]) {
        res <- try(rma(x$yi.f[-i], x$vi.f[-i], mods = cbind(x$X.f[-i, 
            ]), method = x$method, weighted = x$weighted, intercept = FALSE, 
            knha = x$knha, control = x$control, ...), silent = TRUE)
        if (class(res) == "try-error") 
            next
        tau2.del[i] <- res$tau2
        Xi <- matrix(x$X.f[i, ], nrow = 1)
        delpred[i] <- Xi %*% res$b
        vdelpred[i] <- Xi %*% tcrossprod(res$vb, Xi)
    }
    delresid <- x$yi.f - delpred
    delresid[abs(delresid) < 100 * .Machine$double.eps * median(abs(delresid), 
        na.rm = TRUE)] <- 0
    sedelresid <- sqrt(x$vi.f + vdelpred + tau2.del)
    standelres <- delresid/sedelresid
    if (na.act == "na.omit") {
        out <- list(resid = delresid[x$not.na], se = sedelresid[x$not.na], 
            z = standelres[x$not.na])
        out$slab <- x$slab[x$not.na]
    }
    if (na.act == "na.exclude") {
        out <- list(resid = delresid, se = sedelresid, z = standelres)
        out$slab <- x$slab
    }
    if (na.act == "na.fail") 
        stop("Missing values in results.")
    out$digits <- digits
    class(out) <- c("rma.list")
    return(out)
}
summary.rma <-
function (object, digits = object$digits, showfit = TRUE, signif.legend = TRUE, 
    ...) 
{
    if (!is.element("rma", class(object))) 
        stop("Argument 'object' must be an object of class \"rma\".")
    print(x = object, digits = digits, showfit = showfit, signif.legend = signif.legend, 
        ...)
}
transf.exp.int <-
function (x, targs = NULL, ...) 
{
    if (is.null(targs$tau2)) {
        targs$tau2 <- 0
    }
    if (is.null(targs$lower)) {
        targs$lower <- -5 * sqrt(targs$tau2)
    }
    if (is.null(targs$upper)) {
        targs$lower <- +5 * sqrt(targs$tau2)
    }
    toint <- function(zval, x, tau2) {
        exp(zval) * dnorm(zval, mean = x, sd = sqrt(tau2))
    }
    cfunc <- function(x, tau2, lower, upper) {
        integrate(toint, lower = lower, upper = upper, x = x, 
            tau2 = tau2)$value
    }
    z <- sapply(x, FUN = cfunc, tau2 = targs$tau2, lower = targs$lower, 
        upper = targs$upper)
    return(z)
}
transf.ilogit <-
function (x, ...) 
{
    z <- exp(x)/(1 + exp(x))
    return(z)
}
transf.logit <-
function (x, ...) 
{
    z <- log(x/(1 - x))
    return(z)
}
transf.rtoz <-
function (x, ...) 
{
    z <- 1/2 * log((1 + x)/(1 - x))
    return(z)
}
transf.ztor <-
function (x, ...) 
{
    z <- (exp(2 * x) - 1)/(exp(2 * x) + 1)
    z[x == -Inf] <- -1
    z[x == Inf] <- 1
    return(z)
}
transf.ztor.int <-
function (x, targs = NULL, ...) 
{
    if (is.null(targs$tau2)) {
        targs$tau2 <- 0
    }
    if (is.null(targs$lower)) {
        targs$lower <- -5 * sqrt(targs$tau2)
    }
    if (is.null(targs$upper)) {
        targs$lower <- +5 * sqrt(targs$tau2)
    }
    toint <- function(zval, x, tau2) {
        (exp(2 * zval) - 1)/(exp(2 * zval) + 1) * dnorm(zval, 
            mean = x, sd = sqrt(tau2))
    }
    cfunc <- function(x, tau2, lower, upper) {
        integrate(toint, lower = lower, upper = upper, x = x, 
            tau2 = tau2)$value
    }
    z <- sapply(x, FUN = cfunc, tau2 = targs$tau2, lower = targs$lower, 
        upper = targs$upper)
    return(z)
}
trimfill <-
function (x, ...) 
UseMethod("trimfill")
trimfill.rma.uni <-
function (x, estimator = "L0", side = NULL, maxit = 50, verbose = FALSE, 
    ...) 
{
    if (!is.element("rma.uni", class(x))) 
        stop("Argument 'x' must be an object of class \"rma.uni\".")
    if (!x$int.only) 
        stop("Trim-and-fill method only applicable for models without moderators.")
    if (!is.element(estimator, c("L0", "R0"))) 
        stop("Argument 'estimator' should be either 'L0' or 'R0'.")
    yi <- x$yi
    vi <- x$vi
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
    k <- length(yi)
    k0.sav <- -1
    k0 <- 0
    iter <- 0
    while (abs(k0 - k0.sav) > 0) {
        k0.sav <- k0
        iter <- iter + 1
        if (iter > maxit) 
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
            cat("Iteration:", iter, "\tk0 =", k0, "\t  b =", 
                b, "\n")
        if (k0 <= 0) {
            cat("\nEstimated number of missing studies on the", 
                side, "side is zero.\n\n")
            stop
        }
    }
    if (k0 > 0) {
        if (side == "right") {
            yi <- -1 * yi
            yi.c <- -1 * yi.c
        }
        yi.fill <- c(x$yi.f, -1 * yi.c[(k - k0 + 1):k])
        vi.fill <- c(x$vi.f, vi[(k - k0 + 1):k])
        cat("\nEstimated number of missing studies on the", side, 
            "side:", k0, "\n")
        res <- rma(yi.fill, vi.fill, intercept = TRUE, method = x$method, 
            weighted = x$weighted, ...)
        res$fill <- c(rep(0, k), rep(1, k0))
        class(res) <- c("rma.uni.trimfill", class(res))
        res$ids <- c(x$ids, (x$k.f + 1):(x$k.f + k0))
        if (!x$slab.null) {
            res$slab <- c(x$slab, paste("Filled", 1:k0))
            res$slab.null <- FALSE
        }
        else {
            res$slab <- c(paste("Study", x$ids), paste("Filled", 
                1:k0))
            res$slab.null <- FALSE
        }
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
