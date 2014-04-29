.cmicalc <-
function (mi) 
{
    cmi <- gamma(mi/2)/(sqrt(mi/2) * gamma((mi - 1)/2))
    is.na <- is.na(cmi)
    cmi[is.na] <- 1 - 3/(4 * mi[is.na] - 1)
    return(cmi)
}
.dnchg <-
function (parms, ai, bi, ci, di, X.fit, random, verbose = FALSE, 
    digits = 4, dnchgcalc, dnchgprec, intCtrl) 
{
    p <- ncol(X.fit)
    k <- length(ai)
    b <- parms[seq_len(p)]
    tau2 <- ifelse(random, exp(parms[p + 1]), 0)
    mu.i <- X.fit %*% cbind(b)
    lli <- rep(NA, k)
    if (!random) {
        for (i in seq_len(k)) {
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
        for (i in seq_len(k)) {
            res <- try(integrate(.dnchgi, lower = intCtrl$lower, 
                upper = intCtrl$upper, ai = ai[i], bi = bi[i], 
                ci = ci[i], di = di[i], mu.i = mu.i[i], tau2 = tau2, 
                random = random, dnchgcalc = dnchgcalc, dnchgprec = dnchgprec, 
                rel.tol = intCtrl$rel.tol, subdivisions = intCtrl$subdivisions, 
                stop.on.error = FALSE), silent = !verbose)
            if (inherits(res, "try-error")) {
                stop(paste0("Could not integrate over density of non-central hypergeometric distribution in study ", 
                  i, "."))
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
    for (i in seq_len(k)) {
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
            stop(paste0("Could not compute density of non-central hypergeometric distribution in study ", 
                i, "."))
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
    v <- seq_len(k)
    sub <- function(k, v) {
        if (k == 1L) {
            matrix(v, 1, k)
        }
        else {
            X <- NULL
            for (i in seq_len(k)) {
                X <- rbind(X, cbind(v[i], Recall(k - 1, v[-i])))
            }
            X
        }
    }
    return(sub(k, v[seq_len(k)]))
}
.gensigns <-
function (k) 
{
    ncols <- k
    nrows <- 2^k
    out <- matrix(NA, nrow = nrows, ncol = ncols)
    for (i in seq_len(ncols)) {
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
            for (i in seq_len(len.x)) {
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
    sWX <- sqrt(W) %*% X
    res.qrs <- qr.solve(sWX, diag(k))
    return(tcrossprod(res.qrs))
}
.ll.rma.mv <-
function (par = c(sigma2, tau2, rho), Y, M, X, sigma2, tau2, 
    rho, reml, k, p, s.nvals, t.nvals, r.nvals, withS, withG, 
    D.S, Z.G1, Z.G2, struct, g.levels.r, tol, posdefify, verbose, 
    digits, REMLf, sparse) 
{
    sigma2 <- ifelse(is.na(sigma2), exp(par[1:s.nvals]), sigma2)
    tau2 <- ifelse(is.na(tau2), exp(par[(s.nvals + 1):(s.nvals + 
        t.nvals)]), tau2)
    rho <- ifelse(is.na(rho), transf.ztor(par[(s.nvals + t.nvals + 
        1):(s.nvals + t.nvals + r.nvals)]), rho)
    if (withS) {
        for (j in seq_len(s.nvals)) {
            M <- M + sigma2[j] * D.S[[j]]
        }
    }
    if (withG) {
        ncol.Z.G1 <- ncol(Z.G1)
        if (struct == "CS") {
            G <- matrix(rho, nrow = ncol.Z.G1, ncol = ncol.Z.G1)
            diag(G) <- 1
            G <- diag(sqrt(rep(tau2, ncol.Z.G1)), nrow = ncol.Z.G1, 
                ncol = ncol.Z.G1) %*% G %*% diag(sqrt(rep(tau2, 
                ncol.Z.G1)), nrow = ncol.Z.G1, ncol = ncol.Z.G1)
            diag(G) <- tau2
        }
        if (struct == "HCS") {
            G <- matrix(rho, nrow = ncol.Z.G1, ncol = ncol.Z.G1)
            diag(G) <- 1
            G <- diag(sqrt(tau2), nrow = ncol.Z.G1, ncol = ncol.Z.G1) %*% 
                G %*% diag(sqrt(tau2), nrow = ncol.Z.G1, ncol = ncol.Z.G1)
            diag(G) <- tau2
        }
        if (struct == "UN") {
            G <- matrix(NA, nrow = ncol.Z.G1, ncol = ncol.Z.G1)
            G[upper.tri(G)] <- rho
            G[lower.tri(G)] <- t(G)[lower.tri(G)]
            diag(G) <- 1
            G <- diag(sqrt(tau2), nrow = ncol.Z.G1, ncol = ncol.Z.G1) %*% 
                G %*% diag(sqrt(tau2), nrow = ncol.Z.G1, ncol = ncol.Z.G1)
            diag(G) <- tau2
            if (posdefify) 
                G <- as.matrix(nearPD(G)$mat)
        }
        if (struct == "UNHO") {
            G <- matrix(NA, nrow = ncol.Z.G1, ncol = ncol.Z.G1)
            G[upper.tri(G)] <- rho
            G[lower.tri(G)] <- t(G)[lower.tri(G)]
            diag(G) <- 1
            G <- diag(sqrt(rep(tau2, ncol.Z.G1)), nrow = ncol.Z.G1, 
                ncol = ncol.Z.G1) %*% G %*% diag(sqrt(rep(tau2, 
                ncol.Z.G1)), nrow = ncol.Z.G1, ncol = ncol.Z.G1)
            diag(G) <- tau2
            if (posdefify) 
                G <- as.matrix(nearPD(G, keepDiag = TRUE)$mat)
        }
        if (struct == "AR") {
            if (ncol.Z.G1 > 1) {
                G <- toeplitz(ARMAacf(ar = rho, lag.max = ncol.Z.G1 - 
                  1))
            }
            else {
                G <- diag(1)
            }
            G <- diag(sqrt(rep(tau2, ncol.Z.G1)), nrow = ncol.Z.G1, 
                ncol = ncol.Z.G1) %*% G %*% diag(sqrt(rep(tau2, 
                ncol.Z.G1)), nrow = ncol.Z.G1, ncol = ncol.Z.G1)
            diag(G) <- tau2
        }
        if (struct == "HAR") {
            if (ncol.Z.G1 > 1) {
                G <- toeplitz(ARMAacf(ar = rho, lag.max = ncol.Z.G1 - 
                  1))
            }
            else {
                G <- diag(1)
            }
            G <- diag(sqrt(tau2), nrow = ncol.Z.G1, ncol = ncol.Z.G1) %*% 
                G %*% diag(sqrt(tau2), nrow = ncol.Z.G1, ncol = ncol.Z.G1)
            diag(G) <- tau2
        }
        if (any(g.levels.r)) {
            G[g.levels.r, ] <- 0
            G[, g.levels.r] <- 0
        }
        if (sparse) 
            G <- Matrix(G, sparse = TRUE)
        M <- M + (Z.G1 %*% G %*% t(Z.G1)) * tcrossprod(Z.G2)
    }
    if (verbose) {
        L <- try(chol(M), silent = !verbose)
    }
    else {
        L <- suppressWarnings(try(chol(M), silent = !verbose))
    }
    if (inherits(L, "try-error")) {
        llval <- -Inf
    }
    else {
        W <- chol2inv(L)
        U <- chol(W)
        sX <- U %*% X
        sY <- U %*% Y
        b <- solve(crossprod(sX), crossprod(sX, sY))
        RSS.f <- sum(as.vector(sY - sX %*% b)^2)
        if (reml) {
            llval <- -1/2 * (k - p) * log(2 * base::pi) + ifelse(REMLf, 
                1/2 * determinant(crossprod(X), logarithm = TRUE)$modulus, 
                0) - 1/2 * determinant(M, logarithm = TRUE)$modulus - 
                1/2 * determinant(crossprod(X, W) %*% X, logarithm = TRUE)$modulus - 
                1/2 * RSS.f
        }
        else {
            llval <- -1/2 * (k) * log(2 * base::pi) - 1/2 * determinant(M, 
                logarithm = TRUE)$modulus - 1/2 * RSS.f
        }
    }
    if (verbose) {
        if (withS) 
            cat("sigma2 =", ifelse(is.na(sigma2), NA, paste(formatC(sigma2, 
                digits = digits, format = "f", flag = " "), " ", 
                sep = "")), "  ", sep = "")
        if (withG) 
            cat("tau2 =", ifelse(is.na(tau2), NA, paste(formatC(tau2, 
                digits = digits, format = "f", flag = " "), " ", 
                sep = "")), "  ", sep = "")
        if (withG) 
            cat("rho =", ifelse(is.na(rho), NA, paste(formatC(rho, 
                digits = digits, format = "f", flag = " "), " ", 
                sep = "")), "  ", sep = "")
        cat("  ll = ", ifelse(is.na(llval), NA, formatC(llval, 
            digits = digits, format = "f", flag = " ")), sep = "", 
            "\n")
    }
    llval <- -1 * llval
    llval
}
.modfit <-
function (Y, X, V, k, tol = 1e-07) 
{
    eV <- eigen(V, symmetric = TRUE)
    d <- eV$values
    if (any(d <= tol)) {
        stop("Var-cov matrix is not positive definite.")
    }
    A <- diag(1/sqrt(d), nrow = k, ncol = k) %*% t(eV$vectors)
    Ainv <- eV$vectors %*% diag(1/sqrt(d), nrow = k, ncol = k)
    AX <- A %*% X
    res.qrs <- qr.solve(AX, diag(k))
    vb <- tcrossprod(res.qrs)
    b <- qr.solve(AX, A %*% Y)
    res <- list(b = b, vb = vb)
    return(res)
}
.onAttach <-
function (libname, pkgname) 
{
    loadmsg <- "\nLoading 'metafor' package (version 1.9-3). For an overview \nand introduction to the package please type: help(metafor)."
    packageStartupMessage(loadmsg, domain = NULL, appendLF = TRUE)
}
.QE.func <-
function (tau2val, Y, vi, X, k, objective, verbose = FALSE, digits = 4) 
{
    wi <- 1/(vi + tau2val)
    W <- diag(wi, nrow = k, ncol = k)
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
    if (is.element(measure, c("OR", "PETO", "D2OR", "D2ORN", 
        "D2ORL"))) {
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
    if (is.element(measure, c("SMD", "SMDH", "PBIT", "OR2D", 
        "OR2DN", "OR2DL"))) {
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
    if (is.element(measure, c("COR", "UCOR", "RTET", "RBIS"))) {
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
    if (is.element(measure, c("SMCC", "SMCR"))) {
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
    yi.names <- attr(x, "yi.names")
    if (!is.null(yi.names)) {
        for (i in 1:length(yi.names)) {
            if (!is.element(yi.names[i], names(dat))) {
                next
            }
            else {
                eval(parse(text = paste0("attr(dat$", yi.names[i], 
                  ", 'measure') <- attr(x$", yi.names[i], ", 'measure')")))
            }
        }
    }
    attr(dat, "yi.names") <- attr(x, "yi.names")
    attr(dat, "vi.names") <- attr(x, "vi.names")
    attr(dat, "sei.names") <- attr(x, "sei.names")
    attr(dat, "zi.names") <- attr(x, "zi.names")
    attr(dat, "ci.lb.names") <- attr(x, "ci.lb.names")
    attr(dat, "ci.ub.names") <- attr(x, "ci.ub.names")
    if (!is.null(attr(x, "digits"))) 
        attr(dat, "digits") <- attr(x, "digits")
    return(dat)
}
`[.list.rma` <-
function (x, i, ...) 
{
    out <- x
    attr(out, "class") <- NULL
    slab.pos <- which(names(out) == "slab")
    if (!missing(i)) 
        out[seq_len(slab.pos - 1)] <- lapply(out[seq_len(slab.pos - 
            1)], function(r) if (class(r) == "matrix") 
            r[i, ]
        else r[i])
    if (length(out[[1]]) == 0L) 
        return(NULL)
    out$slab <- x$slab[i]
    if (any(is.na(out$slab))) 
        return(NULL)
    out$digits <- x$digits
    out$transf <- x$transf
    out$method <- x$method
    class(out) <- "list.rma"
    return(out)
}
addpoly <-
function (x, ...) 
UseMethod("addpoly")
addpoly.default <-
function (x, vi, sei, ci.lb, ci.ub, rows = -1, level = 95, digits = 2, 
    annotate = TRUE, mlab, transf = FALSE, atransf = FALSE, targs, 
    efac = 1, col, border, cex, ...) 
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
    if (missing(col)) 
        col <- "black"
    if (missing(border)) 
        border <- "black"
    alpha <- ifelse(level > 1, (100 - level)/100, 1 - level)
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
            for (j in seq_len(length(rows.na))) {
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
    ci.bounds[rev.order, ] <- ci.bounds[rev.order, 2:1]
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
    for (i in seq_len(k)) {
        polygon(x = c(ci.lb[i], yi[i], ci.ub[i], yi[i]), y = c(rows[i], 
            rows[i] + (height/100) * cex * efac, rows[i], rows[i] - 
                (height/100) * cex * efac), col = col, border = border, 
            ...)
        if (!is.null(mlab)) {
            text(xlim[1], rows[i], mlab[i], pos = 4, cex = cex, 
                ...)
        }
    }
}
addpoly.rma <-
function (x, row = -2, level = x$level, digits = 2, annotate = TRUE, 
    mlab, transf = FALSE, atransf = FALSE, targs, efac = 1, col, 
    border, cex, ...) 
{
    if (!is.element("rma", class(x))) 
        stop("Argument 'x' must be an object of class \"rma\".")
    if (!x$int.only) 
        stop("Fitted model should not contain moderators.")
    if (missing(mlab)) 
        mlab <- NULL
    if (missing(targs)) 
        targs <- NULL
    if (missing(cex)) 
        cex <- NULL
    if (missing(col)) 
        col <- "black"
    if (missing(border)) 
        border <- "black"
    if (is.null(mlab)) 
        mlab <- ifelse((x$method == "FE"), "FE Model", "RE Model")
    addpoly(x$b, ci.lb = x$ci.lb, ci.ub = x$ci.ub, rows = row, 
        level = level, digits = digits, annotate = annotate, 
        mlab = mlab, transf = transf, atransf = atransf, targs = targs, 
        efac = efac, col = col, border = border, cex = cex, ...)
}
AIC.rma <-
function (object, ..., k = 2, correct = FALSE) 
{
    if (!is.element("rma", class(object))) 
        stop("Argument 'object' must be an object of class \"rma\".")
    if (object$method == "REML") {
        ifelse(correct, object$fit.stats$REML[5], object$fit.stats$REML[3])
    }
    else {
        ifelse(correct, object$fit.stats$ML[5], object$fit.stats$ML[3])
    }
}
anova.rma.uni <-
function (object, object2, btt, digits, ...) 
{
    if (!is.element("rma.uni", class(object))) 
        stop("Argument 'object' must be an object of class \"rma.uni\".")
    if (missing(digits)) 
        digits <- object$digits
    if (missing(object2)) {
        if (missing(btt)) 
            btt <- NULL
        x <- object
        k <- x$k
        p <- x$p
        b <- x$b
        vb <- x$vb
        if (is.null(btt)) {
            if (p > 1) {
                if (x$int.incl) {
                  btt <- seq.int(from = 2, to = p)
                }
                else {
                  btt <- seq_len(p)
                }
            }
            else {
                btt <- 1
            }
        }
        else {
            btt <- btt[(btt >= 1) & (btt <= p)]
            btt <- unique(round(btt))
            if (length(intersect(btt, seq_len(p))) == 0L) {
                stop("Non-existent coefficients specified with 'btt'.")
            }
        }
        bntt <- setdiff(seq_len(p), btt)
        m <- length(btt)
        QM <- c(t(b)[btt] %*% chol2inv(chol(vb[btt, btt])) %*% 
            b[btt])
        if (x$knha || x$robust) {
            QM <- QM/m
            QMp <- pf(QM, df1 = m, df2 = k - p, lower.tail = FALSE)
        }
        else {
            QMp <- pchisq(QM, df = m, lower.tail = FALSE)
        }
        test <- "Wald"
        res <- list(QM = QM, QMp = QMp, btt = btt, k = k, p = p, 
            m = m, knha = x$knha, robust = x$robust, digits = digits, 
            test = test)
        class(res) <- c("anova.rma.uni")
        return(res)
    }
    else {
        if (!is.element("rma.uni", class(object2))) 
            stop("Argument 'object2' must be an object of class \"rma.uni\".")
        m.f <- object
        m.r <- object2
        if (!(identical(c(m.f$yi), c(m.r$yi)) && identical(c(m.f$vi), 
            c(m.r$vi)))) 
            stop("Observed outcomes and/or sampling variances not equal in the full and reduced model.")
        p.f <- m.f$parms
        p.r <- m.r$parms
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
            R2 <- NA
        }
        else {
            R2 <- round(100 * max(0, (m.r$tau2 - m.f$tau2)/m.r$tau2), 
                2)
        }
        test <- "LRT"
        res <- list(fit.stats.f = fit.stats.f, fit.stats.r = fit.stats.r, 
            p.f = p.f, p.r = p.r, LRT = LRT, pval = pval, QE.f = m.f$QE, 
            QE.r = m.r$QE, tau2.f = m.f$tau2, tau2.r = m.r$tau2, 
            R2 = R2, method = m.f$method, digits = digits, test = test)
        class(res) <- c("anova.rma.uni")
        return(res)
    }
}
baujat <-
function (x, ...) 
UseMethod("baujat")
baujat.rma.mh <-
function (x, xlim, ylim, xlab, ylab, cex, grid = TRUE, ...) 
{
    if (!is.element("rma.mh", class(x))) 
        stop("Argument 'x' must be an object of class \"rma.mh\".")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    if (x$k == 1) 
        stop("Stopped because k = 1.")
    delpred <- rep(NA, x$k.f)
    vdelpred <- rep(NA, x$k.f)
    pred.full <- x$X.f %*% x$b
    o.warn <- getOption("warn")
    on.exit(options(warn = o.warn))
    options(warn = -1)
    for (i in seq_len(x$k.f)[x$not.na]) {
        if (is.element(x$measure, c("RR", "OR", "RD"))) {
            res <- try(rma.mh(ai = x$ai.f[-i], bi = x$bi.f[-i], 
                ci = x$ci.f[-i], di = x$di.f[-i], measure = x$measure, 
                add = x$add, to = x$to, drop00 = x$drop00, correct = x$correct), 
                silent = TRUE)
        }
        else {
            res <- try(rma.mh(x1i = x$x1i.f[-i], x2i = x$x2i.f[-i], 
                t1i = x$t1i.f[-i], t2i = x$t2i.f[-i], measure = x$measure, 
                add = x$add, to = x$to, drop00 = x$drop00, correct = x$correct), 
                silent = TRUE)
        }
        if (inherits(res, "try-error")) 
            next
        Xi <- matrix(x$X.f[i, ], nrow = 1)
        delpred[i] <- Xi %*% res$b
        vdelpred[i] <- Xi %*% tcrossprod(res$vb, Xi)
    }
    yhati <- (delpred - pred.full)^2/vdelpred
    options(na.action = "na.pass")
    xhati <- 1/(x$tau2 + x$vi.f) * resid(x)^2
    options(na.action = na.act)
    if (missing(cex)) 
        cex <- 0.8
    if (missing(xlab)) {
        if (x$method == "FE") {
            xlab <- ifelse(x$int.only, "Contribution to Overall Heterogeneity", 
                "Contribution to Residual Heterogeneity")
        }
        else {
            xlab <- "Squared Pearson Residual"
        }
    }
    if (missing(ylab)) 
        ylab <- ifelse(x$int.only, "Influence on Overall Result", 
            "Influence on Fitted Value")
    if (missing(xlim)) 
        xlim <- range(xhati, na.rm = TRUE)
    if (missing(ylim)) 
        ylim <- range(yhati, na.rm = TRUE)
    plot(xhati, yhati, pch = 19, col = "white", xlab = xlab, 
        ylab = ylab, cex = cex, xlim = xlim, ylim = ylim, ...)
    if (grid) 
        grid()
    text(xhati, yhati, x$slab, cex = cex, ...)
    invisible(data.frame(x = xhati, y = yhati))
}
baujat.rma.peto <-
function (x, xlim, ylim, xlab, ylab, cex, grid = TRUE, ...) 
{
    if (!is.element("rma.peto", class(x))) 
        stop("Argument 'x' must be an object of class \"rma.peto\".")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    if (x$k == 1) 
        stop("Stopped because k = 1.")
    delpred <- rep(NA, x$k.f)
    vdelpred <- rep(NA, x$k.f)
    pred.full <- x$X.f %*% x$b
    o.warn <- getOption("warn")
    on.exit(options(warn = o.warn))
    options(warn = -1)
    for (i in seq_len(x$k.f)[x$not.na]) {
        res <- try(rma.peto(ai = x$ai.f[-i], bi = x$bi.f[-i], 
            ci = x$ci.f[-i], di = x$di.f[-i], add = x$add, to = x$to, 
            drop00 = x$drop00), silent = TRUE)
        if (inherits(res, "try-error")) 
            next
        Xi <- matrix(x$X.f[i, ], nrow = 1)
        delpred[i] <- Xi %*% res$b
        vdelpred[i] <- Xi %*% tcrossprod(res$vb, Xi)
    }
    yhati <- (delpred - pred.full)^2/vdelpred
    options(na.action = "na.pass")
    xhati <- 1/(x$tau2 + x$vi.f) * resid(x)^2
    options(na.action = na.act)
    if (missing(cex)) 
        cex <- 0.8
    if (missing(xlab)) {
        if (x$method == "FE") {
            xlab <- ifelse(x$int.only, "Contribution to Overall Heterogeneity", 
                "Contribution to Residual Heterogeneity")
        }
        else {
            xlab <- "Squared Pearson Residual"
        }
    }
    if (missing(ylab)) 
        ylab <- ifelse(x$int.only, "Influence on Overall Result", 
            "Influence on Fitted Value")
    if (missing(xlim)) 
        xlim <- range(xhati, na.rm = TRUE)
    if (missing(ylim)) 
        ylim <- range(yhati, na.rm = TRUE)
    plot(xhati, yhati, pch = 19, col = "white", xlab = xlab, 
        ylab = ylab, cex = cex, xlim = xlim, ylim = ylim, ...)
    if (grid) 
        grid()
    text(xhati, yhati, x$slab, cex = cex, ...)
    invisible(data.frame(x = xhati, y = yhati))
}
baujat.rma.uni <-
function (x, xlim, ylim, xlab, ylab, cex, grid = TRUE, ...) 
{
    if (!is.element("rma.uni", class(x))) 
        stop("Argument 'x' must be an object of class \"rma.uni\".")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    if (x$k == 1) 
        stop("Stopped because k = 1.")
    delpred <- rep(NA, x$k.f)
    vdelpred <- rep(NA, x$k.f)
    pred.full <- x$X.f %*% x$b
    o.warn <- getOption("warn")
    on.exit(options(warn = o.warn))
    options(warn = -1)
    for (i in seq_len(x$k.f)[x$not.na]) {
        res <- try(rma(x$yi.f[-i], x$vi.f[-i], weights = x$weights.f[-i], 
            mods = cbind(x$X.f[-i, ]), method = x$method, weighted = x$weighted, 
            intercept = FALSE, knha = x$knha, control = x$control), 
            silent = TRUE)
        if (inherits(res, "try-error")) 
            next
        Xi <- matrix(x$X.f[i, ], nrow = 1)
        delpred[i] <- Xi %*% res$b
        vdelpred[i] <- Xi %*% tcrossprod(res$vb, Xi)
    }
    yhati <- (delpred - pred.full)^2/vdelpred
    options(na.action = "na.pass")
    xhati <- 1/(x$tau2 + x$vi.f) * resid(x)^2
    options(na.action = na.act)
    if (missing(cex)) 
        cex <- 0.8
    if (missing(xlab)) {
        if (x$method == "FE") {
            xlab <- ifelse(x$int.only, "Contribution to Overall Heterogeneity", 
                "Contribution to Residual Heterogeneity")
        }
        else {
            xlab <- "Squared Pearson Residual"
        }
    }
    if (missing(ylab)) 
        ylab <- ifelse(x$int.only, "Influence on Overall Result", 
            "Influence on Fitted Value")
    if (missing(xlim)) 
        xlim <- range(xhati, na.rm = TRUE)
    if (missing(ylim)) 
        ylim <- range(yhati, na.rm = TRUE)
    plot(xhati, yhati, pch = 19, col = "white", xlab = xlab, 
        ylab = ylab, cex = cex, xlim = xlim, ylim = ylim, ...)
    if (grid) 
        grid()
    text(xhati, yhati, x$slab, cex = cex, ...)
    invisible(data.frame(x = xhati, y = yhati))
}
BIC.rma <-
function (object, ...) 
{
    if (!is.element("rma", class(object))) 
        stop("Argument 'object' must be an object of class \"rma\".")
    if (object$method == "REML") {
        object$fit.stats$REML[4]
    }
    else {
        object$fit.stats$ML[4]
    }
}
bldiag <-
function (...) 
{
    mlist <- list(...)
    if (length(mlist) == 1) 
        mlist <- unlist(mlist, recursive = FALSE)
    csdim <- rbind(c(0, 0), apply(sapply(mlist, dim), 1, cumsum))
    out <- array(0, dim = csdim[length(mlist) + 1, ])
    add1 <- matrix(rep(1:0, 2), ncol = 2)
    for (i in seq(along = mlist)) {
        indx <- apply(csdim[i:(i + 1), ] + add1, 2, function(x) x[1]:x[2])
        if (is.null(dim(indx))) {
            out[indx[[1]], indx[[2]]] <- mlist[[i]]
        }
        else {
            out[indx[, 1], indx[, 2]] <- mlist[[i]]
        }
    }
    return(out)
}
blup <-
function (x, ...) 
UseMethod("blup")
blup.rma.uni <-
function (x, level, digits, transf, targs, ...) 
{
    if (!is.element("rma.uni", class(x))) 
        stop("Argument 'x' must be an object of class \"rma.uni\".")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    if (missing(level)) 
        level <- x$level
    if (missing(digits)) 
        digits <- x$digits
    if (missing(transf)) 
        transf <- FALSE
    if (missing(targs)) 
        targs <- NULL
    alpha <- ifelse(level > 1, (100 - level)/100, 1 - level)
    if (x$knha || x$robust) {
        crit <- qt(alpha/2, df = x$k - x$p, lower.tail = FALSE)
    }
    else {
        crit <- qnorm(alpha/2, lower.tail = FALSE)
    }
    pred <- rep(NA, x$k.f)
    vpred <- rep(NA, x$k.f)
    li <- x$tau2/(x$tau2 + x$vi.f)
    for (i in seq_len(x$k.f)[x$not.na]) {
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
        transf <- TRUE
    }
    pi.bounds <- cbind(pi.lb, pi.ub)
    rev.order <- ifelse(pi.ub < pi.lb, TRUE, FALSE)
    rev.order[is.na(rev.order)] <- FALSE
    pi.bounds[rev.order, ] <- pi.bounds[rev.order, 2:1]
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
    out$transf <- transf
    class(out) <- c("list.rma")
    return(out)
}
cbind.escalc <-
function (..., deparse.level = 1) 
{
    dat <- data.frame(..., check.names = FALSE)
    arguments <- list(...)
    yi.names <- NULL
    vi.names <- NULL
    sei.names <- NULL
    zi.names <- NULL
    ci.lb.names <- NULL
    ci.ub.names <- NULL
    digits <- NULL
    for (arg in arguments) {
        yi.names <- c(attr(arg, "yi.names"), yi.names)
        vi.names <- c(attr(arg, "vi.names"), vi.names)
        sei.names <- c(attr(arg, "sei.names"), sei.names)
        zi.names <- c(attr(arg, "zi.names"), zi.names)
        ci.lb.names <- c(attr(arg, "ci.lb.names"), ci.lb.names)
        ci.ub.names <- c(attr(arg, "ci.ub.names"), ci.ub.names)
        digits <- c(attr(arg, "digits"), digits)
    }
    attr(dat, "yi.names") <- unique(yi.names)
    attr(dat, "vi.names") <- unique(vi.names)
    attr(dat, "sei.names") <- unique(sei.names)
    attr(dat, "zi.names") <- unique(zi.names)
    attr(dat, "ci.lb.names") <- unique(ci.lb.names)
    attr(dat, "ci.ub.names") <- unique(ci.ub.names)
    attr(dat, "digits") <- digits[1]
    class(dat) <- c("escalc", "data.frame")
    return(dat)
}
coef.permutest.rma.uni <-
function (object, ...) 
{
    if (!is.element("permutest.rma.uni", class(object))) 
        stop("Argument 'object' must be an object of class \"permutest.rma.uni\".")
    x <- object
    res.table <- cbind(estimate = x$b, se = x$se, zval = x$zval, 
        pval = x$pval, ci.lb = x$ci.lb, ci.ub = x$ci.ub)
    if (x$knha || x$robust) 
        colnames(res.table)[3] <- "tval"
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
    res.table <- cbind(estimate = as.vector(x$b), se = x$se, 
        zval = x$zval, pval = x$pval, ci.lb = x$ci.lb, ci.ub = x$ci.ub)
    if (x$knha || x$robust) 
        colnames(res.table)[3] <- "tval"
    res.table <- data.frame(res.table)
    return(res.table)
}
confint.rma.glmm <-
function (object, parm, level, digits, ...) 
{
    if (!is.element("rma.glmm", class(object))) 
        stop("Argument 'object' must be an object of class \"rma.glmm\".")
    stop("Method not yet implemented for objects of class \"rma.glmm\". Sorry!")
}
confint.rma.mh <-
function (object, parm, level, digits, ...) 
{
    if (!is.element("rma.mh", class(object))) 
        stop("Argument 'object' must be an object of class \"rma.mh\".")
    x <- object
    if (x$k == 1) 
        stop("Stopped because k = 1.")
    if (missing(level)) 
        level <- x$level
    if (missing(digits)) 
        digits <- x$digits
    alpha <- ifelse(level > 1, (100 - level)/100, 1 - level)
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
function (object, parm, level, digits, ...) 
{
    if (!is.element("rma.peto", class(object))) 
        stop("Argument 'object' must be an object of class \"rma.peto\".")
    x <- object
    if (x$k == 1) 
        stop("Stopped because k = 1.")
    if (missing(level)) 
        level <- x$level
    if (missing(digits)) 
        digits <- x$digits
    alpha <- ifelse(level > 1, (100 - level)/100, 1 - level)
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
function (object, parm, level, fixed = FALSE, random = TRUE, 
    digits, verbose = FALSE, control, ...) 
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
    if (missing(level)) 
        level <- x$level
    if (missing(digits)) 
        digits <- x$digits
    if (missing(control)) 
        control <- list()
    if (!fixed && !random) 
        stop("At least one of the arguments 'fixed' and 'random' must be TRUE.")
    con <- list(tol = .Machine$double.eps^0.25, maxiter = 1000, 
        tau2.min = ifelse(is.null(x$control$tau2.min), 0, x$control$tau2.min), 
        tau2.max = ifelse(is.null(x$control$tau2.max), 50, x$control$tau2.max), 
        verbose = FALSE)
    con[pmatch(names(control), names(con))] <- control
    if (verbose) 
        con$verbose <- verbose
    verbose <- con$verbose
    if (random) {
        if (!x$allvipos) 
            stop("Cannot compute confidence interval for the amount of (residual)\n  heterogeneity with non-positive sampling variances in the data.")
        alpha <- ifelse(level > 1, (100 - level)/100, 1 - level)
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
                  verbose = verbose, digits = digits)$root, silent = TRUE)
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
                verbose = verbose, digits = digits)$root, silent = TRUE)
            if (!is.numeric(tau2.ub)) {
                tau2.ub <- NA
                status.ub <- 0
                conv <- 0
            }
        }
        if (status.lb == 0L) 
            warning("Error in iterative search for the lower bound.")
        if (status.ub == 0L) 
            warning("Error in iterative search for the upper bound.")
        if (conv == 0L) 
            stop("Try increasing tau2.max (via the 'control' argument).")
        wi <- 1/vi
        W <- diag(wi, nrow = k, ncol = k)
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
        res.random <- rbind(`tau^2` = tau2, tau = tau, `I^2(%)` = I2, 
            `H^2` = H2)
        if (x$method == "FE") 
            res.random[, 1] <- NA
        colnames(res.random) <- c("estimate", "ci.lb", "ci.ub")
    }
    if (fixed) {
        alpha <- ifelse(level > 1, (100 - level)/100, 1 - level)
        if (x$knha || x$robust) {
            crit <- qt(alpha/2, df = k - p, lower.tail = FALSE)
        }
        else {
            crit <- qnorm(alpha/2, lower.tail = FALSE)
        }
        ci.lb <- c(x$b - crit * x$se)
        ci.ub <- c(x$b + crit * x$se)
        res.fixed <- cbind(estimate = x$b, ci.lb = ci.lb, ci.ub = ci.ub)
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
    svb <- chol2inv(chol(x$vb))
    o.warn <- getOption("warn")
    on.exit(options(warn = o.warn))
    options(warn = -1)
    for (i in seq_len(x$k.f)[x$not.na]) {
        res <- try(rma(x$yi.f[-i], x$vi.f[-i], weights = x$weights.f[-i], 
            mods = cbind(x$X.f[-i, ]), method = x$method, weighted = x$weighted, 
            intercept = FALSE, knha = x$knha, control = x$control), 
            silent = TRUE)
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
function (x, order, digits, transf, targs, ...) 
{
    if (!is.element("rma.mh", class(x))) 
        stop("Argument 'x' must be an object of class \"rma.mh\".")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    if (missing(order)) 
        order <- NULL
    if (missing(digits)) 
        digits <- x$digits
    if (missing(transf)) 
        transf <- FALSE
    if (missing(targs)) 
        targs <- NULL
    if (is.null(order)) 
        order <- seq_len(x$k.f)
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
    for (i in seq_len(x$k.f)[not.na]) {
        if (is.element(x$measure, c("RR", "OR", "RD"))) {
            res <- try(rma.mh(ai = ai.f[seq_len(i)], bi = bi.f[seq_len(i)], 
                ci = ci.f[seq_len(i)], di = di.f[seq_len(i)], 
                measure = x$measure, add = x$add, to = x$to, 
                drop00 = x$drop00, correct = x$correct), silent = TRUE)
        }
        else {
            res <- try(rma.mh(x1i = x1i.f[seq_len(i)], x2i = x2i.f[seq_len(i)], 
                t1i = t1i.f[seq_len(i)], t2i = t2i.f[seq_len(i)], 
                measure = x$measure, add = x$add, to = x$to, 
                drop00 = x$drop00, correct = x$correct), silent = TRUE)
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
    alpha <- ifelse(x$level > 1, (100 - x$level)/100, 1 - x$level)
    crit <- qnorm(alpha/2, lower.tail = FALSE)
    b[1] <- yi.f[1]
    se[1] <- sqrt(vi.f[1])
    zval[1] <- yi.f[1]/se[1]
    pval[1] <- 2 * pnorm(abs(zval[1]), lower.tail = FALSE)
    ci.lb[1] <- yi.f[1] - crit * se[1]
    ci.ub[1] <- yi.f[1] + crit * se[1]
    QE[1] <- 0
    QEp[1] <- 1
    if (is.logical(transf) && transf && is.element(x$measure, 
        c("OR", "RR", "IRR"))) 
        transf <- exp
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
        transf <- TRUE
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
    out$transf <- transf
    out$slab.null <- x$slab.null
    out$level <- x$level
    out$measure <- x$measure
    out$knha <- x$knha
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
function (x, order, digits, transf, targs, ...) 
{
    if (!is.element("rma.peto", class(x))) 
        stop("Argument 'x' must be an object of class \"rma.peto\".")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    if (missing(order)) 
        order <- NULL
    if (missing(digits)) 
        digits <- x$digits
    if (missing(transf)) 
        transf <- FALSE
    if (missing(targs)) 
        targs <- NULL
    if (is.null(order)) 
        order <- seq_len(x$k.f)
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
    for (i in seq_len(x$k.f)[not.na]) {
        res <- try(rma.peto(ai = ai.f[seq_len(i)], bi = bi.f[seq_len(i)], 
            ci = ci.f[seq_len(i)], di = di.f[seq_len(i)], add = x$add, 
            to = x$to, drop00 = x$drop00), silent = TRUE)
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
    alpha <- ifelse(x$level > 1, (100 - x$level)/100, 1 - x$level)
    crit <- qnorm(alpha/2, lower.tail = FALSE)
    b[1] <- yi.f[1]
    se[1] <- sqrt(vi.f[1])
    zval[1] <- yi.f[1]/se[1]
    pval[1] <- 2 * pnorm(abs(zval[1]), lower.tail = FALSE)
    ci.lb[1] <- yi.f[1] - crit * se[1]
    ci.ub[1] <- yi.f[1] + crit * se[1]
    QE[1] <- 0
    QEp[1] <- 1
    if (is.logical(transf) && transf) 
        transf <- exp
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
        transf <- TRUE
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
    out$transf <- transf
    out$slab.null <- x$slab.null
    out$level <- x$level
    out$measure <- x$measure
    out$knha <- x$knha
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
function (x, order, digits, transf, targs, ...) 
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
    if (missing(digits)) 
        digits <- x$digits
    if (missing(transf)) 
        transf <- FALSE
    if (missing(targs)) 
        targs <- NULL
    if (is.null(order)) 
        order <- seq_len(x$k.f)
    yi.f <- x$yi.f[order]
    vi.f <- x$vi.f[order]
    X.f <- cbind(x$X.f[order, ])
    weights.f <- x$weights.f[order]
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
    for (i in seq_len(x$k.f)[not.na]) {
        res <- try(rma(yi.f[seq_len(i)], vi.f[seq_len(i)], weights = weights.f[seq_len(i)], 
            method = x$method, weighted = x$weighted, intercept = TRUE, 
            knha = x$knha, control = x$control), silent = TRUE)
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
    alpha <- ifelse(x$level > 1, (100 - x$level)/100, 1 - x$level)
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
        transf <- TRUE
    }
    ci.bounds <- cbind(ci.lb, ci.ub)
    rev.order <- ifelse(ci.ub < ci.lb, TRUE, FALSE)
    rev.order[is.na(rev.order)] <- FALSE
    ci.bounds[rev.order, ] <- ci.bounds[rev.order, 2:1]
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
    out$transf <- transf
    out$slab.null <- x$slab.null
    out$level <- x$level
    out$measure <- x$measure
    out$knha <- x$knha
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
    dfbetas <- matrix(NA, nrow = x$k.f, ncol = x$p)
    if (x$weighted && !is.null(x$weights)) {
        A <- diag(x$weights, nrow = x$k, ncol = x$k)
        stXAX <- .invcalc(X = x$X, W = A, k = x$k)
    }
    else {
        stXX <- .invcalc(X = x$X, W = diag(x$k), k = x$k)
    }
    o.warn <- getOption("warn")
    on.exit(options(warn = o.warn))
    options(warn = -1)
    for (i in seq_len(x$k.f)[x$not.na]) {
        res <- try(rma(x$yi.f[-i], x$vi.f[-i], weights = x$weights.f[-i], 
            mods = cbind(x$X.f[-i, ]), method = x$method, weighted = x$weighted, 
            intercept = FALSE, knha = x$knha, control = x$control), 
            silent = TRUE)
        if (inherits(res, "try-error")) 
            next
        tau2.del[i] <- res$tau2
        dfbeta <- x$b - res$b
        if (x$weighted) {
            if (is.null(x$weights)) {
                vb.del <- .invcalc(X = x$X, W = diag(1/(x$vi + 
                  tau2.del[i]), nrow = x$k, ncol = x$k), k = x$k)
            }
            else {
                vb.del <- tcrossprod(stXAX, x$X) %*% A %*% diag(x$vi + 
                  tau2.del[i], nrow = x$k, ncol = x$k) %*% A %*% 
                  x$X %*% stXAX
            }
        }
        else {
            vb.del <- tcrossprod(stXX, x$X) %*% diag(x$vi + tau2.del[i], 
                nrow = x$k, ncol = x$k) %*% x$X %*% stXX
        }
        dfbetas[i, ] <- dfbeta/sqrt(res$s2w * diag(vb.del))
    }
    if (na.act == "na.omit") {
        out <- dfbetas[x$not.na, , drop = FALSE]
        rownames(out) <- x$slab[x$not.na]
    }
    if (na.act == "na.exclude" || na.act == "na.pass") {
        out <- dfbetas
        rownames(out) <- x$slab
    }
    if (na.act == "na.fail" && any(!x$not.na)) 
        stop("Missing values in results.")
    colnames(out) <- rownames(x$b)
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
    vtype = "LS", var.names = c("yi", "vi"), add.measure = FALSE, 
    append = TRUE, replace = TRUE, digits = 4, ...) 
{
    if (!is.element(measure, c("RR", "OR", "PETO", "RD", "AS", 
        "PHI", "YUQ", "YUY", "RTET", "PBIT", "OR2D", "OR2DN", 
        "OR2DL", "IRR", "IRD", "IRSD", "MD", "SMD", "SMDH", "ROM", 
        "RPB", "RBIS", "D2OR", "D2ORN", "D2ORL", "COR", "UCOR", 
        "ZCOR", "PR", "PLN", "PLO", "PAS", "PFT", "IR", "IRLN", 
        "IRS", "IRFT", "MN", "MC", "SMCC", "SMCR", "ARAW", "AHW", 
        "ABT"))) 
        stop("Unknown 'measure' specified.")
    if (!is.element(to, c("all", "only0", "if0all", "none"))) 
        stop("Unknown 'to' argument specified.")
    if (any(!is.element(vtype, c("UB", "LS", "HO", "ST", "CS")), 
        na.rm = TRUE)) 
        stop("Unknown 'vtype' argument specified.")
    if (add.measure) {
        if (length(var.names) == 2) 
            var.names <- c(var.names, "measure")
        if (length(var.names) != 3) 
            stop("Argument var.names must be of length 2 or 3.")
        if (any(var.names != make.names(var.names, unique = TRUE))) {
            var.names <- make.names(var.names, unique = TRUE)
            warning(paste0("Argument 'var.names' does not contain syntactically valid variable names.\n  Variable names adjusted to: var.names = c('", 
                var.names[1], "', '", var.names[2], "', '", var.names[3], 
                "')."))
        }
    }
    else {
        if (length(var.names) == 2) 
            var.names <- var.names[1:2]
        if (length(var.names) != 2) 
            stop("Argument var.names must be of length 2.")
        if (any(var.names != make.names(var.names, unique = TRUE))) {
            var.names <- make.names(var.names, unique = TRUE)
            warning(paste0("Argument 'var.names' does not contain syntactically valid variable names.\n  Variable names adjusted to: var.names = c('", 
                var.names[1], "', '", var.names[2], "')."))
        }
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
        "PHI", "YUQ", "YUY", "RTET", "PBIT", "OR2D", "OR2DN", 
        "OR2DL"))) {
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
        if (!is.null(subset)) {
            ai <- ai[subset]
            bi <- bi[subset]
            ci <- ci[subset]
            di <- di[subset]
        }
        if (length(ai) == 0L || length(bi) == 0L || length(ci) == 
            0L || length(di) == 0L) 
            stop("Cannot compute outcomes. Check that all of the required \n  information is specified via the appropriate arguments.")
        if (!all(length(ai) == c(length(ai), length(bi), length(ci), 
            length(di)))) 
            stop("Supplied data vectors are not all of the same length.")
        if (any(c(ai, bi, ci, di) < 0, na.rm = TRUE)) 
            stop("One or more counts are negative.")
        ni.u <- ai + bi + ci + di
        k <- length(ai)
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
        if (is.element(measure, c("OR", "OR2D", "OR2DN", "OR2DL"))) {
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
            if (length(vtype) == 1L) 
                vtype <- rep(vtype, k)
            vi <- rep(NA, k)
            mnwp1i <- sum(ai, na.rm = TRUE)/sum(n1i, na.rm = TRUE)
            mnwp2i <- sum(ci, na.rm = TRUE)/sum(n2i, na.rm = TRUE)
            for (i in seq_len(k)) {
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
            yi <- (ai/bi)/(ci/di)
            yi <- (yi - 1)/(yi + 1)
            vi <- 1/4 * (1 - yi^2)^2 * (1/ai + 1/bi + 1/ci + 
                1/di)
        }
        if (measure == "YUY") {
            yi <- (ai/bi)/(ci/di)
            yi <- (sqrt(yi) - 1)/(sqrt(yi) + 1)
            vi <- 1/16 * (1 - yi^2)^2 * (1/ai + 1/bi + 1/ci + 
                1/di)
        }
        if (measure == "RTET") {
            if (!require(polycor)) 
                stop("Please install the 'polycor' package to compute this measure.")
            warn.before <- getOption("warn")
            options(warn = 2)
            yi <- rep(NA, k)
            vi <- rep(NA, k)
            for (i in seq_len(k)) {
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
            z1i <- qnorm(p1i)
            z2i <- qnorm(p2i)
            yi <- z1i - z2i
            vi <- 2 * pi * p1i * (1 - p1i) * exp(z1i^2)/n1i + 
                2 * pi * p2i * (1 - p2i) * exp(z2i^2)/n2i
        }
        if (is.element(measure, c("OR2D", "OR2DL"))) {
            yi <- sqrt(3)/pi * yi
            vi <- 3/pi^2 * vi
        }
        if (measure == "OR2DN") {
            yi <- yi/1.65
            vi <- vi/1.65^2
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
        if (!is.null(subset)) {
            x1i <- x1i[subset]
            x2i <- x2i[subset]
            t1i <- t1i[subset]
            t2i <- t2i[subset]
        }
        if (length(x1i) == 0L || length(x2i) == 0L || length(t1i) == 
            0L || length(t2i) == 0L) 
            stop("Cannot compute outcomes. Check that all of the required \n  information is specified via the appropriate arguments.")
        if (!all(length(x1i) == c(length(x1i), length(x2i), length(t1i), 
            length(t2i)))) 
            stop("Supplied data vectors are not all of the same length.")
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
        "RBIS", "D2OR", "D2ORN", "D2ORL"))) {
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
        if (!is.null(subset)) {
            m1i <- m1i[subset]
            m2i <- m2i[subset]
            sd1i <- sd1i[subset]
            sd2i <- sd2i[subset]
            n1i <- n1i[subset]
            n2i <- n2i[subset]
        }
        if (length(m1i) == 0L || length(m2i) == 0L || length(sd1i) == 
            0L || length(sd2i) == 0L || length(n1i) == 0L || 
            length(n2i) == 0L) 
            stop("Cannot compute outcomes. Check that all of the required \n  information is specified via the appropriate arguments.")
        if (!all(length(m1i) == c(length(m1i), length(m2i), length(sd1i), 
            length(sd2i), length(n1i), length(n2i)))) 
            stop("Supplied data vectors are not all of the same length.")
        if (any(c(sd1i, sd2i) < 0, na.rm = TRUE)) 
            stop("One or more standard deviations are negative.")
        if (any(c(n1i, n2i) < 0, na.rm = TRUE)) 
            stop("One or more sample sizes are negative.")
        ni.u <- n1i + n2i
        k <- length(m1i)
        ni <- ni.u
        mi <- ni - 2
        spi <- sqrt(((n1i - 1) * sd1i^2 + (n2i - 1) * sd2i^2)/mi)
        di <- (m1i - m2i)/spi
        if (measure == "MD") {
            yi <- m1i - m2i
            vi <- sd1i^2/n1i + sd2i^2/n2i
        }
        if (measure == "SMD") {
            warn.before <- getOption("warn")
            options(warn = -1)
            cmi <- .cmicalc(mi)
            options(warn = warn.before)
            yi <- cmi * di
            if (length(vtype) == 1L) 
                vtype <- rep(vtype, k)
            vi <- rep(NA, k)
            mnwyi <- sum(ni * yi, na.rm = TRUE)/sum(ni, na.rm = TRUE)
            for (i in seq_len(k)) {
                if (vtype[i] == "UB") 
                  vi[i] <- 1/n1i[i] + 1/n2i[i] + (1 - (mi[i] - 
                    2)/(mi[i] * cmi[i]^2)) * yi[i]^2
                if (vtype[i] == "LS") 
                  vi[i] <- 1/n1i[i] + 1/n2i[i] + yi[i]^2/(2 * 
                    ni[i])
                if (vtype[i] == "HO") 
                  vi[i] <- 1/n1i[i] + 1/n2i[i] + mnwyi^2/(2 * 
                    ni[i])
            }
        }
        if (measure == "SMDH") {
            warn.before <- getOption("warn")
            options(warn = -1)
            cmi <- .cmicalc(mi)
            options(warn = warn.before)
            si <- sqrt((sd1i^2 + sd2i^2)/2)
            yi <- cmi * (m1i - m2i)/si
            vi <- yi^2 * (sd1i^4/(n1i - 1) + sd2i^4/(n2i - 1))/(2 * 
                (sd1i^2 + sd2i^2)^2) + (sd1i^2/(n1i - 1) + sd2i^2/(n2i - 
                1))/((sd1i^2 + sd2i^2)/2)
            vi <- cmi^2 * vi
        }
        if (measure == "ROM") {
            yi <- log(m1i/m2i)
            vi <- sd1i^2/(n1i * m1i^2) + sd2i^2/(n2i * m2i^2)
        }
        if (is.element(measure, c("RPB", "RBIS"))) {
            hi <- mi/n1i + mi/n2i
            yi <- di/sqrt(di^2 + hi)
            if (length(vtype) == 1L) 
                vtype <- rep(vtype, k)
            vi <- rep(NA, k)
            for (i in seq_len(k)) {
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
        if (is.element(measure, c("D2OR", "D2ORL"))) {
            yi <- pi/sqrt(3) * di
            vi <- pi^2/3 * (1/n1i + 1/n2i + di^2/(2 * ni))
        }
        if (measure == "D2ORN") {
            yi <- 1.65 * di
            vi <- 1.65^2 * (1/n1i + 1/n2i + di^2/(2 * ni))
        }
    }
    if (is.element(measure, c("COR", "UCOR", "ZCOR"))) {
        mf.ri <- mf[[match("ri", names(mf))]]
        mf.ni <- mf[[match("ni", names(mf))]]
        ri <- eval(mf.ri, data, enclos = sys.frame(sys.parent()))
        ni <- eval(mf.ni, data, enclos = sys.frame(sys.parent()))
        if (!is.null(subset)) {
            ri <- ri[subset]
            ni <- ni[subset]
        }
        if (length(ri) == 0L || length(ni) == 0L) 
            stop("Cannot compute outcomes. Check that all of the required \n  information is specified via the appropriate arguments.")
        if (length(ri) != length(ni)) 
            stop("Supplied data vectors are not of the same length.")
        if (any(abs(ri) > 1, na.rm = TRUE)) 
            stop("One or more correlations are > 1 or < -1.")
        if (any(ni < 0, na.rm = TRUE)) 
            stop("One or more sample sizes are negative.")
        if (any(abs(ni) <= 4, na.rm = TRUE)) 
            warning("Cannot estimate sampling variance when ni <= 4.")
        ni.u <- ni
        k <- length(ri)
        if (measure == "COR") {
            yi <- ri
        }
        if (measure == "UCOR") {
            yi <- ri + ri * (1 - ri^2)/(2 * (ni - 4))
            yi[ni <= 4] <- NA
        }
        if (is.element(measure, c("COR", "UCOR"))) {
            if (length(vtype) == 1L) 
                vtype <- rep(vtype, k)
            vi <- rep(NA, k)
            mnwyi <- sum(ni * yi, na.rm = TRUE)/sum(ni, na.rm = TRUE)
            for (i in seq_len(k)) {
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
        }
        if (measure == "ZCOR") {
            yi <- 1/2 * log((1 + ri)/(1 - ri))
            vi <- 1/(ni - 3)
        }
        vi[ni <= 4] <- NA
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
        if (!is.null(subset)) {
            xi <- xi[subset]
            mi <- mi[subset]
        }
        if (length(xi) == 0L || length(mi) == 0L) 
            stop("Cannot compute outcomes. Check that all of the required \n  information is specified via the appropriate arguments.")
        if (length(xi) != length(mi)) 
            stop("Supplied data vectors are not all of the same length.")
        if (any(c(xi, mi) < 0, na.rm = TRUE)) 
            stop("One or more counts are negative.")
        ni.u <- xi + mi
        k <- length(xi)
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
            if (length(vtype) == 1L) 
                vtype <- rep(vtype, k)
            vi <- rep(NA, k)
            mnwpri <- sum(xi, na.rm = TRUE)/sum(ni, na.rm = TRUE)
            for (i in seq_len(k)) {
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
            if (length(vtype) == 1L) 
                vtype <- rep(vtype, k)
            vi <- rep(NA, k)
            mnwpri <- sum(xi, na.rm = TRUE)/sum(ni, na.rm = TRUE)
            for (i in seq_len(k)) {
                if (vtype[i] == "LS") 
                  vi[i] <- 1/xi[i] - 1/ni[i]
                if (vtype[i] == "HO") 
                  vi[i] <- 1/(mnwpri * ni[i]) - 1/ni[i]
            }
        }
        if (measure == "PLO") {
            yi <- log(pri/(1 - pri))
            if (length(vtype) == 1L) 
                vtype <- rep(vtype, k)
            vi <- rep(NA, k)
            mnwpri <- sum(xi, na.rm = TRUE)/sum(ni, na.rm = TRUE)
            for (i in seq_len(k)) {
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
        if (!is.null(subset)) {
            xi <- xi[subset]
            ti <- ti[subset]
        }
        if (length(xi) == 0L || length(ti) == 0L) 
            stop("Cannot compute outcomes. Check that all of the required \n  information is specified via the appropriate arguments.")
        if (length(xi) != length(ti)) 
            stop("Supplied data vectors are not all of the same length.")
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
        if (!is.null(subset)) {
            mi <- mi[subset]
            sdi <- sdi[subset]
            ni <- ni[subset]
        }
        if (length(mi) == 0L || length(sdi) == 0L || length(ni) == 
            0L) 
            stop("Cannot compute outcomes. Check that all of the required \n  information is specified via the appropriate arguments.")
        if (!all(length(mi) == c(length(mi), length(sdi), length(ni)))) 
            stop("Supplied data vectors are not all of the same length.")
        if (any(sdi < 0, na.rm = TRUE)) 
            stop("One or more standard deviations are negative.")
        if (any(ni < 0, na.rm = TRUE)) 
            stop("One or more sample sizes are negative.")
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
        if (!is.null(subset)) {
            m1i <- m1i[subset]
            m2i <- m2i[subset]
            sd1i <- sd1i[subset]
            sd2i <- sd2i[subset]
            ni <- ni[subset]
            ri <- ri[subset]
        }
        if (is.element(measure, c("MC", "SMCC"))) {
            if (length(m1i) == 0L || length(m2i) == 0L || length(sd1i) == 
                0L || length(sd2i) == 0L || length(ni) == 0L || 
                length(ri) == 0L) 
                stop("Cannot compute outcomes. Check that all of the required \n  information is specified via the appropriate arguments.")
            if (!all(length(m1i) == c(length(m1i), length(m2i), 
                length(sd1i), length(sd2i), length(ni), length(ri)))) 
                stop("Supplied data vectors are not all of the same length.")
            if (any(c(sd1i, sd2i) < 0, na.rm = TRUE)) 
                stop("One or more standard deviations are negative.")
        }
        else {
            if (length(m1i) == 0L || length(m2i) == 0L || length(sd1i) == 
                0L || length(ni) == 0L || length(ri) == 0L) 
                stop("Cannot compute outcomes. Check that all of the required \n  information is specified via the appropriate arguments.")
            if (!all(length(m1i) == c(length(m1i), length(m2i), 
                length(sd1i), length(ni), length(ri)))) 
                stop("Supplied data vectors are not all of the same length.")
            if (any(sd1i < 0, na.rm = TRUE)) 
                stop("One or more standard deviations are negative.")
        }
        if (any(abs(ri) > 1, na.rm = TRUE)) 
            stop("One or more correlations are > 1 or < -1.")
        if (any(ni < 0, na.rm = TRUE)) 
            stop("One or more sample sizes are negative.")
        ni.u <- ni
        ni <- ni.u
        mi <- ni - 1
        if (measure == "MC") {
            yi <- m1i - m2i
            vi <- (sd1i^2 + sd2i^2 - 2 * ri * sd1i * sd2i)/ni
        }
        if (measure == "SMCC") {
            warn.before <- getOption("warn")
            options(warn = -1)
            cmi <- .cmicalc(mi)
            options(warn = warn.before)
            sddi <- sqrt(sd1i^2 + sd2i^2 - 2 * ri * sd1i * sd2i)
            yi <- cmi * (m1i - m2i)/sddi
            vi <- 1/ni + yi^2/(2 * ni)
        }
        if (measure == "SMCR") {
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
        if (!is.null(subset)) {
            ai <- ai[subset]
            mi <- mi[subset]
            ni <- ni[subset]
        }
        if (length(ai) == 0L || length(mi) == 0L || length(ni) == 
            0L) 
            stop("Cannot compute outcomes. Check that all of the required \n  information is specified via the appropriate arguments.")
        if (!all(length(ai) == c(length(ai), length(mi), length(ni)))) 
            stop("Supplied data vectors are not all of the same length.")
        if (any(ai > 1, na.rm = TRUE)) 
            stop("One or more alphas are > 1.")
        if (any(mi < 2, na.rm = TRUE)) 
            stop("One or more mi's are < 2.")
        if (any(ni < 0, na.rm = TRUE)) 
            stop("One or more sample sizes are negative.")
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
    is.inf <- is.infinite(yi) | is.infinite(vi)
    if (any(is.inf)) {
        warning("Some yi and/or vi values equal to +-Inf. Recoded to NAs.")
        yi[is.inf] <- NA
        vi[is.inf] <- NA
    }
    is.NaN <- is.nan(yi) | is.nan(vi)
    if (any(is.NaN)) {
        yi[is.NaN] <- NA
        vi[is.NaN] <- NA
    }
    vi[vi < 0] <- NA
    if (!is.null(slab)) {
        if (!is.null(subset)) 
            slab <- slab[subset]
        if (any(is.na(slab))) 
            stop("NAs in study labels.")
        if (any(duplicated(slab))) 
            slab <- make.unique(slab)
        if (length(slab) != length(yi)) 
            stop("Study labels not of same length as data.")
        attr(yi, "slab") <- slab
    }
    if (!is.null(subset)) {
        if (!no.data) 
            data <- data[subset, , drop = FALSE]
    }
    attr(yi, "measure") <- measure
    if (!no.data && append) {
        dat <- data.frame(data)
        if (replace) {
            dat[[var.names[1]]] <- yi
            dat[[var.names[2]]] <- vi
            if (add.measure) {
                dat[[var.names[3]]] <- ""
                dat[[var.names[3]]][!is.na(yi)] <- measure
            }
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
            if (add.measure) {
                if (is.element(var.names[3], names(dat))) {
                  is.na.measure <- c(dat[[var.names[3]]] == "") & 
                    !is.na(yi)
                  dat[[var.names[3]]][is.na.measure] <- measure
                }
                else {
                  dat[[var.names[3]]] <- ""
                  dat[[var.names[3]]][!is.na(yi)] <- measure
                }
            }
        }
    }
    else {
        if (add.measure) {
            dat <- data.frame(yi, vi)
            dat$measure <- ""
            dat$measure[!is.na(yi)] <- measure
            names(dat) <- var.names
        }
        else {
            dat <- data.frame(yi, vi)
            names(dat) <- var.names[1:2]
        }
        attr(dat$yi, "ni") <- ni.u
    }
    attr(dat, "digits") <- digits
    attr(dat, "yi.names") <- unique(c(var.names[1], attr(data, 
        "yi.names")))
    attr(dat, "vi.names") <- unique(c(var.names[2], attr(data, 
        "vi.names")))
    attr(dat, "sei.names") <- attr(data, "sei.names")
    attr(dat, "zi.names") <- attr(data, "zi.names")
    attr(dat, "ci.lb.names") <- attr(data, "ci.lb.names")
    attr(dat, "ci.ub.names") <- attr(data, "ci.ub.names")
    class(dat) <- c("escalc", "data.frame")
    return(dat)
}
escalc.formula <-
function (measure, formula, weights, data, add = 1/2, to = "only0", 
    drop00 = FALSE, vtype = "LS", var.names = c("yi", "vi"), 
    digits = 4, ...) 
{
    if (!is.element(measure, c("RR", "OR", "PETO", "RD", "AS", 
        "PHI", "YUQ", "YUY", "RTET", "PBIT", "OR2D", "OR2DN", 
        "OR2DL", "IRR", "IRD", "IRSD", "MD", "SMD", "SMDH", "ROM", 
        "RPB", "RBIS", "D2OR", "D2ORN", "D2ORL", "COR", "UCOR", 
        "ZCOR", "PR", "PLN", "PLO", "PAS", "PFT", "IR", "IRLN", 
        "IRS", "IRFT", "MN", "MC", "SMCC", "SMCR", "ARAW", "AHW", 
        "ABT"))) 
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
        "PHI", "YUQ", "YUY", "RTET", "PBIT", "OR2D", "OR2DN", 
        "OR2DL"))) {
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
        "RBIS", "D2OR", "D2ORN", "D2ORL"))) {
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
        rownames(out) <- c("logLik:", "deviance:", "AIC:", "BIC:", 
            "AICc:")
        colnames(out) <- c("REML")
    }
    else {
        out <- cbind(object$fit.stats$ML)
        rownames(out) <- c("logLik:", "deviance:", "AIC:", "BIC:", 
            "AICc:")
        colnames(out) <- c("ML")
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
function (x, annotate = TRUE, xlim, alim, clim, ylim, at, steps = 5, 
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
    alpha <- ifelse(level > 1, (100 - level)/100, 1 - level)
    yi <- x$estimate
    if (is.null(attr(yi, "measure"))) {
        measure <- "GEN"
    }
    else {
        measure <- attr(yi, "measure")
    }
    vi <- x$se^2
    ci.lb <- x$ci.lb
    ci.ub <- x$ci.ub
    if (length(yi) != length(vi)) 
        stop("Length of yi and vi (or sei) vectors are not the same.")
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
            for (j in seq_len(length(rows.na))) {
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
    ci.bounds[rev.order, ] <- ci.bounds[rev.order, 2:1]
    ci.lb <- ci.bounds[, 1]
    ci.ub <- ci.bounds[, 2]
    if (!missing(clim)) {
        clim <- sort(clim)
        if (length(clim) != 2L) 
            stop("Argument 'clim' must be of length 2.")
        ci.lb[ci.lb < clim[1]] <- clim[1]
        ci.ub[ci.ub > clim[2]] <- clim[2]
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
            at <- seq(from = alim[1], to = alim[2], length.out = steps)
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
    if (missing(xlab)) 
        xlab <- .setxlab(measure, transf.char, atransf.char, 
            gentype = 2)
    mtext(xlab, side = 1, at = min(at) + (max(at) - min(at))/2, 
        line = par("mgp")[1] - 0.5, cex = cex.lab, ...)
    for (i in seq_len(k)) {
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
        for (l in seq_len(NCOL(ilab))) {
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
    res <- list(xlim = par("usr")[1:2], alim = alim, at = at, 
        ylim = ylim, rows = rows, cex = cex, cex.lab = cex.lab, 
        cex.axis = cex.axis)
    invisible(res)
}
forest.default <-
function (x, vi, sei, ci.lb, ci.ub, annotate = TRUE, showweight = FALSE, 
    xlim, alim, clim, ylim, at, steps = 5, level = 95, digits = 2, 
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
    alpha <- ifelse(level > 1, (100 - level)/100, 1 - level)
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
            slab <- paste("Study ", seq_len(k))
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
            for (j in seq_len(length(rows.na))) {
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
    ci.bounds[rev.order, ] <- ci.bounds[rev.order, 2:1]
    ci.lb <- ci.bounds[, 1]
    ci.ub <- ci.bounds[, 2]
    if (!missing(clim)) {
        clim <- sort(clim)
        if (length(clim) != 2L) 
            stop("Argument 'clim' must be of length 2.")
        ci.lb[ci.lb < clim[1]] <- clim[1]
        ci.ub[ci.ub > clim[2]] <- clim[2]
    }
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
            at <- seq(from = alim[1], to = alim[2], length.out = steps)
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
    if (missing(xlab)) 
        xlab <- .setxlab(measure, transf.char, atransf.char, 
            gentype = 1)
    mtext(xlab, side = 1, at = min(at) + (max(at) - min(at))/2, 
        line = par("mgp")[1] - 0.5, cex = cex.lab, ...)
    for (i in seq_len(k)) {
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
        for (l in seq_len(NCOL(ilab))) {
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
    res <- list(xlim = par("usr")[1:2], alim = alim, at = at, 
        ylim = ylim, rows = rows, cex = cex, cex.lab = cex.lab, 
        cex.axis = cex.axis)
    invisible(res)
}
forest.rma <-
function (x, annotate = TRUE, addfit = TRUE, addcred = FALSE, 
    showweight = FALSE, xlim, alim, clim, ylim, at, steps = 5, 
    level = x$level, digits = 2, refline = 0, xlab, slab, mlab, 
    ilab, ilab.xpos, ilab.pos, order, transf = FALSE, atransf = FALSE, 
    targs, rows, efac = 1, pch = 15, psize, col, border, lty, 
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
    if (x$int.only) {
        if (missing(col)) {
            col <- c("black", "gray50")
        }
        else {
            if (length(col) == 1L) 
                col <- c(col, "gray50")
        }
        if (missing(border)) 
            border <- "black"
    }
    else {
        if (missing(col)) 
            col <- "gray"
        if (missing(border)) 
            border <- "gray"
    }
    if (missing(lty)) {
        lty <- c("solid", "dotted")
    }
    else {
        if (length(lty) == 1L) 
            lty <- c(lty, "dotted")
    }
    measure <- x$measure
    if (length(digits) == 1L) 
        digits <- c(digits, digits)
    alpha <- ifelse(level > 1, (100 - level)/100, 1 - level)
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
        if (addcred) {
            pred.ci.lb <- temp$cr.lb
            pred.ci.ub <- temp$cr.ub
        }
        else {
            pred.ci.lb <- temp$ci.lb
            pred.ci.ub <- temp$ci.ub
        }
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
            for (j in seq_len(length(rows.na))) {
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
    ci.bounds[rev.order, ] <- ci.bounds[rev.order, 2:1]
    ci.lb <- ci.bounds[, 1]
    ci.ub <- ci.bounds[, 2]
    pred.ci.bounds <- cbind(pred.ci.lb, pred.ci.ub)
    rev.order <- ifelse(pred.ci.ub < pred.ci.lb, TRUE, FALSE)
    rev.order[is.na(rev.order)] <- FALSE
    pred.ci.bounds[rev.order, ] <- pred.ci.bounds[rev.order, 
        2:1]
    pred.ci.lb <- pred.ci.bounds[, 1]
    pred.ci.ub <- pred.ci.bounds[, 2]
    if (!missing(clim)) {
        clim <- sort(clim)
        if (length(clim) != 2L) 
            stop("Argument 'clim' must be of length 2.")
        ci.lb[ci.lb < clim[1]] <- clim[1]
        ci.ub[ci.ub > clim[2]] <- clim[2]
        pred.ci.lb[pred.ci.lb < clim[1]] <- clim[1]
        pred.ci.ub[pred.ci.ub > clim[2]] <- clim[2]
    }
    if (is.null(psize)) {
        if (is.null(weights)) {
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
        else {
            wi <- weights
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
            at <- seq(from = alim[1], to = alim[2], length.out = steps)
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
        for (i in seq_len(k)) {
            if (is.na(pred[i])) 
                next
            polygon(x = c(max(pred.ci.lb[i], alim[1]), pred[i], 
                min(pred.ci.ub[i], alim[2]), pred[i]), y = c(rows[i], 
                rows[i] + (height/100) * cex * efac, rows[i], 
                rows[i] - (height/100) * cex * efac), col = col, 
                border = border, ...)
        }
    }
    if (addfit && x$int.only) {
        if (is.element("rma.mv", class(x)) && x$withG && is.element(x$struct, 
            c("HCS", "UN", "HAR"))) {
            if (!is.logical(addcred)) {
                temp <- predict(x, level = level, tau2.levels = addcred[1])
                addcred <- TRUE
            }
            else {
                if (addcred) {
                  stop("Need to specify the level of the inner factor via the 'addcred' argument.")
                }
                else {
                  temp <- predict(x, level = level, tau2.levels = 1)
                }
            }
        }
        else {
            temp <- predict(x, level = level)
        }
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
        b.ci.bounds[rev.order, ] <- b.ci.bounds[rev.order, 2:1]
        b.ci.lb <- b.ci.bounds[, 1]
        b.ci.ub <- b.ci.bounds[, 2]
        b.cr.bounds <- cbind(b.cr.lb, b.cr.ub)
        rev.order <- ifelse(b.cr.ub < b.cr.lb, TRUE, FALSE)
        rev.order[is.na(rev.order)] <- FALSE
        b.cr.bounds[rev.order, ] <- b.cr.bounds[rev.order, 2:1]
        b.cr.lb <- b.cr.bounds[, 1]
        b.cr.ub <- b.cr.bounds[, 2]
        if (!missing(clim)) {
            b.ci.lb[b.ci.lb < clim[1]] <- clim[1]
            b.ci.ub[b.ci.ub > clim[2]] <- clim[2]
            b.cr.lb[b.cr.lb < clim[1]] <- clim[1]
            b.cr.ub[b.cr.ub > clim[2]] <- clim[2]
        }
        if (x$method != "FE" && addcred) {
            segments(max(b.cr.lb, alim[1]), -1, min(b.cr.ub, 
                alim[2]), -1, lty = lty[2], col = col[2], ...)
            if (b.cr.lb >= alim[1]) {
                segments(b.cr.lb, -1 - (height/150) * cex * efac, 
                  b.cr.lb, -1 + (height/150) * cex * efac, col = col[2], 
                  ...)
            }
            else {
                polygon(x = c(alim[1], alim[1] + (1.4/100) * 
                  cex * (xlim[2] - xlim[1]), alim[1] + (1.4/100) * 
                  cex * (xlim[2] - xlim[1]), alim[1]), y = c(-1, 
                  -1 + (height/150) * cex * efac, -1 - (height/150) * 
                    cex * efac, -1), col = col[2], border = col[2], 
                  ...)
            }
            if (b.cr.ub <= alim[2]) {
                segments(b.cr.ub, -1 - (height/150) * cex * efac, 
                  b.cr.ub, -1 + (height/150) * cex * efac, col = col[2], 
                  ...)
            }
            else {
                polygon(x = c(alim[2], alim[2] - (1.4/100) * 
                  cex * (xlim[2] - xlim[1]), alim[2] - (1.4/100) * 
                  cex * (xlim[2] - xlim[1]), alim[2]), y = c(-1, 
                  -1 + (height/150) * cex * efac, -1 - (height/150) * 
                    cex * efac, -1), col = col[2], border = col[2], 
                  ...)
            }
        }
        polygon(x = c(b.ci.lb, b, b.ci.ub, b), y = c(-1, -1 + 
            (height/100) * cex * efac, -1, -1 - (height/100) * 
            cex * efac), col = col[1], border = border, ...)
        if (missing(mlab)) 
            mlab <- ifelse((x$method == "FE"), "FE Model", "RE Model")
        text(xlim[1], -1, mlab, pos = 4, cex = cex, ...)
    }
    axis(side = 1, at = at, labels = at.lab, cex.axis = cex.axis, 
        ...)
    if (missing(xlab)) 
        xlab <- .setxlab(measure, transf.char, atransf.char, 
            gentype = 1)
    mtext(xlab, side = 1, at = min(at) + (max(at) - min(at))/2, 
        line = par("mgp")[1] - 0.5, cex = cex.lab, ...)
    for (i in seq_len(k)) {
        if (is.na(yi[i]) || is.na(vi)[i]) 
            next
        segments(max(ci.lb[i], alim[1]), rows[i], min(ci.ub[i], 
            alim[2]), rows[i], lty = lty[1], ...)
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
        for (l in seq_len(NCOL(ilab))) {
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
    res <- list(xlim = par("usr")[1:2], alim = alim, at = at, 
        ylim = ylim, rows = rows, cex = cex, cex.lab = cex.lab, 
        cex.axis = cex.axis)
    invisible(res)
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
    if (is.null(vi)) {
        if (is.null(sei)) {
            stop("Need to specify vi or sei argument.")
        }
        else {
            vi <- sei^2
        }
    }
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
    res <- list(type = type, fsnum = fsnum, alpha = alpha, pval = pval, 
        meanes = meanes, target = target, digits = digits)
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
    na.act <- getOption("na.action")
    yaxis <- match.arg(yaxis, c("sei", "vi", "seinv", "vinv", 
        "ni", "ninv", "sqrtni", "sqrtninv", "lni", "wi"))
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
        if (yaxis == "wi") 
            ylab <- "Weight (in %)"
    }
    if (missing(at)) 
        at <- NULL
    if (missing(targs)) 
        targs <- NULL
    if (is.element(yaxis, c("ni", "ninv", "sqrtni", "sqrtninv", 
        "lni"))) {
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
        if (yaxis == "wi") 
            digits <- c(2, 2)
    }
    else {
        if (length(digits) == 1L) 
            digits <- c(digits, digits)
    }
    if (x$int.only) {
        if (missing(refline)) 
            refline <- x$b
        if (is.element("rma.mv", class(x)) && x$withG && is.element(x$struct, 
            c("HCS", "UN", "HAR"))) {
            if (!is.logical(addtau2)) {
                temp <- predict(x, tau2.levels = addtau2[1])
                if (is.numeric(addtau2)) {
                  tau2 <- x$tau2[addtau2[1]]
                }
                else {
                  tau2 <- x$tau2[pmatch(addtau2[1], x$g.levels.f[[1]])]
                }
            }
            else {
                if (addtau2) {
                  stop("Need to specify the level of the inner factor via the 'addtau2' argument.")
                }
                else {
                  tau2 <- 0
                }
            }
        }
        else {
            tau2 <- ifelse(addtau2, x$tau2, 0)
        }
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
    options(na.action = "na.omit")
    weights <- weights(x)
    options(na.action = na.act)
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
        if (yaxis == "wi") 
            ylim <- c(min(weights), max(weights))
    }
    else {
        if (is.element(yaxis, c("sei", "vi", "ninv", "sqrtninv"))) 
            ylim <- c(max(ylim), min(ylim))
        if (is.element(yaxis, c("seinv", "vinv", "ni", "sqrtni", 
            "lni", "wi"))) 
            ylim <- c(min(ylim), max(ylim))
        if (is.element(yaxis, c("sei", "vi", "ni", "ninv", "sqrtni", 
            "sqrtninv", "lni"))) {
            if (ylim[1] < 0 || ylim[2] < 0) 
                stop("Both limits for the y axis must be >= 0.")
        }
        if (is.element(yaxis, c("seinv", "vinv"))) {
            if (ylim[1] <= 0 || ylim[2] <= 0) 
                stop("Both limits for the y axis must be > 0.")
        }
        if (is.element(yaxis, c("wi"))) {
            if (ylim[1] < 0 || ylim[2] < 0) 
                stop("Both limits for the y axis must be >= 0.")
        }
    }
    if (is.element(yaxis, c("sei", "vi", "seinv", "vinv"))) {
        alpha <- ifelse(level > 1, (100 - level)/100, 1 - level)
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
    if (is.element(yaxis, c("ni", "ninv", "sqrtni", "sqrtninv", 
        "lni", "wi"))) {
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
    axis(side = 2, at = seq(from = ylim[1], to = ylim[2], length.out = steps), 
        labels = formatC(seq(from = ylim[1], to = ylim[2], length.out = steps), 
            digits = digits[2], format = "f"), ...)
    abline(h = seq(from = ylim[1], to = ylim[2], length.out = steps), 
        col = hlines, ...)
    if (is.element(yaxis, c("sei", "vi", "seinv", "vinv"))) {
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
        yi.vals <- seq(from = ylim[1], to = ylim[2], length.out = ci.res)
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
    if (is.element(yaxis, c("sei", "vi", "seinv", "vinv"))) 
        segments(refline, ylim[1], refline, ylim[2], ...)
    if (is.element(yaxis, c("ni", "ninv", "sqrtni", "sqrtninv", 
        "lni", "wi"))) 
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
    if (yaxis == "wi") 
        points(yi, weights, pch = pch, ...)
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
        if (yaxis == "wi") 
            points(yi[x$fill == 1], (weights)[x$fill == 1], pch = pch.fill, 
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
hatvalues.rma.mv <-
function (model, type = "diagonal", ...) 
{
    if (!is.element("rma.mv", class(model))) 
        stop("Argument 'model' must be an object of class \"rma.mv\".")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    type <- match.arg(type, c("diagonal", "matrix"))
    x <- model
    if (is.null(x$W)) {
        W <- chol2inv(chol(x$M))
        H <- as.matrix(x$X %*% x$vb %*% crossprod(x$X, W))
    }
    else {
        A <- x$W
        stXAX <- chol2inv(chol(as.matrix(t(x$X) %*% A %*% x$X)))
        H <- as.matrix(x$X %*% stXAX %*% crossprod(x$X, A))
    }
    if (type == "diagonal") {
        hii <- rep(NA, x$k.f)
        hii[x$not.na] <- as.vector(diag(H))
        hii[hii > 1 - 10 * .Machine$double.eps] <- 1
        names(hii) <- x$slab
        if (na.act == "na.omit") 
            hii <- hii[x$not.na]
        if (na.act == "na.fail" && any(!x$not.na)) 
            stop("Missing values in results.")
        return(hii)
    }
    if (type == "matrix") {
        Hfull <- matrix(NA, nrow = x$k.f, ncol = x$k.f)
        Hfull[x$not.na, x$not.na] <- H
        rownames(Hfull) <- x$slab
        colnames(Hfull) <- x$slab
        if (na.act == "na.omit") 
            Hfull <- Hfull[x$not.na, x$not.na, drop = FALSE]
        if (na.act == "na.fail" && any(!x$not.na)) 
            stop("Missing values in results.")
        return(Hfull)
    }
}
hatvalues.rma.uni <-
function (model, type = "diagonal", ...) 
{
    if (!is.element("rma.uni", class(model))) 
        stop("Argument 'model' must be an object of class \"rma.uni\".")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    type <- match.arg(type, c("diagonal", "matrix"))
    x <- model
    if (x$weighted) {
        if (is.null(x$weights)) {
            wi <- 1/(x$vi + x$tau2)
            W <- diag(wi, nrow = x$k, ncol = x$k)
            H <- x$X %*% x$vb %*% crossprod(x$X, W)
        }
        else {
            A <- diag(x$weights, nrow = x$k, ncol = x$k)
            stXAX <- .invcalc(X = x$X, W = A, k = x$k)
            H <- x$X %*% stXAX %*% crossprod(x$X, A)
        }
    }
    else {
        stXX <- .invcalc(X = x$X, W = diag(x$k), k = x$k)
        H <- x$X %*% tcrossprod(stXX, x$X)
    }
    if (type == "diagonal") {
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
    if (type == "matrix") {
        Hfull <- matrix(NA, nrow = x$k.f, ncol = x$k.f)
        Hfull[x$not.na, x$not.na] <- H
        rownames(Hfull) <- x$slab
        colnames(Hfull) <- x$slab
        if (na.act == "na.omit") 
            Hfull <- Hfull[x$not.na, x$not.na, drop = FALSE]
        if (na.act == "na.fail" && any(!x$not.na)) 
            stop("Missing values in results.")
        return(Hfull)
    }
}
hc <-
function (object, ...) 
UseMethod("hc")
hc.rma.uni <-
function (object, digits, transf, targs, control, ...) 
{
    if (!is.element("rma.uni", class(object))) 
        stop("Argument 'object' must be an object of class \"rma.uni\".")
    x <- object
    if (!x$int.only) 
        stop("Method only applicable for models without moderators.")
    if (missing(digits)) 
        digits <- x$digits
    if (missing(transf)) 
        transf <- FALSE
    if (missing(targs)) 
        targs <- NULL
    yi <- x$yi
    vi <- x$vi
    k <- length(yi)
    if (k == 1) 
        stop("Stopped because k = 1.")
    alpha <- ifelse(x$level > 1, (100 - x$level)/100, 1 - x$level)
    if (missing(control)) 
        control <- list()
    con <- list(tol = .Machine$double.eps^0.25, maxiter = 1000, 
        verbose = FALSE)
    con[pmatch(names(control), names(con))] <- control
    wi <- 1/vi
    W1 <- sum(wi)
    W2 <- sum(wi^2)/W1
    W3 <- sum(wi^3)/W1
    W4 <- sum(wi^4)/W1
    b <- sum(wi * yi)/W1
    Q <- sum(wi * ((yi - b)^2))
    tau2 <- max(0, (Q - (k - 1))/(W1 - W2))
    vb <- (tau2 * W2 + 1)/W1
    se <- sqrt(vb)
    VR <- 1 + tau2 * W2
    SDR <- sqrt(VR)
    EQ <- function(r) (k - 1) + tau2 * (W1 - W2) + (tau2^2) * 
        ((1/VR^2) * (r^2) - 1/VR) * (W3 - W2^2)
    VQ <- function(r) {
        rsq <- r^2
        recipvr2 <- 1/VR^2
        2 * (k - 1) + 4 * tau2 * (W1 - W2) + 2 * tau2^2 * (W1 * 
            W2 - 2 * W3 + W2^2) + 4 * tau2^2 * (recipvr2 * rsq - 
            1/VR) * (W3 - W2^2) + 4 * tau2^3 * (recipvr2 * rsq - 
            1/VR) * (W4 - 2 * W2 * W3 + W2^3) + 2 * tau2^4 * 
            (recipvr2 - 2 * (1/VR^3) * rsq) * (W3 - W2^2)^2
    }
    scale <- function(r) {
        VQ(r)/EQ(r)
    }
    shape <- function(r) {
        EQ(r)^2/VQ(r)
    }
    finv <- function(f) (W1/W2 - 1) * ((f^2) - 1) + (k - 1)
    Eqn <- function(t) {
        integrand <- function(r) {
            pgamma(finv(r/t), scale = scale(SDR * r), shape = shape(SDR * 
                r)) * dnorm(r)
        }
        integral <- integrate(integrand, lower = t, upper = Inf)$value
        val <- integral - alpha/2
        val
    }
    t0 <- try(uniroot(Eqn, lower = 0, upper = 2, tol = con$tol, 
        maxiter = con$maxiter))
    if (inherits(t0, "try-error")) 
        stop("Error in uniroot().")
    t0 <- t0$root
    u0 <- SDR * t0
    ci.lb <- b - u0 * se
    ci.ub <- b + u0 * se
    b.rma <- x$b
    se.rma <- x$se
    ci.lb.rma <- x$ci.lb
    ci.ub.rma <- x$ci.ub
    if (is.function(transf)) {
        if (is.null(targs)) {
            b <- sapply(b, transf)
            b.rma <- sapply(b.rma, transf)
            se <- NA
            se.rma <- NA
            ci.lb <- sapply(ci.lb, transf)
            ci.ub <- sapply(ci.ub, transf)
            ci.lb.rma <- sapply(ci.lb.rma, transf)
            ci.ub.rma <- sapply(ci.ub.rma, transf)
        }
        else {
            b <- sapply(b, transf, targs)
            br.rma <- sapply(b.rma, transf, targs)
            se <- NA
            se.rma <- NA
            ci.lb <- sapply(ci.lb, transf, targs)
            ci.ub <- sapply(ci.ub, transf, targs)
            ci.lb.rma <- sapply(ci.lb.rma, transf, targs)
            ci.ub.rma <- sapply(ci.ub.rma, transf, targs)
        }
    }
    res <- list(b = b, se = se, ci.lb = ci.lb, ci.ub = ci.ub, 
        b.rma = b.rma, se.rma = se.rma, ci.lb.rma = ci.lb.rma, 
        ci.ub.rma = ci.ub.rma, method = "DL", method.rma = x$method, 
        tau2 = tau2, tau2.rma = x$tau2, digits = digits)
    class(res) <- c("hc.rma.uni")
    return(res)
}
head.list.rma <-
function (x, n = 6L, ...) 
{
    stopifnot(length(n) == 1L)
    n <- if (n < 0L) {
        max(length(x[[1]]) + n, 0L)
    }
    else {
        min(n, length(x[[1]]))
    }
    x[seq_len(n), , drop = FALSE]
}
influence.rma.uni <-
function (model, digits, ...) 
{
    if (!is.element("rma.uni", class(model))) 
        stop("Argument 'model' must be an object of class \"rma.uni\".")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    x <- model
    if (missing(digits)) 
        digits <- x$digits
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
    pred.full <- x$X.f %*% x$b
    if (x$weighted) {
        if (is.null(x$weights)) {
            wi <- 1/(x$vi + x$tau2)
            W <- diag(wi, nrow = x$k, ncol = x$k)
            svb <- crossprod(x$X, W) %*% x$X/x$s2w
        }
        else {
            svb <- chol2inv(chol(x$vb))
            A <- diag(x$weights, nrow = x$k, ncol = x$k)
            stXAX <- .invcalc(X = x$X, W = A, k = x$k)
            H <- x$X %*% stXAX %*% t(x$X) %*% A
        }
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
    for (i in seq_len(x$k.f)[x$not.na]) {
        res <- try(rma(x$yi.f[-i], x$vi.f[-i], weights = x$weights.f[-i], 
            mods = cbind(x$X.f[-i, ]), method = x$method, weighted = x$weighted, 
            intercept = FALSE, knha = x$knha, control = x$control), 
            silent = TRUE)
        if (inherits(res, "try-error")) 
            next
        tau2.del[i] <- res$tau2
        QE.del[i] <- res$QE
        Xi <- matrix(x$X.f[i, ], nrow = 1)
        delpred[i] <- Xi %*% res$b
        vdelpred[i] <- Xi %*% tcrossprod(res$vb, Xi)
        if (x$weighted) {
            if (is.null(x$weights)) {
                dffits[i] <- (pred.full[i] - delpred[i])/sqrt(res$s2w * 
                  hii[i] * (tau2.del[i] + x$vi.f[i]))
            }
            else {
                dffits[i] <- (pred.full[i] - delpred[i])/sqrt(res$s2w * 
                  diag(H %*% diag(tau2.del[i] + x$vi, nrow = x$k, 
                    ncol = x$k) %*% t(H)))[i - sum(!x$not.na[1:i])]
            }
        }
        else {
            dffits[i] <- (pred.full[i] - delpred[i])/sqrt(res$s2w * 
                diag(H %*% diag(tau2.del[i] + x$vi, nrow = x$k, 
                  ncol = x$k) %*% t(H)))[i - sum(!x$not.na[1:i])]
        }
        dfbeta <- x$b - res$b
        if (x$weighted) {
            if (is.null(x$weights)) {
                vb.del <- .invcalc(X = x$X, W = diag(1/(x$vi + 
                  tau2.del[i]), nrow = x$k, ncol = x$k), k = x$k)
            }
            else {
                vb.del <- tcrossprod(stXAX, x$X) %*% A %*% diag(x$vi + 
                  tau2.del[i], nrow = x$k, ncol = x$k) %*% A %*% 
                  x$X %*% stXAX
            }
        }
        else {
            vb.del <- tcrossprod(stXX, x$X) %*% diag(x$vi + tau2.del[i], 
                nrow = x$k, ncol = x$k) %*% x$X %*% stXX
        }
        dfbetas[i, ] <- dfbeta/sqrt(res$s2w * diag(vb.del))
        cook.d[i] <- (crossprod(dfbeta, svb) %*% dfbeta)
        cov.r[i] <- det(res$vb)/det(x$vb)
    }
    delresid <- x$yi.f - delpred
    sedelresid <- sqrt(x$vi.f + vdelpred + tau2.del)
    standelres <- delresid/sedelresid
    options(na.action = "na.omit")
    weight[x$not.na] <- weights(x)
    options(na.action = na.act)
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
    rownames(out$inf) <- x$slab
    rownames(out$dfb) <- x$slab
    colnames(out$dfb) <- rownames(x$b)
    colnames(out$inf) <- c("rstudent", "dffits", "cook.d", "cov.r", 
        "tau2.del", "QE.del", "hat", "weight")
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
    order.vec <- order(psize, decreasing = TRUE)
    dat.t$yi <- dat.t$yi[order.vec]
    dat.c$yi <- dat.c$yi[order.vec]
    psize <- psize[order.vec]
    pch <- pch[order.vec]
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
function (x, digits, transf, targs, ...) 
{
    if (!is.element("rma.mh", class(x))) 
        stop("Argument 'x' must be an object of class \"rma.mh\".")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    if (!x$int.only) 
        stop("Method only applicable for models without moderators.")
    if (x$k == 1) 
        stop("Stopped because k = 1.")
    if (missing(digits)) 
        digits <- x$digits
    if (missing(transf)) 
        transf <- FALSE
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
    for (i in seq_len(x$k.f)[x$not.na]) {
        if (is.element(x$measure, c("RR", "OR", "RD"))) {
            res <- try(rma.mh(ai = x$ai.f[-i], bi = x$bi.f[-i], 
                ci = x$ci.f[-i], di = x$di.f[-i], measure = x$measure, 
                add = x$add, to = x$to, drop00 = x$drop00, correct = x$correct), 
                silent = TRUE)
        }
        else {
            res <- try(rma.mh(x1i = x$x1i.f[-i], x2i = x$x2i.f[-i], 
                t1i = x$t1i.f[-i], t2i = x$t2i.f[-i], measure = x$measure, 
                add = x$add, to = x$to, drop00 = x$drop00, correct = x$correct), 
                silent = TRUE)
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
    if (is.logical(transf) && transf && is.element(x$measure, 
        c("OR", "RR", "IRR"))) 
        transf <- exp
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
        transf <- TRUE
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
    out$transf <- transf
    class(out) <- c("list.rma")
    return(out)
}
leave1out.rma.peto <-
function (x, digits, transf, targs, ...) 
{
    if (!is.element("rma.peto", class(x))) 
        stop("Argument 'x' must be an object of class \"rma.peto\".")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    if (!x$int.only) 
        stop("Method only applicable for models without moderators.")
    if (x$k == 1) 
        stop("Stopped because k = 1.")
    if (missing(digits)) 
        digits <- x$digits
    if (missing(transf)) 
        transf <- FALSE
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
    for (i in seq_len(x$k.f)[x$not.na]) {
        res <- try(rma.peto(ai = x$ai.f[-i], bi = x$bi.f[-i], 
            ci = x$ci.f[-i], di = x$di.f[-i], add = x$add, to = x$to, 
            drop00 = x$drop00), silent = TRUE)
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
    if (is.logical(transf) && transf) 
        transf <- exp
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
        transf <- TRUE
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
    out$transf <- transf
    class(out) <- c("list.rma")
    return(out)
}
leave1out.rma.uni <-
function (x, digits, transf, targs, ...) 
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
    if (missing(digits)) 
        digits <- x$digits
    if (missing(transf)) 
        transf <- FALSE
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
    for (i in seq_len(x$k.f)[x$not.na]) {
        res <- try(rma(x$yi.f[-i], x$vi.f[-i], weights = x$weights.f[-i], 
            method = x$method, weighted = x$weighted, intercept = TRUE, 
            knha = x$knha, control = x$control), silent = TRUE)
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
        transf <- TRUE
    }
    ci.bounds <- cbind(ci.lb, ci.ub)
    rev.order <- ifelse(ci.ub < ci.lb, TRUE, FALSE)
    rev.order[is.na(rev.order)] <- FALSE
    ci.bounds[rev.order, ] <- ci.bounds[rev.order, 2:1]
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
    out$transf <- transf
    class(out) <- c("list.rma")
    return(out)
}
llplot <-
function (measure = "OR", ai, bi, ci, di, n1i, n2i, data, subset, 
    drop00 = TRUE, xvals = 1000, xlim, ylim, xlab, ylab, scale = TRUE, 
    lty, lwd, col, level = 99.99, refline = 0, ...) 
{
    if (!is.element(measure, c("OR"))) 
        stop("Currently only measure=\"OR\" can be specified.")
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
    alpha <- ifelse(level > 1, (100 - level)/100, 1 - level)
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
    k <- length(ai)
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
        ai <- ai[subset]
        bi <- bi[subset]
        ci <- ci[subset]
        di <- di[subset]
    }
    dat <- escalc(measure = "OR", ai = ai, bi = bi, ci = ci, 
        di = di, drop00 = drop00)
    yi <- dat$yi
    vi <- dat$vi
    ids <- seq_len(k)
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
    if (!is.null(subset)) {
        ids <- ids[subset]
        lty <- lty[subset]
        lwd <- lwd[subset]
        col <- col[subset]
        id0 <- id0[subset]
        id00 <- id00[subset]
    }
    k <- length(ai)
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
        lwd <- seq(from = 4, to = 0.2, length = k)[rank(vi)]
    if (is.null(col)) 
        col <- paste0("gray", round(seq(from = 0, to = 80, length = k))[rank(vi)])
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
    for (i in (1:k)[order(1/vi)]) {
        lines(logORs, lls[i, ], lty = lty[i], lwd = lwd[i], col = col[i], 
            ...)
    }
    if (is.numeric(refline)) 
        abline(v = refline, lty = "solid", lwd = 2, ...)
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
    digits, ...) 
{
    if (!is.element("rma.uni", class(x))) 
        stop("Argument 'x' must be an object of class \"rma.uni\".")
    if (missing(digits)) 
        digits <- x$digits
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
    if (x$int.only) {
        zval.perm <- try(rep(NA, iter), silent = TRUE)
        if (inherits(zval.perm, "try-error")) 
            stop("Number of iterations requested too large.")
        QM.perm <- try(rep(NA, iter), silent = TRUE)
        if (inherits(QM.perm, "try-error")) 
            stop("Number of iterations requested too large.")
        if (progbar) 
            pbar <- txtProgressBar(min = 0, max = iter, style = 3)
        if (exact) {
            signmat <- .gensigns(x$k)
            for (i in seq_len(iter)) {
                res <- try(rma(signmat[i, ] * x$yi, x$vi, weights = x$weights, 
                  method = x$method, weighted = x$weighted, intercept = TRUE, 
                  knha = x$knha, control = x$control, btt = 1), 
                  silent = FALSE)
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
                res <- try(rma(signs * x$yi, x$vi, weights = x$weights, 
                  method = x$method, weighted = x$weighted, intercept = TRUE, 
                  knha = x$knha, control = x$control, btt = 1), 
                  silent = FALSE)
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
        zval.perm <- suppressWarnings(try(matrix(NA, nrow = iter, 
            ncol = x$p), silent = TRUE))
        if (inherits(zval.perm, "try-error")) 
            stop("Number of iterations requested too large.")
        QM.perm <- try(rep(NA, iter), silent = TRUE)
        if (inherits(QM.perm, "try-error")) 
            stop("Number of iterations requested too large.")
        if (progbar) 
            pbar <- txtProgressBar(min = 0, max = iter, style = 3)
        if (exact) {
            permmat <- .genuperms(indices)
            for (i in seq_len(iter)) {
                res <- try(rma(x$yi, x$vi, weights = x$weights, 
                  mods = cbind(X[permmat[i, ], ]), method = x$method, 
                  weighted = x$weighted, intercept = FALSE, knha = x$knha, 
                  control = x$control, btt = x$btt), silent = FALSE)
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
                res <- try(rma(x$yi, x$vi, weights = x$weights, 
                  mods = cbind(X[sample(x$k), ]), method = x$method, 
                  weighted = x$weighted, intercept = FALSE, knha = x$knha, 
                  control = x$control, btt = x$btt), silent = FALSE)
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
        for (j in seq_len(x$p)) {
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
            which.dfb <- seq_len(x$p)
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
    ids <- switch(slab.style, `1` = x$ids, `2` = rownames(x$inf), 
        `3` = abbreviate(rownames(x$inf), ...))
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
        for (i in seq_len(length(which.inf))) {
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
                laxis(side = 1, at = seq_len(len.ids), labels = lab.ids, 
                  xlab = "", las = las, ...)
                labline(h = 0, lty = "dashed", ...)
                labline(h = c(qnorm(0.025), qnorm(0.975)), lty = "dotted", 
                  ...)
                if (na.act == "na.exclude" || na.act == "na.pass") 
                  llines(seq_len(len.ids)[not.na], zi[not.na], 
                    col = col.na, ...)
                llines(seq_len(len.ids), zi, ...)
                lpoints(seq_len(len.ids), zi, pch = pch, bg = bg, 
                  ...)
                lpoints(seq_len(len.ids)[ids.infl], zi[ids.infl], 
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
                laxis(side = 1, at = seq_len(len.ids), labels = lab.ids, 
                  xlab = "", las = las, ...)
                labline(h = 0, lty = "dashed", ...)
                labline(h = 3 * sqrt(x$p/(x$k - x$p)), lty = "dotted", 
                  ...)
                labline(h = -3 * sqrt(x$p/(x$k - x$p)), lty = "dotted", 
                  ...)
                if (na.act == "na.exclude" || na.act == "na.pass") 
                  llines(seq_len(len.ids)[not.na], zi[not.na], 
                    col = col.na, ...)
                llines(seq_len(len.ids), zi, ...)
                lpoints(seq_len(len.ids), zi, pch = pch, bg = bg, 
                  ...)
                lpoints(seq_len(len.ids)[ids.infl], zi[ids.infl], 
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
                laxis(side = 1, at = seq_len(len.ids), labels = lab.ids, 
                  xlab = "", las = las, ...)
                labline(h = qchisq(0.5, df = x$p), lty = "dotted", 
                  ...)
                if (na.act == "na.exclude" || na.act == "na.pass") 
                  llines(seq_len(len.ids)[not.na], zi[not.na], 
                    col = col.na, ...)
                llines(seq_len(len.ids), zi, ...)
                lpoints(seq_len(len.ids), zi, pch = pch, bg = bg, 
                  ...)
                lpoints(seq_len(len.ids)[ids.infl], zi[ids.infl], 
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
                laxis(side = 1, at = seq_len(len.ids), labels = lab.ids, 
                  xlab = "", las = las, ...)
                labline(h = 1, lty = "dashed", ...)
                if (na.act == "na.exclude" || na.act == "na.pass") 
                  llines(seq_len(len.ids)[not.na], zi[not.na], 
                    col = col.na, ...)
                llines(seq_len(len.ids), zi, ...)
                lpoints(seq_len(len.ids), zi, pch = pch, bg = bg, 
                  ...)
                lpoints(seq_len(len.ids)[ids.infl], zi[ids.infl], 
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
                laxis(side = 1, at = seq_len(len.ids), labels = lab.ids, 
                  xlab = "", las = las, ...)
                labline(h = x$tau2, lty = "dashed", ...)
                if (na.act == "na.exclude" || na.act == "na.pass") 
                  llines(seq_len(len.ids)[not.na], zi[not.na], 
                    col = col.na, ...)
                llines(seq_len(len.ids), zi, ...)
                lpoints(seq_len(len.ids), zi, pch = pch, bg = bg, 
                  ...)
                lpoints(seq_len(len.ids)[ids.infl], zi[ids.infl], 
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
                laxis(side = 1, at = seq_len(len.ids), labels = lab.ids, 
                  xlab = "", las = las, ...)
                labline(h = x$QE, lty = "dashed", ...)
                if (na.act == "na.exclude" || na.act == "na.pass") 
                  llines(seq_len(len.ids)[not.na], zi[not.na], 
                    col = col.na, ...)
                llines(seq_len(len.ids), zi, ...)
                lpoints(seq_len(len.ids), zi, pch = pch, bg = bg, 
                  ...)
                lpoints(seq_len(len.ids)[ids.infl], zi[ids.infl], 
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
                laxis(side = 1, at = seq_len(len.ids), labels = lab.ids, 
                  xlab = "", las = las, ...)
                labline(h = x$p/x$k, lty = "dashed", ...)
                labline(h = 3 * x$p/x$k, lty = "dotted", ...)
                if (na.act == "na.exclude" || na.act == "na.pass") 
                  llines(seq_len(len.ids)[not.na], zi[not.na], 
                    col = col.na, ...)
                llines(seq_len(len.ids), zi, ...)
                lpoints(seq_len(len.ids), zi, pch = pch, bg = bg, 
                  ...)
                lpoints(seq_len(len.ids)[ids.infl], zi[ids.infl], 
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
                laxis(side = 1, at = seq_len(len.ids), labels = lab.ids, 
                  xlab = "", las = las, ...)
                labline(h = 100/x$k, lty = "dashed", ...)
                if (na.act == "na.exclude" || na.act == "na.pass") 
                  llines(seq_len(len.ids)[not.na], zi[not.na], 
                    col = col.na, ...)
                llines(seq_len(len.ids), zi, ...)
                lpoints(seq_len(len.ids), zi, pch = pch, bg = bg, 
                  ...)
                lpoints(seq_len(len.ids)[ids.infl], zi[ids.infl], 
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
        for (i in seq_len(length(which.dfb))) {
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
                colnames(x$dfb)[which.dfb[i]]), xlab = "", ylab = "", 
                las = las, ...)
            laxis(side = 1, at = seq_len(len.ids), labels = lab.ids, 
                xlab = "", las = las, ...)
            labline(h = 0, lty = "dashed", ...)
            labline(h = 1, lty = "dotted", ...)
            labline(h = -1, lty = "dotted", ...)
            if (na.act == "na.exclude" || na.act == "na.pass") 
                llines(seq_len(len.ids)[not.na], zi[not.na], 
                  col = col.na, ...)
            llines(seq_len(len.ids), zi, ...)
            lpoints(seq_len(len.ids), zi, pch = pch, bg = bg, 
                ...)
            lpoints(seq_len(len.ids)[ids.infl], zi[ids.infl], 
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
        lines(seq_len(k)[not.na], z[not.na], col = "lightgray", 
            ...)
        lines(seq_len(k), z, ...)
        points(seq_len(k), z, pch = 21, bg = "black", ...)
        axis(side = 1, at = seq_len(k), labels = ids, ...)
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
        lines(seq_len(k)[not.na], z[not.na], col = "lightgray", 
            ...)
        lines(seq_len(k), z, ...)
        points(seq_len(k), z, pch = 21, bg = "black", ...)
        axis(side = 1, at = seq_len(k), labels = ids, ...)
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
            lines(seq_len(k)[not.na], z[not.na], col = "lightgray", 
                ...)
            lines(seq_len(k), z, ...)
            points(seq_len(k), z, pch = 21, bg = "black", ...)
            axis(side = 1, at = seq_len(k), labels = ids, ...)
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
            lines(seq_len(k)[not.na], z[not.na], col = "lightgray", 
                ...)
            lines(seq_len(k), z, ...)
            points(seq_len(k), z, pch = 21, bg = "black", ...)
            axis(side = 1, at = seq_len(k), labels = ids, ...)
            abline(h = 0, lty = "dashed", ...)
            abline(h = c(qnorm(0.025), qnorm(0.975)), lty = "dotted", 
                ...)
            title("Standardized Residuals", ...)
        }
    }
    invisible()
}
predict.rma <-
function (object, newmods, intercept, tau2.levels, addx = FALSE, 
    level, digits, transf, targs, ...) 
{
    if (!is.element("rma", class(object))) 
        stop("Argument 'object' must be an object of class \"rma\".")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    x <- object
    if (missing(newmods)) 
        newmods <- NULL
    if (missing(intercept)) 
        intercept <- x$intercept
    if (missing(tau2.levels)) 
        tau2.levels <- NULL
    if (missing(level)) 
        level <- x$level
    if (missing(digits)) 
        digits <- x$digits
    if (missing(transf)) 
        transf <- FALSE
    if (missing(targs)) 
        targs <- NULL
    alpha <- ifelse(level > 1, (100 - level)/100, 1 - level)
    if (x$knha || x$robust) {
        crit <- qt(alpha/2, df = x$k - x$p, lower.tail = FALSE)
    }
    else {
        crit <- qnorm(alpha/2, lower.tail = FALSE)
    }
    if (x$int.only && !is.null(newmods)) 
        stop("Cannot specify new moderator values for models without moderators.")
    if (is.null(newmods)) {
        if (!is.element("rma.mv", class(object))) {
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
            if (x$int.only) {
                if (!x$withG || is.element(x$struct, c("CS", 
                  "AR", "UNHO"))) {
                  k.new <- 1
                  X.new <- cbind(1)
                }
                else {
                  if (is.null(tau2.levels)) {
                    k.new <- x$g.nlevels.f[1]
                    X.new <- cbind(rep(1, k.new))
                    tau2.levels <- levels(x$mf.g.f$inner)
                  }
                  else {
                    k.new <- length(tau2.levels)
                    X.new <- cbind(rep(1, k.new))
                  }
                }
            }
            else {
                k.new <- x$k.f
                X.new <- x$X.f
                tau2.levels <- as.character(x$mf.g.f$inner)
            }
        }
    }
    else {
        if ((!x$int.incl && x$p == 1L) || (x$int.incl && x$p == 
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
                k.new <- nrow(newmods)
                X.new <- cbind(newmods)
            }
        }
        if (x$int.incl) {
            if (intercept) {
                X.new <- cbind(intrcpt = rep(1, k.new), X.new)
            }
            else {
                X.new <- cbind(intrcpt = rep(0, k.new), X.new)
            }
        }
        if (ncol(X.new) != x$p) 
            stop("Dimensions of 'newmods' do not match dimensions of the model.")
    }
    if (is.element("rma.mv", class(object)) && x$withG && is.element(x$struct, 
        c("HCS", "UN", "HAR"))) {
        if (is.null(tau2.levels)) 
            stop("Need to specify 'tau2.levels' argument.")
        if (!is.numeric(tau2.levels) && any(is.na(pmatch(tau2.levels, 
            x$g.levels.f[[1]], duplicates.ok = TRUE)))) 
            stop("Non-existing levels specified via 'tau2.levels' argument.")
        if (is.numeric(tau2.levels)) {
            tau2.levels <- round(tau2.levels)
            if (any(tau2.levels < 1) || any(tau2.levels > x$g.nlevels.f[1])) 
                stop("Non-existing tau^2 values specified via 'tau2.levels' argument.")
        }
        if (length(tau2.levels) == 1L) 
            tau2.levels <- rep(tau2.levels, k.new)
        if (length(tau2.levels) != k.new) 
            stop("Length of 'tau2.levels' does not match number of predicted values.")
    }
    pred <- rep(NA, k.new)
    vpred <- rep(NA, k.new)
    for (i in seq_len(k.new)) {
        Xi.new <- X.new[i, , drop = FALSE]
        pred[i] <- Xi.new %*% x$b
        vpred[i] <- Xi.new %*% tcrossprod(x$vb, Xi.new)
    }
    se <- sqrt(vpred)
    ci.lb <- pred - crit * se
    ci.ub <- pred + crit * se
    if (!is.element("rma.mv", class(object))) {
        cr.lb <- pred - crit * sqrt(vpred + x$tau2)
        cr.ub <- pred + crit * sqrt(vpred + x$tau2)
    }
    else {
        if (x$withG) {
            if (is.element(x$struct, c("CS", "AR", "UNHO"))) {
                cr.lb <- pred - crit * sqrt(vpred + sum(x$sigma2) + 
                  x$tau2)
                cr.ub <- pred + crit * sqrt(vpred + sum(x$sigma2) + 
                  x$tau2)
            }
            else {
                if (is.numeric(tau2.levels)) {
                  cr.lb <- pred - crit * sqrt(vpred + sum(x$sigma2) + 
                    x$tau2[tau2.levels])
                  cr.ub <- pred + crit * sqrt(vpred + sum(x$sigma2) + 
                    x$tau2[tau2.levels])
                  tau2.levels <- x$g.levels.f[[1]][tau2.levels]
                }
                else {
                  cr.lb <- pred - crit * sqrt(vpred + sum(x$sigma2) + 
                    x$tau2[pmatch(tau2.levels, x$g.levels.f[[1]], 
                      duplicates.ok = TRUE)])
                  cr.ub <- pred + crit * sqrt(vpred + sum(x$sigma2) + 
                    x$tau2[pmatch(tau2.levels, x$g.levels.f[[1]], 
                      duplicates.ok = TRUE)])
                  tau2.levels <- x$g.levels.f[[1]][pmatch(tau2.levels, 
                    x$g.levels.f[[1]], duplicates.ok = TRUE)]
                }
            }
        }
        else {
            cr.lb <- pred - crit * sqrt(vpred + sum(x$sigma2))
            cr.ub <- pred + crit * sqrt(vpred + sum(x$sigma2))
        }
    }
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
        transf <- TRUE
    }
    ci.bounds <- cbind(ci.lb, ci.ub)
    rev.order <- ifelse(ci.ub < ci.lb, TRUE, FALSE)
    rev.order[is.na(rev.order)] <- FALSE
    ci.bounds[rev.order, ] <- ci.bounds[rev.order, 2:1]
    ci.lb <- ci.bounds[, 1]
    ci.ub <- ci.bounds[, 2]
    cr.bounds <- cbind(cr.lb, cr.ub)
    rev.order <- ifelse(cr.ub < cr.lb, TRUE, FALSE)
    rev.order[is.na(rev.order)] <- FALSE
    cr.bounds[rev.order, ] <- cr.bounds[rev.order, 2:1]
    cr.lb <- cr.bounds[, 1]
    cr.ub <- cr.bounds[, 2]
    if (is.null(newmods) && !x$int.only) {
        slab <- x$slab
    }
    else {
        slab <- seq_len(k.new)
    }
    if (x$int.only) {
        if (k.new == 1L) {
            slab <- ""
        }
        else {
            slab <- seq_len(k.new)
        }
    }
    if (na.act == "na.omit") {
        not.na <- !is.na(pred)
    }
    else {
        not.na <- rep(TRUE, k.new)
    }
    out <- list(pred = pred[not.na], se = se[not.na], ci.lb = ci.lb[not.na], 
        ci.ub = ci.ub[not.na], cr.lb = cr.lb[not.na], cr.ub = cr.ub[not.na])
    if (is.element("rma.mv", class(object)) && x$withG && is.element(x$struct, 
        c("HCS", "UN", "HAR"))) 
        out$tau2.levels <- tau2.levels
    if (addx) 
        out$X <- matrix(X.new[not.na, ], ncol = x$p)
    out$slab <- slab[not.na]
    if (na.act == "na.fail" && any(!x$not.na)) 
        stop("Missing values in results.")
    if (addx) 
        colnames(out$X) <- colnames(x$X)
    if (x$method == "FE") {
        out$cr.lb <- NULL
        out$cr.ub <- NULL
    }
    out$digits <- digits
    out$method <- x$method
    out$transf <- transf
    class(out) <- c("list.rma")
    return(out)
}
print.anova.rma.uni <-
function (x, digits, ...) 
{
    if (class(x) != "anova.rma.uni") 
        stop("Argument 'x' must be an object of class \"anova.rma.uni\".")
    if (missing(digits)) 
        digits <- x$digits
    cutoff <- paste(c(".", rep(0, digits - 1), 1), collapse = "")
    ncutoff <- as.numeric(cutoff)
    if (x$test == "Wald") {
        QMp <- x$QMp
        if (QMp > ncutoff) {
            QMp <- paste0("= ", formatC(QMp, digits = digits, 
                format = "f"))
        }
        else {
            QMp <- paste0("< ", cutoff)
        }
        cat("\n")
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
    else {
        res.table <- rbind(c(x$p.f, x$fit.stats.f[3], x$fit.stats.f[4], 
            x$fit.stats.f[5], x$fit.stats.f[1], NA, NA, x$QE.f, 
            x$tau2.f, NA), c(x$p.r, x$fit.stats.r[3], x$fit.stats.r[4], 
            x$fit.stats.r[5], x$fit.stats.r[1], x$LRT, x$pval, 
            x$QE.r, x$tau2.r, NA))
        res.table[, seq.int(from = 2, to = 10)] <- formatC(res.table[, 
            seq.int(from = 2, to = 10)], digits = digits, format = "f")
        colnames(res.table) <- c("df", "AIC", "BIC", "AICc", 
            "logLik", "LRT", "pval", "QE", "tau^2", "R^2")
        rownames(res.table) <- c("Full", "Reduced")
        pval <- x$pval
        if (pval > ncutoff) {
            res.table[2, 7] <- formatC(pval, digits = digits, 
                format = "f")
        }
        else {
            res.table[2, 7] <- paste0("<", cutoff)
        }
        res.table[1, c(6, 7)] <- ""
        res.table[1, 10] <- ""
        res.table[2, 10] <- paste0(x$R2, "%")
        if (x$method == "FE") {
            res.table <- res.table[, seq_len(8)]
        }
        print(res.table, quote = FALSE, right = TRUE)
    }
    invisible()
}
print.confint.rma <-
function (x, digits, ...) 
{
    if (!is.element("confint.rma", class(x))) 
        stop("Argument 'x' must be an object of class \"confint.rma\".")
    if (missing(digits)) 
        digits <- x$digits
    cat("\n")
    if (names(x)[1] == "fixed") {
        res.fixed <- formatC(x$fixed, digits = digits, format = "f")
        print(res.fixed, quote = FALSE, right = TRUE)
    }
    if (names(x)[1] == "random" || names(x)[2] == "random") {
        if (names(x)[1] == "fixed") 
            cat("\n")
        res.random <- formatC(x$random, digits = digits, format = "f")
        print(res.random, quote = FALSE, right = TRUE)
        if (is.na(x$random[1, 2]) && is.na(x$random[1, 3])) 
            message("\nThe upper and lower CI bounds for tau^2 both fall below ", 
                x$tau2.min, ".\nThe CIs are therefore equal to the null/empty set.", 
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
    if (missing(digits)) 
        digits <- attr(x, "digits")
    if (is.null(digits)) 
        digits <- 4
    yi.pos <- na.omit(match(attr(x, "yi.names"), names(x)))
    vi.pos <- na.omit(match(attr(x, "vi.names"), names(x)))
    sei.pos <- na.omit(match(attr(x, "sei.names"), names(x)))
    zi.pos <- na.omit(match(attr(x, "zi.names"), names(x)))
    ci.lb.pos <- na.omit(match(attr(x, "ci.lb.names"), names(x)))
    ci.ub.pos <- na.omit(match(attr(x, "ci.ub.names"), names(x)))
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
    print(x, ...)
}
print.fsn <-
function (x, digits, ...) 
{
    if (class(x) != "fsn") 
        stop("Argument 'x' must be an object of class \"fsn\".")
    if (missing(digits)) 
        digits <- x$digits
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
            pval <- paste0("<", cutoff)
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
            pval <- paste0("<", cutoff)
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
print.hc.rma.uni <-
function (x, digits, ...) 
{
    if (!is.element("hc.rma.uni", class(x))) 
        stop("Argument 'x' must be an object of class \"hc.rma.uni\".")
    if (missing(digits)) 
        digits <- x$digits
    cutoff <- paste(c(".", rep(0, digits - 1), 1), collapse = "")
    ncutoff <- as.numeric(cutoff)
    res.table <- data.frame(method = c(x$method.rma, x$method), 
        tau2 = formatC(c(x$tau2.rma, x$tau2), digits = digits, 
            format = "f"), estimate = formatC(c(x$b.rma, x$b), 
            digits = digits, format = "f"), se = c(ifelse(is.na(x$se.rma), 
            NA, formatC(x$se.rma, digits = digits, format = "f")), 
            ifelse(is.na(x$se), NA, formatC(x$se, digits = digits, 
                format = "f"))), ci.lb = formatC(c(x$ci.lb.rma, 
            x$ci.lb), digits = digits, format = "f"), ci.ub = formatC(c(x$ci.ub.rma, 
            x$ci.ub), digits = digits, format = "f"), stringsAsFactors = FALSE)
    rownames(res.table) <- c("rma", "hc")
    print(res.table, quote = FALSE, right = TRUE)
    invisible(res.table)
}
print.infl.rma.uni <-
function (x, digits, ...) 
{
    if (class(x) != "infl.rma.uni") 
        stop("Argument 'x' must be an object of class \"infl.rma.uni\".")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    if (missing(digits)) 
        digits <- x$digits
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
function (x, digits, ...) 
{
    if (!is.element("list.rma", class(x))) 
        stop("Argument 'x' must be an object of class \"list.rma\".")
    if (missing(digits)) 
        digits <- x$digits
    attr(x, "class") <- NULL
    slab.pos <- which(names(x) == "slab")
    out <- x[seq_len(slab.pos - 1)]
    out <- data.frame(out, row.names = x$slab)
    if (nrow(out) == 0L) 
        stop("All values are NA.", call. = FALSE)
    transf.true <- 0
    if (exists("transf", where = x, inherits = FALSE) && x$transf) {
        transf.true <- 1
        out$se <- NULL
    }
    if (exists("method", where = x, inherits = FALSE)) {
        min.pos <- slab.pos - is.element("tau2.levels", names(x)) - 
            is.element("X", names(x)) - transf.true
    }
    else {
        min.pos <- slab.pos - transf.true
    }
    out[, seq_len(min.pos - 1)] <- apply(out[, seq_len(min.pos - 
        1), drop = FALSE], 2, formatC, digits = digits, format = "f")
    print(out, quote = FALSE, right = TRUE)
}
print.permutest.rma.uni <-
function (x, digits, signif.stars = getOption("show.signif.stars"), 
    signif.legend = signif.stars, ...) 
{
    if (class(x) != "permutest.rma.uni") 
        stop("Argument 'x' must be an object of class \"permutest.rma.uni\".")
    if (missing(digits)) 
        digits <- x$digits
    cutoff <- paste(c(".", rep(0, digits - 1), 1), collapse = "")
    ncutoff <- as.numeric(cutoff)
    QMp <- x$QMp
    if (QMp > ncutoff) {
        QMp <- paste("=", formatC(QMp, digits = digits, format = "f"))
    }
    else {
        QMp <- paste0("< ", cutoff)
    }
    cat("\n")
    if (!x$int.only) {
        cat("Test of Moderators (coefficient(s) ", paste0(x$btt, 
            collapse = ","), "): \n")
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
    res.table <- cbind(estimate = x$b, se = x$se, zval = x$zval, 
        `pval*` = x$pval, ci.lb = x$ci.lb, ci.ub = x$ci.ub)
    if (x$knha || x$robust) 
        colnames(res.table)[3] <- "tval"
    signif <- symnum(x$pval, corr = FALSE, na = FALSE, cutpoints = c(0, 
        0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", 
        "*", ".", " "))
    if (signif.stars) {
        res.table <- cbind(formatC(res.table, digits = digits, 
            format = "f"), signif)
        colnames(res.table)[7] <- ""
    }
    else {
        res.table <- formatC(res.table, digits = digits, format = "f")
    }
    res.table[x$pval > ncutoff, 4] <- formatC(x$pval[x$pval > 
        ncutoff], digits = digits, format = "f")
    res.table[x$pval < ncutoff, 4] <- paste0("<", cutoff)
    cat("Model Results:")
    cat("\n\n")
    print(res.table, quote = FALSE, right = TRUE, print.gap = 2)
    cat("\n")
    if (signif.legend) 
        cat("---\nSignif. codes: ", attr(signif, "legend"), "\n\n")
    invisible()
}
print.ranktest.rma <-
function (x, digits, ...) 
{
    if (class(x) != "ranktest.rma") 
        stop("Argument 'x' must be an object of class \"ranktest.rma\".")
    if (missing(digits)) 
        digits <- x$digits
    cutoff <- paste(c(".", rep(0, digits - 1), 1), collapse = "")
    ncutoff <- as.numeric(cutoff)
    pval <- x$pval
    if (pval > ncutoff) {
        pval <- paste("=", formatC(pval, digits = digits, format = "f"))
    }
    else {
        pval <- paste0("< ", cutoff)
    }
    cat("\n")
    cat("Rank Correlation Test for Funnel Plot Asymmetry\n\n")
    cat("Kendall's tau = ", formatC(x$tau, digits = digits, format = "f"), 
        ", p ", pval, "\n\n", sep = "")
    invisible()
}
print.regtest.rma <-
function (x, digits, ret.fit, ...) 
{
    if (class(x) != "regtest.rma") 
        stop("Argument 'x' must be an object of class \"regtest.rma\".")
    if (missing(digits)) 
        digits <- x$digits
    if (missing(ret.fit)) 
        ret.fit <- x$ret.fit
    cutoff <- paste(c(".", rep(0, digits - 1), 1), collapse = "")
    ncutoff <- as.numeric(cutoff)
    pval <- x$pval
    if (pval > ncutoff) {
        pval <- paste("=", formatC(pval, digits = digits, format = "f"))
    }
    else {
        pval <- paste0("< ", cutoff)
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
function (x, digits, showfit = FALSE, signif.stars = getOption("show.signif.stars"), 
    signif.legend = signif.stars, ...) 
{
    if (!is.element("rma.glmm", class(x))) 
        stop("Argument 'x' must be an object of class \"rma.glmm\".")
    if (missing(digits)) 
        digits <- x$digits
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
        names(fs) <- c("logLik", "deviance", "AIC", "BIC", "AICc")
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
                ifelse(is.na(x$se.tau2), "", paste0(" (SE = ", 
                  formatC(x$se.tau2, digits = digits, format = "f"), 
                  ")")), "\n", sep = "")
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
                ifelse(is.na(x$se.tau2), "", paste0(" (SE = ", 
                  formatC(x$se.tau2, digits = digits, format = "f"), 
                  ")")), "\n", sep = "")
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
            QEp.Wld <- paste0("< ", cutoff)
        }
        if (QEp.LRT > ncutoff) {
            QEp.LRT <- paste("=", formatC(QEp.LRT, digits = digits, 
                format = "f"))
        }
        else {
            QEp.LRT <- paste0("< ", cutoff)
        }
        QE.Wld <- formatC(round(x$QE.Wld, digits = digits), digits = digits, 
            format = "f")
        QE.LRT <- formatC(round(x$QE.LRT, digits = digits), digits = digits, 
            format = "f")
        if (nchar(QE.Wld) > nchar(QE.LRT)) 
            QE.LRT <- paste0(paste(rep(" ", nchar(QE.Wld) - nchar(QE.LRT)), 
                collapse = ""), QE.LRT)
        if (nchar(QE.LRT) > nchar(QE.Wld)) 
            QE.Wld <- paste0(paste(rep(" ", nchar(QE.LRT) - nchar(QE.Wld)), 
                collapse = ""), QE.Wld)
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
        QMp <- paste0("< ", cutoff)
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
        res.table <- c(estimate = x$b, se = x$se, zval = x$zval, 
            pval = x$pval, ci.lb = x$ci.lb, ci.ub = x$ci.ub)
        if (x$knha || x$robust) 
            names(res.table)[3] <- "tval"
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
        res.table[4][x$pval < ncutoff] <- paste0("<", cutoff)
    }
    else {
        res.table <- cbind(estimate = x$b, se = x$se, zval = x$zval, 
            pval = x$pval, ci.lb = x$ci.lb, ci.ub = x$ci.ub)
        if (x$knha || x$robust) 
            colnames(res.table)[3] <- "tval"
        signif <- symnum(x$pval, corr = FALSE, na = FALSE, cutpoints = c(0, 
            0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", 
            "*", ".", " "))
        if (signif.stars) {
            res.table <- cbind(formatC(res.table, digits = digits, 
                format = "f"), signif)
            colnames(res.table)[7] <- ""
        }
        else {
            res.table <- formatC(res.table, digits = digits, 
                format = "f")
        }
        res.table[x$pval > ncutoff, 4] <- formatC(x$pval[x$pval > 
            ncutoff], digits = digits, format = "f")
        res.table[x$pval < ncutoff, 4] <- paste0("<", cutoff)
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
function (x, digits, showfit = FALSE, ...) 
{
    if (!is.element("rma.mh", class(x))) 
        stop("Argument 'x' must be an object of class \"rma.mh\".")
    if (missing(digits)) 
        digits <- x$digits
    cutoff <- paste(c(".", rep(0, digits - 1), 1), collapse = "")
    ncutoff <- as.numeric(cutoff)
    cat("\n")
    cat("Fixed-Effects Model (k = ", x$k, ")", sep = "")
    if (showfit) {
        cat("\n")
        fs <- c(formatC(x$fit.stats$ML, digits = digits, format = "f"))
        names(fs) <- c("logLik", "deviance", "AIC", "BIC", "AICc")
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
            QEp <- paste0("< ", cutoff)
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
            res.table[4][x$pval < ncutoff] <- paste0("<", cutoff)
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
            if (is.na(x$MH)) {
                cat("Cochran-Mantel-Haenszel Test:     test value not computable for these data \n", 
                  sep = "")
            }
            else {
                pval <- x$MHp
                if (pval > ncutoff) {
                  pval <- paste("=", formatC(pval, digits = digits, 
                    format = "f"))
                }
                else {
                  pval <- paste0("< ", cutoff)
                }
                cat("Cochran-Mantel-Haenszel Test:     CMH = ", 
                  formatC(x$MH, digits, format = "f"), ", df = 1,", 
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
                  pval <- paste0("< ", cutoff)
                }
                cat("Tarone's Test for Heterogeneity:  X^2 = ", 
                  formatC(x$TA, digits, format = "f"), ", df = ", 
                  x$k.pos - 1, ", p-val ", pval, "\n\n", sep = "")
            }
        }
        if (x$measure == "IRR") {
            if (is.na(x$MH)) {
                cat("Mantel-Haenszel Test:     test value not computable for these data \n", 
                  sep = "")
            }
            else {
                pval <- x$MHp
                if (pval > ncutoff) {
                  pval <- paste("=", formatC(pval, digits = digits, 
                    format = "f"))
                }
                else {
                  pval <- paste0("< ", cutoff)
                }
                cat("Mantel-Haenszel Test: MH = ", formatC(x$MH, 
                  digits, format = "f"), ", df = 1, p-val ", 
                  pval, "\n\n", sep = "")
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
            res.table[4][x$pval < ncutoff] <- paste0("<", cutoff)
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
print.rma.mv <-
function (x, digits, showfit = FALSE, signif.stars = getOption("show.signif.stars"), 
    signif.legend = signif.stars, ...) 
{
    if (!is.element("rma.mv", class(x))) 
        stop("Argument 'x' must be an object of class \"rma.mv\".")
    if (missing(digits)) 
        digits <- x$digits
    cutoff <- paste(c(".", rep(0, digits - 1), 1), collapse = "")
    ncutoff <- as.numeric(cutoff)
    cat("\n")
    cat("Multivariate Meta-Analysis Model (k = ", x$k, "; ", 
        sep = "")
    cat("method: ", x$method, ")", sep = "")
    if (showfit) {
        cat("\n")
        if (x$method == "REML") {
            fs <- c(formatC(round(x$fit.stats$REML, digits = digits), 
                digits = digits, format = "f"))
        }
        else {
            fs <- c(formatC(round(x$fit.stats$ML, digits = digits), 
                digits = digits, format = "f"))
        }
        names(fs) <- c("logLik", "Deviance", "AIC", "BIC", "AICc")
        cat("\n")
        print(fs, quote = FALSE, print.gap = 2)
        cat("\n")
    }
    else {
        cat("\n\n")
    }
    sigma2 <- formatC(x$sigma2, digits = digits, format = "f")
    tau2 <- formatC(x$tau2, digits = digits, format = "f")
    rho <- formatC(x$rho, digits = digits, format = "f")
    sigma <- formatC(sqrt(x$sigma2), digits = digits, format = "f")
    tau <- formatC(sqrt(x$tau2), digits = digits, format = "f")
    cat("Variance Components: ")
    right <- TRUE
    if (!x$withS && !x$withG) {
        cat("none\n\n")
    }
    else {
        cat("\n\n")
        if (x$withS) {
            vc <- cbind(estim = sigma2, sqrt = sigma, nlvls = x$s.nlevels, 
                fixed = ifelse(x$vc.fix$sigma2, "yes", "no"), 
                factor = x$s.names, R = ifelse(x$Rfix, "yes", 
                  "no"))
            if (!x$withR) 
                vc <- vc[, -6, drop = FALSE]
            if (length(x$sigma2) == 1) {
                rownames(vc) <- "sigma^2  "
            }
            else {
                rownames(vc) <- paste("sigma^2.", 1:length(x$sigma2), 
                  sep = "")
            }
            print(vc, quote = FALSE, right = right, print.gap = 2)
            cat("\n")
        }
        if (x$withG) {
            mng <- max(nchar(x$g.names))
            cat("outer factor: ", paste0(x$g.names[2], paste(rep(" ", 
                max(0, mng - nchar(x$g.names[2]))), collapse = ""), 
                collapse = ""), " (nlvls = ", x$g.nlevels[2], 
                ")\n", sep = "")
            cat("inner factor: ", paste0(x$g.names[1], paste(rep(" ", 
                max(0, mng - nchar(x$g.names[1]))), collapse = ""), 
                collapse = ""), " (nlvls = ", x$g.nlevels.f[1], 
                ")\n", sep = "")
            cat("\n")
            if (is.element(x$struct, c("CS", "AR"))) {
                vc <- cbind(tau2, tau, ifelse(x$vc.fix$tau2, 
                  "yes", "no"))
                vc <- rbind(vc, c(rho, "", ifelse(x$vc.fix$rho, 
                  "yes", "no")))
                colnames(vc) <- c("estim", "sqrt", "fixed")
                rownames(vc) <- c("tau^2    ", "rho")
                print(vc, quote = FALSE, right = right, print.gap = 2)
            }
            if (is.element(x$struct, c("HCS", "HAR"))) {
                vc <- cbind(tau2, tau, x$g.levels.k, ifelse(x$vc.fix$tau2, 
                  "yes", "no"), x$g.levels.f[[1]])
                vc <- rbind(vc, c(rho, "", "", ifelse(x$vc.fix$rho, 
                  "yes", "no"), ""))
                colnames(vc) <- c("estim", "sqrt", "k.lvl", "fixed", 
                  "level")
                if (length(x$tau2) == 1) {
                  rownames(vc) <- c("tau^2   ", "rho")
                }
                else {
                  rownames(vc) <- c(paste("tau^2.", 1:length(x$tau2), 
                    "  ", sep = ""), "rho")
                }
                print(vc, quote = FALSE, right = right, print.gap = 2)
            }
            if (is.element(x$struct, c("UN", "UNHO"))) {
                if (x$struct == "UN") {
                  vc <- cbind(tau2, tau, x$g.levels.k, ifelse(x$vc.fix$tau2, 
                    "yes", "no"), x$g.levels.f[[1]])
                }
                else {
                  vc <- cbind(rep(tau2, length(x$g.levels.k)), 
                    rep(tau, length(x$g.levels.k)), x$g.levels.k, 
                    ifelse(rep(x$vc.fix$tau2, length(x$g.levels.k)), 
                      "yes", "no"), x$g.levels.f[[1]])
                }
                colnames(vc) <- c("estim", "sqrt", "k.lvl", "fixed", 
                  "level")
                if (length(x$g.levels.k) == 1) {
                  rownames(vc) <- c("tau^2")
                }
                else {
                  rownames(vc) <- paste("tau^2.", 1:length(x$g.levels.k), 
                    "  ", sep = "")
                }
                print(vc, quote = FALSE, right = right, print.gap = 2)
                cat("\n")
                if (length(x$rho) == 1) {
                  G <- matrix(NA, nrow = 2, ncol = 2)
                }
                else {
                  G <- matrix(NA, nrow = x$g.nlevels.f[1], ncol = x$g.nlevels.f[1])
                }
                G[upper.tri(G)] <- rho
                G[lower.tri(G)] <- t(G)[lower.tri(G)]
                diag(G) <- 1
                if (length(x$rho) == 1) {
                  G.info <- matrix(NA, nrow = 2, ncol = 2)
                }
                else {
                  G.info <- matrix(NA, nrow = x$g.nlevels.f[1], 
                    ncol = x$g.nlevels.f[1])
                }
                G.info[upper.tri(G.info)] <- x$g.levels.comb.k
                G.info[lower.tri(G.info)] <- t(G.info)[lower.tri(G.info)]
                G.info[upper.tri(G.info)] <- ifelse(x$vc.fix$rho, 
                  "yes", "no")
                diag(G.info) <- "-"
                vc <- cbind(G, "", G.info)
                colnames(vc) <- c(paste("rho.", abbreviate(x$g.levels.f[[1]]), 
                  sep = ""), "", abbreviate(x$g.levels.f[[1]]))
                rownames(vc) <- x$g.levels.f[[1]]
                print(vc, quote = FALSE, right = right, print.gap = 2)
            }
            cat("\n")
        }
    }
    if (!is.na(x$QE)) {
        QEp <- x$QEp
        if (QEp > ncutoff) {
            QEp <- paste("=", formatC(QEp, digits = digits, format = "f"))
        }
        else {
            QEp <- paste0("< ", cutoff)
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
        QMp <- paste0("< ", cutoff)
    }
    if (x$p > 1) {
        cat("Test of Moderators (coefficient(s) ", paste(x$btt, 
            collapse = ","), "): \n", sep = "")
        cat("QM(df = ", x$m, ") = ", formatC(x$QM, digits = digits, 
            format = "f"), ", p-val ", QMp, "\n\n", sep = "")
    }
    if (x$int.only) {
        res.table <- c(estimate = x$b, se = x$se, zval = x$zval, 
            pval = x$pval, ci.lb = x$ci.lb, ci.ub = x$ci.ub)
        if (x$knha || x$robust) 
            names(res.table)[3] <- "tval"
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
        res.table[4][x$pval < ncutoff] <- paste0("<", cutoff)
    }
    else {
        res.table <- cbind(estimate = x$b, se = x$se, zval = x$zval, 
            pval = x$pval, ci.lb = x$ci.lb, ci.ub = x$ci.ub)
        if (x$knha || x$robust) 
            colnames(res.table)[3] <- "tval"
        signif <- symnum(x$pval, corr = FALSE, na = FALSE, cutpoints = c(0, 
            0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", 
            "*", ".", " "))
        if (signif.stars) {
            res.table <- cbind(formatC(res.table, digits = digits, 
                format = "f"), signif)
            colnames(res.table)[7] <- ""
        }
        else {
            res.table <- formatC(res.table, digits = digits, 
                format = "f")
        }
        res.table[x$pval > ncutoff, 4] <- formatC(x$pval[x$pval > 
            ncutoff], digits = digits, format = "f")
        res.table[x$pval < ncutoff, 4] <- paste0("<", cutoff)
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
print.rma.peto <-
function (x, digits, showfit = FALSE, ...) 
{
    if (!is.element("rma.peto", class(x))) 
        stop("Argument 'x' must be an object of class \"rma.peto\".")
    if (missing(digits)) 
        digits <- x$digits
    cutoff <- paste(c(".", rep(0, digits - 1), 1), collapse = "")
    ncutoff <- as.numeric(cutoff)
    cat("\n")
    cat("Fixed-Effects Model (k = ", x$k, ")", sep = "")
    if (showfit) {
        cat("\n")
        fs <- c(formatC(x$fit.stats$ML, digits = digits, format = "f"))
        names(fs) <- c("logLik", "deviance", "AIC", "BIC", "AICc")
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
            QEp <- paste0("< ", cutoff)
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
        res.table[4][x$pval < ncutoff] <- paste0("<", cutoff)
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
function (x, digits, showfit = FALSE, signif.stars = getOption("show.signif.stars"), 
    signif.legend = signif.stars, ...) 
{
    if (!is.element("rma.uni", class(x))) 
        stop("Argument 'x' must be an object of class \"rma.uni\".")
    if (missing(digits)) 
        digits <- x$digits
    cutoff <- paste(c(".", rep(0, digits - 1), 1), collapse = "")
    ncutoff <- as.numeric(cutoff)
    if (is.element("rma.uni.trimfill", class(x))) {
        cat("\n")
        cat("Estimated number of missing studies on the ", x$side, 
            " side: ", x$k0, " (SE = ", ifelse(is.na(x$se.k0), 
                NA, formatC(x$se.k0, digits = digits, format = "f")), 
            ")\n", sep = "")
        if (x$k0.est == "R0") {
            p.k0 <- x$p.k0
            if (p.k0 > ncutoff) {
                p.k0 <- paste("=", formatC(p.k0, digits = digits, 
                  format = "f"))
            }
            else {
                p.k0 <- paste0("< ", cutoff)
            }
            cat("Test of H0: no missing studies on the ", x$side, 
                " side: p-val ", p.k0, "\n", sep = "")
        }
    }
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
        }
        else {
            fs <- c(formatC(round(x$fit.stats$ML, digits = digits), 
                digits = digits, format = "f"))
        }
        names(fs) <- c("logLik", "deviance", "AIC", "BIC", "AICc")
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
                ifelse(is.na(x$se.tau2), "", paste0(" (SE = ", 
                  formatC(x$se.tau2, digits = digits, format = "f"), 
                  ")")), "\n", sep = "")
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
                ifelse(is.na(x$se.tau2), "", paste0(" (SE = ", 
                  formatC(x$se.tau2, digits = digits, format = "f"), 
                  ")")), "\n", sep = "")
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
            if (!is.null(x$R2)) 
                cat("\nR^2 (amount of heterogeneity accounted for):            ", 
                  ifelse(is.na(x$R2), NA, formatC(x$R2, digits = 2, 
                    format = "f")), "%", sep = "")
        }
        cat("\n\n")
    }
    if (!is.na(x$QE)) {
        QEp <- x$QEp
        if (QEp > ncutoff) {
            QEp <- paste("=", formatC(QEp, digits = digits, format = "f"))
        }
        else {
            QEp <- paste0("< ", cutoff)
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
        QMp <- paste0("< ", cutoff)
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
        res.table <- c(estimate = x$b, se = x$se, zval = x$zval, 
            pval = x$pval, ci.lb = x$ci.lb, ci.ub = x$ci.ub)
        if (x$knha || x$robust) 
            names(res.table)[3] <- "tval"
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
        res.table[4][x$pval < ncutoff] <- paste0("<", cutoff)
    }
    else {
        res.table <- cbind(estimate = x$b, se = x$se, zval = x$zval, 
            pval = x$pval, ci.lb = x$ci.lb, ci.ub = x$ci.ub)
        if (x$knha || x$robust) 
            colnames(res.table)[3] <- "tval"
        signif <- symnum(x$pval, corr = FALSE, na = FALSE, cutpoints = c(0, 
            0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", 
            "*", ".", " "))
        if (signif.stars) {
            res.table <- cbind(formatC(res.table, digits = digits, 
                format = "f"), signif)
            colnames(res.table)[7] <- ""
        }
        else {
            res.table <- formatC(res.table, digits = digits, 
                format = "f")
        }
        res.table[x$pval > ncutoff, 4] <- formatC(x$pval[x$pval > 
            ncutoff], digits = digits, format = "f")
        res.table[x$pval < ncutoff, 4] <- paste0("<", cutoff)
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
print.summary.rma <-
function (x, digits, showfit = TRUE, signif.stars = getOption("show.signif.stars"), 
    signif.legend = signif.stars, ...) 
{
    if (!is.element("summary.rma", class(x))) 
        stop("Argument 'x' must be an object of class \"summary.rma\".")
    if (missing(digits)) 
        digits <- x$digits
    class(x) <- class(x)[-1]
    print(x, digits = digits, showfit = showfit, signif.stars = signif.stars, 
        signif.legend = signif.legend, ...)
    invisible()
}
profile.rma.mv <-
function (fitted, sigma2, tau2, rho, xlim, ylim, steps = 20, 
    startmethod = "init", progbar = TRUE, plot = TRUE, pch = 19, 
    ...) 
{
    if (!is.element("rma.mv", class(fitted))) 
        stop("Argument 'fitted' must be an object of class \"rma.mv\".")
    if (steps < 2) 
        stop("Argument 'steps' must be >= 2.")
    x <- fitted
    if (missing(sigma2) && missing(tau2) && missing(rho)) 
        stop("Must specify one of the arguments 'sigma2', 'tau2', or 'rho'.")
    if (sum(!missing(sigma2), !missing(tau2), !missing(rho)) > 
        1L) 
        stop("Must specify only one of the arguments 'sigma2', 'tau2', or 'rho'.")
    if (!missing(sigma2) && (all(is.na(x$vc.fix$sigma2)) || all(x$vc.fix$sigma2))) 
        stop("Model does not contain any (estimated) 'sigma2' components.")
    if (!missing(tau2) && (all(is.na(x$vc.fix$tau2)) || all(x$vc.fix$tau2))) 
        stop("Model does not contain any (estimated) 'tau2' components.")
    if (!missing(rho) && c(all(is.na(x$vc.fix$rho)) || all(x$vc.fix$rho))) 
        stop("Model does not contain any (estimated) 'rho' components.")
    if (!missing(sigma2) && (length(sigma2) > 1L)) 
        stop("Can only specify one 'sigma2' component.")
    if (!missing(tau2) && (length(tau2) > 1L)) 
        stop("Can only specify one 'tau2' component.")
    if (!missing(rho) && (length(rho) > 1L)) 
        stop("Can only specify one 'rho' component.")
    if (!missing(sigma2) && is.logical(sigma2)) 
        stop("Must specify the number for the 'sigma2' component.")
    if (!missing(tau2) && is.logical(tau2)) 
        stop("Must specify the number for the 'tau2' component.")
    if (!missing(rho) && is.logical(rho)) 
        stop("Must specify the number for the 'rho' component.")
    if (!missing(sigma2) && (sigma2 > length(x$vc.fix$sigma2) || 
        sigma2 <= 0)) 
        stop("No such 'sigma2' component in the model.")
    if (!missing(tau2) && (tau2 > length(x$vc.fix$tau2) || tau2 <= 
        0)) 
        stop("No such 'tau2' component in the model.")
    if (!missing(rho) && (rho > length(x$vc.fix$rho) || rho <= 
        0)) 
        stop("No such 'rho' component in the model.")
    if (!missing(sigma2) && x$vc.fix$sigma2[sigma2]) 
        stop("Specified 'sigma2' component was fixed.")
    if (!missing(tau2) && x$vc.fix$tau2[tau2]) 
        stop("Specified 'tau2' component was fixed.")
    if (!missing(rho) && x$vc.fix$rho[rho]) 
        stop("Specified 'rho' component was fixed.")
    if (!missing(sigma2)) {
        vc <- x$sigma2[sigma2]
        comp <- "sigma2"
    }
    if (!missing(tau2)) {
        vc <- x$tau2[tau2]
        comp <- "tau2"
    }
    if (!missing(rho)) {
        vc <- x$rho[rho]
        comp <- "rho"
    }
    if (missing(xlim)) {
        if (comp != "rho") {
            vc.lb <- max(0, vc/4)
            vc.ub <- max(0.1, vc * 4)
        }
        else {
            vc.lb <- max(-0.99999, vc - 0.5)
            vc.ub <- min(+0.99999, vc + 0.5)
        }
        if (is.na(vc.lb) || is.na(vc.ub)) 
            stop("Cannot set 'xlim' automatically. Please set this argument manually.")
        xlim <- c(vc.lb, vc.ub)
    }
    else {
        if (length(xlim) != 2L) 
            stop("Argument 'xlim' should be a vector of length 2.")
        xlim <- sort(xlim)
    }
    vcs <- seq(xlim[1], xlim[2], length = steps)
    if (length(vcs) <= 1) 
        stop("Cannot set 'xlim' automatically. Please set this argument manually.")
    lls <- rep(NA, length(vcs))
    b <- matrix(NA, nrow = length(vcs), ncol = x$p)
    ci.lb <- matrix(NA, nrow = length(vcs), ncol = x$p)
    ci.ub <- matrix(NA, nrow = length(vcs), ncol = x$p)
    x.control <- x$control
    if (length(x.control) > 0) {
        con.pos.sigma2.init <- pmatch(names(x.control), "sigma2.init")
        con.pos.tau2.init <- pmatch(names(x.control), "tau2.init")
        con.pos.rho.init <- pmatch(names(x.control), "rho.init")
    }
    else {
        con.pos.sigma2.init <- NA
        con.pos.tau2.init <- NA
        con.pos.rho.init <- NA
    }
    if (progbar) 
        pbar <- txtProgressBar(min = 0, max = steps, style = 3)
    for (i in 1:length(vcs)) {
        sigma2.arg <- ifelse(x$vc.fix$sigma2, x$sigma2, NA)
        tau2.arg <- ifelse(x$vc.fix$tau2, x$tau2, NA)
        rho.arg <- ifelse(x$vc.fix$rho, x$rho, NA)
        if (comp == "sigma2") 
            sigma2.arg[sigma2] <- vcs[i]
        if (comp == "tau2") 
            tau2.arg[tau2] <- vcs[i]
        if (comp == "rho") 
            rho.arg[rho] <- vcs[i]
        if (startmethod == "prev" && i > 1 && !inherits(res, 
            "try-error")) {
            if (is.na(con.pos.sigma2.init)) {
                x.control$sigma2.init <- res$sigma2
            }
            else {
                x.control[[con.pos.sigma2.init]] <- res$sigma2
            }
            if (is.na(con.pos.tau2.init)) {
                x.control$tau2.init <- res$tau2
            }
            else {
                x.control[[con.pos.tau2.init]] <- res$tau2
            }
            if (is.na(con.pos.rho.init)) {
                x.control$rho.init <- res$rho
            }
            else {
                x.control[[con.pos.rho.init]] <- res$rho
            }
            res <- try(rma.mv(x$yi, x$V, x$W, mods = x$X, random = x$random, 
                struct = x$struct, intercept = FALSE, method = x$method, 
                tdist = x$knha, level = x$level, R = x$R, Rscale = x$Rscale, 
                data = x$mf.r, sigma2 = sigma2.arg, tau2 = tau2.arg, 
                rho = rho.arg, control = x.control), silent = TRUE)
        }
        else {
            res <- try(rma.mv(x$yi, x$V, x$W, mods = x$X, random = x$random, 
                struct = x$struct, intercept = FALSE, method = x$method, 
                tdist = x$knha, level = x$level, R = x$R, Rscale = x$Rscale, 
                data = x$mf.r, sigma2 = sigma2.arg, tau2 = tau2.arg, 
                rho = rho.arg, control = x$control), silent = TRUE)
        }
        if (inherits(res, "try-error")) 
            next
        lls[i] <- c(logLik(res))
        b[i, ] <- c(res$b)
        ci.lb[i, ] <- c(res$ci.lb)
        ci.ub[i, ] <- c(res$ci.ub)
        if (progbar) 
            setTxtProgressBar(pbar, i)
    }
    if (progbar) 
        close(pbar)
    b <- data.frame(b)
    ci.lb <- data.frame(ci.lb)
    ci.ub <- data.frame(ci.ub)
    names(b) <- rownames(x$b)
    names(ci.lb) <- rownames(x$b)
    names(ci.ub) <- rownames(x$b)
    res <- list(vc = vcs, ll = lls, b = b, ci.lb = ci.lb, ci.ub = ci.ub)
    names(res)[1] <- switch(comp, sigma2 = "sigma2", tau2 = "tau2", 
        rho = "rho")
    class(res) <- c("profile.rma.mv")
    if (plot) {
        if (missing(ylim)) {
            ylim <- range(lls, na.rm = TRUE)
            ylim[1] <- ylim[1] - 0.1
            ylim[2] <- ylim[2] + 0.1
        }
        if (comp == "sigma2") 
            xlab <- bquote(sigma[.(sigma2)]^2 ~ "Value")
        if (comp == "tau2") 
            xlab <- bquote(tau[.(tau2)]^2 ~ "Value")
        if (comp == "rho") 
            xlab <- bquote(rho[.(rho)] ~ "Value")
        plot(vcs, lls, type = "o", xlab = xlab, ylab = paste(ifelse(x$method == 
            "REML", "Restricted", ""), " Log-Likelihood", sep = ""), 
            bty = "l", pch = pch, ylim = ylim, ...)
        abline(v = vc, lty = "dotted")
        abline(h = logLik(x), lty = "dotted")
        if (comp == "sigma2") 
            title(bquote("Profile Plot for" ~ sigma[.(sigma2)]^2))
        if (comp == "tau2") 
            title(bquote("Profile Plot for" ~ tau[.(tau2)]^2))
        if (comp == "rho") 
            title(bquote("Profile Plot for" ~ rho[.(rho)]))
    }
    invisible(res)
}
profile.rma.uni <-
function (fitted, xlim, ylim, steps = 20, progbar = TRUE, plot = TRUE, 
    pch = 19, ...) 
{
    if (!is.element("rma.uni", class(fitted))) 
        stop("Argument 'fitted' must be an object of class \"rma.uni\".")
    if (steps < 2) 
        stop("Argument 'steps' must be >= 2.")
    x <- fitted
    if (missing(xlim)) {
        tau2.ci <- suppressWarnings(try(confint(x), silent = TRUE))
        if (inherits(tau2.ci, "try-error")) {
            tau2.lb <- NA
            tau2.ub <- NA
        }
        else {
            tau2.lb <- min(x$tau2, tau2.ci$random[1, 2])
            tau2.ub <- max(x$tau2, tau2.ci$random[1, 3])
        }
        if (is.na(tau2.lb) || is.na(tau2.ub)) {
            tau2.lb <- max(0, x$tau2 - 1.96 * x$se.tau2)
            tau2.ub <- max(0, x$tau2 + 1.96 * x$se.tau2)
        }
        if (is.na(tau2.lb) || is.na(tau2.ub)) {
            tau2.lb <- max(0, x$tau2/4)
            tau2.ub <- max(0.1, x$tau2 * 4)
        }
        if (is.na(tau2.lb) || is.na(tau2.ub)) 
            stop("Cannot set 'xlim' automatically. Please set this argument manually.")
        xlim <- c(tau2.lb, tau2.ub)
    }
    else {
        if (length(xlim) != 2L) 
            stop("Argument 'xlim' should be a vector of length 2.")
        xlim <- sort(xlim)
    }
    tau2s <- seq(xlim[1], xlim[2], length = steps)
    if (length(tau2s) <= 1) 
        stop("Cannot set 'xlim' automatically. Please set this argument manually.")
    lls <- rep(NA, length(tau2s))
    b <- matrix(NA, nrow = length(tau2s), ncol = x$p)
    ci.lb <- matrix(NA, nrow = length(tau2s), ncol = x$p)
    ci.ub <- matrix(NA, nrow = length(tau2s), ncol = x$p)
    if (progbar) 
        pbar <- txtProgressBar(min = 0, max = steps, style = 3)
    for (i in 1:length(tau2s)) {
        res <- try(rma(x$yi, x$vi, weights = x$weights, mods = x$X, 
            method = x$method, weighted = x$weighted, intercept = FALSE, 
            knha = x$knha, level = x$level, control = x$control, 
            tau2 = tau2s[i]), silent = TRUE)
        if (inherits(res, "try-error")) 
            next
        lls[i] <- c(logLik(res))
        b[i, ] <- c(res$b)
        ci.lb[i, ] <- c(res$ci.lb)
        ci.ub[i, ] <- c(res$ci.ub)
        if (progbar) 
            setTxtProgressBar(pbar, i)
    }
    if (progbar) 
        close(pbar)
    b <- data.frame(b)
    ci.lb <- data.frame(ci.lb)
    ci.ub <- data.frame(ci.ub)
    names(b) <- rownames(res$b)
    names(ci.lb) <- rownames(res$b)
    names(ci.ub) <- rownames(res$b)
    res <- list(tau2 = tau2s, ll = lls, b = b, ci.lb = ci.lb, 
        ci.ub = ci.ub)
    class(res) <- c("profile.rma.uni")
    if (plot) {
        if (missing(ylim)) {
            ylim <- range(lls, na.rm = TRUE)
            ylim[1] <- ylim[1] - 0.1
            ylim[2] <- ylim[2] + 0.1
        }
        xlab <- expression(paste(tau^2, " Value"))
        plot(tau2s, lls, type = "o", xlab = xlab, ylab = paste(ifelse(x$method == 
            "REML", "Restricted", ""), " Log-Likelihood", sep = ""), 
            bty = "l", pch = pch, ylim = ylim, ...)
        abline(v = x$tau2, lty = "dotted")
        abline(h = logLik(x), lty = "dotted")
        title(expression(paste("Profile Plot for ", tau^2)))
    }
    invisible(res)
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
        for (i in seq_len(x$k)) {
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
            for (i in seq_len(x$k)) {
                text(pos.x[i], pos.y[i], slab[i], pos = ifelse(pos.x[i] >= 
                  0, 2, 4), offset = offset, ...)
            }
        }
    }
    invisible(sav)
}
qqnorm.rma.mv <-
function (y, ...) 
{
    if (!is.element("rma.mv", class(y))) 
        stop("Argument 'y' must be an object of class \"rma.mv\".")
    stop("Method not yet implemented for objects of class \"rma.mv\". Sorry!")
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
        for (i in seq_len(x$k)) {
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
            for (i in seq_len(x$k)) {
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
    na.act <- getOption("na.action")
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
        alpha <- ifelse(level > 1, (100 - level)/100, 1 - level)
        dat <- matrix(rnorm(x$k * reps), nrow = x$k, ncol = reps)
        options(na.action = "na.omit")
        H <- hatvalues(x, type = "matrix")
        options(na.action = na.act)
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
        for (i in seq_len(x$k)) {
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
            for (i in seq_len(x$k)) {
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
        tau2 <- 1/mean(1/x$tau2)
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
    alpha <- ifelse(level > 1, (100 - level)/100, 1 - level)
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
    for (i in seq_len(length(atyis))) {
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
    for (i in seq_len(length(atyis))) {
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
    for (i in seq_len(length(atyis))) {
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
    for (i in seq_len(length(atyis))) {
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
    for (i in seq_len(length(atyis))) {
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
    res <- list(tau = tau, pval = pval, digits = x$digits)
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
    weights <- x$weights
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
        fit <- rma(yi, vi, weights = weights, mods = X, method = x$method, 
            weighted = x$weighted, intercept = FALSE, knha = x$knha, 
            control = x$control, ...)
        zval <- fit$zval[p + 1]
        pval <- fit$pval[p + 1]
        dfs <- ifelse(x$knha || x$robust, fit$k - fit$p, NA)
    }
    else {
        if (ncol(X) >= x$k) 
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
    res <- list(model = model, predictor = predictor, zval = zval, 
        pval = pval, dfs = dfs, method = x$method, digits = x$digits, 
        ret.fit = ret.fit, fit = fit)
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
        "AS", "PHI", "YUQ", "YUY", "RTET", "PBIT", "OR2D", "OR2DN", 
        "OR2DL", "IRR", "IRD", "IRSD", "MD", "SMD", "SMDH", "ROM", 
        "RPB", "RBIS", "D2OR", "D2ORN", "D2ORL", "COR", "UCOR", 
        "ZCOR", "PR", "PLN", "PLO", "PAS", "PFT", "IR", "IRLN", 
        "IRS", "IRFT", "MN", "MC", "SMCC", "SMCR", "ARAW", "AHW", 
        "ABT"))) 
        stop("Unknown 'measure' specified.")
    if (!is.element(method, c("FE", "HS", "HE", "DL", "SJ", "ML", 
        "REML", "EB", "DLIT", "SJIT", "PM"))) 
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
    very.verbose <- ifelse(!is.logical(verbose) && verbose > 
        1, TRUE, FALSE)
    if (very.verbose) 
        message("Extracting/computing yi/vi values ...")
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
    mf.weights <- mf[[match("weights", names(mf))]]
    mf.slab <- mf[[match("slab", names(mf))]]
    mf.subset <- mf[[match("subset", names(mf))]]
    mf.mods <- mf[[match("mods", names(mf))]]
    yi <- eval(mf.yi, data, enclos = sys.frame(sys.parent()))
    weights <- eval(mf.weights, data, enclos = sys.frame(sys.parent()))
    slab <- eval(mf.slab, data, enclos = sys.frame(sys.parent()))
    subset <- eval(mf.subset, data, enclos = sys.frame(sys.parent()))
    mods <- eval(mf.mods, data, enclos = sys.frame(sys.parent()))
    ai <- bi <- ci <- di <- x1i <- x2i <- t1i <- t2i <- NA
    is.formula <- FALSE
    if (!is.null(yi)) {
        if (class(yi) == "formula") {
            options(na.action = "na.pass")
            mods <- model.matrix(yi, data = data)
            attr(mods, "assign") <- NULL
            yi <- model.response(model.frame(yi, data = data))
            options(na.action = na.act)
            names(yi) <- NULL
            intercept <- FALSE
            is.formula <- TRUE
        }
        if (is.matrix(yi)) 
            yi <- as.vector(yi)
        k <- length(yi)
        if (measure == "GEN") {
            if (!is.null(attr(yi, "measure"))) 
                measure <- attr(yi, "measure")
        }
        attr(yi, "measure") <- measure
        mf.vi <- mf[[match("vi", names(mf))]]
        mf.sei <- mf[[match("sei", names(mf))]]
        mf.ni <- mf[[match("ni", names(mf))]]
        vi <- eval(mf.vi, data, enclos = sys.frame(sys.parent()))
        sei <- eval(mf.sei, data, enclos = sys.frame(sys.parent()))
        ni <- eval(mf.ni, data, enclos = sys.frame(sys.parent()))
        if (is.null(vi)) {
            if (is.null(sei)) {
                stop("Need to specify vi or sei argument.")
            }
            else {
                vi <- sei^2
            }
        }
        if (is.matrix(vi)) 
            vi <- as.vector(vi)
        if (length(vi) != k) 
            stop("Length of yi and vi (or sei) vectors are not the same.")
        if (is.null(ni) && !is.null(attr(yi, "ni"))) 
            ni <- attr(yi, "ni")
        if (!is.null(ni) && (length(ni) != k)) 
            ni <- NULL
        if (!is.null(ni)) 
            attr(yi, "ni") <- ni
        if (is.null(slab) & !is.null(attr(yi, "slab"))) 
            slab <- attr(yi, "slab")
        if (!is.null(subset)) {
            yi <- yi[subset]
            vi <- vi[subset]
            ni <- ni[subset]
            attr(yi, "measure") <- measure
            attr(yi, "ni") <- ni
        }
    }
    else {
        if (is.element(measure, c("RR", "OR", "PETO", "RD", "AS", 
            "PHI", "YUQ", "YUY", "RTET", "PBIT", "OR2D", "OR2DN", 
            "OR2DL"))) {
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
            k <- length(ai)
            if (!is.null(subset)) {
                ai <- ai[subset]
                bi <- bi[subset]
                ci <- ci[subset]
                di <- di[subset]
            }
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
            k <- length(x1i)
            if (!is.null(subset)) {
                x1i <- x1i[subset]
                x2i <- x2i[subset]
                t1i <- t1i[subset]
                t2i <- t2i[subset]
            }
            dat <- escalc(measure, x1i = x1i, x2i = x2i, t1i = t1i, 
                t2i = t2i, add = add, to = to, drop00 = drop00, 
                vtype = vtype)
        }
        if (is.element(measure, c("MD", "SMD", "SMDH", "ROM", 
            "RPB", "RBIS", "D2OR", "D2ORN", "D2ORL"))) {
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
            k <- length(m1i)
            if (!is.null(subset)) {
                m1i <- m1i[subset]
                m2i <- m2i[subset]
                sd1i <- sd1i[subset]
                sd2i <- sd2i[subset]
                n1i <- n1i[subset]
                n2i <- n2i[subset]
            }
            dat <- escalc(measure, m1i = m1i, m2i = m2i, sd1i = sd1i, 
                sd2i = sd2i, n1i = n1i, n2i = n2i, vtype = vtype)
        }
        if (is.element(measure, c("COR", "UCOR", "ZCOR"))) {
            mf.ri <- mf[[match("ri", names(mf))]]
            mf.ni <- mf[[match("ni", names(mf))]]
            ri <- eval(mf.ri, data, enclos = sys.frame(sys.parent()))
            ni <- eval(mf.ni, data, enclos = sys.frame(sys.parent()))
            k <- length(ri)
            if (!is.null(subset)) {
                ri <- ri[subset]
                ni <- ni[subset]
            }
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
            k <- length(xi)
            if (!is.null(subset)) {
                xi <- xi[subset]
                mi <- mi[subset]
            }
            dat <- escalc(measure, xi = xi, mi = mi, add = add, 
                to = to, vtype = vtype)
        }
        if (is.element(measure, c("IR", "IRLN", "IRS", "IRFT"))) {
            mf.xi <- mf[[match("xi", names(mf))]]
            mf.ti <- mf[[match("ti", names(mf))]]
            xi <- eval(mf.xi, data, enclos = sys.frame(sys.parent()))
            ti <- eval(mf.ti, data, enclos = sys.frame(sys.parent()))
            k <- length(xi)
            if (!is.null(subset)) {
                xi <- xi[subset]
                ti <- ti[subset]
            }
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
            k <- length(mi)
            if (!is.null(subset)) {
                mi <- mi[subset]
                sdi <- sdi[subset]
                ni <- ni[subset]
            }
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
            k <- length(m1i)
            if (!is.null(subset)) {
                m1i <- m1i[subset]
                m2i <- m2i[subset]
                sd1i <- sd1i[subset]
                sd2i <- sd2i[subset]
                ni <- ni[subset]
                ri <- ri[subset]
            }
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
            k <- length(ai)
            if (!is.null(subset)) {
                ai <- ai[subset]
                mi <- mi[subset]
                ni <- ni[subset]
            }
            dat <- escalc(measure, ai = ai, mi = mi, ni = ni, 
                vtype = vtype)
        }
        if (is.element(measure, "GEN")) 
            stop("Specify the desired outcome measure via the 'measure' argument.")
        yi <- dat$yi
        vi <- dat$vi
        ni <- attr(yi, "ni")
    }
    if (length(weights) == 1) 
        weights <- rep(weights, k)
    if (!is.null(weights) && (length(weights) != k)) 
        stop("Length of yi and weights vectors are not the same.")
    if (!is.null(subset)) 
        weights <- weights[subset]
    if (very.verbose) 
        message("Creating model matrix ...")
    if (class(mods) == "formula") {
        options(na.action = "na.pass")
        mods <- model.matrix(mods, data = data)
        attr(mods, "assign") <- NULL
        options(na.action = na.act)
        intercept <- FALSE
        is.formula <- TRUE
    }
    if (is.vector(mods)) 
        mods <- cbind(mods)
    if (is.data.frame(mods)) 
        mods <- as.matrix(mods)
    if (is.character(mods)) 
        stop("Model matrix contains character variables.")
    if (!is.null(mods) && (nrow(mods) != k)) 
        stop("Number of rows of the model matrix does not match length of yi argument.")
    ids <- seq_len(k)
    if (very.verbose) 
        message("Generating/extracting study labels ...")
    if (is.null(slab)) {
        slab.null <- TRUE
        slab <- ids
    }
    else {
        if (any(is.na(slab))) 
            stop("NAs in study labels.")
        if (any(duplicated(slab))) 
            slab <- make.unique(slab)
        if (length(slab) != k) 
            stop("Study labels not of same length as data.")
        slab.null <- FALSE
    }
    if (!is.null(subset)) {
        if (very.verbose) 
            message("Subsetting ...")
        mods <- mods[subset, , drop = FALSE]
        slab <- slab[subset]
        ids <- ids[subset]
    }
    k <- length(yi)
    if (any(vi <= 0, na.rm = TRUE)) {
        allvipos <- FALSE
        warning("There are outcomes with non-positive sampling variances.")
        vi.neg <- vi < 0
        if (any(vi.neg, na.rm = TRUE)) {
            vi[vi.neg] <- 0
            warning("Negative sampling variances constrained to zero.")
        }
    }
    else {
        allvipos <- TRUE
    }
    if (any(weights < 0, na.rm = TRUE)) 
        stop("Negative weights not allowed.")
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
    weights.f <- weights
    ni.f <- ni
    mods.f <- mods
    k.f <- k
    YVMW.na <- is.na(cbind(yi, vi, mods, weights))
    if (any(YVMW.na)) {
        if (very.verbose) 
            message("Handling NAs ...")
        not.na <- apply(YVMW.na, MARGIN = 1, sum) == 0L
        if (na.act == "na.omit" || na.act == "na.exclude" || 
            na.act == "na.pass") {
            yi <- yi[not.na]
            vi <- vi[not.na]
            weights <- weights[not.na]
            ni <- ni[not.na]
            mods <- mods[not.na, , drop = FALSE]
            k <- length(yi)
            warning("Studies with NAs omitted from model fitting.")
            attr(yi, "measure") <- measure
            attr(yi, "ni") <- ni
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
    if (is.null(mods) && !intercept) {
        warning("Must either include an intercept and/or moderators in model.\n  Coerced intercept into the model.")
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
                btt <- seq.int(from = 2, to = p)
            }
            else {
                btt <- seq_len(p)
            }
        }
        else {
            btt <- 1
        }
    }
    else {
        btt <- btt[(btt >= 1) & (btt <= p)]
        btt <- unique(round(btt))
        if (length(btt) == 0L) 
            stop("Non-existent coefficients specified with 'btt'.")
    }
    bntt <- setdiff(seq_len(p), btt)
    m <- length(btt)
    con <- list(tau2.init = NULL, tau2.min = 0, tau2.max = 50, 
        threshold = 10^-5, maxiter = 100, stepadj = 1, REMLf = TRUE, 
        verbose = FALSE, tol = 1e-07)
    con.pos <- pmatch(names(control), names(con))
    con[c(na.omit(con.pos))] <- control[!is.na(con.pos)]
    if (verbose) 
        con$verbose <- verbose
    verbose <- con$verbose
    if (any(eigen(crossprod(X), symmetric = TRUE, only.values = TRUE)$values <= 
        con$tol)) 
        stop("Model matrix not of full rank. Cannot fit model.")
    iter <- 0
    se.tau2 <- I2 <- H2 <- QE <- QEp <- NA
    s2w <- 1
    alpha <- ifelse(level > 1, (100 - level)/100, 1 - level)
    robust <- FALSE
    Y <- as.matrix(yi)
    if (is.numeric(tau2)) {
        tau2.fix <- TRUE
        tau2.val <- tau2
    }
    else {
        tau2.fix <- FALSE
        tau2.val <- NA
    }
    if (very.verbose && !tau2.fix) 
        message("Estimating tau^2 value ...")
    if (method == "HS") {
        if (!allvipos) 
            stop("HS estimator cannot be used with non-positive sampling variances.")
        wi <- 1/vi
        W <- diag(wi, nrow = k, ncol = k)
        stXWX <- .invcalc(X = X, W = W, k = k)
        P <- W - W %*% X %*% stXWX %*% crossprod(X, W)
        RSS <- crossprod(Y, P) %*% Y
        tau2 <- ifelse(tau2.fix, tau2.val, (RSS - k)/sum(wi))
    }
    if (is.element(method, c("HE", "ML", "REML", "EB"))) {
        stXX <- .invcalc(X = X, W = diag(k), k = k)
        P <- diag(k) - X %*% tcrossprod(stXX, X)
        RSS <- crossprod(Y, P) %*% Y
        V <- diag(vi, nrow = k, ncol = k)
        PV <- P %*% V
        trPV <- .tr(PV)
        tau2 <- ifelse(tau2.fix, tau2.val, (RSS - trPV)/(k - 
            p))
    }
    if (method == "DL") {
        if (!allvipos) 
            stop("DL estimator cannot be used with non-positive sampling variances.")
        wi <- 1/vi
        W <- diag(wi, nrow = k, ncol = k)
        stXWX <- .invcalc(X = X, W = W, k = k)
        P <- W - W %*% X %*% stXWX %*% crossprod(X, W)
        RSS <- crossprod(Y, P) %*% Y
        trP <- .tr(P)
        tau2 <- ifelse(tau2.fix, tau2.val, (RSS - (k - p))/trP)
    }
    if (method == "SJ") {
        if (is.null(con$tau2.init)) {
            tau2.0 <- c(var(yi) * (k - 1)/k)
        }
        else {
            tau2.0 <- con$tau2.init
        }
        wi <- 1/(vi + tau2.0)
        W <- diag(wi, nrow = k, ncol = k)
        stXWX <- .invcalc(X = X, W = W, k = k)
        P <- W - W %*% X %*% stXWX %*% crossprod(X, W)
        RSS <- crossprod(Y, P) %*% Y
        V <- diag(vi, nrow = k, ncol = k)
        PV <- P %*% V
        tau2 <- ifelse(tau2.fix, tau2.val, tau2.0 * RSS/(k - 
            p))
    }
    if (method == "DLIT") {
        conv <- 1
        change <- con$threshold + 1
        if (is.null(con$tau2.init)) {
            tau2 <- 0
        }
        else {
            tau2 <- con$tau2.init
        }
        while (change > con$threshold) {
            if (verbose) 
                cat("Iteration", iter, "\ttau^2 =", round(tau2, 
                  digits), "\n")
            iter <- iter + 1
            tau2.old <- tau2
            wi <- 1/(vi + tau2)
            W <- diag(wi, nrow = k, ncol = k)
            stXWX <- .invcalc(X = X, W = W, k = k)
            P <- W - W %*% X %*% stXWX %*% crossprod(X, W)
            RSS <- crossprod(Y, P) %*% Y
            trP <- .tr(P)
            tau2 <- ifelse(tau2.fix, tau2.val, (RSS - (k - p))/trP)
            tau2[tau2 < con$tau2.min] <- con$tau2.min
            change <- abs(tau2.old - tau2)
            if (iter > con$maxiter) {
                conv <- 0
                break
            }
        }
        if (conv == 0L) 
            stop("Algorithm did not converge.")
    }
    if (method == "SJIT") {
        conv <- 1
        change <- con$threshold + 1
        if (is.null(con$tau2.init)) {
            tau2 <- var(yi) * (k - 1)/k
            tau2.0 <- tau2
        }
        else {
            tau2 <- con$tau2.init
            tau2.0 <- tau2
        }
        while (change > con$threshold) {
            if (verbose) 
                cat("Iteration", iter, "\ttau^2 =", round(tau2, 
                  digits), "\n")
            iter <- iter + 1
            tau2.old <- tau2
            wi <- 1/(vi + tau2)
            W <- diag(wi, nrow = k, ncol = k)
            stXWX <- .invcalc(X = X, W = W, k = k)
            P <- W - W %*% X %*% stXWX %*% crossprod(X, W)
            RSS <- crossprod(Y, P) %*% Y
            V <- diag(vi, nrow = k, ncol = k)
            PV <- P %*% V
            tau2 <- ifelse(tau2.fix, tau2.val, tau2 * RSS/(k - 
                p))
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
        if (.QE.func(con$tau2.min, Y = Y, vi = vi, X = X, k = k, 
            objective = k - p) < 0) {
            tau2 <- con$tau2.min
        }
        else {
            tau2 <- ifelse(tau2.fix, tau2.val, try(uniroot(.QE.func, 
                interval = c(con$tau2.min, con$tau2.max), tol = con$threshold, 
                Y = Y, vi = vi, X = X, k = k, objective = k - 
                  p)$root, silent = TRUE))
            if (!is.numeric(tau2)) 
                stop("Error in iterative search for tau2. Try increasing tau2.max or switch to another 'method'.")
            wi <- 1/(vi + tau2)
            W <- diag(wi, nrow = k, ncol = k)
            stXWX <- .invcalc(X = X, W = W, k = k)
            P <- W - W %*% X %*% stXWX %*% crossprod(X, W)
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
            if (verbose) 
                cat("Iteration", iter, "\ttau^2 =", round(tau2, 
                  digits), "\n")
            iter <- iter + 1
            tau2.old <- tau2
            wi <- 1/(vi + tau2)
            W <- diag(wi, nrow = k, ncol = k)
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
                adj <- (crossprod(Y, P) %*% Y * k/(k - p) - k)/sum(wi)
            }
            adj <- adj * con$stepadj
            while (tau2 + adj < con$tau2.min) {
                adj <- adj/2
            }
            tau2 <- ifelse(tau2.fix, tau2.val, tau2 + adj)
            change <- abs(tau2.old - tau2)
            if (iter > con$maxiter) {
                conv <- 0
                break
            }
        }
        if (conv == 0L) 
            stop("Fisher scoring algorithm did not converge. Try increasing the number of iterations,\n  adjusting the threshold, or use a different estimator for tau^2.")
    }
    tau2 <- max(con$tau2.min, c(tau2))
    if (verbose && is.element(method, c("ML", "REML", "EB"))) 
        cat("Fisher scoring algorithm converged after", iter, 
            "iterations.\n")
    if (method == "HS") {
        se.tau2 <- sqrt(1/sum(wi)^2 * (2 * (k - p) + 4 * max(tau2, 
            0) * .tr(P) + 2 * max(tau2, 0)^2 * sum(P * P)))
    }
    if (method == "HE") {
        se.tau2 <- sqrt(1/(k - p)^2 * (2 * sum(PV * t(PV)) + 
            4 * max(tau2, 0) * trPV + 2 * max(tau2, 0)^2 * (k - 
            p)))
    }
    if (method == "DL" || method == "DLIT") {
        se.tau2 <- sqrt(1/trP^2 * (2 * (k - p) + 4 * max(tau2, 
            0) * trP + 2 * max(tau2, 0)^2 * sum(P * P)))
    }
    if (method == "SJ") {
        se.tau2 <- sqrt(tau2.0^2/(k - p)^2 * (2 * sum(PV * t(PV)) + 
            4 * max(tau2, 0) * sum(PV * P) + 2 * max(tau2, 0)^2 * 
            sum(P * P)))
    }
    if (method == "ML") {
        se.tau2 <- sqrt(2/sum(wi^2))
    }
    if (method == "REML") {
        se.tau2 <- sqrt(2/sum(P * P))
    }
    if (method == "EB" || method == "PM" || method == "SJIT") {
        V <- diag(vi, nrow = k, ncol = k)
        PV <- P %*% V
        se.tau2 <- sqrt((k/(k - p))^2/sum(wi)^2 * (2 * sum(PV * 
            t(PV)) + 4 * max(tau2, 0) * sum(PV * P) + 2 * max(tau2, 
            0)^2 * sum(P * P)))
    }
    if (method == "FE") {
        tau2 <- 0
        if (!allvipos && weighted) 
            stop("Weighted estimation cannot be used with a fixed-effects\n  model when there are non-positive sampling variances.")
    }
    if (very.verbose) 
        message("Model fitting ...")
    wi <- 1/(vi + tau2)
    W <- diag(wi, nrow = k, ncol = k)
    if (weighted) {
        if (is.null(weights)) {
            stXWX <- .invcalc(X = X, W = W, k = k)
            b <- stXWX %*% crossprod(X, W) %*% Y
            vb <- stXWX
            RSS.f <- sum(wi * (yi - X %*% b)^2)
        }
        else {
            A <- diag(weights, nrow = k, ncol = k)
            stXAX <- .invcalc(X = X, W = A, k = k)
            b <- stXAX %*% crossprod(X, A) %*% Y
            vb <- stXAX %*% t(X) %*% A %*% diag(vi + tau2, nrow = k, 
                ncol = k) %*% A %*% X %*% stXAX
            RSS.f <- sum(wi * (yi - X %*% b)^2)
        }
        if (robust) {
            ei <- c(yi - X %*% b)
            vb <- vb %*% t(X) %*% W %*% diag(ei^2, nrow = k, 
                ncol = k) %*% W %*% X %*% vb
            vb <- vb * k/(k - p)
        }
        if (knha) {
            if (RSS.f <= .Machine$double.eps) {
                s2w <- 1
            }
            else {
                s2w <- RSS.f/(k - p)
            }
            vb <- s2w * vb
            if (method == "FE") 
                warning("Knapp & Hartung (2003) method is not meant to be used in the context of fixed-effects models.")
        }
        QM <- c(t(b)[btt] %*% chol2inv(chol(vb[btt, btt])) %*% 
            b[btt])
    }
    else {
        stXX <- .invcalc(X = X, W = diag(k), k = k)
        b <- stXX %*% crossprod(X, Y)
        vb <- tcrossprod(stXX, X) %*% diag(vi + tau2, nrow = k, 
            ncol = k) %*% X %*% stXX
        RSS.f <- sum(wi * (yi - X %*% b)^2)
        if (robust) {
            ei <- c(Y - X %*% b)
            vb <- stXX %*% t(X) %*% diag(ei^2, nrow = k, ncol = k) %*% 
                X %*% stXX
            vb <- vb * k/(k - p)
        }
        if (knha) {
            stXWX <- .invcalc(X = X, W = W, k = k)
            b.knha <- stXWX %*% crossprod(X, W) %*% Y
            RSS.knha <- sum(wi * (yi - X %*% b.knha)^2)
            if (RSS.knha <= .Machine$double.eps) {
                s2w <- 1
            }
            else {
                s2w <- RSS.knha/(k - p)
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
    if (very.verbose) 
        message("Heterogeneity testing ...")
    if (allvipos) {
        wi <- 1/vi
        W.FE <- diag(wi, nrow = k, ncol = k)
        stXWX <- .invcalc(X = X, W = W.FE, k = k)
        P <- W.FE - W.FE %*% X %*% stXWX %*% crossprod(X, W.FE)
        QE <- max(0, c(crossprod(Y, P) %*% Y))
        QEp <- ifelse(k - p >= 1, pchisq(QE, df = k - p, lower.tail = FALSE), 
            1)
        vi.avg <- (k - p)/.tr(P)
        I2 <- 100 * tau2/(vi.avg + tau2)
        H2 <- tau2/vi.avg + 1
    }
    if (!int.only && int.incl && method != "FE") {
        if (very.verbose) {
            message("Fitting RE model for R^2 computation ...")
            res.RE <- try(rma(yi, vi, weights = weights, method = method, 
                weighted = weighted, knha = knha, verbose = ifelse(verbose, 
                  TRUE, FALSE), control = con, digits = digits), 
                silent = FALSE)
        }
        else {
            res.RE <- suppressWarnings(try(rma(yi, vi, weights = weights, 
                method = method, weighted = weighted, knha = knha, 
                verbose = ifelse(verbose, TRUE, FALSE), control = con, 
                digits = digits), silent = FALSE))
        }
        if (!inherits(res.RE, "try-error")) {
            tau2.RE <- res.RE$tau2
            if (identical(tau2.RE, 0)) {
                R2 <- NA
            }
            else {
                R2 <- round(max(0, 100 * (tau2.RE - tau2)/tau2.RE), 
                  2)
            }
        }
        else {
            R2 <- NA
        }
    }
    else {
        R2 <- NULL
    }
    if (very.verbose) 
        message("Computing fit statistics and log likelihood ...")
    parms <- p + ifelse(method == "FE" || tau2.fix, 0, 1)
    ll.ML <- -1/2 * (k) * log(2 * base::pi) - 1/2 * sum(log(vi + 
        tau2)) - 1/2 * RSS.f
    ll.REML <- -1/2 * (k - p) * log(2 * base::pi) + ifelse(con$REMLf, 
        1/2 * determinant(crossprod(X), logarithm = TRUE)$modulus, 
        0) - 1/2 * sum(log(vi + tau2)) - 1/2 * determinant(crossprod(X, 
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
    AICc.ML <- -2 * ll.ML + 2 * parms * max(k, parms + 2)/(max(k, 
        parms + 2) - parms - 1)
    dev.REML <- -2 * (ll.REML - 0)
    AIC.REML <- -2 * ll.REML + 2 * parms
    BIC.REML <- -2 * ll.REML + parms * log(k - p)
    AICc.REML <- -2 * ll.REML + 2 * parms * max(k - p, parms + 
        2)/(max(k - p, parms + 2) - parms - 1)
    fit.stats <- matrix(c(ll.ML, dev.ML, AIC.ML, BIC.ML, AICc.ML, 
        ll.REML, dev.REML, AIC.REML, BIC.REML, AICc.REML), ncol = 2, 
        byrow = FALSE)
    dimnames(fit.stats) <- list(c("ll", "dev", "AIC", "BIC", 
        "AICc"), c("ML", "REML"))
    fit.stats <- data.frame(fit.stats)
    if (very.verbose) 
        message("Preparing output ...")
    p.eff <- p
    k.eff <- k
    res <- list(b = b, se = se, zval = zval, pval = pval, ci.lb = ci.lb, 
        ci.ub = ci.ub, vb = vb, tau2 = tau2, se.tau2 = se.tau2, 
        k = k, k.f = k.f, k.eff = k.eff, p = p, p.eff = p.eff, 
        parms = parms, m = m, QE = QE, QEp = QEp, QM = QM, QMp = QMp, 
        I2 = I2, H2 = H2, R2 = R2, int.only = int.only, int.incl = int.incl, 
        allvipos = allvipos, yi = yi, vi = vi, X = X, weights = weights, 
        yi.f = yi.f, vi.f = vi.f, X.f = X.f, weights.f = weights.f, 
        ai.f = ai.f, bi.f = bi.f, ci.f = ci.f, di.f = di.f, x1i.f = x1i.f, 
        x2i.f = x2i.f, t1i.f = t1i.f, t2i.f = t2i.f, ni = ni, 
        ni.f = ni.f, ids = ids, not.na = not.na, slab = slab, 
        slab.null = slab.null, measure = measure, method = method, 
        weighted = weighted, knha = knha, robust = robust, s2w = s2w, 
        btt = btt, intercept = intercept, digits = digits, level = level, 
        control = control, verbose = verbose, add = add, to = to, 
        drop00 = drop00, fit.stats = fit.stats)
    class(res) <- c("rma.uni", "rma")
    return(res)
}
rma.glmm <-
function (ai, bi, ci, di, n1i, n2i, x1i, x2i, t1i, t2i, xi, mi, 
    ti, ni, mods, measure, intercept = TRUE, data, slab, subset, 
    add = 1/2, to = "only0", drop00 = TRUE, vtype = "LS", model = "UM.FS", 
    method = "ML", tdist = FALSE, level = 95, digits = 4, btt, 
    nAGQ = 7, verbose = FALSE, control) 
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
    if (model == "CM.AL" && measure == "IRR") 
        model <- "CM.EL"
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    if (missing(btt)) 
        btt <- NULL
    if (missing(control)) 
        control <- list()
    if (is.element(measure, c("OR", "IRR")) && model == "UM.RS" && 
        method == "ML" && nAGQ > 1) {
        warning("Currently not possible to fit RE/ME model='UM.RS' with nAGQ > 1. nAGQ automatically set to 1.")
        nAGQ <- 1
    }
    knha <- tdist
    very.verbose <- ifelse(!is.logical(verbose) && verbose > 
        1, TRUE, FALSE)
    if (very.verbose) 
        message("Extracting data and computing yi/vi values ...")
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
        k <- length(ai)
        if (!is.null(subset)) {
            ai <- ai[subset]
            bi <- bi[subset]
            ci <- ci[subset]
            di <- di[subset]
        }
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
        k <- length(x1i)
        if (!is.null(subset)) {
            x1i <- x1i[subset]
            x2i <- x2i[subset]
            t1i <- t1i[subset]
            t2i <- t2i[subset]
        }
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
        k <- length(xi)
        if (!is.null(subset)) {
            xi <- xi[subset]
            mi <- mi[subset]
        }
        dat <- escalc(measure, xi = xi, mi = mi, add = add, to = to, 
            vtype = vtype)
    }
    if (is.element(measure, c("IRLN"))) {
        mf.xi <- mf[[match("xi", names(mf))]]
        mf.ti <- mf[[match("ti", names(mf))]]
        xi <- eval(mf.xi, data, enclos = sys.frame(sys.parent()))
        ti <- eval(mf.ti, data, enclos = sys.frame(sys.parent()))
        k <- length(xi)
        if (!is.null(subset)) {
            xi <- xi[subset]
            ti <- ti[subset]
        }
        dat <- escalc(measure, xi = xi, ti = ti, add = add, to = to, 
            vtype = vtype)
    }
    yi <- dat$yi
    vi <- dat$vi
    ni <- attr(yi, "ni")
    ids <- seq_len(k)
    if (very.verbose) 
        message("Creating model matrix ...")
    is.formula <- class(mods) == "formula"
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
    if (is.character(mods)) 
        stop("Model matrix contains character variables.")
    if (!is.null(mods) && (nrow(mods) != k)) 
        stop("Number of rows of the model matrix does not match length of yi argument.")
    if (very.verbose) 
        message("Generating/extracting study labels ...")
    if (is.null(slab)) {
        if (!is.null(attr(yi, "slab"))) {
            slab.null <- FALSE
            slab <- attr(yi, "slab")
        }
        else {
            slab.null <- TRUE
            slab <- ids
        }
    }
    else {
        if (any(is.na(slab))) 
            stop("NAs in study labels.")
        if (any(duplicated(slab))) 
            slab <- make.unique(slab)
        if (length(slab) != k) 
            stop("Study labels not of same length as data.")
        slab.null <- FALSE
    }
    if (!is.null(subset)) {
        if (very.verbose) 
            message("Subsetting ...")
        mods <- mods[subset, , drop = FALSE]
        slab <- slab[subset]
        ids <- ids[subset]
    }
    k <- length(yi)
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
            if (very.verbose) 
                message("Handling NAs in table data ...")
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
            if (very.verbose) 
                message("Handling NAs in table data ...")
            not.na <- apply(x1ix2it1it2imods.na, MARGIN = 1, 
                sum) == 0L
            if (na.act == "na.omit" || na.act == "na.exclude" || 
                na.act == "na.pass") {
                x1i <- x1i[not.na]
                x2i <- x2i[not.na]
                t1i <- t1i[not.na]
                t2i <- t2i[not.na]
                mods <- mods[not.na, , drop = FALSE]
                k <- length(x1i)
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
            if (very.verbose) 
                message("Handling NAs in table data ...")
            not.na <- apply(ximimods.na, MARGIN = 1, sum) == 
                0L
            if (na.act == "na.omit" || na.act == "na.exclude" || 
                na.act == "na.pass") {
                xi <- xi[not.na]
                mi <- mi[not.na]
                mods <- mods[not.na, , drop = FALSE]
                k <- length(xi)
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
            if (very.verbose) 
                message("Handling NAs in table data ...")
            not.na <- apply(xitimods.na, MARGIN = 1, sum) == 
                0L
            if (na.act == "na.omit" || na.act == "na.exclude" || 
                na.act == "na.pass") {
                xi <- xi[not.na]
                ti <- ti[not.na]
                mods <- mods[not.na, , drop = FALSE]
                k <- length(xi)
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
    yivi.na <- is.na(cbind(yi, vi, mods.f))
    if (any(yivi.na)) {
        if (very.verbose) 
            message("Handling NAs in yi/vi ...")
        not.na.yivi <- apply(yivi.na, MARGIN = 1, sum) == 0L
        if (na.act == "na.omit" || na.act == "na.exclude" || 
            na.act == "na.pass") {
            yi <- yi[not.na.yivi]
            ni <- ni[not.na.yivi]
            vi <- vi[not.na.yivi]
            warning("Some yi/vi values are NA.")
            attr(yi, "measure") <- measure
            attr(yi, "ni") <- ni
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
        warning("Must either include an intercept and/or moderators in model.\n  Coerced intercept into the model.")
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
    if (method == "FE" && p > k) 
        stop("Number of parameters to be estimated is larger than the number of observations.")
    if (method != "FE" && (p + 1) > k) 
        stop("Number of parameters to be estimated is larger than the number of observations.")
    if (is.null(btt)) {
        if (p > 1) {
            if (int.incl) {
                btt <- seq.int(from = 2, to = p)
            }
            else {
                btt <- seq_len(p)
            }
        }
        else {
            btt <- 1
        }
    }
    else {
        btt <- btt[(btt >= 1) & (btt <= p)]
        btt <- unique(round(btt))
        if (length(btt) == 0L) 
            stop("Non-existent coefficients specified with 'btt'.")
    }
    m <- length(btt)
    con <- list(verbose = FALSE, optimizer = "optim", optmethod = "BFGS", 
        scale = TRUE, tol = 1e-07, dnchgcalc = "dFNCHypergeo", 
        dnchgprec = 1e-10)
    con.pos <- pmatch(names(control), names(con))
    con[c(na.omit(con.pos))] <- control[!is.na(con.pos)]
    if (verbose) 
        con$verbose <- verbose
    verbose <- con$verbose
    silent <- !verbose
    pos.optCtrl <- pmatch(names(control), "optCtrl", nomatch = 0)
    if (sum(pos.optCtrl) > 0) {
        optCtrl <- control[[which(pos.optCtrl == 1)]]
    }
    else {
        optCtrl <- list()
    }
    if (con$optimizer == "optim") {
        con.pos <- pmatch(names(optCtrl), "REPORT", nomatch = 0)
        if (sum(con.pos) > 0) {
            optCtrl[which(con.pos == 1)] <- 1
            names(optCtrl)[which(con.pos == 1)] <- "REPORT"
        }
        else {
            optCtrl$REPORT <- 1
        }
        optCtrl$trace <- con$verbose
    }
    if (con$optimizer == "nlminb") 
        optCtrl$trace <- ifelse(con$verbose > 0, 1, 0)
    if (is.element(con$optimizer, c("uobyqa", "newuoa", "bobyqa"))) 
        optCtrl$iprint <- ifelse(con$verbose > 0, 3, 0)
    pos.clogitCtrl <- pmatch(names(control), "clogitCtrl", nomatch = 0)
    if (sum(pos.clogitCtrl) > 0) {
        clogitCtrl <- control[[which(pos.clogitCtrl == 1)]]
    }
    else {
        clogitCtrl <- list()
    }
    pos.clogisticCtrl <- pmatch(names(control), "clogisticCtrl", 
        nomatch = 0)
    if (sum(pos.clogisticCtrl) > 0) {
        clogisticCtrl <- control[[which(pos.clogisticCtrl == 
            1)]]
    }
    else {
        clogisticCtrl <- list()
    }
    pos.glmCtrl <- pmatch(names(control), "glmCtrl", nomatch = 0)
    if (sum(pos.glmCtrl) > 0) {
        glmCtrl <- control[[which(pos.glmCtrl == 1)]]
    }
    else {
        glmCtrl <- list()
    }
    glmCtrl$trace <- ifelse(con$verbose > 0, TRUE, FALSE)
    pos.glmerCtrl <- pmatch(names(control), "glmerCtrl", nomatch = 0)
    if (sum(pos.glmerCtrl) > 0) {
        glmerCtrl <- control[[which(pos.glmerCtrl == 1)]]
    }
    else {
        glmerCtrl <- list()
    }
    pos.intCtrl <- pmatch(names(control), "intCtrl", nomatch = 0)
    if (sum(pos.intCtrl) > 0) {
        intCtrl <- control[[which(pos.intCtrl == 1)]]
    }
    else {
        intCtrl <- list()
    }
    con.pos <- pmatch(names(intCtrl), "lower", nomatch = 0)
    if (sum(con.pos) > 0) {
        intCtrl[which(con.pos == 1)] <- -Inf
        names(intCtrl)[which(con.pos == 1)] <- "lower"
    }
    else {
        intCtrl$lower <- -Inf
    }
    con.pos <- pmatch(names(intCtrl), "upper", nomatch = 0)
    if (sum(con.pos) > 0) {
        intCtrl[which(con.pos == 1)] <- Inf
        names(intCtrl)[which(con.pos == 1)] <- "upper"
    }
    else {
        intCtrl$upper <- Inf
    }
    con.pos <- pmatch(names(intCtrl), "subdivisions", nomatch = 0)
    if (sum(con.pos) > 0) {
        intCtrl[which(con.pos == 1)] <- 100L
        names(intCtrl)[which(con.pos == 1)] <- "subdivisions"
    }
    else {
        intCtrl$subdivisions <- 100L
    }
    con.pos <- pmatch(names(intCtrl), "rel.tol", nomatch = 0)
    if (sum(con.pos) > 0) {
        intCtrl[which(con.pos == 1)] <- .Machine$double.eps^0.25
        names(intCtrl)[which(con.pos == 1)] <- "rel.tol"
    }
    else {
        intCtrl$rel.tol <- .Machine$double.eps^0.25
    }
    if (!is.element(con$optimizer, c("optim", "nlminb", "uobyqa", 
        "newuoa", "bobyqa", "clogit", "clogistic"))) 
        stop("Unknown optimizer specified.")
    if (con$dnchgcalc != "dnoncenhypergeom" && con$dnchgcalc != 
        "dFNCHypergeo") 
        stop("Unknown dnchgcalc method specified.")
    if (is.element(con$optimizer, c("clogit", "clogistic")) && 
        method == "ML") 
        stop("Cannot use 'clogit' or 'clogistic' with method='ML'.")
    if (is.element(measure, c("OR", "IRR"))) {
        if ((model == "UM.FS" && method == "ML") || (model == 
            "UM.RS") || (model == "CM.AL" && method == "ML") || 
            (model == "CM.EL" && method == "ML")) {
            if (!require(lme4)) 
                stop("Please install the 'lme4' package to fit this model.")
        }
    }
    if (is.element(measure, c("PLO", "IRLN")) && method == "ML") {
        if (!require(lme4)) 
            stop("Please install the 'lme4' package to fit this model.")
    }
    if (is.element(measure, c("OR", "IRR"))) {
        if ((model == "UM.FS" && method == "ML") || (model == 
            "UM.RS") || (model == "CM.AL" && method == "ML") || 
            (model == "CM.EL" && method == "ML")) {
            if (!require(lme4)) 
                stop("Please install the 'lme4' package to fit this model.")
        }
    }
    if (is.element(measure, c("PLO", "IRLN")) && method == "ML") {
        if (!require(lme4)) 
            stop("Please install the 'lme4' package to fit this model.")
    }
    if (measure == "OR" && model == "CM.EL") {
        if (is.element(con$optimizer, c("uobyqa", "newuoa", "bobyqa"))) {
            if (!require(minqa)) 
                stop("Please install the 'minqa' package to fit this model.")
            minqa <- get(con$optimizer)
            con$optimizer <- "minqa"
        }
        if (con$optimizer == "optim" || con$optimizer == "nlminb" || 
            con$optimizer == "minqa") {
            if (!require(numDeriv)) 
                stop("Please install the 'numDeriv' package to fit this model.")
            if (con$dnchgcalc == "dFNCHypergeo") {
                if (!require(BiasedUrn)) 
                  stop("Please install the 'BiasedUrn' package to fit this model.")
            }
        }
        if (con$optimizer == "clogit") {
            if (!require(survival)) 
                stop("Please install the 'survival' package to fit this model.")
        }
        if (con$optimizer == "clogistic") {
            if (!require(Epi)) 
                stop("Please install the 'Epi' package to fit this model.")
        }
    }
    if (any(eigen(crossprod(X), symmetric = TRUE, only.values = TRUE)$values <= 
        con$tol)) 
        stop("Model matrix not of full rank. Cannot fit model.")
    se.tau2 <- I2 <- H2 <- QE <- QEp <- NA
    alpha <- ifelse(level > 1, (100 - level)/100, 1 - level)
    robust <- FALSE
    if (!int.only && int.incl && con$scale) {
        Xsave <- X
        meanX <- apply(X[, 2:p, drop = FALSE], 2, mean)
        sdX <- apply(X[, 2:p, drop = FALSE], 2, sd)
        is.d <- apply(X, 2, function(x) all(sapply(x, identical, 
            0) | sapply(x, identical, 1)))
        X[, !is.d] <- apply(X[, !is.d, drop = FALSE], 2, scale)
    }
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
            study <- factor(rep(seq_len(k), each = 2))
            const <- cbind(rep(1, 2 * k))
            X.fit <- X[rep(seq(k), each = 2), , drop = FALSE]
            X.fit <- cbind(group1 * X.fit[, , drop = FALSE])
            row.names(X.fit) <- seq_len(2 * k)
            if (model == "UM.FS") {
                if (verbose) 
                  message("Fitting FE model ...")
                if (k > 1) {
                  res.FE <- try(glm(dat.grp ~ -1 + X.fit + study, 
                    offset = dat.off, family = dat.fam, control = glmCtrl), 
                    silent = silent)
                }
                else {
                  res.FE <- try(glm(dat.grp ~ -1 + X.fit + const, 
                    offset = dat.off, family = dat.fam, control = glmCtrl), 
                    silent = silent)
                }
                if (inherits(res.FE, "try-error")) 
                  stop("Cannot fit FE model.")
                if (verbose) 
                  message("Fitting saturated model ...")
                if (k > 1) {
                  X.QE <- model.matrix(~-1 + X.fit + study + 
                    study:group1)
                  res.QE <- try(glm(dat.grp ~ -1 + X.QE, offset = dat.off, 
                    family = dat.fam, control = glmCtrl), silent = silent)
                }
                else {
                  res.QE <- res.FE
                }
                if (inherits(res.QE, "try-error")) 
                  stop("Cannot fit saturated model.")
                ll.FE <- c(logLik(res.FE))
                ll.QE <- c(logLik(res.QE))
                b2.QE <- cbind(na.omit(coef(res.QE)[-seq_len(k + 
                  p)]))
                vb2.QE <- vcov(res.QE)[-seq_len(k + p), -seq_len(k + 
                  p), drop = FALSE]
                if (method == "ML") {
                  if (verbose) 
                    message("Fitting ML model ...")
                  res.ML <- try(glmer(dat.grp ~ -1 + X.fit + 
                    study + (group12 - 1 | study), family = dat.fam, 
                    nAGQ = nAGQ, verbose = verbose, control = do.call(glmerControl, 
                      glmerCtrl)), silent = silent)
                  if (inherits(res.ML, "try-error")) 
                    stop("Cannot fit ML model.")
                  ll.ML <- ll.QE - 1/2 * deviance(res.ML)
                }
                if (method == "FE") {
                  b <- cbind(coef(res.FE)[seq_len(p)])
                  vb <- vcov(res.FE)[seq_len(p), seq_len(p), 
                    drop = FALSE]
                  tau2 <- 0
                  sigma2 <- NA
                  parms <- p + k
                  p.eff <- p + k
                  k.eff <- 2 * k
                }
                if (method == "ML") {
                  b <- cbind(fixef(res.ML)[seq_len(p)])
                  vb <- as.matrix(vcov(res.ML))[seq_len(p), seq_len(p), 
                    drop = FALSE]
                  tau2 <- VarCorr(res.ML)[[1]][1]
                  sigma2 <- NA
                  parms <- p + k + 1
                  p.eff <- p + k
                  k.eff <- 2 * k
                }
            }
            if (model == "UM.RS") {
                if (verbose) 
                  message("Fitting FE model ...")
                res.FE <- try(glmer(dat.grp ~ -1 + X.fit + const + 
                  (1 | study), offset = dat.off, family = dat.fam, 
                  nAGQ = nAGQ, verbose = verbose, control = do.call(glmerControl, 
                    glmerCtrl)), silent = silent)
                if (inherits(res.FE, "try-error")) 
                  stop("Cannot fit FE model.")
                if (verbose) 
                  message("Fitting saturated model ...")
                if (k > 1) {
                  X.QE <- model.matrix(~-1 + X.fit + const + 
                    study:group1)
                  res.QE <- try(glm(dat.grp ~ -1 + X.QE, offset = dat.off, 
                    family = dat.fam, control = glmCtrl), silent = silent)
                  X.QE <- X.QE[, !is.na(coef(res.QE)), drop = FALSE]
                  res.QE <- try(glmer(dat.grp ~ -1 + X.QE + (1 | 
                    study), offset = dat.off, family = dat.fam, 
                    start = c(sqrt(VarCorr(res.FE)[[1]][1])), 
                    nAGQ = nAGQ, verbose = verbose, control = do.call(glmerControl, 
                      glmerCtrl)), silent = silent)
                }
                else {
                  res.QE <- res.FE
                }
                if (inherits(res.QE, "try-error")) 
                  stop("Cannot fit saturated model.")
                ll.FE <- c(logLik(res.FE))
                ll.QE <- c(logLik(res.QE))
                b2.QE <- cbind(fixef(res.QE)[-seq_len(p + 1)])
                vb2.QE <- as.matrix(vcov(res.QE))[-seq_len(p + 
                  1), -seq_len(p + 1), drop = FALSE]
                if (method == "ML") {
                  if (verbose) 
                    message("Fitting ML model ...")
                  res.ML <- try(glmer(dat.grp ~ -1 + X.fit + 
                    const + (1 | study) + (group12 - 1 | study), 
                    offset = dat.off, family = dat.fam, nAGQ = nAGQ, 
                    verbose = verbose, control = do.call(glmerControl, 
                      glmerCtrl)), silent = silent)
                  if (inherits(res.ML, "try-error")) 
                    stop("Cannot fit ML model.")
                  ll.ML <- c(logLik(res.ML))
                }
                if (method == "FE") {
                  b <- cbind(fixef(res.FE)[seq_len(p)])
                  vb <- as.matrix(vcov(res.FE))[seq_len(p), seq_len(p), 
                    drop = FALSE]
                  tau2 <- 0
                  sigma2 <- VarCorr(res.FE)[[1]][1]
                  parms <- p + 1 + 1
                  p.eff <- p + 1
                  k.eff <- 2 * k
                }
                if (method == "ML") {
                  b <- cbind(fixef(res.ML)[seq_len(p)])
                  vb <- as.matrix(vcov(res.ML))[seq_len(p), seq_len(p), 
                    drop = FALSE]
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
            study <- factor(seq_len(k))
            X.fit <- X
            if (verbose) 
                message("Fitting FE model ...")
            res.FE <- try(glm(dat.grp ~ -1 + X.fit, offset = dat.off, 
                family = binomial, control = glmCtrl), silent = silent)
            if (inherits(res.FE, "try-error")) 
                stop("Cannot fit FE model.")
            if (verbose) 
                message("Fitting saturated model ...")
            if (k > 1) {
                X.QE <- model.matrix(~-1 + X.fit + study)
                res.QE <- try(glm(dat.grp ~ -1 + X.QE, offset = dat.off, 
                  family = binomial, control = glmCtrl), silent = silent)
            }
            else {
                res.QE <- res.FE
            }
            if (inherits(res.QE, "try-error")) 
                stop("Cannot fit saturated model.")
            ll.FE <- c(logLik(res.FE))
            ll.QE <- c(logLik(res.QE))
            b2.QE <- cbind(na.omit(coef(res.QE)[-seq_len(p)]))
            vb2.QE <- vcov(res.QE)[-seq_len(p), -seq_len(p), 
                drop = FALSE]
            if (method == "ML") {
                if (verbose) 
                  message("Fitting ML model ...")
                if (verbose) {
                  res.ML <- try(glmer(dat.grp ~ -1 + X.fit + 
                    (1 | study), offset = dat.off, family = binomial, 
                    nAGQ = nAGQ, verbose = verbose, control = do.call(glmerControl, 
                      glmerCtrl)), silent = silent)
                }
                else {
                  res.ML <- suppressMessages(try(glmer(dat.grp ~ 
                    -1 + X.fit + (1 | study), offset = dat.off, 
                    family = binomial, nAGQ = nAGQ, verbose = verbose, 
                    control = do.call(glmerControl, glmerCtrl)), 
                    silent = silent))
                }
                if (inherits(res.ML, "try-error")) 
                  stop("Cannot fit ML model.")
                ll.ML <- ll.QE - 1/2 * deviance(res.ML)
            }
            if (method == "FE") {
                b <- cbind(coef(res.FE)[seq_len(p)])
                vb <- vcov(res.FE)[seq_len(p), seq_len(p), drop = FALSE]
                tau2 <- 0
                sigma2 <- NA
                parms <- p
                p.eff <- p
                k.eff <- k
            }
            if (method == "ML") {
                b <- cbind(fixef(res.ML)[seq_len(p)])
                vb <- as.matrix(vcov(res.ML))[seq_len(p), seq_len(p), 
                  drop = FALSE]
                tau2 <- VarCorr(res.ML)[[1]][1]
                sigma2 <- NA
                parms <- p + 1
                p.eff <- p
                k.eff <- k
            }
        }
        if (measure == "OR" && model == "CM.EL") {
            if (verbose) 
                message("Fitting FE model ...")
            if (con$optimizer == "optim" || con$optimizer == 
                "nlminb" || con$optimizer == "minqa") {
                if (con$optimizer == "optim") {
                  res.FE <- try(optim(par = c(coef(res.FE)[seq_len(p)], 
                    0), .dnchg, method = con$optmethod, hessian = TRUE, 
                    ai = ai, bi = bi, ci = ci, di = di, X.fit = X.fit, 
                    random = FALSE, verbose = verbose, digits = digits, 
                    dnchgcalc = con$dnchgcalc, dnchgprec = con$dnchgprec, 
                    control = optCtrl), silent = silent)
                }
                if (con$optimizer == "nlminb") {
                  res.FE <- try(nlminb(start = c(coef(res.FE)[seq_len(p)], 
                    0), .dnchg, ai = ai, bi = bi, ci = ci, di = di, 
                    X.fit = X.fit, random = FALSE, verbose = verbose, 
                    digits = digits, dnchgcalc = con$dnchgcalc, 
                    dnchgprec = con$dnchgprec, control = optCtrl), 
                    silent = silent)
                }
                if (con$optimizer == "minqa") {
                  res.FE <- try(minqa(par = c(coef(res.FE)[seq_len(p)], 
                    0), .dnchg, ai = ai, bi = bi, ci = ci, di = di, 
                    X.fit = X.fit, random = FALSE, verbose = verbose, 
                    digits = digits, dnchgcalc = con$dnchgcalc, 
                    dnchgprec = con$dnchgprec, control = optCtrl), 
                    silent = silent)
                }
                if (con$optimizer == "optim" || con$optimizer == 
                  "nlminb") {
                  if (inherits(res.FE, "try-error") || res.FE$convergence != 
                    0) 
                    stop("Cannot fit FE model.")
                }
                if (con$optimizer == "minqa") {
                  if (inherits(res.FE, "try-error") || res.FE$ierr != 
                    0) 
                    stop("Cannot fit FE model.")
                }
                if (very.verbose) 
                  message("Computing Hessian ...")
                h.FE <- hessian(.dnchg, x = res.FE$par, ai = ai, 
                  bi = bi, ci = ci, di = di, X.fit = X.fit, random = FALSE, 
                  verbose = verbose, digits = digits, dnchgcalc = con$dnchgcalc, 
                  dnchgprec = con$dnchgprec)
                if (verbose) 
                  message("Fitting saturated model ...")
                if (k > 1) {
                  is.aliased <- is.na(coef(res.QE))
                  X.QE <- X.QE[, !is.aliased, drop = FALSE]
                  if (con$optimizer == "optim") {
                    res.QE <- try(optim(par = c(coef(res.QE)[!is.aliased], 
                      0), .dnchg, method = con$optmethod, hessian = TRUE, 
                      ai = ai, bi = bi, ci = ci, di = di, X.fit = X.QE, 
                      random = FALSE, verbose = verbose, digits = digits, 
                      dnchgcalc = con$dnchgcalc, dnchgprec = con$dnchgprec, 
                      control = optCtrl), silent = silent)
                  }
                  if (con$optimizer == "nlminb") {
                    res.QE <- try(nlminb(start = c(coef(res.QE)[!is.aliased], 
                      0), .dnchg, ai = ai, bi = bi, ci = ci, 
                      di = di, X.fit = X.QE, random = FALSE, 
                      verbose = verbose, digits = digits, dnchgcalc = con$dnchgcalc, 
                      dnchgprec = con$dnchgprec, control = optCtrl), 
                      silent = silent)
                  }
                  if (con$optimizer == "minqa") {
                    res.QE <- try(minqa(par = c(coef(res.QE)[!is.aliased], 
                      0), .dnchg, ai = ai, bi = bi, ci = ci, 
                      di = di, X.fit = X.QE, random = FALSE, 
                      verbose = verbose, digits = digits, dnchgcalc = con$dnchgcalc, 
                      dnchgprec = con$dnchgprec, control = optCtrl), 
                      silent = silent)
                  }
                  if (con$optimizer == "optim" || con$optimizer == 
                    "nlminb") {
                    if (inherits(res.QE, "try-error") || res.QE$convergence != 
                      0) 
                      stop("Cannot fit saturated model.")
                  }
                  if (con$optimizer == "minqa") {
                    if (inherits(res.QE, "try-error") || res.QE$ierr != 
                      0) 
                      stop("Cannot fit saturated model.")
                  }
                  if (very.verbose) 
                    message("Computing Hessian ...")
                  h.QE <- hessian(.dnchg, x = res.QE$par, ai = ai, 
                    bi = bi, ci = ci, di = di, X.fit = X.QE, 
                    random = FALSE, verbose = verbose, digits = digits, 
                    dnchgcalc = con$dnchgcalc, dnchgprec = con$dnchgprec)
                }
                else {
                  res.QE <- res.FE
                  h.QE <- h.FE
                }
                if (con$optimizer == "optim") {
                  ll.FE <- -1 * res.FE$value
                  ll.QE <- -1 * res.QE$value
                }
                if (con$optimizer == "nlminb") {
                  ll.FE <- -1 * res.FE$objective
                  ll.QE <- -1 * res.QE$objective
                }
                if (con$optimizer == "minqa") {
                  ll.FE <- -1 * res.FE$fval
                  ll.QE <- -1 * res.QE$fval
                }
                b2.QE <- res.QE$par
                hessian <- h.QE
                p.QE <- length(b2.QE)
                b2.QE <- b2.QE[-p.QE]
                hessian <- hessian[-p.QE, -p.QE, drop = FALSE]
                p.QE <- length(b2.QE)
                is.0 <- apply(hessian == 0L, 2, sum) == p.QE
                b2.QE <- b2.QE[!is.0]
                hessian <- hessian[!is.0, !is.0, drop = FALSE]
                b2.QE <- cbind(b2.QE[-seq_len(p)])
                h.A <- hessian[seq_len(p), seq_len(p), drop = FALSE]
                h.B <- hessian[seq_len(p), -seq_len(p), drop = FALSE]
                h.C <- hessian[-seq_len(p), seq_len(p), drop = FALSE]
                h.D <- hessian[-seq_len(p), -seq_len(p), drop = FALSE]
                chol.h.A <- try(chol(h.A), silent = silent)
                if (class(chol.h.A) == "try-error") 
                  stop("Cannot invert Hessian for saturated model.")
                Ivb2.QE <- h.D - h.C %*% chol2inv(chol.h.A) %*% 
                  h.B
                QE.Wld <- c(t(b2.QE) %*% Ivb2.QE %*% b2.QE)
            }
            if (con$optimizer == "clogit" || con$optimizer == 
                "clogistic") {
                event <- unlist(lapply(seq_len(k), function(i) c(rep.int(1, 
                  ai[i]), rep.int(0, bi[i]), rep.int(1, ci[i]), 
                  rep.int(0, di[i]))))
                group1 <- unlist(lapply(seq_len(k), function(i) c(rep.int(1, 
                  ai[i]), rep.int(1, bi[i]), rep.int(0, ci[i]), 
                  rep.int(0, di[i]))))
                study.l <- factor(rep(seq_len(k), times = ni))
                X.fit.l <- X[rep(seq_len(k), times = ni), , drop = FALSE]
                X.fit.l <- cbind(group1 * X.fit.l)
                const <- rep(1, length(event))
                if (k > 1) {
                  if (con$optimizer == "clogit") {
                    args.clogit <- clogitCtrl
                    args.clogit$formula <- event ~ X.fit.l + 
                      strata(study.l)
                    res.FE <- try(do.call(clogit, args.clogit), 
                      silent = silent)
                  }
                  if (con$optimizer == "clogistic") {
                    args.clogistic <- clogisticCtrl
                    args.clogistic$formula <- event ~ X.fit.l
                    args.clogistic$strata <- study.l
                    res.FE <- try(do.call(clogistic, args.clogistic), 
                      silent = silent)
                  }
                }
                else {
                  stop(paste("Cannot use '", con$optimizer, "' optimizer when k=1.", 
                    sep = ""))
                }
                if (inherits(res.FE, "try-error")) 
                  stop("Cannot fit FE model.")
                if (verbose) 
                  message("Fitting saturated model ...")
                X.QE.l <- model.matrix(~-1 + X.fit.l + study.l:group1)
                X.QE.l <- X.QE.l[, !is.na(coef(res.QE)), drop = FALSE]
                X.QE <- X.QE[, !is.na(coef(res.QE)), drop = FALSE]
                if (con$optimizer == "clogit") {
                  args.clogit <- clogitCtrl
                  args.clogit$formula <- event ~ X.QE.l + strata(study.l)
                  if (verbose) {
                    res.QE <- try(do.call(clogit, args.clogit), 
                      silent = silent)
                  }
                  else {
                    res.QE <- suppressWarnings(try(do.call(clogit, 
                      args.clogit), silent = silent))
                  }
                }
                if (con$optimizer == "clogistic") {
                  args.clogistic <- clogisticCtrl
                  args.clogistic$formula <- event ~ X.QE.l
                  args.clogistic$strata <- study.l
                  res.QE <- try(do.call(clogistic, args.clogistic), 
                    silent = silent)
                }
                if (inherits(res.QE, "try-error")) 
                  stop("Cannot fit saturated model.")
                ll.FE <- -1 * .dnchg(c(cbind(coef(res.FE)), 0), 
                  ai = ai, bi = bi, ci = ci, di = di, X.fit = X.fit, 
                  random = FALSE, verbose = verbose, digits = digits, 
                  dnchgcalc = con$dnchgcalc, dnchgprec = con$dnchgprec)
                ll.QE <- -1 * .dnchg(c(cbind(coef(res.QE)), 0), 
                  ai = ai, bi = bi, ci = ci, di = di, X.fit = X.QE, 
                  random = FALSE, verbose = verbose, digits = digits, 
                  dnchgcalc = con$dnchgcalc, dnchgprec = con$dnchgprec)
                b2.QE <- cbind(coef(res.QE)[-seq_len(p)])
                vb2.QE <- vcov(res.QE)[-seq_len(p), -seq_len(p), 
                  drop = FALSE]
            }
            if (method == "ML") {
                if (verbose) 
                  message("Fitting ML model ...")
                if (con$optimizer == "optim") {
                  res.ML <- try(optim(par = c(b, log(tau2 + 0.001)), 
                    .dnchg, method = con$optmethod, hessian = TRUE, 
                    ai = ai, bi = bi, ci = ci, di = di, X.fit = X.fit, 
                    random = TRUE, verbose = verbose, digits = digits, 
                    dnchgcalc = con$dnchgcalc, dnchgprec = con$dnchgprec, 
                    intCtrl = intCtrl, control = optCtrl), silent = silent)
                }
                if (con$optimizer == "nlminb") {
                  res.ML <- try(nlminb(start = c(b, log(tau2 + 
                    0.001)), .dnchg, ai = ai, bi = bi, ci = ci, 
                    di = di, X.fit = X.fit, random = TRUE, verbose = verbose, 
                    digits = digits, dnchgcalc = con$dnchgcalc, 
                    dnchgprec = con$dnchgprec, intCtrl = intCtrl, 
                    control = optCtrl), silent = silent)
                }
                if (con$optimizer == "minqa") {
                  res.ML <- try(minqa(par = c(b, log(tau2 + 0.001)), 
                    .dnchg, ai = ai, bi = bi, ci = ci, di = di, 
                    X.fit = X.fit, random = TRUE, verbose = verbose, 
                    digits = digits, dnchgcalc = con$dnchgcalc, 
                    dnchgprec = con$dnchgprec, intCtrl = intCtrl, 
                    control = optCtrl), silent = silent)
                }
                if (con$optimizer == "optim" || con$optimizer == 
                  "nlminb") {
                  if (inherits(res.ML, "try-error") || res.ML$convergence != 
                    0) 
                    stop("Cannot fit ML model.")
                }
                if (con$optimizer == "minqa") {
                  if (inherits(res.ML, "try-error") || res.ML$ierr != 
                    0) 
                    stop("Cannot fit ML model.")
                }
                if (very.verbose) 
                  message("Computing Hessian ...")
                h.ML <- hessian(.dnchg, x = res.ML$par, method.args = list(r = 8), 
                  ai = ai, bi = bi, ci = ci, di = di, X.fit = X.fit, 
                  random = TRUE, verbose = verbose, digits = digits, 
                  dnchgcalc = con$dnchgcalc, dnchgprec = con$dnchgprec, 
                  intCtrl = intCtrl)
                if (con$optimizer == "optim") {
                  ll.ML <- -1 * res.ML$value
                }
                if (con$optimizer == "nlminb") {
                  ll.ML <- -1 * res.ML$objective
                }
                if (con$optimizer == "minqa") {
                  ll.ML <- -1 * res.ML$fval
                }
            }
            if (method == "FE") {
                if (con$optimizer == "optim" || con$optimizer == 
                  "nlminb" || con$optimizer == "minqa") {
                  b <- cbind(res.FE$par[seq_len(p)])
                  chol.h <- try(chol(h.FE[seq_len(p), seq_len(p)]), 
                    silent = silent)
                  if (class(chol.h) == "try-error") 
                    stop("Cannot invert Hessian for FE model.")
                  vb <- chol2inv(chol.h)
                }
                if (con$optimizer == "clogit" || con$optimizer == 
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
                b <- cbind(res.ML$par[seq_len(p)])
                chol.h <- try(chol(h.ML), silent = silent)
                if (class(chol.h) == "try-error") 
                  stop("Cannot invert Hessian for ML model.")
                vb.f <- chol2inv(chol.h)
                vb <- vb.f[seq_len(p), seq_len(p), drop = FALSE]
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
        study <- factor(seq_len(k))
        X.fit <- X
        if (verbose) 
            message("Fitting FE model ...")
        res.FE <- try(glm(dat.grp ~ -1 + X.fit, offset = dat.off, 
            family = dat.fam, control = glmCtrl), silent = silent)
        if (inherits(res.FE, "try-error")) 
            stop("Cannot fit FE model.")
        if (verbose) 
            message("Fitting saturated model ...")
        if (k > 1) {
            X.QE <- model.matrix(~-1 + X.fit + study)
            if (verbose) {
                res.QE <- try(glm(dat.grp ~ -1 + X.QE, offset = dat.off, 
                  family = dat.fam, control = glmCtrl), silent = silent)
            }
            else {
                res.QE <- suppressWarnings(try(glm(dat.grp ~ 
                  -1 + X.QE, offset = dat.off, family = dat.fam, 
                  control = glmCtrl), silent = silent))
            }
        }
        else {
            res.QE <- res.FE
        }
        if (inherits(res.QE, "try-error")) 
            stop("Cannot fit saturated model.")
        ll.FE <- c(logLik(res.FE))
        ll.QE <- c(logLik(res.QE))
        b2.QE <- cbind(na.omit(coef(res.QE)[-seq_len(p)]))
        vb2.QE <- vcov(res.QE)[-seq_len(p), -seq_len(p), drop = FALSE]
        if (method == "ML") {
            if (verbose) 
                message("Fitting ML model ...")
            if (verbose) {
                res.ML <- try(glmer(dat.grp ~ -1 + X.fit + (1 | 
                  study), offset = dat.off, family = dat.fam, 
                  nAGQ = nAGQ, verbose = verbose, control = do.call(glmerControl, 
                    glmerCtrl)), silent = silent)
            }
            else {
                res.ML <- suppressMessages(try(glmer(dat.grp ~ 
                  -1 + X.fit + (1 | study), offset = dat.off, 
                  family = dat.fam, nAGQ = nAGQ, verbose = verbose, 
                  control = do.call(glmerControl, glmerCtrl)), 
                  silent = silent))
            }
            if (inherits(res.ML, "try-error")) 
                stop("Cannot fit ML model.")
            ll.ML <- ll.QE - 1/2 * deviance(res.ML)
        }
        if (method == "FE") {
            b <- cbind(coef(res.FE)[seq_len(p)])
            vb <- vcov(res.FE)[seq_len(p), seq_len(p), drop = FALSE]
            tau2 <- 0
            sigma2 <- NA
            parms <- p
            p.eff <- p
            k.eff <- k
        }
        if (method == "ML") {
            b <- cbind(fixef(res.ML)[seq_len(p)])
            vb <- as.matrix(vcov(res.ML))[seq_len(p), seq_len(p), 
                drop = FALSE]
            tau2 <- VarCorr(res.ML)[[1]][1]
            sigma2 <- NA
            parms <- p + 1
            p.eff <- p
            k.eff <- k
        }
    }
    rownames(vb) <- colnames(vb) <- rownames(b) <- colnames(X)
    if (very.verbose) 
        message("Heterogeneity testing ...")
    if (measure != "OR" || model != "CM.EL" || con$optimizer != 
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
    W <- diag(wi, nrow = k.I2, ncol = k.I2)
    stXWX <- .invcalc(X = X.I2, W = W, k = k.I2)
    P <- W - W %*% X.I2 %*% stXWX %*% crossprod(X.I2, W)
    vi.avg <- (k.I2 - p)/.tr(P)
    I2 <- 100 * tau2/(vi.avg + tau2)
    H2 <- tau2/vi.avg + 1
    chol.h <- try(chol(vb[btt, btt]), silent = silent)
    if (class(chol.h) == "try-error") 
        stop("Cannot invert Hessian for QM test.")
    QM <- c(t(b)[btt] %*% chol2inv(chol.h) %*% b[btt])
    if (!int.only && int.incl && con$scale) {
        mX <- rbind(c(1, -1 * ifelse(is.d[-1], 0, meanX/sdX)), 
            cbind(0, diag(ifelse(is.d[-1], 1, 1/sdX), nrow = length(is.d) - 
                1, ncol = length(is.d) - 1)))
        b <- mX %*% b
        vb <- mX %*% vb %*% t(mX)
        X <- Xsave
        rownames(vb) <- colnames(vb) <- rownames(b) <- colnames(X)
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
    if (very.verbose) 
        message("Computing fit statistics and log likelihood ...")
    ll.ML <- ifelse(method == "FE", ll.FE, ll.ML)
    ll.REML <- NA
    dev.ML <- -2 * (ll.ML - ll.QE)
    AIC.ML <- -2 * ll.ML + 2 * parms
    BIC.ML <- -2 * ll.ML + parms * log(k.eff)
    AICc.ML <- -2 * ll.ML + 2 * parms * max(k.eff, parms + 2)/(max(k.eff, 
        parms + 2) - parms - 1)
    dev.REML <- NA
    AIC.REML <- NA
    BIC.REML <- NA
    AICc.REML <- NA
    fit.stats <- matrix(c(ll.ML, dev.ML, AIC.ML, AICc.ML, BIC.ML, 
        ll.REML, dev.REML, AIC.REML, BIC.REML, AICc.REML), ncol = 2, 
        byrow = FALSE)
    dimnames(fit.stats) <- list(c("ll", "dev", "AIC", "BIC", 
        "AICc"), c("ML", "REML"))
    fit.stats <- data.frame(fit.stats)
    if (very.verbose) 
        message("Preparing output ...")
    weighted <- TRUE
    robust <- FALSE
    res <- list(b = b, se = se, zval = zval, pval = pval, ci.lb = ci.lb, 
        ci.ub = ci.ub, vb = vb, tau2 = tau2, se.tau2 = se.tau2, 
        sigma2 = sigma2, k = k, k.f = k.f, k.yi = k.yi, k.eff = k.eff, 
        p = p, p.eff = p.eff, parms = parms, m = m, QE.Wld = QE.Wld, 
        QEp.Wld = QEp.Wld, QE.LRT = QE.LRT, QEp.LRT = QEp.LRT, 
        QE.df = QE.df, QM = QM, QMp = QMp, I2 = I2, H2 = H2, 
        int.only = int.only, int.incl = int.incl, yi = yi, vi = vi, 
        X = X, yi.f = yi.f, vi.f = vi.f, X.f = X.f, ai = ai, 
        bi = bi, ci = ci, di = di, ai.f = ai.f, bi.f = bi.f, 
        ci.f = ci.f, di.f = di.f, x1i = x1i, x2i = x2i, t1i = t1i, 
        t2i = t2i, x1i.f = x1i.f, x2i.f = x2i.f, t1i.f = t1i.f, 
        t2i.f = t2i.f, xi = xi, mi = mi, ti = ti, xi.f = xi.f, 
        mi.f = mi.f, ti.f = ti.f, ni = ni, ni.f = ni.f, ids = ids, 
        not.na = not.na, not.na.yivi = not.na.yivi, slab = slab, 
        slab.null = slab.null, measure = measure, method = method, 
        model = model, weighted = weighted, knha = knha, robust = robust, 
        btt = btt, intercept = intercept, digits = digits, level = level, 
        control = control, verbose = verbose, add = add, to = to, 
        drop00 = drop00, fit.stats = fit.stats)
    class(res) <- c("rma.glmm", "rma")
    return(res)
}
rma.mh <-
function (ai, bi, ci, di, n1i, n2i, x1i, x2i, t1i, t2i, measure = "OR", 
    data, slab, subset, add = 1/2, to = "only0", drop00 = TRUE, 
    correct = TRUE, level = 95, digits = 4, verbose = FALSE) 
{
    if (!is.element(measure, c("OR", "RR", "RD", "IRR", "IRD"))) 
        stop("Mantel-Haenszel method can only be used with measures OR, RR, RD, IRR, and IRD.")
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
    if (verbose) 
        message("Extracting data and computing yi/vi values ...")
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
        ids <- seq_len(k)
        if (verbose) 
            message("Generating/extracting study labels ...")
        if (is.null(slab)) {
            slab.null <- TRUE
            slab <- ids
        }
        else {
            if (any(is.na(slab))) 
                stop("NAs in study labels.")
            if (any(duplicated(slab))) 
                slab <- make.unique(slab)
            if (length(slab) != k) 
                stop("Study labels not of same length as data.")
            slab.null <- FALSE
        }
        if (!is.null(subset)) {
            if (verbose) 
                message("Subsetting ...")
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
            if (verbose) 
                message("Handling NAs in table data ...")
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
            if (verbose) 
                message("Handling NAs in yi/vi ...")
            not.na.yivi <- apply(yivi.na, MARGIN = 1, sum) == 
                0L
            if (na.act == "na.omit" || na.act == "na.exclude" || 
                na.act == "na.pass") {
                yi <- yi[not.na.yivi]
                vi <- vi[not.na.yivi]
                ni <- ni[not.na.yivi]
                warning("Some yi/vi values are NA.")
                attr(yi, "measure") <- measure
                attr(yi, "ni") <- ni
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
        n1i <- ai + bi
        n2i <- ci + di
        Ni <- ai + bi + ci + di
    }
    if (is.element(measure, c("IRR", "IRD"))) {
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
        ids <- seq_len(k)
        if (verbose) 
            message("Generating/extracting study labels ...")
        if (is.null(slab)) {
            slab.null <- TRUE
            slab <- ids
        }
        else {
            if (any(is.na(slab))) 
                stop("NAs in study labels.")
            if (any(duplicated(slab))) 
                slab <- make.unique(slab)
            if (length(slab) != k) 
                stop("Study labels not of same length as data.")
            slab.null <- FALSE
        }
        if (!is.null(subset)) {
            if (verbose) 
                message("Subsetting ...")
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
            if (verbose) 
                message("Handling NAs in table data ...")
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
            if (verbose) 
                message("Handling NAs in yi/vi ...")
            not.na.yivi <- apply(yivi.na, MARGIN = 1, sum) == 
                0L
            if (na.act == "na.omit" || na.act == "na.exclude" || 
                na.act == "na.pass") {
                yi <- yi[not.na.yivi]
                vi <- vi[not.na.yivi]
                ni <- ni[not.na.yivi]
                warning("Some yi/vi values are NA.")
                attr(yi, "measure") <- measure
                attr(yi, "ni") <- ni
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
        Ti <- t1i + t2i
    }
    alpha <- ifelse(level > 1, (100 - level)/100, 1 - level)
    CO <- COp <- MH <- MHp <- BD <- BDp <- TA <- TAp <- k.pos <- NA
    if (verbose) 
        message("Model fitting ...")
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
            MH <- NA
            MHp <- NA
        }
        else {
            CO <- sum(ai - (n1i/Ni) * xt)^2/sum((n1i/Ni) * (n2i/Ni) * 
                (xt * (yt/Ni)))
            COp <- pchisq(CO, df = 1, lower.tail = FALSE)
            MH <- (abs(sum(ai - (n1i/Ni) * xt)) - ifelse(correct, 
                0.5, 0))^2/sum((n1i/Ni) * (n2i/Ni) * (xt * (yt/(Ni - 
                1))))
            MHp <- pchisq(MH, df = 1, lower.tail = FALSE)
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
        xt <- x1i + x2i
        if (identical(sum(xt), 0)) {
            MH <- NA
            MHp <- NA
        }
        else {
            MH <- (abs(sum(x1i - xt * (t1i/Ti))) - ifelse(correct, 
                0.5, 0))^2/sum(xt * (t1i/Ti) * (t2i/Ti))
            MHp <- pchisq(MH, df = 1, lower.tail = FALSE)
        }
    }
    if (measure == "IRD") {
        b <- sum((x1i * t2i - x2i * t1i)/Ti)/sum((t1i/Ti) * t2i)
        se <- sqrt(sum(((t1i/Ti) * t2i)^2 * (x1i/t1i^2 + x2i/t2i^2)))/sum((t1i/Ti) * 
            t2i)
        zval <- b/se
        pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)
        ci.lb <- b - qnorm(alpha/2, lower.tail = FALSE) * se
        ci.ub <- b + qnorm(alpha/2, lower.tail = FALSE) * se
        names(b) <- "intrcpt"
        vb <- matrix(se^2, dimnames = list("intrcpt", "intrcpt"))
    }
    if (verbose) 
        message("Heterogeneity testing ...")
    wi <- 1/vi
    QE <- sum(wi * (yi - b)^2)
    if (k.yi - 1 >= 1) {
        QEp <- pchisq(QE, df = k.yi - 1, lower.tail = FALSE)
    }
    else {
        QEp <- 1
    }
    if (verbose) 
        message("Computing fit statistics and log likelihood ...")
    ll.ML <- -1/2 * (k.yi) * log(2 * base::pi) - 1/2 * sum(log(vi)) - 
        1/2 * QE
    ll.REML <- -1/2 * (k.yi - 1) * log(2 * base::pi) + 1/2 * 
        log(k.yi) - 1/2 * sum(log(vi)) - 1/2 * log(sum(wi)) - 
        1/2 * QE
    dev.ML <- -2 * (ll.ML - sum(dnorm(yi, mean = yi, sd = sqrt(vi), 
        log = TRUE)))
    AIC.ML <- -2 * ll.ML + 2
    BIC.ML <- -2 * ll.ML + log(k.yi)
    AICc.ML <- -2 * ll.ML + 2 * max(k.yi, 3)/(max(k.yi, 3) - 
        2)
    dev.REML <- -2 * (ll.REML - 0)
    AIC.REML <- -2 * ll.REML + 2
    BIC.REML <- -2 * ll.REML + log(k.yi - 1)
    AICc.REML <- -2 * ll.REML + 2 * max(k.yi - 1, 3)/(max(k.yi - 
        1, 3) - 2)
    fit.stats <- matrix(c(ll.ML, dev.ML, AIC.ML, BIC.ML, AICc.ML, 
        ll.REML, dev.REML, AIC.REML, BIC.REML, AICc.REML), ncol = 2, 
        byrow = FALSE)
    dimnames(fit.stats) <- list(c("ll", "dev", "AIC", "BIC", 
        "AICc"), c("ML", "REML"))
    fit.stats <- data.frame(fit.stats)
    if (verbose) 
        message("Preparing output ...")
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
    res <- list(b = b, se = se, zval = zval, pval = pval, ci.lb = ci.lb, 
        ci.ub = ci.ub, vb = vb, tau2 = tau2, k = k, k.f = k.f, 
        k.yi = k.yi, k.pos = k.pos, k.eff = k.eff, p = p, parms = parms, 
        QE = QE, QEp = QEp, CO = CO, COp = COp, MH = MH, MHp = MHp, 
        BD = BD, BDp = BDp, TA = TA, TAp = TAp, int.only = int.only, 
        yi = yi, vi = vi, yi.f = yi.f, vi.f = vi.f, X.f = X.f, 
        ai = ai, bi = bi, ci = ci, di = di, ai.f = ai.f, bi.f = bi.f, 
        ci.f = ci.f, di.f = di.f, x1i = x1i, x2i = x2i, t1i = t1i, 
        t2i = t2i, x1i.f = x1i.f, x2i.f = x2i.f, t1i.f = t1i.f, 
        t2i.f = t2i.f, ni = ni, ni.f = ni.f, ids = ids, not.na = not.na, 
        not.na.yivi = not.na.yivi, slab = slab, slab.null = slab.null, 
        measure = measure, method = method, weighted = weighted, 
        knha = knha, robust = robust, intercept = intercept, 
        digits = digits, level = level, add = add, to = to, drop00 = drop00, 
        correct = correct, fit.stats = fit.stats)
    class(res) <- c("rma.mh", "rma")
    return(res)
}
rma.mv <-
function (yi, V, W, mods, random, struct = "CS", intercept = TRUE, 
    data, slab, subset, method = "REML", tdist = FALSE, level = 95, 
    digits = 4, btt, R, Rscale = "cor", sigma2, tau2, rho, sparse = FALSE, 
    verbose = FALSE, control) 
{
    if (!is.element(method, c("FE", "ML", "REML"))) 
        stop("Unknown 'method' specified.")
    if (!is.element(struct, c("CS", "HCS", "UN", "AR", "HAR", 
        "UNHO"))) 
        stop("Unknown 'struct' specified.")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    if (missing(random)) 
        random <- NULL
    if (missing(btt)) 
        btt <- NULL
    if (missing(R)) 
        R <- NULL
    if (missing(sigma2)) 
        sigma2 <- NULL
    if (missing(tau2)) 
        tau2 <- NULL
    if (missing(rho)) 
        rho <- NULL
    if (missing(control)) 
        control <- list()
    knha <- tdist
    if (is.character(Rscale)) 
        Rscale <- match.arg(Rscale, c("none", "cor", "cor0", 
            "cov0"))
    if (is.logical(Rscale)) 
        Rscale <- ifelse(Rscale, "cor", "none")
    if (is.numeric(Rscale)) {
        Rscale <- round(Rscale)
        if (Rscale > 3 | Rscale < 0) 
            stop("Unknown 'Rscale' value specified.")
        Rscale <- switch(as.character(Rscale), `0` = "none", 
            `1` = "cor", `2` = "cor0", `3` = "cov0")
    }
    very.verbose <- ifelse(!is.logical(verbose) && verbose > 
        1, TRUE, FALSE)
    if (very.verbose) 
        message("Extracting yi/V values ...")
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
    mf.V <- mf[[match("V", names(mf))]]
    mf.W <- mf[[match("W", names(mf))]]
    mf.ni <- mf[[match("ni", names(mf))]]
    mf.slab <- mf[[match("slab", names(mf))]]
    mf.subset <- mf[[match("subset", names(mf))]]
    mf.mods <- mf[[match("mods", names(mf))]]
    yi <- eval(mf.yi, data, enclos = sys.frame(sys.parent()))
    V <- eval(mf.V, data, enclos = sys.frame(sys.parent()))
    W <- eval(mf.W, data, enclos = sys.frame(sys.parent()))
    ni <- eval(mf.ni, data, enclos = sys.frame(sys.parent()))
    slab <- eval(mf.slab, data, enclos = sys.frame(sys.parent()))
    subset <- eval(mf.subset, data, enclos = sys.frame(sys.parent()))
    mods <- eval(mf.mods, data, enclos = sys.frame(sys.parent()))
    is.formula <- FALSE
    if (class(yi) == "formula") {
        options(na.action = "na.pass")
        mods <- model.matrix(yi, data = data)
        attr(mods, "assign") <- NULL
        yi <- model.response(model.frame(yi, data = data))
        options(na.action = na.act)
        names(yi) <- NULL
        intercept <- FALSE
        is.formula <- TRUE
    }
    if (is.matrix(yi)) 
        yi <- as.vector(yi)
    k <- length(yi)
    measure <- "GEN"
    if (!is.null(attr(yi, "measure"))) 
        measure <- attr(yi, "measure")
    attr(yi, "measure") <- measure
    if (is.null(V)) 
        stop("Need to specify V argument.")
    if (is.list(V)) {
        rows <- sapply(V, nrow)
        cols <- sapply(V, ncol)
        if (any(rows != cols)) 
            stop("List elements in V must be square matrices.")
        if (sparse) {
            V <- bdiag(V)
        }
        else {
            V <- bldiag(V)
        }
    }
    if (is.vector(V) || nrow(V) == 1L || ncol(V) == 1L) 
        V <- diag(as.vector(V), nrow = k, ncol = k)
    if (is.data.frame(V)) 
        V <- as.matrix(V)
    V <- unname(V)
    if (dim(V)[1] != dim(V)[2]) 
        stop("V must be a square matrix.")
    if (!isSymmetric(V)) 
        stop("V must be a symmetric matrix.")
    if (dim(V)[1] != k) 
        stop("Length of yi and length/dimensions of V are not the same.")
    if (sparse) 
        V <- Matrix(V, sparse = TRUE)
    if (!is.null(W)) {
        if (is.vector(W) || nrow(W) == 1L || ncol(W) == 1L) {
            W <- as.vector(W)
            if (length(W) == 1L) 
                W <- rep(W, k)
            A <- diag(W, nrow = length(W), ncol = length(W))
        }
        else {
            A <- W
        }
        if (is.data.frame(A)) 
            A <- as.matrix(A)
        A <- unname(A)
        if (dim(A)[1] != dim(A)[2]) 
            stop("W must be a square matrix.")
        if (!isSymmetric(A)) 
            stop("W must be a symmetric matrix.")
        if (dim(A)[1] != k) 
            stop("Length of yi and length/dimensions of W are not the same.")
        if (sparse) 
            A <- Matrix(A, sparse = TRUE)
    }
    else {
        A <- NULL
    }
    if (is.null(ni) && !is.null(attr(yi, "ni"))) 
        ni <- attr(yi, "ni")
    if (!is.null(ni) && (length(ni) != k)) 
        ni <- NULL
    if (very.verbose) 
        message("Creating model matrix ...")
    if (class(mods) == "formula") {
        options(na.action = "na.pass")
        mods <- model.matrix(mods, data = data)
        attr(mods, "assign") <- NULL
        options(na.action = na.act)
        intercept <- FALSE
        is.formula <- TRUE
    }
    if (is.vector(mods)) 
        mods <- cbind(mods)
    if (is.data.frame(mods)) 
        mods <- as.matrix(mods)
    if (is.character(mods)) 
        stop("Model matrix contains character variables.")
    if (!is.null(mods) && (nrow(mods) != k)) 
        stop("Number of rows of the model matrix does not match length of yi argument.")
    if (method != "FE" && !is.null(random)) {
        if (very.verbose) 
            message("Processing 'random' argument ...")
        if (!is.list(random)) 
            random <- list(random)
        mf.r <- lapply(random, get_all_vars, data = data)
        mf.r.ncols <- sapply(mf.r, ncol)
        if (any(mf.r.ncols > 2)) 
            stop("No more than two elements allowed in each formula of the 'random' argument.")
        if (sum(mf.r.ncols == 2) > 1) 
            stop("Only one formula with two elements allowed in the 'random' argument.")
        if (any(mf.r.ncols == 2)) {
            if (length(mf.r) > 1) {
                mf.g <- mf.r[which(mf.r.ncols == 2)]
                mf.s <- mf.r[-which(mf.r.ncols == 2)]
            }
            else {
                mf.g <- mf.r
                mf.s <- NULL
            }
        }
        else {
            mf.g <- NULL
            mf.s <- mf.r
        }
        withS <- !is.null(mf.s)
        withG <- !is.null(mf.g)
        mf.r.nrows <- sapply(mf.r, nrow)
        if (any(mf.r.nrows != k)) 
            stop("Length of variables specified via the 'random' argument does not match length of the data.")
    }
    else {
        mf.s <- NULL
        mf.g <- NULL
        mf.r <- NULL
        withS <- FALSE
        withG <- FALSE
    }
    s.names <- sapply(mf.s, names)
    g.names <- sapply(mf.g, names)
    ids <- seq_len(k)
    if (very.verbose) 
        message("Generating/extracting study labels ...")
    if (is.null(slab)) {
        if (!is.null(attr(yi, "slab"))) {
            slab.null <- FALSE
            slab <- attr(yi, "slab")
        }
        else {
            slab.null <- TRUE
            slab <- ids
        }
    }
    else {
        if (any(is.na(slab))) 
            stop("NAs in study labels.")
        if (any(duplicated(slab))) 
            slab <- make.unique(slab)
        if (length(slab) != k) 
            stop("Study labels not of same length as data.")
        slab.null <- FALSE
    }
    if (!is.null(subset)) {
        if (very.verbose) 
            message("Subsetting ...")
        yi <- yi[subset]
        V <- V[subset, subset, drop = FALSE]
        A <- A[subset, subset, drop = FALSE]
        ni <- ni[subset]
        mods <- mods[subset, , drop = FALSE]
        slab <- slab[subset]
        mf.s <- lapply(mf.s, function(x) x[subset, , drop = FALSE])
        mf.g <- lapply(mf.g, function(x) x[subset, , drop = FALSE])
        mf.r <- lapply(mf.r, function(x) x[subset, , drop = FALSE])
        ids <- ids[subset]
        k <- length(yi)
        attr(yi, "measure") <- measure
        attr(yi, "ni") <- ni
    }
    vi <- diag(V)
    if (any(vi <= 0, na.rm = TRUE)) {
        allvipos <- FALSE
        warning("There are outcomes with non-positive sampling variances.")
        vi.neg <- vi < 0
        if (any(vi.neg, na.rm = TRUE)) {
            V[vi.neg, , drop = FALSE] <- 0
            V[, vi.neg, drop = FALSE] <- 0
            vi[vi.neg] <- 0
            warning("Negative sampling variances constrained to zero.")
        }
    }
    else {
        allvipos <- TRUE
    }
    yi.f <- yi
    vi.f <- vi
    V.f <- V
    W.f <- A
    ni.f <- ni
    mods.f <- mods
    mf.g.f <- mf.g
    k.f <- k
    if (withS) {
        if (very.verbose) 
            message(paste0("Processing '", paste(as.character(random[mf.r.ncols == 
                1]), collapse = ", "), "' term(s) ..."))
        mf.s <- lapply(mf.s, function(x) factor(x[[1]]))
        if (any(sapply(lapply(mf.s, is.na), any))) 
            stop("No NAs allowed in variables specified in the 'random' argument.")
        s.nvals <- length(mf.s)
        if (is.null(sigma2)) 
            sigma2 <- rep(NA, s.nvals)
        if (length(sigma2) == 1) 
            sigma2 <- rep(sigma2, s.nvals)
        if (length(sigma2) != s.nvals) 
            stop(paste("Length of 'sigma2' argument (", length(sigma2), 
                ") does not match actual number of variance components (", 
                s.nvals, ").", sep = ""))
        if (any(sigma2 < 0, na.rm = TRUE)) 
            stop("Specified value(s) of 'sigma2' must be non-negative.")
        s.nlevels <- sapply(mf.s, nlevels)
        s.levels <- lapply(mf.s, levels)
        if (is.null(R)) {
            withR <- FALSE
            Rfix <- rep(FALSE, s.nvals)
        }
        else {
            if (very.verbose) 
                message("Processing 'R' argument ...")
            withR <- TRUE
            if (is.data.frame(R) || !is.list(R)) 
                R <- list(R)
            if (is.null(names(R)) || any(nchar(names(R)) == 0)) 
                stop("Argument 'R' must be a *named* list.")
            R <- R[!sapply(R, is.null)]
            R <- lapply(R, as.matrix)
            R <- R[s.names]
            names(R) <- s.names
            Rfix <- !sapply(R, is.null)
            if (any(Rfix)) {
                if (any(sapply(R[Rfix], function(x) dim(x)[1] != 
                  dim(x)[2]))) 
                  stop("Elements of 'R' must be square matrices.")
                if (any(sapply(R[Rfix], function(x) !isSymmetric(unname(x))))) 
                  stop("Elements of 'R' must be symmetric matrices.")
                for (j in 1:length(R)) {
                  if (!Rfix[j]) 
                    next
                  if (is.null(rownames(R[[j]]))) 
                    rownames(R[[j]]) <- colnames(R[[j]])
                  if (is.null(colnames(R[[j]]))) 
                    colnames(R[[j]]) <- rownames(R[[j]])
                  if (is.null(colnames(R[[j]]))) 
                    stop("Elements of 'R' must have dimension names.")
                }
                R[Rfix] <- lapply(R[Rfix], unique, MARGIN = 1)
                R[Rfix] <- lapply(R[Rfix], unique, MARGIN = 2)
                if (any(sapply(R[Rfix], function(x) length(colnames(x)) != 
                  length(unique(colnames(x)))))) 
                  stop("Each element of 'R' must have unique dimension names.")
                for (j in 1:length(R)) {
                  if (!Rfix[j]) 
                    next
                  if (any(is.na(R[[j]]))) 
                    stop("No missing values allowed in matrix specified via 'R'.")
                  if (any(!is.element(s.levels[[j]], colnames(R[[j]])))) 
                    stop(paste0("There are levels in '", s.names[j], 
                      "' for which there are no rows/columns in the corresponding 'R' matrix."))
                  if (any(!is.element(colnames(R[[j]]), s.levels[[j]]))) 
                    warning(paste0("There are rows/columns in the 'R' matrix for '", 
                      s.names[j], "' for which there are no data."))
                }
            }
            else {
                withR <- FALSE
                Rfix <- rep(FALSE, s.nvals)
                R <- NULL
            }
        }
    }
    else {
        s.nvals <- 1
        sigma2 <- 0
        s.nlevels <- NULL
        s.levels <- NULL
        withR <- FALSE
        Rfix <- FALSE
        R <- NULL
    }
    if (withG) {
        if (very.verbose) 
            message(paste0("Processing '", as.character(random[mf.r.ncols == 
                2]), "' term ..."))
        mf.g.f <- data.frame(inner = factor(mf.g.f[[1]][[1]]), 
            outer = factor(mf.g.f[[1]][[2]]))
        mf.g <- data.frame(inner = factor(mf.g[[1]][[1]]), outer = factor(mf.g[[1]][[2]]))
        if (any(is.na(mf.g))) 
            stop("No NAs allowed in variables specified in the 'random' argument.")
        g.nlevels <- c(nlevels(mf.g[[1]]), nlevels(mf.g[[2]]))
        g.levels <- list(levels(mf.g[[1]]), levels(mf.g[[2]]))
        if (struct == "CS") {
            t.nvals <- 1
            r.nvals <- 1
        }
        if (struct == "HCS") {
            t.nvals <- g.nlevels[1]
            r.nvals <- 1
        }
        if (struct == "UN") {
            t.nvals <- g.nlevels[1]
            r.nvals <- ifelse(g.nlevels[1] > 1, g.nlevels[1] * 
                (g.nlevels[1] - 1)/2, 1)
        }
        if (struct == "AR") {
            t.nvals <- 1
            r.nvals <- 1
        }
        if (struct == "HAR") {
            t.nvals <- g.nlevels[1]
            r.nvals <- 1
        }
        if (struct == "UNHO") {
            t.nvals <- 1
            r.nvals <- ifelse(g.nlevels[1] > 1, g.nlevels[1] * 
                (g.nlevels[1] - 1)/2, 1)
        }
        if (is.null(tau2)) 
            tau2 <- rep(NA, t.nvals)
        if (is.null(rho)) 
            rho <- rep(NA, r.nvals)
        if (length(tau2) == 1 && is.element(struct, c("HCS", 
            "UN", "HAR"))) 
            tau2 <- rep(tau2, t.nvals)
        if (length(rho) == 1 && is.element(struct, c("UN", "UNHO"))) 
            rho <- rep(rho, r.nvals)
        if (length(tau2) != t.nvals) 
            stop(paste("Length of 'tau2' argument (", length(tau2), 
                ") does not match actual number of variance components (", 
                t.nvals, ").", sep = ""))
        if (length(rho) != r.nvals) 
            stop(paste("Length of 'rho' argument (", length(rho), 
                ") does not match actual number of correlations (", 
                r.nvals, ").", sep = ""))
        if (any(tau2 < 0, na.rm = TRUE)) 
            stop("Specified value(s) of 'tau2' must be non-negative.")
        if (any(rho > 1 | rho < -1, na.rm = TRUE)) 
            stop("Specified value(s) of 'rho' must be in [-1,1].")
        if (g.nlevels[1] == 1) {
            Z.G1 <- cbind(rep(1, k))
        }
        else {
            if (sparse) {
                Z.G1 <- sparse.model.matrix(~mf.g[[1]] - 1)
            }
            else {
                Z.G1 <- model.matrix(~mf.g[[1]] - 1)
            }
        }
        if (g.nlevels[2] == 1) {
            Z.G2 <- cbind(rep(1, k))
        }
        else {
            if (sparse) {
                Z.G2 <- sparse.model.matrix(~mf.g[[2]] - 1)
            }
            else {
                Z.G2 <- model.matrix(~mf.g[[2]] - 1)
            }
        }
    }
    else {
        t.nvals <- 1
        r.nvals <- 1
        tau2 <- 0
        rho <- 0
        Z.G1 <- NULL
        Z.G2 <- NULL
        g.nlevels <- NULL
        g.levels <- NULL
    }
    Vlt <- V
    Vlt[upper.tri(Vlt)] <- 0
    Vlt[lower.tri(V, diag = TRUE)] <- V[lower.tri(V, diag = TRUE)]
    if (is.null(A)) {
        has.na <- is.na(cbind(yi, mods)) | apply(is.na(Vlt), 
            MARGIN = 1, any)
    }
    else {
        has.na <- is.na(cbind(yi, mods)) | apply(is.na(Vlt), 
            MARGIN = 1, any) | apply(is.na(A), MARGIN = 1, any)
    }
    if (any(has.na)) {
        if (very.verbose) 
            message("Handling NAs ...")
        not.na <- apply(has.na, MARGIN = 1, sum) == 0L
        if (na.act == "na.omit" || na.act == "na.exclude" || 
            na.act == "na.pass") {
            yi <- yi[not.na]
            V <- V[not.na, not.na, drop = FALSE]
            A <- A[not.na, not.na, drop = FALSE]
            vi <- vi[not.na]
            ni <- ni[not.na]
            mods <- mods[not.na, , drop = FALSE]
            mf.s <- lapply(mf.s, function(x) x[not.na])
            mf.g <- mf.g[not.na, , drop = FALSE]
            mf.r <- lapply(mf.r, function(x) x[not.na, , drop = FALSE])
            Z.G1 <- Z.G1[not.na, , drop = FALSE]
            Z.G2 <- Z.G2[not.na, , drop = FALSE]
            k <- length(yi)
            warning("Rows with NAs omitted from model fitting.")
            attr(yi, "measure") <- measure
            attr(yi, "ni") <- ni
        }
        if (na.act == "na.fail") 
            stop("Missing values in data.")
    }
    else {
        not.na <- rep(TRUE, k)
    }
    if (k <= 1) 
        stop("Processing terminated since k <= 1.")
    if (is.null(mods) && !intercept) {
        warning("Must either include an intercept and/or moderators in model.\n  Coerced intercept into the model.")
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
    if (is.null(btt)) {
        if (p > 1) {
            if (int.incl) {
                btt <- seq.int(from = 2, to = p)
            }
            else {
                btt <- seq_len(p)
            }
        }
        else {
            btt <- 1
        }
    }
    else {
        btt <- btt[(btt >= 1) & (btt <= p)]
        btt <- unique(round(btt))
        if (length(btt) == 0L) 
            stop("Non-existent coefficients specified with 'btt'.")
    }
    bntt <- setdiff(seq_len(p), btt)
    m <- length(btt)
    if (withS) {
        mf.s <- lapply(mf.s, factor)
        s.nlevels <- sapply(mf.s, nlevels)
        s.levels <- lapply(mf.s, levels)
        if (any(is.na(sigma2) & s.nlevels == 1)) {
            sigma2[is.na(sigma2) & s.nlevels == 1] <- 0
            warning("Single-level factor(s) found in 'random' argument. Corresponding sigma2 value(s) fixed to 0.")
        }
        Z.S <- vector(mode = "list", length = s.nvals)
        for (j in 1:s.nvals) {
            if (s.nlevels[j] == 1) {
                Z.S[[j]] <- cbind(rep(1, k))
            }
            else {
                if (sparse) {
                  Z.S[[j]] <- sparse.model.matrix(~mf.s[[j]] - 
                    1)
                }
                else {
                  Z.S[[j]] <- model.matrix(~mf.s[[j]] - 1)
                }
            }
        }
    }
    else {
        Z.S <- NULL
    }
    if (withR) {
        for (j in 1:length(R)) {
            if (!Rfix[j]) 
                next
            R[[j]] <- R[[j]][s.levels[[j]], s.levels[[j]]]
        }
        if (Rscale == "cor" || Rscale == "cor0") 
            R[Rfix] <- lapply(R[Rfix], cov2cor)
        if (Rscale == "cor0") 
            R[Rfix] <- lapply(R[Rfix], function(x) (x - min(x))/(1 - 
                min(x)))
        if (Rscale == "cov0") 
            R[Rfix] <- lapply(R[Rfix], function(x) (x - min(x)))
    }
    if (withS) {
        D.S <- vector(mode = "list", length = s.nvals)
        for (j in seq_len(s.nvals)) {
            if (Rfix[j]) {
                if (sparse) {
                  D.S[[j]] <- Z.S[[j]] %*% Matrix(R[[j]], sparse = TRUE) %*% 
                    t(Z.S[[j]])
                }
                else {
                  D.S[[j]] <- Z.S[[j]] %*% R[[j]] %*% t(Z.S[[j]])
                }
            }
            else {
                D.S[[j]] <- tcrossprod(Z.S[[j]])
            }
        }
    }
    if (withG) {
        g.nlevels.f <- g.nlevels
        g.levels.f <- g.levels
        mf.g <- data.frame(inner = factor(mf.g[[1]]), outer = factor(mf.g[[2]]))
        g.nlevels <- c(nlevels(mf.g[[1]]), nlevels(mf.g[[2]]))
        g.levels <- list(levels(mf.g[[1]]), levels(mf.g[[2]]))
        g.levels.r <- !is.element(g.levels.f[[1]], g.levels[[1]])
        if (any(g.levels.r)) 
            warning("One or more levels of inner factor removed due to NAs.")
        if (g.nlevels[1] == 1 && is.element(struct, c("CS", "HCS", 
            "AR", "HAR")) && is.na(rho)) {
            rho <- 0
            warning("Inner factor has only a single level, so fixed value of 'rho' to 0.")
        }
        g.levels.k <- table(factor(mf.g[[1]], levels = g.levels.f[[1]]))
        if (is.element(struct, c("HCS", "UN", "HAR"))) {
            if (any(is.na(tau2) & g.levels.k == 1)) {
                tau2[is.na(tau2) & g.levels.k == 1] <- 0
                warning("Inner factor has k=1 for one or more levels. Corresponding tau2 value(s) fixed to 0.")
            }
        }
        g.levels.comb.k <- crossprod(Z.G2, Z.G1)
        g.levels.comb.k <- split(g.levels.comb.k, 1:nrow(g.levels.comb.k))
        if (all(unlist(lapply(g.levels.comb.k, sum)) == 1)) {
            if (is.element(struct, c("CS", "HCS", "AR", "HAR")) && 
                is.na(rho)) {
                rho <- 0
                warning("Each level of the outer factor contains only a single level of the inner factor, so fixed value of 'rho' to 0.")
            }
        }
        g.levels.comb.k <- lapply(g.levels.comb.k, function(x) outer(x, 
            x, FUN = "&"))
        g.levels.comb.k <- lapply(g.levels.comb.k, function(x) ifelse(x, 
            1, 0))
        g.levels.comb.k <- Reduce("+", g.levels.comb.k)
        g.levels.comb.k <- g.levels.comb.k[upper.tri(g.levels.comb.k)]
        if (is.element(struct, c("UN", "UNHO")) && any(g.levels.comb.k == 
            0 & is.na(rho))) {
            rho[g.levels.comb.k == 0] <- 0
            warning("Some combinations of the levels of the inner factor never occurred. Corresponding 'rho' value(s) fixed to 0.")
        }
        if (is.element(struct, c("UN", "UNHO")) && g.nlevels.f[1] == 
            1 && is.na(rho)) {
            rho <- 0
            warning("Inner factor has only a single level, so fixed value of 'rho' to 0.")
        }
        if (struct == "CS") {
            G <- matrix(rho, nrow = g.nlevels.f[1], ncol = g.nlevels.f[1])
            diag(G) <- 1
            G <- diag(sqrt(rep(tau2, g.nlevels.f[1])), nrow = g.nlevels.f[1], 
                ncol = g.nlevels.f[1]) %*% G %*% diag(sqrt(rep(tau2, 
                g.nlevels.f[1])), nrow = g.nlevels.f[1], ncol = g.nlevels.f[1])
            diag(G) <- tau2
        }
        if (struct == "HCS") {
            G <- matrix(rho, nrow = g.nlevels.f[1], ncol = g.nlevels.f[1])
            diag(G) <- 1
            G <- diag(sqrt(tau2), nrow = g.nlevels.f[1], ncol = g.nlevels.f[1]) %*% 
                G %*% diag(sqrt(tau2), nrow = g.nlevels.f[1], 
                ncol = g.nlevels.f[1])
            diag(G) <- tau2
        }
        if (struct == "UN") {
            G <- matrix(NA, nrow = g.nlevels.f[1], ncol = g.nlevels.f[1])
            G[upper.tri(G)] <- rho
            G[lower.tri(G)] <- t(G)[lower.tri(G)]
            diag(G) <- 1
            G <- diag(sqrt(tau2), nrow = g.nlevels.f[1], ncol = g.nlevels.f[1]) %*% 
                G %*% diag(sqrt(tau2), nrow = g.nlevels.f[1], 
                ncol = g.nlevels.f[1])
            diag(G) <- tau2
        }
        if (struct == "UNHO") {
            G <- matrix(NA, nrow = g.nlevels.f[1], ncol = g.nlevels.f[1])
            G[upper.tri(G)] <- rho
            G[lower.tri(G)] <- t(G)[lower.tri(G)]
            diag(G) <- 1
            G <- diag(sqrt(rep(tau2, g.nlevels.f[1])), nrow = g.nlevels.f[1], 
                ncol = g.nlevels.f[1]) %*% G %*% diag(sqrt(rep(tau2, 
                g.nlevels.f[1])), nrow = g.nlevels.f[1], ncol = g.nlevels.f[1])
            diag(G) <- tau2
        }
        if (struct == "AR") {
            if (is.na(rho)) {
                G <- matrix(NA, nrow = g.nlevels.f[1], ncol = g.nlevels.f[1])
            }
            else {
                if (g.nlevels.f[1] > 1) {
                  G <- toeplitz(ARMAacf(ar = rho, lag.max = g.nlevels.f[1] - 
                    1))
                }
                else {
                  G <- diag(1)
                }
            }
            G <- diag(sqrt(rep(tau2, g.nlevels.f[1])), nrow = g.nlevels.f[1], 
                ncol = g.nlevels.f[1]) %*% G %*% diag(sqrt(rep(tau2, 
                g.nlevels.f[1])), nrow = g.nlevels.f[1], ncol = g.nlevels.f[1])
            diag(G) <- tau2
        }
        if (struct == "HAR") {
            if (is.na(rho)) {
                G <- matrix(NA, nrow = g.nlevels.f[1], ncol = g.nlevels.f[1])
            }
            else {
                if (g.nlevels.f[1] > 1) {
                  G <- toeplitz(ARMAacf(ar = rho, lag.max = g.nlevels.f[1] - 
                    1))
                }
                else {
                  G <- diag(1)
                }
            }
            G <- diag(sqrt(tau2), nrow = g.nlevels.f[1], ncol = g.nlevels.f[1]) %*% 
                G %*% diag(sqrt(tau2), nrow = g.nlevels.f[1], 
                ncol = g.nlevels.f[1])
            diag(G) <- tau2
        }
        if (any(g.levels.r) && is.element(struct, c("CS", "AR"))) {
            G[g.levels.r, ] <- 0
            G[, g.levels.r] <- 0
        }
        if (any(g.levels.r) && is.element(struct, c("HCS", "HAR"))) {
            G[g.levels.r, ] <- 0
            G[, g.levels.r] <- 0
            tau2[g.levels.r] <- 0
            warning("Fixed 'tau2' to 0 for removed level(s).")
        }
        if (any(g.levels.r) && struct == "UN") {
            G[g.levels.r, ] <- 0
            G[, g.levels.r] <- 0
            tau2[g.levels.r] <- 0
            rho <- G[upper.tri(G)]
            warning("Fixed 'tau2' and corresponding 'rho' value(s) to 0 for removed level(s).")
        }
        if (any(g.levels.r) && struct == "UNHO") {
            G[g.levels.r, ] <- 0
            G[, g.levels.r] <- 0
            diag(G) <- tau2
            rho <- G[upper.tri(G)]
            warning("Fixed 'rho' value(s) to 0 corresponding to removed level(s).")
        }
        if (g.nlevels.f[1] == 2) {
            if (is.element(struct, c("CS", "AR", "UNHO")) && 
                !is.na(tau2) && tau2 == 0) 
                rho <- 0
            if (is.element(struct, c("HCS", "UN", "HAR")) && 
                ((!is.na(tau2[1]) && tau2[1] == 0) || (!is.na(tau2[2]) && 
                  tau2[2] == 0))) 
                rho <- 0
        }
    }
    else {
        G <- NULL
        g.levels.f <- NULL
        g.levels.r <- NULL
        g.levels.k <- NULL
        g.levels.comb.k <- NULL
        g.nlevels.f <- NULL
    }
    Y <- as.matrix(yi)
    if (very.verbose) 
        message("Extracting/computing initial values ...")
    if (verbose) {
        L.FE <- try(chol(V), silent = TRUE)
    }
    else {
        L.FE <- suppressWarnings(try(chol(V), silent = !verbose))
    }
    if (inherits(L.FE, "try-error")) {
        sigma2.init <- rep(0.001, s.nvals)
        tau2.init <- rep(0.001, t.nvals)
        rho.init <- rep(0.5, r.nvals)
    }
    else {
        W <- chol2inv(L.FE)
        U <- chol(W)
        sX <- U %*% X
        sY <- U %*% Y
        b.FE <- solve(crossprod(sX), crossprod(sX, sY))
        total <- max(0.001 * (s.nvals + t.nvals), var(as.vector(Y - 
            X %*% b.FE)) - 1/mean(1/diag(V)))
        sigma2.init <- rep(total/(s.nvals + t.nvals), s.nvals)
        tau2.init <- rep(total/(s.nvals + t.nvals), t.nvals)
        rho.init <- rep(0.5, r.nvals)
        QE <- sum(as.vector(sY - sX %*% b.FE)^2)
    }
    con <- list(verbose = FALSE, optimizer = "optim", optmethod = "Nelder-Mead", 
        sigma2.init = sigma2.init, tau2.init = tau2.init, rho.init = rho.init, 
        REMLf = TRUE, tol = 1e-07, posdefify = FALSE)
    con.pos <- pmatch(names(control), names(con))
    con[c(na.omit(con.pos))] <- control[!is.na(con.pos)]
    if (verbose) 
        con$verbose <- verbose
    verbose <- con$verbose
    silent <- !verbose
    if (withS && any(con$sigma2.init <= 0)) 
        stop("Values of 'sigma2.init' must be positive.")
    if (withG && any(con$tau2.init <= 0)) 
        stop("Values of 'tau2.init' must be positive.")
    if (withG && any(con$rho.init <= -1 | con$rho.init >= 1)) 
        stop("Values of 'rho.init' must be in (-1,1).")
    con <- list(verbose = verbose, optimizer = con$optimizer, 
        optmethod = con$optmethod, sigma2.init = log(con$sigma2.init), 
        tau2.init = log(con$tau2.init), rho.init = transf.rtoz(con$rho.init), 
        REMLf = con$REMLf, tol = con$tol, posdefify = con$posdefify)
    optimizer <- match.arg(con$optimizer, c("optim", "nlminb", 
        "uobyqa", "newuoa", "bobyqa"))
    optmethod <- con$optmethod
    tol <- con$tol
    posdefify <- con$posdefify
    optcontrol <- control[is.na(con.pos)]
    if (length(optcontrol) == 0) 
        optcontrol <- list()
    reml <- ifelse(method == "REML", TRUE, FALSE)
    if (is.element(optimizer, c("uobyqa", "newuoa", "bobyqa"))) {
        if (!require(minqa)) 
            stop("Please install the 'minqa' package to use this optimizer.")
    }
    if (withS) {
        if (length(con$sigma2.init) != s.nvals) 
            stop(paste("Length of 'sigma2.init' argument (", 
                length(con$sigma2.init), ") does not match actual number of variance components (", 
                s.nvals, ").", sep = ""))
    }
    else {
        con$sigma2.init <- -Inf
    }
    if (withG) {
        if (length(con$tau2.init) != t.nvals) 
            stop(paste("Length of 'tau2.init' argument (", length(con$tau2.init), 
                ") does not match actual number of variance components (", 
                t.nvals, ").", sep = ""))
    }
    else {
        con$tau2.init <- -Inf
    }
    if (withG) {
        if (length(con$rho.init) != r.nvals) 
            stop(paste("Length of 'rho.init' argument (", length(con$rho.init), 
                ") does not match actual number of correlations (", 
                r.nvals, ").", sep = ""))
    }
    else {
        con$rho.init <- 0
    }
    if (any(eigen(crossprod(X), symmetric = TRUE, only.values = TRUE)$values <= 
        tol)) 
        stop("Model matrix not of full rank. Cannot fit model.")
    if (withS) {
        sigma2.fix <- !is.na(sigma2)
    }
    else {
        sigma2.fix <- NA
    }
    if (withG) {
        tau2.fix <- !is.na(tau2)
        rho.fix <- !is.na(rho)
    }
    else {
        tau2.fix <- NA
        rho.fix <- NA
    }
    vc.fix <- list(sigma2 = sigma2.fix, tau2 = tau2.fix, rho = rho.fix)
    if (verbose) {
        cat("\nVariance Components in Model:")
        if (!withS && !withG) {
            cat(" none\n\n")
        }
        else {
            cat("\n\n")
            vcs <- rbind(c(sigma2 = round(exp(con$sigma2.init), 
                digits = digits), tau2 = round(exp(con$tau2.init), 
                digits = digits), rho = round(transf.ztor(con$rho.init), 
                digits = digits)), round(c(sigma2, tau2, rho), 
                digits = digits))
            vcs <- data.frame(vcs)
            rownames(vcs) <- c("initial", "specified")
            vcs <- rbind(included = ifelse(c(rep(withS, s.nvals), 
                rep(withG, t.nvals), rep(withG, r.nvals)), "Yes", 
                "No"), fixed = unlist(vc.fix), vcs)
            print(vcs)
            cat("\n")
        }
    }
    alpha <- ifelse(level > 1, (100 - level)/100, 1 - level)
    robust <- FALSE
    if (very.verbose) 
        message("Model fitting ...")
    if (optimizer == "optim") 
        par.arg <- "par"
    if (optimizer == "nlminb") 
        par.arg <- "start"
    if (is.element(optimizer, c("uobyqa", "newuoa", "bobyqa"))) 
        par.arg <- "par"
    optcall <- paste(optimizer, "(", par.arg, "=c(ifelse(is.na(sigma2), con$sigma2.init, sigma2),\n                                                 ifelse(is.na(tau2), con$tau2.init, tau2),\n                                                 ifelse(is.na(rho), con$rho.init, rho)), .ll.rma.mv, ", 
        ifelse(optimizer == "optim", "method=optmethod, ", ""), 
        "Y=Y, M=V, X=X, sigma2=sigma2, tau2=tau2, rho=rho, reml=reml,\n                                                 k=k, p=p,\n                                                 s.nvals=s.nvals, t.nvals=t.nvals, r.nvals=r.nvals,\n                                                 withS=withS, withG=withG, D.S=D.S, Z.G1=Z.G1, Z.G2=Z.G2,\n                                                 struct=struct, g.levels.r=g.levels.r,\n                                                 tol=tol, posdefify=posdefify, verbose=verbose, digits=digits,\n                                                 REMLf=con$REMLf, sparse=sparse, control=optcontrol)", 
        sep = "")
    opt.res <- try(eval(parse(text = optcall)), silent = silent)
    if (inherits(opt.res, "try-error")) 
        stop("Error during optimization.")
    if (is.element(optimizer, c("optim", "nlminb")) && opt.res$convergence != 
        0) 
        stop(paste0("Optimizer (", optimizer, ") did not achieve convergence (convergence = ", 
            opt.res$convergence, ")."))
    if (is.element(optimizer, c("uobyqa", "newuoa", "bobyqa")) && 
        opt.res$ierr != 0) 
        stop(paste0("Optimizer (", optimizer, ") did not achieve convergence (ierr = ", 
            opt.res$ierr, ")."))
    if (optimizer == "optim") 
        ll <- -1 * c(opt.res$value)
    if (optimizer == "nlminb") 
        ll <- -1 * c(opt.res$objective)
    if (is.element(optimizer, c("uobyqa", "newuoa", "bobyqa"))) 
        ll <- -1 * c(opt.res$fval)
    if (verbose) {
        cat("\n")
        print(opt.res)
    }
    if (p < k) {
        sigma2 <- ifelse(is.na(sigma2), exp(opt.res$par[1:s.nvals]), 
            sigma2)
        tau2 <- ifelse(is.na(tau2), exp(opt.res$par[(s.nvals + 
            1):(s.nvals + t.nvals)]), tau2)
        rho <- ifelse(is.na(rho), transf.ztor(opt.res$par[(s.nvals + 
            t.nvals + 1):(s.nvals + t.nvals + r.nvals)]), rho)
        sigma2 <- ifelse(sigma2 <= .Machine$double.eps * 10, 
            0, sigma2)
        tau2 <- ifelse(tau2 <= .Machine$double.eps * 10, 0, tau2)
    }
    else {
        sigma2[is.na(sigma2)] <- 0
        tau2[is.na(tau2)] <- 0
        rho[is.na(rho)] <- 0
    }
    M <- V
    if (withS) {
        for (j in seq_len(s.nvals)) {
            M <- M + sigma2[j] * D.S[[j]]
        }
    }
    if (withG) {
        ncol.Z.G1 <- ncol(Z.G1)
        if (struct == "CS") {
            G <- matrix(rho, nrow = ncol.Z.G1, ncol = ncol.Z.G1)
            diag(G) <- 1
            G <- diag(sqrt(rep(tau2, ncol.Z.G1)), nrow = ncol.Z.G1, 
                ncol = ncol.Z.G1) %*% G %*% diag(sqrt(rep(tau2, 
                ncol.Z.G1)), nrow = ncol.Z.G1, ncol = ncol.Z.G1)
            diag(G) <- tau2
        }
        if (struct == "HCS") {
            G <- matrix(rho, nrow = ncol.Z.G1, ncol = ncol.Z.G1)
            diag(G) <- 1
            G <- diag(sqrt(tau2), nrow = ncol.Z.G1, ncol = ncol.Z.G1) %*% 
                G %*% diag(sqrt(tau2), nrow = ncol.Z.G1, ncol = ncol.Z.G1)
            diag(G) <- tau2
        }
        if (struct == "UN") {
            G <- matrix(NA, nrow = ncol.Z.G1, ncol = ncol.Z.G1)
            G[upper.tri(G)] <- rho
            G[lower.tri(G)] <- t(G)[lower.tri(G)]
            diag(G) <- 1
            G <- diag(sqrt(tau2), nrow = ncol.Z.G1, ncol = ncol.Z.G1) %*% 
                G %*% diag(sqrt(tau2), nrow = ncol.Z.G1, ncol = ncol.Z.G1)
            diag(G) <- tau2
            if (posdefify) {
                G <- as.matrix(nearPD(G)$mat)
                tau2 <- diag(G)
                rho <- cov2cor(G)[upper.tri(G)]
            }
        }
        if (struct == "UNHO") {
            G <- matrix(NA, nrow = ncol.Z.G1, ncol = ncol.Z.G1)
            G[upper.tri(G)] <- rho
            G[lower.tri(G)] <- t(G)[lower.tri(G)]
            diag(G) <- 1
            G <- diag(sqrt(rep(tau2, ncol.Z.G1)), nrow = ncol.Z.G1, 
                ncol = ncol.Z.G1) %*% G %*% diag(sqrt(rep(tau2, 
                ncol.Z.G1)), nrow = ncol.Z.G1, ncol = ncol.Z.G1)
            diag(G) <- tau2
            if (posdefify) {
                G <- as.matrix(nearPD(G, keepDiag = TRUE)$mat)
                tau2 <- G[1, 1]
                rho <- cov2cor(G)[upper.tri(G)]
            }
        }
        if (struct == "AR") {
            if (ncol.Z.G1 > 1) {
                G <- toeplitz(ARMAacf(ar = rho, lag.max = ncol.Z.G1 - 
                  1))
            }
            else {
                G <- diag(1)
            }
            G <- diag(sqrt(rep(tau2, ncol.Z.G1)), nrow = ncol.Z.G1, 
                ncol = ncol.Z.G1) %*% G %*% diag(sqrt(rep(tau2, 
                ncol.Z.G1)), nrow = ncol.Z.G1, ncol = ncol.Z.G1)
            diag(G) <- tau2
        }
        if (struct == "HAR") {
            if (ncol.Z.G1 > 1) {
                G <- toeplitz(ARMAacf(ar = rho, lag.max = ncol.Z.G1 - 
                  1))
            }
            else {
                G <- diag(1)
            }
            G <- diag(sqrt(tau2), nrow = ncol.Z.G1, ncol = ncol.Z.G1) %*% 
                G %*% diag(sqrt(tau2), nrow = ncol.Z.G1, ncol = ncol.Z.G1)
            diag(G) <- tau2
        }
        if (any(g.levels.r)) {
            G[g.levels.r, ] <- 0
            G[, g.levels.r] <- 0
        }
        if (sparse) 
            G <- Matrix(G, sparse = TRUE)
        M <- M + (Z.G1 %*% G %*% t(Z.G1)) * tcrossprod(Z.G2)
        colnames(G) <- rownames(G) <- g.levels.f[[1]]
    }
    if (verbose) {
        L <- try(chol(M), silent = !verbose)
    }
    else {
        L <- suppressWarnings(try(chol(M), silent = !verbose))
    }
    if (inherits(L, "try-error")) {
        stop("Final variance-covariance matrix not positive definite.")
    }
    else {
        W <- chol2inv(L)
        U <- chol(W)
        sX <- U %*% X
    }
    if (is.null(A)) {
        sY <- U %*% Y
        vb <- matrix(solve(crossprod(sX)), nrow = p, ncol = p)
        b <- matrix(vb %*% crossprod(sX, sY), ncol = 1)
        rownames(b) <- colnames(X)
        RSS.f <- sum(as.vector(sY - sX %*% b)^2)
    }
    else {
        stXAX <- chol2inv(chol(as.matrix(t(X) %*% A %*% X)))
        b <- matrix(stXAX %*% crossprod(X, A) %*% Y, ncol = 1)
        vb <- matrix(stXAX %*% t(X) %*% A %*% M %*% A %*% X %*% 
            stXAX, nrow = p, ncol = p)
        rownames(b) <- colnames(X)
        RSS.f <- as.vector(t(Y - X %*% b) %*% W %*% (Y - X %*% 
            b))
    }
    se <- sqrt(diag(vb))
    names(se) <- NULL
    zval <- c(b/se)
    pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)
    QM <- as.vector(t(b)[btt] %*% chol2inv(chol(vb[btt, btt])) %*% 
        b[btt])
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
    if (very.verbose) 
        message("Heterogeneity testing ...")
    QE.df <- k - p
    if (QE.df > 0L) {
        if (inherits(L.FE, "try-error")) {
            QE <- NA
            QEp <- NA
        }
        else {
            QEp <- pchisq(QE, df = QE.df, lower.tail = FALSE)
        }
    }
    else {
        QE <- 0
        QEp <- 1
    }
    ll.QE <- -1/2 * (k) * log(2 * base::pi) - 1/2 * determinant(V, 
        logarithm = TRUE)$modulus
    if (very.verbose) 
        message("Computing fit statistics and log likelihood ...")
    parms <- p + ifelse(withS, sum(ifelse(sigma2.fix, 0, 1)), 
        0) + ifelse(withG, sum(ifelse(tau2.fix, 0, 1)), 0) + 
        ifelse(withG, sum(ifelse(rho.fix, 0, 1)), 0)
    ll.ML <- -1/2 * (k) * log(2 * base::pi) - 1/2 * determinant(M, 
        logarithm = TRUE)$modulus - 1/2 * RSS.f
    ll.REML <- -1/2 * (k - p) * log(2 * base::pi) + ifelse(con$REMLf, 
        1/2 * determinant(crossprod(X), logarithm = TRUE)$modulus, 
        0) - 1/2 * determinant(M, logarithm = TRUE)$modulus - 
        1/2 * determinant(crossprod(X, W) %*% X, logarithm = TRUE)$modulus - 
        1/2 * RSS.f
    dev.ML <- -2 * (ll.ML - ll.QE)
    AIC.ML <- -2 * ll.ML + 2 * parms
    BIC.ML <- -2 * ll.ML + parms * log(k)
    AICc.ML <- -2 * ll.ML + 2 * parms * max(k, parms + 2)/(max(k, 
        parms + 2) - parms - 1)
    dev.REML <- -2 * (ll.REML - 0)
    AIC.REML <- -2 * ll.REML + 2 * parms
    BIC.REML <- -2 * ll.REML + parms * log(k - p)
    AICc.REML <- -2 * ll.REML + 2 * parms * max(k - p, parms + 
        2)/(max(k - p, parms + 2) - parms - 1)
    fit.stats <- matrix(c(ll.ML, dev.ML, AIC.ML, BIC.ML, AICc.ML, 
        ll.REML, dev.REML, AIC.REML, BIC.REML, AICc.REML), ncol = 2, 
        byrow = FALSE)
    dimnames(fit.stats) <- list(c("ll", "dev", "AIC", "BIC", 
        "AICc"), c("ML", "REML"))
    fit.stats <- data.frame(fit.stats)
    if (very.verbose) 
        message("Preparing output ...")
    p.eff <- p
    k.eff <- k
    weighted <- TRUE
    robust <- FALSE
    res <- list(b = b, se = se, zval = zval, pval = pval, ci.lb = ci.lb, 
        ci.ub = ci.ub, vb = vb, sigma2 = sigma2, tau2 = tau2, 
        rho = rho, k = k, k.f = k.f, k.eff = k.eff, p = p, p.eff = p.eff, 
        parms = parms, m = m, QE = QE, QEp = QEp, QM = QM, QMp = QMp, 
        int.only = int.only, int.incl = int.incl, allvipos = allvipos, 
        yi = yi, vi = vi, V = V, W = A, X = X, yi.f = yi.f, vi.f = vi.f, 
        V.f = V.f, X.f = X.f, ni = ni, ni.f = ni.f, M = M, G = G, 
        ids = ids, not.na = not.na, slab = slab, slab.null = slab.null, 
        measure = measure, method = method, weighted = weighted, 
        knha = knha, robust = robust, btt = btt, intercept = intercept, 
        digits = digits, level = level, control = control, fit.stats = fit.stats, 
        vc.fix = vc.fix, withS = withS, withG = withG, withR = withR, 
        s.nvals = s.nvals, t.nvals = t.nvals, r.nvals = r.nvals, 
        s.names = s.names, g.names = g.names, s.nlevels = s.nlevels, 
        g.nlevels.f = g.nlevels.f, g.nlevels = g.nlevels, g.levels.f = g.levels.f, 
        g.levels.k = g.levels.k, g.levels.comb.k = g.levels.comb.k, 
        struct = struct, Rfix = Rfix, R = R, Rscale = Rscale, 
        mf.r = mf.r, mf.g.f = mf.g.f, random = random)
    class(res) <- c("rma.mv", "rma")
    return(res)
}
rma.peto <-
function (ai, bi, ci, di, n1i, n2i, data, slab, subset, add = 1/2, 
    to = "only0", drop00 = TRUE, level = 95, digits = 4, verbose = FALSE) 
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
    measure <- "PETO"
    if (verbose) 
        message("Extracting data and computing yi/vi values ...")
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
    ids <- seq_len(k)
    if (verbose) 
        message("Generating/extracting study labels ...")
    if (is.null(slab)) {
        slab.null <- TRUE
        slab <- ids
    }
    else {
        if (any(is.na(slab))) 
            stop("NAs in study labels.")
        if (any(duplicated(slab))) 
            slab <- make.unique(slab)
        if (length(slab) != k) 
            stop("Study labels not of same length as data.")
        slab.null <- FALSE
    }
    if (!is.null(subset)) {
        if (verbose) 
            message("Subsetting ...")
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
        if (verbose) 
            message("Handling NAs in table data ...")
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
        if (verbose) 
            message("Handling NAs in yi/vi ...")
        not.na.yivi <- apply(yivi.na, MARGIN = 1, sum) == 0L
        if (na.act == "na.omit" || na.act == "na.exclude" || 
            na.act == "na.pass") {
            yi <- yi[not.na.yivi]
            vi <- vi[not.na.yivi]
            ni <- ni[not.na.yivi]
            warning("Some yi/vi values are NA.")
            attr(yi, "measure") <- measure
            attr(yi, "ni") <- ni
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
    n1i <- ai + bi
    n2i <- ci + di
    Ni <- ai + bi + ci + di
    alpha <- ifelse(level > 1, (100 - level)/100, 1 - level)
    if (verbose) 
        message("Model fitting ...")
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
    if (verbose) 
        message("Heterogeneity testing ...")
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
    if (verbose) 
        message("Computing fit statistics and log likelihood ...")
    ll.ML <- -1/2 * (k.yi) * log(2 * base::pi) - 1/2 * sum(log(vi)) - 
        1/2 * RSS
    ll.REML <- -1/2 * (k.yi - 1) * log(2 * base::pi) + 1/2 * 
        log(k.yi) - 1/2 * sum(log(vi)) - 1/2 * log(sum(wi)) - 
        1/2 * RSS
    dev.ML <- -2 * (ll.ML - sum(dnorm(yi, mean = yi, sd = sqrt(vi), 
        log = TRUE)))
    AIC.ML <- -2 * ll.ML + 2
    BIC.ML <- -2 * ll.ML + log(k.yi)
    AICc.ML <- -2 * ll.ML + 2 * max(k.yi, 3)/(max(k.yi, 3) - 
        2)
    dev.REML <- -2 * (ll.REML - 0)
    AIC.REML <- -2 * ll.REML + 2
    BIC.REML <- -2 * ll.REML + log(k.yi - 1)
    AICc.REML <- -2 * ll.REML + 2 * max(k.yi - 1, 3)/(max(k.yi - 
        1, 3) - 2)
    fit.stats <- matrix(c(ll.ML, dev.ML, AIC.ML, BIC.ML, AICc.ML, 
        ll.REML, dev.REML, AIC.REML, BIC.REML, AICc.REML), ncol = 2, 
        byrow = FALSE)
    dimnames(fit.stats) <- list(c("ll", "dev", "AIC", "BIC", 
        "AICc"), c("ML", "REML"))
    fit.stats <- data.frame(fit.stats)
    if (verbose) 
        message("Preparing output ...")
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
    res <- list(b = b, se = se, zval = zval, pval = pval, ci.lb = ci.lb, 
        ci.ub = ci.ub, vb = vb, tau2 = tau2, k = k, k.f = k.f, 
        k.yi = k.yi, k.pos = k.pos, k.eff = k.eff, p = p, parms = parms, 
        QE = QE, QEp = QEp, int.only = int.only, yi = yi, vi = vi, 
        yi.f = yi.f, vi.f = vi.f, X.f = X.f, ai = ai, bi = bi, 
        ci = ci, di = di, ai.f = ai.f, bi.f = bi.f, ci.f = ci.f, 
        di.f = di.f, ni = ni, ni.f = ni.f, ids = ids, not.na = not.na, 
        not.na.yivi = not.na.yivi, slab = slab, slab.null = slab.null, 
        measure = measure, method = method, weighted = weighted, 
        knha = knha, robust = robust, intercept = intercept, 
        digits = digits, level = level, add = add, to = to, drop00 = drop00, 
        fit.stats = fit.stats)
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
        "AS", "PHI", "YUQ", "YUY", "RTET", "PBIT", "OR2D", "OR2DN", 
        "OR2DL", "IRR", "IRD", "IRSD", "MD", "SMD", "SMDH", "ROM", 
        "RPB", "RBIS", "D2OR", "D2ORN", "D2ORL", "COR", "UCOR", 
        "ZCOR", "PR", "PLN", "PLO", "PAS", "PFT", "IR", "IRLN", 
        "IRS", "IRFT", "MN", "MC", "SMCC", "SMCR", "ARAW", "AHW", 
        "ABT"))) 
        stop("Unknown 'measure' specified.")
    if (!is.element(method, c("FE", "HS", "HE", "DL", "SJ", "ML", 
        "REML", "EB", "DLIT", "SJIT", "PM"))) 
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
    very.verbose <- ifelse(!is.logical(verbose) && verbose > 
        1, TRUE, FALSE)
    if (very.verbose) 
        message("Extracting/computing yi/vi values ...")
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
    mf.weights <- mf[[match("weights", names(mf))]]
    mf.slab <- mf[[match("slab", names(mf))]]
    mf.subset <- mf[[match("subset", names(mf))]]
    mf.mods <- mf[[match("mods", names(mf))]]
    yi <- eval(mf.yi, data, enclos = sys.frame(sys.parent()))
    weights <- eval(mf.weights, data, enclos = sys.frame(sys.parent()))
    slab <- eval(mf.slab, data, enclos = sys.frame(sys.parent()))
    subset <- eval(mf.subset, data, enclos = sys.frame(sys.parent()))
    mods <- eval(mf.mods, data, enclos = sys.frame(sys.parent()))
    ai <- bi <- ci <- di <- x1i <- x2i <- t1i <- t2i <- NA
    is.formula <- FALSE
    if (!is.null(yi)) {
        if (class(yi) == "formula") {
            options(na.action = "na.pass")
            mods <- model.matrix(yi, data = data)
            attr(mods, "assign") <- NULL
            yi <- model.response(model.frame(yi, data = data))
            options(na.action = na.act)
            names(yi) <- NULL
            intercept <- FALSE
            is.formula <- TRUE
        }
        if (is.matrix(yi)) 
            yi <- as.vector(yi)
        k <- length(yi)
        if (measure == "GEN") {
            if (!is.null(attr(yi, "measure"))) 
                measure <- attr(yi, "measure")
        }
        attr(yi, "measure") <- measure
        mf.vi <- mf[[match("vi", names(mf))]]
        mf.sei <- mf[[match("sei", names(mf))]]
        mf.ni <- mf[[match("ni", names(mf))]]
        vi <- eval(mf.vi, data, enclos = sys.frame(sys.parent()))
        sei <- eval(mf.sei, data, enclos = sys.frame(sys.parent()))
        ni <- eval(mf.ni, data, enclos = sys.frame(sys.parent()))
        if (is.null(vi)) {
            if (is.null(sei)) {
                stop("Need to specify vi or sei argument.")
            }
            else {
                vi <- sei^2
            }
        }
        if (is.matrix(vi)) 
            vi <- as.vector(vi)
        if (length(vi) != k) 
            stop("Length of yi and vi (or sei) vectors are not the same.")
        if (is.null(ni) && !is.null(attr(yi, "ni"))) 
            ni <- attr(yi, "ni")
        if (!is.null(ni) && (length(ni) != k)) 
            ni <- NULL
        if (!is.null(ni)) 
            attr(yi, "ni") <- ni
        if (is.null(slab) & !is.null(attr(yi, "slab"))) 
            slab <- attr(yi, "slab")
        if (!is.null(subset)) {
            yi <- yi[subset]
            vi <- vi[subset]
            ni <- ni[subset]
            attr(yi, "measure") <- measure
            attr(yi, "ni") <- ni
        }
    }
    else {
        if (is.element(measure, c("RR", "OR", "PETO", "RD", "AS", 
            "PHI", "YUQ", "YUY", "RTET", "PBIT", "OR2D", "OR2DN", 
            "OR2DL"))) {
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
            k <- length(ai)
            if (!is.null(subset)) {
                ai <- ai[subset]
                bi <- bi[subset]
                ci <- ci[subset]
                di <- di[subset]
            }
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
            k <- length(x1i)
            if (!is.null(subset)) {
                x1i <- x1i[subset]
                x2i <- x2i[subset]
                t1i <- t1i[subset]
                t2i <- t2i[subset]
            }
            dat <- escalc(measure, x1i = x1i, x2i = x2i, t1i = t1i, 
                t2i = t2i, add = add, to = to, drop00 = drop00, 
                vtype = vtype)
        }
        if (is.element(measure, c("MD", "SMD", "SMDH", "ROM", 
            "RPB", "RBIS", "D2OR", "D2ORN", "D2ORL"))) {
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
            k <- length(m1i)
            if (!is.null(subset)) {
                m1i <- m1i[subset]
                m2i <- m2i[subset]
                sd1i <- sd1i[subset]
                sd2i <- sd2i[subset]
                n1i <- n1i[subset]
                n2i <- n2i[subset]
            }
            dat <- escalc(measure, m1i = m1i, m2i = m2i, sd1i = sd1i, 
                sd2i = sd2i, n1i = n1i, n2i = n2i, vtype = vtype)
        }
        if (is.element(measure, c("COR", "UCOR", "ZCOR"))) {
            mf.ri <- mf[[match("ri", names(mf))]]
            mf.ni <- mf[[match("ni", names(mf))]]
            ri <- eval(mf.ri, data, enclos = sys.frame(sys.parent()))
            ni <- eval(mf.ni, data, enclos = sys.frame(sys.parent()))
            k <- length(ri)
            if (!is.null(subset)) {
                ri <- ri[subset]
                ni <- ni[subset]
            }
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
            k <- length(xi)
            if (!is.null(subset)) {
                xi <- xi[subset]
                mi <- mi[subset]
            }
            dat <- escalc(measure, xi = xi, mi = mi, add = add, 
                to = to, vtype = vtype)
        }
        if (is.element(measure, c("IR", "IRLN", "IRS", "IRFT"))) {
            mf.xi <- mf[[match("xi", names(mf))]]
            mf.ti <- mf[[match("ti", names(mf))]]
            xi <- eval(mf.xi, data, enclos = sys.frame(sys.parent()))
            ti <- eval(mf.ti, data, enclos = sys.frame(sys.parent()))
            k <- length(xi)
            if (!is.null(subset)) {
                xi <- xi[subset]
                ti <- ti[subset]
            }
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
            k <- length(mi)
            if (!is.null(subset)) {
                mi <- mi[subset]
                sdi <- sdi[subset]
                ni <- ni[subset]
            }
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
            k <- length(m1i)
            if (!is.null(subset)) {
                m1i <- m1i[subset]
                m2i <- m2i[subset]
                sd1i <- sd1i[subset]
                sd2i <- sd2i[subset]
                ni <- ni[subset]
                ri <- ri[subset]
            }
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
            k <- length(ai)
            if (!is.null(subset)) {
                ai <- ai[subset]
                mi <- mi[subset]
                ni <- ni[subset]
            }
            dat <- escalc(measure, ai = ai, mi = mi, ni = ni, 
                vtype = vtype)
        }
        if (is.element(measure, "GEN")) 
            stop("Specify the desired outcome measure via the 'measure' argument.")
        yi <- dat$yi
        vi <- dat$vi
        ni <- attr(yi, "ni")
    }
    if (length(weights) == 1) 
        weights <- rep(weights, k)
    if (!is.null(weights) && (length(weights) != k)) 
        stop("Length of yi and weights vectors are not the same.")
    if (!is.null(subset)) 
        weights <- weights[subset]
    if (very.verbose) 
        message("Creating model matrix ...")
    if (class(mods) == "formula") {
        options(na.action = "na.pass")
        mods <- model.matrix(mods, data = data)
        attr(mods, "assign") <- NULL
        options(na.action = na.act)
        intercept <- FALSE
        is.formula <- TRUE
    }
    if (is.vector(mods)) 
        mods <- cbind(mods)
    if (is.data.frame(mods)) 
        mods <- as.matrix(mods)
    if (is.character(mods)) 
        stop("Model matrix contains character variables.")
    if (!is.null(mods) && (nrow(mods) != k)) 
        stop("Number of rows of the model matrix does not match length of yi argument.")
    ids <- seq_len(k)
    if (very.verbose) 
        message("Generating/extracting study labels ...")
    if (is.null(slab)) {
        slab.null <- TRUE
        slab <- ids
    }
    else {
        if (any(is.na(slab))) 
            stop("NAs in study labels.")
        if (any(duplicated(slab))) 
            slab <- make.unique(slab)
        if (length(slab) != k) 
            stop("Study labels not of same length as data.")
        slab.null <- FALSE
    }
    if (!is.null(subset)) {
        if (very.verbose) 
            message("Subsetting ...")
        mods <- mods[subset, , drop = FALSE]
        slab <- slab[subset]
        ids <- ids[subset]
    }
    k <- length(yi)
    if (any(vi <= 0, na.rm = TRUE)) {
        allvipos <- FALSE
        warning("There are outcomes with non-positive sampling variances.")
        vi.neg <- vi < 0
        if (any(vi.neg, na.rm = TRUE)) {
            vi[vi.neg] <- 0
            warning("Negative sampling variances constrained to zero.")
        }
    }
    else {
        allvipos <- TRUE
    }
    if (any(weights < 0, na.rm = TRUE)) 
        stop("Negative weights not allowed.")
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
    weights.f <- weights
    ni.f <- ni
    mods.f <- mods
    k.f <- k
    YVMW.na <- is.na(cbind(yi, vi, mods, weights))
    if (any(YVMW.na)) {
        if (very.verbose) 
            message("Handling NAs ...")
        not.na <- apply(YVMW.na, MARGIN = 1, sum) == 0L
        if (na.act == "na.omit" || na.act == "na.exclude" || 
            na.act == "na.pass") {
            yi <- yi[not.na]
            vi <- vi[not.na]
            weights <- weights[not.na]
            ni <- ni[not.na]
            mods <- mods[not.na, , drop = FALSE]
            k <- length(yi)
            warning("Studies with NAs omitted from model fitting.")
            attr(yi, "measure") <- measure
            attr(yi, "ni") <- ni
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
    if (is.null(mods) && !intercept) {
        warning("Must either include an intercept and/or moderators in model.\n  Coerced intercept into the model.")
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
                btt <- seq.int(from = 2, to = p)
            }
            else {
                btt <- seq_len(p)
            }
        }
        else {
            btt <- 1
        }
    }
    else {
        btt <- btt[(btt >= 1) & (btt <= p)]
        btt <- unique(round(btt))
        if (length(btt) == 0L) 
            stop("Non-existent coefficients specified with 'btt'.")
    }
    bntt <- setdiff(seq_len(p), btt)
    m <- length(btt)
    con <- list(tau2.init = NULL, tau2.min = 0, tau2.max = 50, 
        threshold = 10^-5, maxiter = 100, stepadj = 1, REMLf = TRUE, 
        verbose = FALSE, tol = 1e-07)
    con.pos <- pmatch(names(control), names(con))
    con[c(na.omit(con.pos))] <- control[!is.na(con.pos)]
    if (verbose) 
        con$verbose <- verbose
    verbose <- con$verbose
    if (any(eigen(crossprod(X), symmetric = TRUE, only.values = TRUE)$values <= 
        con$tol)) 
        stop("Model matrix not of full rank. Cannot fit model.")
    iter <- 0
    se.tau2 <- I2 <- H2 <- QE <- QEp <- NA
    s2w <- 1
    alpha <- ifelse(level > 1, (100 - level)/100, 1 - level)
    robust <- FALSE
    Y <- as.matrix(yi)
    if (is.numeric(tau2)) {
        tau2.fix <- TRUE
        tau2.val <- tau2
    }
    else {
        tau2.fix <- FALSE
        tau2.val <- NA
    }
    if (very.verbose && !tau2.fix) 
        message("Estimating tau^2 value ...")
    if (method == "HS") {
        if (!allvipos) 
            stop("HS estimator cannot be used with non-positive sampling variances.")
        wi <- 1/vi
        W <- diag(wi, nrow = k, ncol = k)
        stXWX <- .invcalc(X = X, W = W, k = k)
        P <- W - W %*% X %*% stXWX %*% crossprod(X, W)
        RSS <- crossprod(Y, P) %*% Y
        tau2 <- ifelse(tau2.fix, tau2.val, (RSS - k)/sum(wi))
    }
    if (is.element(method, c("HE", "ML", "REML", "EB"))) {
        stXX <- .invcalc(X = X, W = diag(k), k = k)
        P <- diag(k) - X %*% tcrossprod(stXX, X)
        RSS <- crossprod(Y, P) %*% Y
        V <- diag(vi, nrow = k, ncol = k)
        PV <- P %*% V
        trPV <- .tr(PV)
        tau2 <- ifelse(tau2.fix, tau2.val, (RSS - trPV)/(k - 
            p))
    }
    if (method == "DL") {
        if (!allvipos) 
            stop("DL estimator cannot be used with non-positive sampling variances.")
        wi <- 1/vi
        W <- diag(wi, nrow = k, ncol = k)
        stXWX <- .invcalc(X = X, W = W, k = k)
        P <- W - W %*% X %*% stXWX %*% crossprod(X, W)
        RSS <- crossprod(Y, P) %*% Y
        trP <- .tr(P)
        tau2 <- ifelse(tau2.fix, tau2.val, (RSS - (k - p))/trP)
    }
    if (method == "SJ") {
        if (is.null(con$tau2.init)) {
            tau2.0 <- c(var(yi) * (k - 1)/k)
        }
        else {
            tau2.0 <- con$tau2.init
        }
        wi <- 1/(vi + tau2.0)
        W <- diag(wi, nrow = k, ncol = k)
        stXWX <- .invcalc(X = X, W = W, k = k)
        P <- W - W %*% X %*% stXWX %*% crossprod(X, W)
        RSS <- crossprod(Y, P) %*% Y
        V <- diag(vi, nrow = k, ncol = k)
        PV <- P %*% V
        tau2 <- ifelse(tau2.fix, tau2.val, tau2.0 * RSS/(k - 
            p))
    }
    if (method == "DLIT") {
        conv <- 1
        change <- con$threshold + 1
        if (is.null(con$tau2.init)) {
            tau2 <- 0
        }
        else {
            tau2 <- con$tau2.init
        }
        while (change > con$threshold) {
            if (verbose) 
                cat("Iteration", iter, "\ttau^2 =", round(tau2, 
                  digits), "\n")
            iter <- iter + 1
            tau2.old <- tau2
            wi <- 1/(vi + tau2)
            W <- diag(wi, nrow = k, ncol = k)
            stXWX <- .invcalc(X = X, W = W, k = k)
            P <- W - W %*% X %*% stXWX %*% crossprod(X, W)
            RSS <- crossprod(Y, P) %*% Y
            trP <- .tr(P)
            tau2 <- ifelse(tau2.fix, tau2.val, (RSS - (k - p))/trP)
            tau2[tau2 < con$tau2.min] <- con$tau2.min
            change <- abs(tau2.old - tau2)
            if (iter > con$maxiter) {
                conv <- 0
                break
            }
        }
        if (conv == 0L) 
            stop("Algorithm did not converge.")
    }
    if (method == "SJIT") {
        conv <- 1
        change <- con$threshold + 1
        if (is.null(con$tau2.init)) {
            tau2 <- var(yi) * (k - 1)/k
            tau2.0 <- tau2
        }
        else {
            tau2 <- con$tau2.init
            tau2.0 <- tau2
        }
        while (change > con$threshold) {
            if (verbose) 
                cat("Iteration", iter, "\ttau^2 =", round(tau2, 
                  digits), "\n")
            iter <- iter + 1
            tau2.old <- tau2
            wi <- 1/(vi + tau2)
            W <- diag(wi, nrow = k, ncol = k)
            stXWX <- .invcalc(X = X, W = W, k = k)
            P <- W - W %*% X %*% stXWX %*% crossprod(X, W)
            RSS <- crossprod(Y, P) %*% Y
            V <- diag(vi, nrow = k, ncol = k)
            PV <- P %*% V
            tau2 <- ifelse(tau2.fix, tau2.val, tau2 * RSS/(k - 
                p))
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
        if (.QE.func(con$tau2.min, Y = Y, vi = vi, X = X, k = k, 
            objective = k - p) < 0) {
            tau2 <- con$tau2.min
        }
        else {
            tau2 <- ifelse(tau2.fix, tau2.val, try(uniroot(.QE.func, 
                interval = c(con$tau2.min, con$tau2.max), tol = con$threshold, 
                Y = Y, vi = vi, X = X, k = k, objective = k - 
                  p)$root, silent = TRUE))
            if (!is.numeric(tau2)) 
                stop("Error in iterative search for tau2. Try increasing tau2.max or switch to another 'method'.")
            wi <- 1/(vi + tau2)
            W <- diag(wi, nrow = k, ncol = k)
            stXWX <- .invcalc(X = X, W = W, k = k)
            P <- W - W %*% X %*% stXWX %*% crossprod(X, W)
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
            if (verbose) 
                cat("Iteration", iter, "\ttau^2 =", round(tau2, 
                  digits), "\n")
            iter <- iter + 1
            tau2.old <- tau2
            wi <- 1/(vi + tau2)
            W <- diag(wi, nrow = k, ncol = k)
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
                adj <- (crossprod(Y, P) %*% Y * k/(k - p) - k)/sum(wi)
            }
            adj <- adj * con$stepadj
            while (tau2 + adj < con$tau2.min) {
                adj <- adj/2
            }
            tau2 <- ifelse(tau2.fix, tau2.val, tau2 + adj)
            change <- abs(tau2.old - tau2)
            if (iter > con$maxiter) {
                conv <- 0
                break
            }
        }
        if (conv == 0L) 
            stop("Fisher scoring algorithm did not converge. Try increasing the number of iterations,\n  adjusting the threshold, or use a different estimator for tau^2.")
    }
    tau2 <- max(con$tau2.min, c(tau2))
    if (verbose && is.element(method, c("ML", "REML", "EB"))) 
        cat("Fisher scoring algorithm converged after", iter, 
            "iterations.\n")
    if (method == "HS") {
        se.tau2 <- sqrt(1/sum(wi)^2 * (2 * (k - p) + 4 * max(tau2, 
            0) * .tr(P) + 2 * max(tau2, 0)^2 * sum(P * P)))
    }
    if (method == "HE") {
        se.tau2 <- sqrt(1/(k - p)^2 * (2 * sum(PV * t(PV)) + 
            4 * max(tau2, 0) * trPV + 2 * max(tau2, 0)^2 * (k - 
            p)))
    }
    if (method == "DL" || method == "DLIT") {
        se.tau2 <- sqrt(1/trP^2 * (2 * (k - p) + 4 * max(tau2, 
            0) * trP + 2 * max(tau2, 0)^2 * sum(P * P)))
    }
    if (method == "SJ") {
        se.tau2 <- sqrt(tau2.0^2/(k - p)^2 * (2 * sum(PV * t(PV)) + 
            4 * max(tau2, 0) * sum(PV * P) + 2 * max(tau2, 0)^2 * 
            sum(P * P)))
    }
    if (method == "ML") {
        se.tau2 <- sqrt(2/sum(wi^2))
    }
    if (method == "REML") {
        se.tau2 <- sqrt(2/sum(P * P))
    }
    if (method == "EB" || method == "PM" || method == "SJIT") {
        V <- diag(vi, nrow = k, ncol = k)
        PV <- P %*% V
        se.tau2 <- sqrt((k/(k - p))^2/sum(wi)^2 * (2 * sum(PV * 
            t(PV)) + 4 * max(tau2, 0) * sum(PV * P) + 2 * max(tau2, 
            0)^2 * sum(P * P)))
    }
    if (method == "FE") {
        tau2 <- 0
        if (!allvipos && weighted) 
            stop("Weighted estimation cannot be used with a fixed-effects\n  model when there are non-positive sampling variances.")
    }
    if (very.verbose) 
        message("Model fitting ...")
    wi <- 1/(vi + tau2)
    W <- diag(wi, nrow = k, ncol = k)
    if (weighted) {
        if (is.null(weights)) {
            stXWX <- .invcalc(X = X, W = W, k = k)
            b <- stXWX %*% crossprod(X, W) %*% Y
            vb <- stXWX
            RSS.f <- sum(wi * (yi - X %*% b)^2)
        }
        else {
            A <- diag(weights, nrow = k, ncol = k)
            stXAX <- .invcalc(X = X, W = A, k = k)
            b <- stXAX %*% crossprod(X, A) %*% Y
            vb <- stXAX %*% t(X) %*% A %*% diag(vi + tau2, nrow = k, 
                ncol = k) %*% A %*% X %*% stXAX
            RSS.f <- sum(wi * (yi - X %*% b)^2)
        }
        if (robust) {
            ei <- c(yi - X %*% b)
            vb <- vb %*% t(X) %*% W %*% diag(ei^2, nrow = k, 
                ncol = k) %*% W %*% X %*% vb
            vb <- vb * k/(k - p)
        }
        if (knha) {
            if (RSS.f <= .Machine$double.eps) {
                s2w <- 1
            }
            else {
                s2w <- RSS.f/(k - p)
            }
            vb <- s2w * vb
            if (method == "FE") 
                warning("Knapp & Hartung (2003) method is not meant to be used in the context of fixed-effects models.")
        }
        QM <- c(t(b)[btt] %*% chol2inv(chol(vb[btt, btt])) %*% 
            b[btt])
    }
    else {
        stXX <- .invcalc(X = X, W = diag(k), k = k)
        b <- stXX %*% crossprod(X, Y)
        vb <- tcrossprod(stXX, X) %*% diag(vi + tau2, nrow = k, 
            ncol = k) %*% X %*% stXX
        RSS.f <- sum(wi * (yi - X %*% b)^2)
        if (robust) {
            ei <- c(Y - X %*% b)
            vb <- stXX %*% t(X) %*% diag(ei^2, nrow = k, ncol = k) %*% 
                X %*% stXX
            vb <- vb * k/(k - p)
        }
        if (knha) {
            stXWX <- .invcalc(X = X, W = W, k = k)
            b.knha <- stXWX %*% crossprod(X, W) %*% Y
            RSS.knha <- sum(wi * (yi - X %*% b.knha)^2)
            if (RSS.knha <= .Machine$double.eps) {
                s2w <- 1
            }
            else {
                s2w <- RSS.knha/(k - p)
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
    if (very.verbose) 
        message("Heterogeneity testing ...")
    if (allvipos) {
        wi <- 1/vi
        W.FE <- diag(wi, nrow = k, ncol = k)
        stXWX <- .invcalc(X = X, W = W.FE, k = k)
        P <- W.FE - W.FE %*% X %*% stXWX %*% crossprod(X, W.FE)
        QE <- max(0, c(crossprod(Y, P) %*% Y))
        QEp <- ifelse(k - p >= 1, pchisq(QE, df = k - p, lower.tail = FALSE), 
            1)
        vi.avg <- (k - p)/.tr(P)
        I2 <- 100 * tau2/(vi.avg + tau2)
        H2 <- tau2/vi.avg + 1
    }
    if (!int.only && int.incl && method != "FE") {
        if (very.verbose) {
            message("Fitting RE model for R^2 computation ...")
            res.RE <- try(rma(yi, vi, weights = weights, method = method, 
                weighted = weighted, knha = knha, verbose = ifelse(verbose, 
                  TRUE, FALSE), control = con, digits = digits), 
                silent = FALSE)
        }
        else {
            res.RE <- suppressWarnings(try(rma(yi, vi, weights = weights, 
                method = method, weighted = weighted, knha = knha, 
                verbose = ifelse(verbose, TRUE, FALSE), control = con, 
                digits = digits), silent = FALSE))
        }
        if (!inherits(res.RE, "try-error")) {
            tau2.RE <- res.RE$tau2
            if (identical(tau2.RE, 0)) {
                R2 <- NA
            }
            else {
                R2 <- round(max(0, 100 * (tau2.RE - tau2)/tau2.RE), 
                  2)
            }
        }
        else {
            R2 <- NA
        }
    }
    else {
        R2 <- NULL
    }
    if (very.verbose) 
        message("Computing fit statistics and log likelihood ...")
    parms <- p + ifelse(method == "FE" || tau2.fix, 0, 1)
    ll.ML <- -1/2 * (k) * log(2 * base::pi) - 1/2 * sum(log(vi + 
        tau2)) - 1/2 * RSS.f
    ll.REML <- -1/2 * (k - p) * log(2 * base::pi) + ifelse(con$REMLf, 
        1/2 * determinant(crossprod(X), logarithm = TRUE)$modulus, 
        0) - 1/2 * sum(log(vi + tau2)) - 1/2 * determinant(crossprod(X, 
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
    AICc.ML <- -2 * ll.ML + 2 * parms * max(k, parms + 2)/(max(k, 
        parms + 2) - parms - 1)
    dev.REML <- -2 * (ll.REML - 0)
    AIC.REML <- -2 * ll.REML + 2 * parms
    BIC.REML <- -2 * ll.REML + parms * log(k - p)
    AICc.REML <- -2 * ll.REML + 2 * parms * max(k - p, parms + 
        2)/(max(k - p, parms + 2) - parms - 1)
    fit.stats <- matrix(c(ll.ML, dev.ML, AIC.ML, BIC.ML, AICc.ML, 
        ll.REML, dev.REML, AIC.REML, BIC.REML, AICc.REML), ncol = 2, 
        byrow = FALSE)
    dimnames(fit.stats) <- list(c("ll", "dev", "AIC", "BIC", 
        "AICc"), c("ML", "REML"))
    fit.stats <- data.frame(fit.stats)
    if (very.verbose) 
        message("Preparing output ...")
    p.eff <- p
    k.eff <- k
    res <- list(b = b, se = se, zval = zval, pval = pval, ci.lb = ci.lb, 
        ci.ub = ci.ub, vb = vb, tau2 = tau2, se.tau2 = se.tau2, 
        k = k, k.f = k.f, k.eff = k.eff, p = p, p.eff = p.eff, 
        parms = parms, m = m, QE = QE, QEp = QEp, QM = QM, QMp = QMp, 
        I2 = I2, H2 = H2, R2 = R2, int.only = int.only, int.incl = int.incl, 
        allvipos = allvipos, yi = yi, vi = vi, X = X, weights = weights, 
        yi.f = yi.f, vi.f = vi.f, X.f = X.f, weights.f = weights.f, 
        ai.f = ai.f, bi.f = bi.f, ci.f = ci.f, di.f = di.f, x1i.f = x1i.f, 
        x2i.f = x2i.f, t1i.f = t1i.f, t2i.f = t2i.f, ni = ni, 
        ni.f = ni.f, ids = ids, not.na = not.na, slab = slab, 
        slab.null = slab.null, measure = measure, method = method, 
        weighted = weighted, knha = knha, robust = robust, s2w = s2w, 
        btt = btt, intercept = intercept, digits = digits, level = level, 
        control = control, verbose = verbose, add = add, to = to, 
        drop00 = drop00, fit.stats = fit.stats)
    class(res) <- c("rma.uni", "rma")
    return(res)
}
rstandard.rma.mh <-
function (model, digits, ...) 
{
    if (!is.element("rma.mh", class(model))) 
        stop("Argument 'model' must be an object of class \"rma.mh\".")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    x <- model
    if (missing(digits)) 
        digits <- x$digits
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
rstandard.rma.mv <-
function (model, digits, ...) 
{
    if (!is.element("rma.mv", class(model))) 
        stop("Argument 'model' must be an object of class \"rma.mv\".")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    x <- model
    if (missing(digits)) 
        digits <- x$digits
    M <- x$M
    options(na.action = "na.omit")
    H <- hatvalues(x, type = "matrix")
    options(na.action = na.act)
    ImH <- diag(x$k) - H
    ei <- ImH %*% cbind(x$yi)
    ei[abs(ei) < 100 * .Machine$double.eps] <- 0
    ve <- ImH %*% tcrossprod(M, ImH)
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
rstandard.rma.peto <-
function (model, digits, ...) 
{
    if (!is.element("rma.peto", class(model))) 
        stop("Argument 'model' must be an object of class \"rma.peto\".")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    x <- model
    if (missing(digits)) 
        digits <- x$digits
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
function (model, digits, ...) 
{
    if (!is.element("rma.uni", class(model))) 
        stop("Argument 'model' must be an object of class \"rma.uni\".")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    x <- model
    if (missing(digits)) 
        digits <- x$digits
    M <- diag(x$vi + x$tau2, nrow = x$k, ncol = x$k)
    options(na.action = "na.omit")
    H <- hatvalues(x, type = "matrix")
    options(na.action = na.act)
    ImH <- diag(x$k) - H
    ei <- ImH %*% cbind(x$yi)
    ei[abs(ei) < 100 * .Machine$double.eps] <- 0
    ve <- ImH %*% tcrossprod(M, ImH)
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
function (model, digits, ...) 
{
    if (!is.element("rma.mh", class(model))) 
        stop("Argument 'model' must be an object of class \"rma.mh\".")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    x <- model
    if (missing(digits)) 
        digits <- x$digits
    delpred <- rep(NA, x$k.f)
    vdelpred <- rep(NA, x$k.f)
    o.warn <- getOption("warn")
    on.exit(options(warn = o.warn))
    options(warn = -1)
    for (i in seq_len(x$k.f)[x$not.na]) {
        if (is.element(x$measure, c("RR", "OR", "RD"))) {
            res <- try(rma.mh(ai = x$ai.f[-i], bi = x$bi.f[-i], 
                ci = x$ci.f[-i], di = x$di.f[-i], measure = x$measure, 
                add = x$add, to = x$to, drop00 = x$drop00, correct = x$correct), 
                silent = TRUE)
        }
        else {
            res <- try(rma.mh(x1i = x$x1i.f[-i], x2i = x$x2i.f[-i], 
                t1i = x$t1i.f[-i], t2i = x$t2i.f[-i], measure = x$measure, 
                add = x$add, to = x$to, drop00 = x$drop00, correct = x$correct), 
                silent = TRUE)
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
function (model, digits, ...) 
{
    if (!is.element("rma.peto", class(model))) 
        stop("Argument 'model' must be an object of class \"rma.peto\".")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    x <- model
    if (missing(digits)) 
        digits <- x$digits
    delpred <- rep(NA, x$k.f)
    vdelpred <- rep(NA, x$k.f)
    o.warn <- getOption("warn")
    on.exit(options(warn = o.warn))
    options(warn = -1)
    for (i in seq_len(x$k.f)[x$not.na]) {
        res <- try(rma.peto(ai = x$ai.f[-i], bi = x$bi.f[-i], 
            ci = x$ci.f[-i], di = x$di.f[-i], add = x$add, to = x$to, 
            drop00 = x$drop00), silent = TRUE)
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
function (model, digits, ...) 
{
    if (!is.element("rma.uni", class(model))) 
        stop("Argument 'model' must be an object of class \"rma.uni\".")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    x <- model
    if (missing(digits)) 
        digits <- x$digits
    tau2.del <- rep(NA, x$k.f)
    delpred <- rep(NA, x$k.f)
    vdelpred <- rep(NA, x$k.f)
    o.warn <- getOption("warn")
    on.exit(options(warn = o.warn))
    options(warn = -1)
    for (i in seq_len(x$k.f)[x$not.na]) {
        res <- try(rma(x$yi.f[-i], x$vi.f[-i], weights = x$weights.f[-i], 
            mods = cbind(x$X.f[-i, ]), method = x$method, weighted = x$weighted, 
            intercept = FALSE, knha = x$knha, control = x$control), 
            silent = TRUE)
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
    var.names, H0 = 0, append = TRUE, replace = TRUE, level = 95, 
    digits, transf = FALSE, ...) 
{
    if (!is.element("escalc", class(object))) 
        stop("Argument 'object' must be an object of class \"escalc\".")
    x <- object
    alpha <- ifelse(level > 1, (100 - level)/100, 1 - level)
    crit <- qnorm(alpha/2, lower.tail = FALSE)
    if (length(out.names) != 4) 
        stop("Argument out.names must be of length 4.")
    if (any(out.names != make.names(out.names, unique = TRUE))) {
        out.names <- make.names(out.names, unique = TRUE)
        warning(paste0("Argument 'out.names' does not contain syntactically valid variable names.\n  Variable names adjusted to: out.names = c('", 
            out.names[1], "', '", out.names[2], "', '", out.names[3], 
            "', '", out.names[4], "')."))
    }
    if (missing(var.names)) {
        if (!is.null(attr(x, "yi.names"))) {
            yi.name <- attr(x, "yi.names")[1]
        }
        else {
            if (!is.element("yi", names(x))) 
                stop("Cannot determine name of the 'yi' variable.")
            yi.name <- "yi"
        }
        if (!is.null(attr(x, "vi.names"))) {
            vi.name <- attr(x, "vi.names")[1]
        }
        else {
            if (!is.element("vi", names(x))) 
                stop("Cannot determine name of the 'vi' variable.")
            vi.name <- "vi"
        }
    }
    else {
        if (length(var.names) != 2) 
            stop("Argument 'var.names' must be of length 2.")
        if (any(var.names != make.names(var.names, unique = TRUE))) {
            var.names <- make.names(var.names, unique = TRUE)
            warning(paste0("Argument 'var.names' does not contain syntactically valid variable names.\n  Variable names adjusted to: var.names = c('", 
                var.names[1], "', '", var.names[2], "')."))
        }
        yi.name <- var.names[1]
        vi.name <- var.names[2]
    }
    yi <- x[[yi.name]]
    vi <- x[[vi.name]]
    if (is.null(yi) || is.null(vi)) 
        stop(paste0("Cannot find variables '", yi.name, "' and/or '", 
            vi.name, "' in the data frame."))
    k <- length(yi)
    if (length(H0) == 1) 
        H0 <- rep(H0, k)
    sei <- sqrt(vi)
    zi <- (yi - H0)/sei
    if (is.function(transf)) {
        ci.lb <- mapply(transf, yi - crit * sei, ...)
        ci.ub <- mapply(transf, yi + crit * sei, ...)
        yi <- mapply(transf, yi, ...)
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
    if (!missing(digits)) {
        attr(dat, "digits") <- digits
    }
    else {
        attr(dat, "digits") <- attr(x, "digits")
    }
    if (is.null(attr(dat, "digits"))) 
        attr(dat, "digits") <- 4
    if (!missing(var.names)) {
        attr(dat, "yi.names") <- unique(c(var.names[1], attr(object, 
            "yi.names")))
    }
    else {
        attr(dat, "yi.names") <- unique(c(yi.name, attr(object, 
            "yi.names")))
    }
    if (!missing(var.names)) {
        attr(dat, "vi.names") <- unique(c(var.names[2], attr(object, 
            "vi.names")))
    }
    else {
        attr(dat, "vi.names") <- unique(c(vi.name, attr(object, 
            "vi.names")))
    }
    attr(dat, "sei.names") <- unique(c(out.names[1], attr(object, 
        "sei.names")))
    attr(dat, "zi.names") <- unique(c(out.names[2], attr(object, 
        "zi.names")))
    attr(dat, "ci.lb.names") <- unique(c(out.names[3], attr(object, 
        "ci.lb.names")))
    attr(dat, "ci.ub.names") <- unique(c(out.names[4], attr(object, 
        "ci.ub.names")))
    class(dat) <- c("escalc", "data.frame")
    return(dat)
}
summary.rma <-
function (object, digits, showfit = TRUE, ...) 
{
    if (!is.element("rma", class(object))) 
        stop("Argument 'object' must be an object of class \"rma\".")
    if (missing(digits)) 
        digits <- object$digits
    class(object) <- c("summary.rma", class(object))
    return(object)
}
tail.list.rma <-
function (x, n = 6L, ...) 
{
    stopifnot(length(n) == 1L)
    nrx <- length(x[[1]])
    n <- if (n < 0L) {
        max(nrx + n, 0L)
    }
    else {
        min(n, nrx)
    }
    x[seq.int(to = nrx, length.out = n), , drop = FALSE]
}
to.long <-
function (measure, ai, bi, ci, di, n1i, n2i, x1i, x2i, t1i, t2i, 
    m1i, m2i, sd1i, sd2i, xi, mi, ri, ti, sdi, ni, data, slab, 
    subset, add = 1/2, to = "none", drop00 = FALSE, vlong = FALSE, 
    append = TRUE, var.names) 
{
    if (!is.element(measure, c("RR", "OR", "PETO", "RD", "AS", 
        "PHI", "YUQ", "YUY", "RTET", "PBIT", "OR2D", "OR2DN", 
        "OR2DL", "IRR", "IRD", "IRSD", "MD", "SMD", "SMDH", "ROM", 
        "RPB", "RBIS", "D2OR", "D2ORN", "D2ORL", "COR", "UCOR", 
        "ZCOR", "PR", "PLN", "PLO", "PAS", "PFT", "IR", "IRLN", 
        "IRS", "IRFT", "MN", "MC", "SMCC", "SMCR", "ARAW", "AHW", 
        "ABT"))) 
        stop("Unknown 'measure' specified.")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    if (!is.element(to, c("all", "only0", "if0all", "none"))) 
        stop("Unknown 'to' argument specified.")
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
        "PHI", "YUQ", "YUY", "RTET", "PBIT", "OR2D", "OR2DN", 
        "OR2DL"))) {
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
        k <- length(ai)
        if (!is.null(subset)) {
            ai <- ai[subset]
            bi <- bi[subset]
            ci <- ci[subset]
            di <- di[subset]
        }
        if (length(ai) == 0L || length(bi) == 0L || length(ci) == 
            0L || length(di) == 0L) 
            stop("Cannot compute outcomes. Check that all of the required \n  information is specified via the appropriate arguments.")
        if (!all(length(ai) == c(length(ai), length(bi), length(ci), 
            length(di)))) 
            stop("Supplied data vectors are not all of the same length.")
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
        k <- length(x1i)
        if (!is.null(subset)) {
            x1i <- x1i[subset]
            x2i <- x2i[subset]
            t1i <- t1i[subset]
            t2i <- t2i[subset]
        }
        if (length(x1i) == 0L || length(x2i) == 0L || length(t1i) == 
            0L || length(t2i) == 0L) 
            stop("Cannot compute outcomes. Check that all of the required \n  information is specified via the appropriate arguments.")
        if (!all(length(x1i) == c(length(x1i), length(x2i), length(t1i), 
            length(t2i)))) 
            stop("Supplied data vectors are not all of the same length.")
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
    }
    if (is.element(measure, c("MD", "SMD", "SMDH", "ROM", "RPB", 
        "RBIS", "D2OR", "D2ORN", "D2ORL"))) {
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
        k <- length(m1i)
        if (!is.null(subset)) {
            m1i <- m1i[subset]
            m2i <- m2i[subset]
            sd1i <- sd1i[subset]
            sd2i <- sd2i[subset]
            n1i <- n1i[subset]
            n2i <- n2i[subset]
        }
        if (length(m1i) == 0L || length(m2i) == 0L || length(sd1i) == 
            0L || length(sd2i) == 0L || length(n1i) == 0L || 
            length(n2i) == 0L) 
            stop("Cannot compute outcomes. Check that all of the required \n  information is specified via the appropriate arguments.")
        if (!all(length(m1i) == c(length(m1i), length(m2i), length(sd1i), 
            length(sd2i), length(n1i), length(n2i)))) 
            stop("Supplied data vectors are not all of the same length.")
        if (any(c(sd1i, sd2i) < 0, na.rm = TRUE)) 
            stop("One or more standard deviations are negative.")
        if (any(c(n1i, n2i) < 0, na.rm = TRUE)) 
            stop("One or more sample sizes are negative.")
        ni.u <- n1i + n2i
    }
    if (is.element(measure, c("COR", "UCOR", "ZCOR"))) {
        mf.ri <- mf[[match("ri", names(mf))]]
        mf.ni <- mf[[match("ni", names(mf))]]
        ri <- eval(mf.ri, data, enclos = sys.frame(sys.parent()))
        ni <- eval(mf.ni, data, enclos = sys.frame(sys.parent()))
        k <- length(ri)
        if (!is.null(subset)) {
            ri <- ri[subset]
            ni <- ni[subset]
        }
        if (length(ri) == 0L || length(ni) == 0L) 
            stop("Cannot compute outcomes. Check that all of the required \n  information is specified via the appropriate arguments.")
        if (length(ri) != length(ni)) 
            stop("Supplied data vectors are not of the same length.")
        if (any(abs(ri) > 1, na.rm = TRUE)) 
            stop("One or more correlations are > 1 or < -1.")
        if (any(ni < 0, na.rm = TRUE)) 
            stop("One or more sample sizes are negative.")
        ni.u <- ni
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
        k <- length(xi)
        if (!is.null(subset)) {
            xi <- xi[subset]
            mi <- mi[subset]
        }
        if (length(xi) == 0L || length(mi) == 0L) 
            stop("Cannot compute outcomes. Check that all of the required \n  information is specified via the appropriate arguments.")
        if (length(xi) != length(mi)) 
            stop("Supplied data vectors are not all of the same length.")
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
    }
    if (is.element(measure, c("IR", "IRLN", "IRS", "IRFT"))) {
        mf.xi <- mf[[match("xi", names(mf))]]
        mf.ti <- mf[[match("ti", names(mf))]]
        xi <- eval(mf.xi, data, enclos = sys.frame(sys.parent()))
        ti <- eval(mf.ti, data, enclos = sys.frame(sys.parent()))
        k <- length(xi)
        if (!is.null(subset)) {
            xi <- xi[subset]
            ti <- ti[subset]
        }
        if (length(xi) == 0L || length(ti) == 0L) 
            stop("Cannot compute outcomes. Check that all of the required \n  information is specified via the appropriate arguments.")
        if (length(xi) != length(ti)) 
            stop("Supplied data vectors are not all of the same length.")
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
    }
    if (is.element(measure, c("MN"))) {
        mf.mi <- mf[[match("mi", names(mf))]]
        mf.sdi <- mf[[match("sdi", names(mf))]]
        mf.ni <- mf[[match("ni", names(mf))]]
        mi <- eval(mf.mi, data, enclos = sys.frame(sys.parent()))
        sdi <- eval(mf.sdi, data, enclos = sys.frame(sys.parent()))
        ni <- eval(mf.ni, data, enclos = sys.frame(sys.parent()))
        k <- length(mi)
        if (!is.null(subset)) {
            mi <- mi[subset]
            sdi <- sdi[subset]
            ni <- ni[subset]
        }
        if (length(mi) == 0L || length(sdi) == 0L || length(ni) == 
            0L) 
            stop("Cannot compute outcomes. Check that all of the required \n  information is specified via the appropriate arguments.")
        if (!all(length(mi) == c(length(mi), length(sdi), length(ni)))) 
            stop("Supplied data vectors are not all of the same length.")
        if (any(sdi < 0, na.rm = TRUE)) 
            stop("One or more standard deviations are negative.")
        if (any(ni < 0, na.rm = TRUE)) 
            stop("One or more sample sizes are negative.")
        ni.u <- ni
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
        k <- length(m1i)
        if (!is.null(subset)) {
            m1i <- m1i[subset]
            m2i <- m2i[subset]
            sd1i <- sd1i[subset]
            sd2i <- sd2i[subset]
            ni <- ni[subset]
            ri <- ri[subset]
        }
        if (is.element(measure, c("MC", "SMCC"))) {
            if (length(m1i) == 0L || length(m2i) == 0L || length(sd1i) == 
                0L || length(sd2i) == 0L || length(ni) == 0L || 
                length(ri) == 0L) 
                stop("Cannot compute outcomes. Check that all of the required \n  information is specified via the appropriate arguments.")
            if (!all(length(m1i) == c(length(m1i), length(m2i), 
                length(sd1i), length(sd2i), length(ni), length(ri)))) 
                stop("Supplied data vectors are not all of the same length.")
            if (any(c(sd1i, sd2i) < 0, na.rm = TRUE)) 
                stop("One or more standard deviations are negative.")
        }
        else {
            if (length(m1i) == 0L || length(m2i) == 0L || length(sd1i) == 
                0L || length(ni) == 0L || length(ri) == 0L) 
                stop("Cannot compute outcomes. Check that all of the required \n  information is specified via the appropriate arguments.")
            if (!all(length(m1i) == c(length(m1i), length(m2i), 
                length(sd1i), length(ni), length(ri)))) 
                stop("Supplied data vectors are not all of the same length.")
            if (any(sd1i < 0, na.rm = TRUE)) 
                stop("One or more standard deviations are negative.")
        }
        if (any(abs(ri) > 1, na.rm = TRUE)) 
            stop("One or more correlations are > 1 or < -1.")
        if (any(ni < 0, na.rm = TRUE)) 
            stop("One or more sample sizes are negative.")
        ni.u <- ni
    }
    if (is.element(measure, c("ARAW", "AHW", "ABT"))) {
        mf.ai <- mf[[match("ai", names(mf))]]
        mf.mi <- mf[[match("mi", names(mf))]]
        mf.ni <- mf[[match("ni", names(mf))]]
        ai <- eval(mf.ai, data, enclos = sys.frame(sys.parent()))
        mi <- eval(mf.mi, data, enclos = sys.frame(sys.parent()))
        ni <- eval(mf.ni, data, enclos = sys.frame(sys.parent()))
        k <- length(ai)
        if (!is.null(subset)) {
            ai <- ai[subset]
            mi <- mi[subset]
            ni <- ni[subset]
        }
        if (length(ai) == 0L || length(mi) == 0L || length(ni) == 
            0L) 
            stop("Cannot compute outcomes. Check that all of the required \n  information is specified via the appropriate arguments.")
        if (!all(length(ai) == c(length(ai), length(mi), length(ni)))) 
            stop("Supplied data vectors are not all of the same length.")
        if (any(ai > 1, na.rm = TRUE)) 
            stop("One or more alphas are > 1.")
        if (any(mi < 2, na.rm = TRUE)) 
            stop("One or more mi's are < 2.")
        if (any(ni < 0, na.rm = TRUE)) 
            stop("One or more sample sizes are negative.")
        ni.u <- ni
    }
    if (is.null(slab)) {
        slab <- seq_len(k)
    }
    else {
        if (any(is.na(slab))) 
            stop("NAs in study labels.")
        if (any(duplicated(slab))) 
            slab <- make.unique(slab)
        if (length(slab) != k) 
            stop("Study labels not of same length as data.")
    }
    if (!is.null(subset)) {
        slab <- slab[subset]
        if (!no.data) 
            data <- data[subset, ]
    }
    if (is.element(measure, c("RR", "OR", "RD", "AS", "PETO", 
        "PHI", "YUQ", "YUY", "RTET", "PBIT", "OR2D", "OR2DN", 
        "OR2DL"))) {
        aibicidi.na <- is.na(cbind(ai, bi, ci, di))
        if (any(aibicidi.na)) {
            not.na <- apply(aibicidi.na, MARGIN = 1, sum) == 
                0L
            if (na.act == "na.omit") {
                ai <- ai[not.na]
                bi <- bi[not.na]
                ci <- ci[not.na]
                di <- di[not.na]
                slab <- slab[not.na]
                if (!no.data) 
                  data <- data[not.na, ]
                warning("Tables with NAs omitted.")
            }
            if (na.act == "na.fail") 
                stop("Missing values in tables.")
        }
        k <- length(ai)
        if (k < 1) 
            stop("Processing terminated since k = 0.")
        if (vlong) {
            dat <- matrix(NA, nrow = 4 * k, ncol = 4)
            dat[, 1] <- rep(slab, each = 4)
            dat[, 2] <- rep(c(1, 1, 2, 2), k)
            dat[, 3] <- rep(c(1, 2, 1, 2), k)
            dat[, 4] <- c(rbind(ai, bi, ci, di))
            if (missing(var.names)) {
                colnames(dat) <- c("study", "group", "outcome", 
                  "freq")
            }
            else {
                if (length(var.names) != 4) 
                  stop("Variable names not of length 4.")
                colnames(dat) <- var.names
            }
            dat <- data.frame(dat)
            dat[, 1] <- factor(dat[, 1])
            dat[, 2] <- factor(dat[, 2])
            dat[, 3] <- factor(dat[, 3])
            if (!no.data && append) 
                dat <- data.frame(data[rep(1:k, each = 4), ], 
                  dat)
        }
        else {
            dat <- matrix(NA, nrow = 2 * k, ncol = 4)
            dat[, 1] <- rep(slab, each = 2)
            dat[, 2] <- rep(c(1, 2), k)
            dat[, 3] <- c(rbind(ai, ci))
            dat[, 4] <- c(rbind(bi, di))
            if (missing(var.names)) {
                colnames(dat) <- c("study", "group", "out1", 
                  "out2")
            }
            else {
                if (length(var.names) != 4) 
                  stop("Variable names not of length 4.")
                colnames(dat) <- var.names
            }
            dat <- data.frame(dat)
            dat[, 1] <- factor(dat[, 1])
            dat[, 2] <- factor(dat[, 2])
            if (!no.data && append) 
                dat <- data.frame(data[rep(1:k, each = 2), ], 
                  dat)
        }
    }
    if (is.element(measure, c("IRR", "IRD", "IRSD"))) {
        x1ix2it1it2i.na <- is.na(cbind(x1i, x2i, t1i, t2i))
        if (any(x1ix2it1it2i.na)) {
            not.na <- apply(x1ix2it1it2i.na, MARGIN = 1, sum) == 
                0L
            if (na.act == "na.omit") {
                x1i <- x1i[not.na]
                x2i <- x2i[not.na]
                t1i <- t1i[not.na]
                t2i <- t2i[not.na]
                slab <- slab[not.na]
                if (!no.data) 
                  data <- data[not.na, ]
                warning("Tables with NAs omitted.")
            }
            if (na.act == "na.fail") 
                stop("Missing values in tables.")
        }
        k <- length(x1i)
        if (k < 1) 
            stop("Processing terminated since k = 0.")
        dat <- matrix(NA, nrow = 2 * k, ncol = 4)
        dat[, 1] <- rep(slab, each = 2)
        dat[, 2] <- rep(c(1, 2), k)
        dat[, 3] <- c(rbind(x1i, x2i))
        dat[, 4] <- c(rbind(t1i, t2i))
        if (missing(var.names)) {
            colnames(dat) <- c("study", "group", "events", "ptime")
        }
        else {
            if (length(var.names) != 4) 
                stop("Variable names not of length 4.")
            colnames(dat) <- var.names
        }
        dat <- data.frame(dat)
        dat[, 1] <- factor(dat[, 1])
        dat[, 2] <- factor(dat[, 2])
        if (!no.data && append) 
            dat <- data.frame(data[rep(1:k, each = 2), ], dat)
    }
    if (is.element(measure, c("MD", "SMD", "SMDH", "ROM", "RPB", 
        "RBIS", "D2OR", "D2ORN", "D2ORL"))) {
        m1im2isd1isd2in1in2i.na <- is.na(cbind(m1i, m2i, sd1i, 
            sd2i, n1i, n2i))
        if (any(m1im2isd1isd2in1in2i.na)) {
            not.na <- apply(m1im2isd1isd2in1in2i.na, MARGIN = 1, 
                sum) == 0L
            if (na.act == "na.omit") {
                m1i <- m1i[not.na]
                m2i <- m2i[not.na]
                sd1i <- sd1i[not.na]
                sd2i <- sd2i[not.na]
                n1i <- n1i[not.na]
                n2i <- n2i[not.na]
                slab <- slab[not.na]
                if (!no.data) 
                  data <- data[not.na, ]
                warning("Tables with NAs omitted.")
            }
            if (na.act == "na.fail") 
                stop("Missing values in tables.")
        }
        k <- length(m1i)
        if (k < 1) 
            stop("Processing terminated since k = 0.")
        dat <- matrix(NA, nrow = 2 * k, ncol = 5)
        dat[, 1] <- rep(slab, each = 2)
        dat[, 2] <- rep(c(1, 2), k)
        dat[, 3] <- c(rbind(m1i, m2i))
        dat[, 4] <- c(rbind(sd1i, sd2i))
        dat[, 5] <- c(rbind(n1i, n2i))
        if (missing(var.names)) {
            colnames(dat) <- c("study", "group", "mean", "sd", 
                "n")
        }
        else {
            if (length(var.names) != 5) 
                stop("Variable names not of length 5.")
            colnames(dat) <- var.names
        }
        dat <- data.frame(dat)
        dat[, 1] <- factor(dat[, 1])
        dat[, 2] <- factor(dat[, 2])
        if (!no.data && append) 
            dat <- data.frame(data[rep(1:k, each = 2), ], dat)
    }
    if (is.element(measure, c("COR", "UCOR", "ZCOR"))) {
        rini.na <- is.na(cbind(ri, ni))
        if (any(rini.na)) {
            not.na <- apply(rini.na, MARGIN = 1, sum) == 0L
            if (na.act == "na.omit") {
                ri <- ri[not.na]
                ni <- ni[not.na]
                slab <- slab[not.na]
                if (!no.data) 
                  data <- data[not.na, ]
                warning("Tables with NAs omitted.")
            }
            if (na.act == "na.fail") 
                stop("Missing values in tables.")
        }
        k <- length(ri)
        if (k < 1) 
            stop("Processing terminated since k = 0.")
        dat <- matrix(NA, nrow = k, ncol = 3)
        dat[, 1] <- slab
        dat[, 2] <- ri
        dat[, 3] <- ni
        if (missing(var.names)) {
            colnames(dat) <- c("study", "r", "n")
        }
        else {
            if (length(var.names) != 3) 
                stop("Variable names not of length 3.")
            colnames(dat) <- var.names
        }
        dat <- data.frame(dat)
        dat[, 1] <- factor(dat[, 1])
        if (!no.data && append) 
            dat <- data.frame(data, dat)
    }
    if (is.element(measure, c("PR", "PLN", "PLO", "PAS", "PFT"))) {
        ximi.na <- is.na(cbind(xi, mi))
        if (any(ximi.na)) {
            not.na <- apply(ximi.na, MARGIN = 1, sum) == 0L
            if (na.act == "na.omit") {
                xi <- xi[not.na]
                mi <- mi[not.na]
                slab <- slab[not.na]
                if (!no.data) 
                  data <- data[not.na, ]
                warning("Tables with NAs omitted.")
            }
            if (na.act == "na.fail") 
                stop("Missing values in tables.")
        }
        k <- length(xi)
        if (k < 1) 
            stop("Processing terminated since k = 0.")
        if (vlong) {
            dat <- matrix(NA, nrow = 2 * k, ncol = 3)
            dat[, 1] <- rep(slab, each = 2)
            dat[, 2] <- rep(c(1, 2), k)
            dat[, 3] <- c(rbind(xi, mi))
            if (missing(var.names)) {
                colnames(dat) <- c("study", "outcome", "freq")
            }
            else {
                if (length(var.names) != 3) 
                  stop("Variable names not of length 3.")
                colnames(dat) <- var.names
            }
            dat <- data.frame(dat)
            dat[, 1] <- factor(dat[, 1])
            dat[, 2] <- factor(dat[, 2])
            if (!no.data && append) 
                dat <- data.frame(data[rep(1:k, each = 2), ], 
                  dat)
        }
        else {
            dat <- matrix(NA, nrow = k, ncol = 3)
            dat[, 1] <- slab
            dat[, 2] <- xi
            dat[, 3] <- mi
            if (missing(var.names)) {
                colnames(dat) <- c("study", "out1", "out2")
            }
            else {
                if (length(var.names) != 3) 
                  stop("Variable names not of length 3.")
                colnames(dat) <- var.names
            }
            dat <- data.frame(dat)
            if (!no.data && append) 
                dat <- data.frame(data, dat)
        }
    }
    if (is.element(measure, c("IR", "IRLN", "IRS", "IRFT"))) {
        xiti.na <- is.na(cbind(xi, ti))
        if (any(xiti.na)) {
            not.na <- apply(xiti.na, MARGIN = 1, sum) == 0L
            if (na.act == "na.omit") {
                xi <- xi[not.na]
                ti <- ti[not.na]
                slab <- slab[not.na]
                if (!no.data) 
                  data <- data[not.na, ]
                warning("Tables with NAs omitted.")
            }
            if (na.act == "na.fail") 
                stop("Missing values in tables.")
        }
        k <- length(xi)
        if (k < 1) 
            stop("Processing terminated since k = 0.")
        dat <- matrix(NA, nrow = k, ncol = 3)
        dat[, 1] <- slab
        dat[, 2] <- xi
        dat[, 3] <- ti
        if (missing(var.names)) {
            colnames(dat) <- c("study", "events", "ptime")
        }
        else {
            if (length(var.names) != 3) 
                stop("Variable names not of length 3.")
            colnames(dat) <- var.names
        }
        dat <- data.frame(dat)
        dat[, 1] <- factor(dat[, 1])
        if (!no.data && append) 
            dat <- data.frame(data, dat)
    }
    if (is.element(measure, c("MN"))) {
        misdini.na <- is.na(cbind(mi, sdi, ni))
        if (any(misdini.na)) {
            not.na <- apply(misdini.na, MARGIN = 1, sum) == 0L
            if (na.act == "na.omit") {
                mi <- mi[not.na]
                sdi <- sdi[not.na]
                ni <- ni[not.na]
                slab <- slab[not.na]
                if (!no.data) 
                  data <- data[not.na, ]
                warning("Tables with NAs omitted.")
            }
            if (na.act == "na.fail") 
                stop("Missing values in tables.")
        }
        k <- length(mi)
        if (k < 1) 
            stop("Processing terminated since k = 0.")
        dat <- matrix(NA, nrow = k, ncol = 4)
        dat[, 1] <- slab
        dat[, 2] <- mi
        dat[, 3] <- sdi
        dat[, 4] <- ni
        if (missing(var.names)) {
            colnames(dat) <- c("study", "mean", "sd", "n")
        }
        else {
            if (length(var.names) != 4) 
                stop("Variable names not of length 4.")
            colnames(dat) <- var.names
        }
        dat <- data.frame(dat)
        dat[, 1] <- factor(dat[, 1])
        if (!no.data && append) 
            dat <- data.frame(data, dat)
    }
    if (is.element(measure, c("MC", "SMCC", "SMCR"))) {
        if (is.element(measure, c("MC", "SMCC"))) {
            m1im2isdiniri.na <- is.na(cbind(m1i, m2i, sd1i, sd2i, 
                ni, ri))
        }
        else {
            m1im2isdiniri.na <- is.na(cbind(m1i, m2i, sd1i, ni, 
                ri))
        }
        if (any(m1im2isdiniri.na)) {
            not.na <- apply(m1im2isdiniri.na, MARGIN = 1, sum) == 
                0L
            if (na.act == "na.omit") {
                m1i <- m1i[not.na]
                m2i <- m2i[not.na]
                sd1i <- sd1i[not.na]
                if (is.element(measure, c("MC", "SMCC"))) 
                  sd2i <- sd2i[not.na]
                ni <- ni[not.na]
                ri <- ri[not.na]
                slab <- slab[not.na]
                if (!no.data) 
                  data <- data[not.na, ]
                warning("Tables with NAs omitted.")
            }
            if (na.act == "na.fail") 
                stop("Missing values in tables.")
        }
        k <- length(m1i)
        if (k < 1) 
            stop("Processing terminated since k = 0.")
        if (is.element(measure, c("MC", "SMCC"))) {
            dat <- matrix(NA, nrow = k, ncol = 7)
            dat[, 1] <- slab
            dat[, 2] <- m1i
            dat[, 3] <- m2i
            dat[, 4] <- sd1i
            dat[, 5] <- sd2i
            dat[, 6] <- ni
            dat[, 7] <- ri
            if (missing(var.names)) {
                colnames(dat) <- c("study", "mean1", "mean2", 
                  "sd1", "sd2", "n", "r")
            }
            else {
                if (length(var.names) != 7) 
                  stop("Variable names not of length 7.")
                colnames(dat) <- var.names
            }
            dat <- data.frame(dat)
            dat[, 1] <- factor(dat[, 1])
            if (!no.data && append) 
                dat <- data.frame(data, dat)
        }
        else {
            dat <- matrix(NA, nrow = k, ncol = 6)
            dat[, 1] <- slab
            dat[, 2] <- m1i
            dat[, 3] <- m2i
            dat[, 4] <- sd1i
            dat[, 5] <- ni
            dat[, 6] <- ri
            if (missing(var.names)) {
                colnames(dat) <- c("study", "mean1", "mean2", 
                  "sd1", "n", "r")
            }
            else {
                if (length(var.names) != 6) 
                  stop("Variable names not of length 6.")
                colnames(dat) <- var.names
            }
            dat <- data.frame(dat)
            dat[, 1] <- factor(dat[, 1])
            if (!no.data && append) 
                dat <- data.frame(data, dat)
        }
    }
    if (is.element(measure, c("ARAW", "AHW", "ABT"))) {
        aimini.na <- is.na(cbind(ai, mi, ni))
        if (any(aimini.na)) {
            not.na <- apply(aimini.na, MARGIN = 1, sum) == 0L
            if (na.act == "na.omit") {
                ai <- ai[not.na]
                mi <- mi[not.na]
                ni <- ni[not.na]
                slab <- slab[not.na]
                if (!no.data) 
                  data <- data[not.na, ]
                warning("Tables with NAs omitted.")
            }
            if (na.act == "na.fail") 
                stop("Missing values in tables.")
        }
        k <- length(ai)
        if (k < 1) 
            stop("Processing terminated since k = 0.")
        dat <- matrix(NA, nrow = k, ncol = 4)
        dat[, 1] <- slab
        dat[, 2] <- ai
        dat[, 3] <- mi
        dat[, 4] <- ni
        if (missing(var.names)) {
            colnames(dat) <- c("study", "alpha", "m", "n")
        }
        else {
            if (length(var.names) != 4) 
                stop("Variable names not of length 4.")
            colnames(dat) <- var.names
        }
        dat <- data.frame(dat)
        dat[, 1] <- factor(dat[, 1])
        if (!no.data && append) 
            dat <- data.frame(data, dat)
    }
    rownames(dat) <- 1:nrow(dat)
    return(dat)
}
to.table <-
function (measure, ai, bi, ci, di, n1i, n2i, x1i, x2i, t1i, t2i, 
    m1i, m2i, sd1i, sd2i, xi, mi, ri, ti, sdi, ni, data, slab, 
    subset, add = 1/2, to = "none", drop00 = FALSE, rows, cols) 
{
    if (!is.element(measure, c("RR", "OR", "PETO", "RD", "AS", 
        "PHI", "YUQ", "YUY", "RTET", "PBIT", "OR2D", "OR2DN", 
        "OR2DL", "IRR", "IRD", "IRSD", "MD", "SMD", "SMDH", "ROM", 
        "RPB", "RBIS", "D2OR", "D2ORN", "D2ORL", "COR", "UCOR", 
        "ZCOR", "PR", "PLN", "PLO", "PAS", "PFT", "IR", "IRLN", 
        "IRS", "IRFT", "MN", "MC", "SMCC", "SMCR", "ARAW", "AHW", 
        "ABT"))) 
        stop("Unknown 'measure' specified.")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    if (!is.element(to, c("all", "only0", "if0all", "none"))) 
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
    if (is.element(measure, c("RR", "OR", "RD", "AS", "PETO", 
        "PHI", "YUQ", "YUY", "RTET", "PBIT", "OR2D", "OR2DN", 
        "OR2DL"))) {
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
        k <- length(ai)
        if (!is.null(subset)) {
            ai <- ai[subset]
            bi <- bi[subset]
            ci <- ci[subset]
            di <- di[subset]
        }
        if (length(ai) == 0L || length(bi) == 0L || length(ci) == 
            0L || length(di) == 0L) 
            stop("Cannot compute outcomes. Check that all of the required \n  information is specified via the appropriate arguments.")
        if (!all(length(ai) == c(length(ai), length(bi), length(ci), 
            length(di)))) 
            stop("Supplied data vectors are not all of the same length.")
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
        k <- length(x1i)
        if (!is.null(subset)) {
            x1i <- x1i[subset]
            x2i <- x2i[subset]
            t1i <- t1i[subset]
            t2i <- t2i[subset]
        }
        if (length(x1i) == 0L || length(x2i) == 0L || length(t1i) == 
            0L || length(t2i) == 0L) 
            stop("Cannot compute outcomes. Check that all of the required \n  information is specified via the appropriate arguments.")
        if (!all(length(x1i) == c(length(x1i), length(x2i), length(t1i), 
            length(t2i)))) 
            stop("Supplied data vectors are not all of the same length.")
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
    }
    if (is.element(measure, c("MD", "SMD", "SMDH", "ROM", "RPB", 
        "RBIS", "D2OR", "D2ORN", "D2ORL"))) {
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
        k <- length(m1i)
        if (!is.null(subset)) {
            m1i <- m1i[subset]
            m2i <- m2i[subset]
            sd1i <- sd1i[subset]
            sd2i <- sd2i[subset]
            n1i <- n1i[subset]
            n2i <- n2i[subset]
        }
        if (length(m1i) == 0L || length(m2i) == 0L || length(sd1i) == 
            0L || length(sd2i) == 0L || length(n1i) == 0L || 
            length(n2i) == 0L) 
            stop("Cannot compute outcomes. Check that all of the required \n  information is specified via the appropriate arguments.")
        if (!all(length(m1i) == c(length(m1i), length(m2i), length(sd1i), 
            length(sd2i), length(n1i), length(n2i)))) 
            stop("Supplied data vectors are not all of the same length.")
        if (any(c(sd1i, sd2i) < 0, na.rm = TRUE)) 
            stop("One or more standard deviations are negative.")
        if (any(c(n1i, n2i) < 0, na.rm = TRUE)) 
            stop("One or more sample sizes are negative.")
        ni.u <- n1i + n2i
    }
    if (is.element(measure, c("COR", "UCOR", "ZCOR"))) {
        mf.ri <- mf[[match("ri", names(mf))]]
        mf.ni <- mf[[match("ni", names(mf))]]
        ri <- eval(mf.ri, data, enclos = sys.frame(sys.parent()))
        ni <- eval(mf.ni, data, enclos = sys.frame(sys.parent()))
        k <- length(ri)
        if (!is.null(subset)) {
            ri <- ri[subset]
            ni <- ni[subset]
        }
        if (length(ri) == 0L || length(ni) == 0L) 
            stop("Cannot compute outcomes. Check that all of the required \n  information is specified via the appropriate arguments.")
        if (length(ri) != length(ni)) 
            stop("Supplied data vectors are not of the same length.")
        if (any(abs(ri) > 1, na.rm = TRUE)) 
            stop("One or more correlations are > 1 or < -1.")
        if (any(ni < 0, na.rm = TRUE)) 
            stop("One or more sample sizes are negative.")
        ni.u <- ni
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
        k <- length(xi)
        if (!is.null(subset)) {
            xi <- xi[subset]
            mi <- mi[subset]
        }
        if (length(xi) == 0L || length(mi) == 0L) 
            stop("Cannot compute outcomes. Check that all of the required \n  information is specified via the appropriate arguments.")
        if (length(xi) != length(mi)) 
            stop("Supplied data vectors are not all of the same length.")
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
    }
    if (is.element(measure, c("IR", "IRLN", "IRS", "IRFT"))) {
        mf.xi <- mf[[match("xi", names(mf))]]
        mf.ti <- mf[[match("ti", names(mf))]]
        xi <- eval(mf.xi, data, enclos = sys.frame(sys.parent()))
        ti <- eval(mf.ti, data, enclos = sys.frame(sys.parent()))
        k <- length(xi)
        if (!is.null(subset)) {
            xi <- xi[subset]
            ti <- ti[subset]
        }
        if (length(xi) == 0L || length(ti) == 0L) 
            stop("Cannot compute outcomes. Check that all of the required \n  information is specified via the appropriate arguments.")
        if (length(xi) != length(ti)) 
            stop("Supplied data vectors are not all of the same length.")
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
    }
    if (is.element(measure, c("MN"))) {
        mf.mi <- mf[[match("mi", names(mf))]]
        mf.sdi <- mf[[match("sdi", names(mf))]]
        mf.ni <- mf[[match("ni", names(mf))]]
        mi <- eval(mf.mi, data, enclos = sys.frame(sys.parent()))
        sdi <- eval(mf.sdi, data, enclos = sys.frame(sys.parent()))
        ni <- eval(mf.ni, data, enclos = sys.frame(sys.parent()))
        k <- length(mi)
        if (!is.null(subset)) {
            mi <- mi[subset]
            sdi <- sdi[subset]
            ni <- ni[subset]
        }
        if (length(mi) == 0L || length(sdi) == 0L || length(ni) == 
            0L) 
            stop("Cannot compute outcomes. Check that all of the required \n  information is specified via the appropriate arguments.")
        if (!all(length(mi) == c(length(mi), length(sdi), length(ni)))) 
            stop("Supplied data vectors are not all of the same length.")
        if (any(sdi < 0, na.rm = TRUE)) 
            stop("One or more standard deviations are negative.")
        if (any(ni < 0, na.rm = TRUE)) 
            stop("One or more sample sizes are negative.")
        ni.u <- ni
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
        k <- length(m1i)
        if (!is.null(subset)) {
            m1i <- m1i[subset]
            m2i <- m2i[subset]
            sd1i <- sd1i[subset]
            sd2i <- sd2i[subset]
            ni <- ni[subset]
            ri <- ri[subset]
        }
        if (is.element(measure, c("MC", "SMCC"))) {
            if (length(m1i) == 0L || length(m2i) == 0L || length(sd1i) == 
                0L || length(sd2i) == 0L || length(ni) == 0L || 
                length(ri) == 0L) 
                stop("Cannot compute outcomes. Check that all of the required \n  information is specified via the appropriate arguments.")
            if (!all(length(m1i) == c(length(m1i), length(m2i), 
                length(sd1i), length(sd2i), length(ni), length(ri)))) 
                stop("Supplied data vectors are not all of the same length.")
            if (any(c(sd1i, sd2i) < 0, na.rm = TRUE)) 
                stop("One or more standard deviations are negative.")
        }
        else {
            if (length(m1i) == 0L || length(m2i) == 0L || length(sd1i) == 
                0L || length(ni) == 0L || length(ri) == 0L) 
                stop("Cannot compute outcomes. Check that all of the required \n  information is specified via the appropriate arguments.")
            if (!all(length(m1i) == c(length(m1i), length(m2i), 
                length(sd1i), length(ni), length(ri)))) 
                stop("Supplied data vectors are not all of the same length.")
            if (any(sd1i < 0, na.rm = TRUE)) 
                stop("One or more standard deviations are negative.")
        }
        if (any(abs(ri) > 1, na.rm = TRUE)) 
            stop("One or more correlations are > 1 or < -1.")
        if (any(ni < 0, na.rm = TRUE)) 
            stop("One or more sample sizes are negative.")
        ni.u <- ni
    }
    if (is.element(measure, c("ARAW", "AHW", "ABT"))) {
        mf.ai <- mf[[match("ai", names(mf))]]
        mf.mi <- mf[[match("mi", names(mf))]]
        mf.ni <- mf[[match("ni", names(mf))]]
        ai <- eval(mf.ai, data, enclos = sys.frame(sys.parent()))
        mi <- eval(mf.mi, data, enclos = sys.frame(sys.parent()))
        ni <- eval(mf.ni, data, enclos = sys.frame(sys.parent()))
        k <- length(ai)
        if (!is.null(subset)) {
            ai <- ai[subset]
            mi <- mi[subset]
            ni <- ni[subset]
        }
        if (length(ai) == 0L || length(mi) == 0L || length(ni) == 
            0L) 
            stop("Cannot compute outcomes. Check that all of the required \n  information is specified via the appropriate arguments.")
        if (!all(length(ai) == c(length(ai), length(mi), length(ni)))) 
            stop("Supplied data vectors are not all of the same length.")
        if (any(ai > 1, na.rm = TRUE)) 
            stop("One or more alphas are > 1.")
        if (any(mi < 2, na.rm = TRUE)) 
            stop("One or more mi's are < 2.")
        if (any(ni < 0, na.rm = TRUE)) 
            stop("One or more sample sizes are negative.")
        ni.u <- ni
    }
    if (is.null(slab)) {
        slab <- seq_len(k)
    }
    else {
        if (any(is.na(slab))) 
            stop("NAs in study labels.")
        if (any(duplicated(slab))) 
            slab <- make.unique(slab)
        if (length(slab) != k) 
            stop("Study labels not of same length as data.")
    }
    if (!is.null(subset)) 
        slab <- slab[subset]
    if (is.element(measure, c("RR", "OR", "RD", "AS", "PETO", 
        "PHI", "YUQ", "YUY", "RTET", "PBIT", "OR2D", "OR2DN", 
        "OR2DL"))) {
        aibicidi.na <- is.na(cbind(ai, bi, ci, di))
        if (any(aibicidi.na)) {
            not.na <- apply(aibicidi.na, MARGIN = 1, sum) == 
                0L
            if (na.act == "na.omit") {
                ai <- ai[not.na]
                bi <- bi[not.na]
                ci <- ci[not.na]
                di <- di[not.na]
                slab <- slab[not.na]
                warning("Tables with NAs omitted.")
            }
            if (na.act == "na.fail") 
                stop("Missing values in tables.")
        }
        k <- length(ai)
        if (k < 1) 
            stop("Processing terminated since k = 0.")
        if (missing(rows)) {
            rows <- c("Grp1", "Grp2")
        }
        else {
            if (length(rows) != 2) 
                stop("Group names not of length 2.")
        }
        if (missing(cols)) {
            cols <- c("Out1", "Out2")
        }
        else {
            if (length(cols) != 2) 
                stop("Outcome names not of length 2.")
        }
        dat <- array(NA, dim = c(2, 2, k), dimnames = list(rows, 
            cols, slab))
        for (i in 1:k) {
            tab.i <- rbind(c(ai[i], bi[i]), c(ci[i], di[i]))
            dat[, , i] <- tab.i
        }
    }
    if (is.element(measure, c("IRR", "IRD", "IRSD"))) {
        x1ix2it1it2i.na <- is.na(cbind(x1i, x2i, t1i, t2i))
        if (any(x1ix2it1it2i.na)) {
            not.na <- apply(x1ix2it1it2i.na, MARGIN = 1, sum) == 
                0L
            if (na.act == "na.omit") {
                x1i <- x1i[not.na]
                x2i <- x2i[not.na]
                t1i <- t1i[not.na]
                t2i <- t2i[not.na]
                slab <- slab[not.na]
                warning("Tables with NAs omitted.")
            }
            if (na.act == "na.fail") 
                stop("Missing values in tables.")
        }
        k <- length(x1i)
        if (k < 1) 
            stop("Processing terminated since k = 0.")
        if (missing(rows)) {
            rows <- c("Grp1", "Grp2")
        }
        else {
            if (length(rows) != 2) 
                stop("Group names not of length 2.")
        }
        if (missing(cols)) {
            cols <- c("Events", "Person-Time")
        }
        else {
            if (length(cols) != 2) 
                stop("Outcome names not of length 2.")
        }
        dat <- array(NA, dim = c(2, 2, k), dimnames = list(rows, 
            cols, slab))
        for (i in 1:k) {
            tab.i <- rbind(c(x1i[i], t1i[i]), c(x2i[i], t2i[i]))
            dat[, , i] <- tab.i
        }
    }
    if (is.element(measure, c("MD", "SMD", "SMDH", "ROM", "RPB", 
        "RBIS", "D2OR", "D2ORN", "D2ORL"))) {
        m1im2isd1isd2in1in2i.na <- is.na(cbind(m1i, m2i, sd1i, 
            sd2i, n1i, n2i))
        if (any(m1im2isd1isd2in1in2i.na)) {
            not.na <- apply(m1im2isd1isd2in1in2i.na, MARGIN = 1, 
                sum) == 0L
            if (na.act == "na.omit") {
                m1i <- m1i[not.na]
                m2i <- m2i[not.na]
                sd1i <- sd1i[not.na]
                sd2i <- sd2i[not.na]
                n1i <- n1i[not.na]
                n2i <- n2i[not.na]
                slab <- slab[not.na]
                warning("Tables with NAs omitted.")
            }
            if (na.act == "na.fail") 
                stop("Missing values in tables.")
        }
        k <- length(m1i)
        if (k < 1) 
            stop("Processing terminated since k = 0.")
        if (missing(rows)) {
            rows <- c("Grp1", "Grp2")
        }
        else {
            if (length(rows) != 2) 
                stop("Group names not of length 2.")
        }
        if (missing(cols)) {
            cols <- c("Mean", "SD", "n")
        }
        else {
            if (length(cols) != 3) 
                stop("Outcome names not of length 3.")
        }
        dat <- array(NA, dim = c(2, 3, k), dimnames = list(rows, 
            cols, slab))
        for (i in 1:k) {
            tab.i <- rbind(c(m1i[i], sd1i[i], n1i[i]), c(m2i[i], 
                sd2i[i], n2i[i]))
            dat[, , i] <- tab.i
        }
    }
    if (is.element(measure, c("COR", "UCOR", "ZCOR"))) {
        rini.na <- is.na(cbind(ri, ni))
        if (any(rini.na)) {
            not.na <- apply(rini.na, MARGIN = 1, sum) == 0L
            if (na.act == "na.omit") {
                ri <- ri[not.na]
                ni <- ni[not.na]
                slab <- slab[not.na]
                warning("Tables with NAs omitted.")
            }
            if (na.act == "na.fail") 
                stop("Missing values in tables.")
        }
        k <- length(ri)
        if (k < 1) 
            stop("Processing terminated since k = 0.")
        if (missing(rows)) {
            rows <- c("Grp")
        }
        else {
            if (length(rows) != 1) 
                stop("Group names not of length 1.")
        }
        if (missing(cols)) {
            cols <- c("r", "n")
        }
        else {
            if (length(cols) != 2) 
                stop("Outcome names not of length 2.")
        }
        dat <- array(NA, dim = c(1, 2, k), dimnames = list(rows, 
            cols, slab))
        for (i in 1:k) {
            tab.i <- c(ri[i], ni[i])
            dat[, , i] <- tab.i
        }
    }
    if (is.element(measure, c("PR", "PLN", "PLO", "PAS", "PFT"))) {
        ximi.na <- is.na(cbind(xi, mi))
        if (any(ximi.na)) {
            not.na <- apply(ximi.na, MARGIN = 1, sum) == 0L
            if (na.act == "na.omit") {
                xi <- xi[not.na]
                mi <- mi[not.na]
                slab <- slab[not.na]
                warning("Tables with NAs omitted.")
            }
            if (na.act == "na.fail") 
                stop("Missing values in tables.")
        }
        k <- length(xi)
        if (k < 1) 
            stop("Processing terminated since k = 0.")
        if (missing(rows)) {
            rows <- c("Grp")
        }
        else {
            if (length(rows) != 1) 
                stop("Group names not of length 1.")
        }
        if (missing(cols)) {
            cols <- c("Out1", "Out2")
        }
        else {
            if (length(cols) != 2) 
                stop("Outcome names not of length 2.")
        }
        dat <- array(NA, dim = c(1, 2, k), dimnames = list(rows, 
            cols, slab))
        for (i in 1:k) {
            tab.i <- c(xi[i], mi[i])
            dat[, , i] <- tab.i
        }
    }
    if (is.element(measure, c("IR", "IRLN", "IRS", "IRFT"))) {
        xiti.na <- is.na(cbind(xi, ti))
        if (any(xiti.na)) {
            not.na <- apply(xiti.na, MARGIN = 1, sum) == 0L
            if (na.act == "na.omit") {
                xi <- xi[not.na]
                ti <- ti[not.na]
                slab <- slab[not.na]
                warning("Tables with NAs omitted.")
            }
            if (na.act == "na.fail") 
                stop("Missing values in tables.")
        }
        k <- length(xi)
        if (k < 1) 
            stop("Processing terminated since k = 0.")
        if (missing(rows)) {
            rows <- c("Grp")
        }
        else {
            if (length(rows) != 1) 
                stop("Group names not of length 1.")
        }
        if (missing(cols)) {
            cols <- c("Events", "Person-Time")
        }
        else {
            if (length(cols) != 2) 
                stop("Outcome names not of length 2.")
        }
        dat <- array(NA, dim = c(1, 2, k), dimnames = list(rows, 
            cols, slab))
        for (i in 1:k) {
            tab.i <- c(xi[i], ti[i])
            dat[, , i] <- tab.i
        }
    }
    if (is.element(measure, c("MN"))) {
        misdini.na <- is.na(cbind(mi, sdi, ni))
        if (any(misdini.na)) {
            not.na <- apply(misdini.na, MARGIN = 1, sum) == 0L
            if (na.act == "na.omit") {
                mi <- mi[not.na]
                sdi <- sdi[not.na]
                ni <- ni[not.na]
                slab <- slab[not.na]
                warning("Tables with NAs omitted.")
            }
            if (na.act == "na.fail") 
                stop("Missing values in tables.")
        }
        k <- length(mi)
        if (k < 1) 
            stop("Processing terminated since k = 0.")
        if (missing(rows)) {
            rows <- c("Grp")
        }
        else {
            if (length(rows) != 1) 
                stop("Group names not of length 1.")
        }
        if (missing(cols)) {
            cols <- c("Mean", "SD", "n")
        }
        else {
            if (length(cols) != 3) 
                stop("Outcome names not of length 3.")
        }
        dat <- array(NA, dim = c(1, 3, k), dimnames = list(rows, 
            cols, slab))
        for (i in 1:k) {
            tab.i <- c(mi[i], sdi[i], ni[i])
            dat[, , i] <- tab.i
        }
    }
    if (is.element(measure, c("MC", "SMCC", "SMCR"))) {
        if (is.element(measure, c("MC", "SMCC"))) {
            m1im2isdiniri.na <- is.na(cbind(m1i, m2i, sd1i, sd2i, 
                ni, ri))
        }
        else {
            m1im2isdiniri.na <- is.na(cbind(m1i, m2i, sd1i, ni, 
                ri))
        }
        if (any(m1im2isdiniri.na)) {
            not.na <- apply(m1im2isdiniri.na, MARGIN = 1, sum) == 
                0L
            if (na.act == "na.omit") {
                m1i <- m1i[not.na]
                m2i <- m2i[not.na]
                sd1i <- sd1i[not.na]
                if (is.element(measure, c("MC", "SMCC"))) 
                  sd2i <- sd2i[not.na]
                ni <- ni[not.na]
                ri <- ri[not.na]
                slab <- slab[not.na]
                warning("Tables with NAs omitted.")
            }
            if (na.act == "na.fail") 
                stop("Missing values in tables.")
        }
        k <- length(m1i)
        if (k < 1) 
            stop("Processing terminated since k = 0.")
        if (missing(rows)) {
            rows <- c("Grp")
        }
        else {
            if (length(rows) != 1) 
                stop("Group names not of length 1.")
        }
        if (is.element(measure, c("MC", "SMCC"))) {
            if (missing(cols)) {
                cols <- c("Mean1", "Mean2", "SD1", "SD2", "n", 
                  "r")
            }
            else {
                if (length(cols) != 6) 
                  stop("Outcome names not of length 6.")
            }
        }
        else {
            if (missing(cols)) {
                cols <- c("Mean1", "Mean2", "SD1", "n", "r")
            }
            else {
                if (length(cols) != 5) 
                  stop("Outcome names not of length 5.")
            }
        }
        if (is.element(measure, c("MC", "SMCC"))) {
            dat <- array(NA, dim = c(1, 6, k), dimnames = list(rows, 
                cols, slab))
            for (i in 1:k) {
                tab.i <- c(m1i[i], m2i[i], sd1i[i], sd2i[i], 
                  ni[i], ri[i])
                dat[, , i] <- tab.i
            }
        }
        else {
            dat <- array(NA, dim = c(1, 5, k), dimnames = list(rows, 
                cols, slab))
            for (i in 1:k) {
                tab.i <- c(m1i[i], m2i[i], sd1i[i], ni[i], ri[i])
                dat[, , i] <- tab.i
            }
        }
    }
    if (is.element(measure, c("ARAW", "AHW", "ABT"))) {
        aimini.na <- is.na(cbind(ai, mi, ni))
        if (any(aimini.na)) {
            not.na <- apply(aimini.na, MARGIN = 1, sum) == 0L
            if (na.act == "na.omit") {
                ai <- ai[not.na]
                mi <- mi[not.na]
                ni <- ni[not.na]
                slab <- slab[not.na]
                warning("Tables with NAs omitted.")
            }
            if (na.act == "na.fail") 
                stop("Missing values in tables.")
        }
        k <- length(ai)
        if (k < 1) 
            stop("Processing terminated since k = 0.")
        if (missing(rows)) {
            rows <- c("Grp")
        }
        else {
            if (length(rows) != 1) 
                stop("Group names not of length 1.")
        }
        if (missing(cols)) {
            cols <- c("alpha", "m", "n")
        }
        else {
            if (length(cols) != 3) 
                stop("Outcome names not of length 3.")
        }
        dat <- array(NA, dim = c(1, 3, k), dimnames = list(rows, 
            cols, slab))
        for (i in 1:k) {
            tab.i <- c(ai[i], mi[i], ni[i])
            dat[, , i] <- tab.i
        }
    }
    return(dat)
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
    o.warn <- getOption("warn")
    on.exit(options(warn = o.warn))
    options(warn = -1)
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
    estimator <- match.arg(estimator, c("L0", "R0", "Q0"))
    if (x$k == 1) 
        stop("Stopped because k = 1.")
    yi <- x$yi
    vi <- x$vi
    weights <- x$weights
    ni <- x$ni
    if (is.null(side)) {
        res <- rma(yi, vi, weights = weights, mods = sqrt(vi), 
            intercept = TRUE, method = x$method, weighted = x$weighted, 
            ...)
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
    weights <- weights[idix]
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
        weights.t <- weights[1:(k - k0)]
        res <- rma(yi.t, vi.t, weights = weights.t, intercept = TRUE, 
            method = x$method, weighted = x$weighted, ...)
        b <- c(res$b)
        yi.c <- yi - b
        yi.c.r <- rank(abs(yi.c), ties.method = "first")
        yi.c.r.s <- sign(yi.c) * yi.c.r
        if (estimator == "R0") {
            k0 <- (k - max(-1 * yi.c.r.s[yi.c.r.s < 0])) - 1
            se.k0 <- sqrt(2 * max(0, k0) + 2)
        }
        if (estimator == "L0") {
            Sr <- sum(yi.c.r.s[yi.c.r.s > 0])
            k0 <- (4 * Sr - k * (k + 1))/(2 * k - 1)
            varSr <- 1/24 * (k * (k + 1) * (2 * k + 1) + 10 * 
                k0^3 + 27 * k0^2 + 17 * k0 - 18 * k * k0^2 - 
                18 * k * k0 + 6 * k^2 * k0)
            se.k0 <- 4 * sqrt(varSr)/(2 * k - 1)
        }
        if (estimator == "Q0") {
            Sr <- sum(yi.c.r.s[yi.c.r.s > 0])
            k0 <- k - 1/2 - sqrt(2 * k^2 - 4 * Sr + 1/4)
            varSr <- 1/24 * (k * (k + 1) * (2 * k + 1) + 10 * 
                k0^3 + 27 * k0^2 + 17 * k0 - 18 * k * k0^2 - 
                18 * k * k0 + 6 * k^2 * k0)
            se.k0 <- 2 * sqrt(varSr)/sqrt((k - 1/2)^2 - k0 * 
                (2 * k - k0 - 1))
        }
        k0 <- max(0, k0)
        k0 <- round(k0)
        se.k0 <- max(0, se.k0)
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
        weights.fill <- c(x$weights.f, weights[(k - k0 + 1):k])
        ni.fill <- c(x$ni.f, ni[(k - k0 + 1):k])
        attr(yi.fill, "measure") <- x$measure
        res <- rma(yi.fill, vi.fill, weights = weights.fill, 
            ni = ni.fill, intercept = TRUE, method = x$method, 
            weighted = x$weighted, ...)
        res$fill <- c(rep(0, k), rep(1, k0))
        res$ids <- c(x$ids, (x$k.f + 1):(x$k.f + k0))
        if (x$slab.null) {
            res$slab <- c(paste("Study", x$ids), paste("Filled", 
                seq_len(k0)))
            res$slab.null <- FALSE
        }
        else {
            res$slab <- c(x$slab, paste("Filled", seq_len(k0)))
            res$slab.null <- FALSE
        }
    }
    else {
        res <- x
        res$fill <- rep(0, k)
    }
    res$k0 <- k0
    res$se.k0 <- se.k0
    res$side <- side
    res$k0.est <- estimator
    if (estimator == "R0") {
        m <- -1:(k0 - 1)
        res$p.k0 <- 1 - sum(choose(0 + m + 1, m + 1) * 0.5^(0 + 
            m + 2))
    }
    else {
        res$p.k0 <- NA
    }
    class(res) <- c("rma.uni.trimfill", class(res))
    return(res)
}
vcov.rma <-
function (object, ...) 
{
    if (!is.element("rma", class(object))) 
        stop("Argument 'object' must be an object of class \"rma\".")
    return(object$vb)
}
weights.rma.mh <-
function (object, type = "diagonal", ...) 
{
    if (!is.element("rma.mh", class(object))) 
        stop("Argument 'object' must be an object of class \"rma.mh\".")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    type <- match.arg(type, c("diagonal", "matrix"))
    x <- object
    if (is.element(x$measure, c("RR", "OR", "RD"))) {
        Ni <- x$ai + x$bi + x$ci + x$di
    }
    else {
        Ti <- x$t1i + x$t2i
    }
    if (x$measure == "OR") 
        wi <- (x$bi/Ni) * x$ci
    if (x$measure == "RR") 
        wi <- (x$ci/Ni) * (x$ai + x$bi)
    if (x$measure == "RD") 
        wi <- ((x$ai + x$bi)/Ni) * (x$ci + x$di)
    if (x$measure == "IRR") 
        wi <- (x$x2i/Ti) * x$t1i
    if (x$measure == "IRD") 
        wi <- (x$t1i/Ti) * x$t2i
    if (type == "diagonal") {
        weight <- rep(NA, x$k.f)
        weight[x$not.na] <- wi/sum(wi) * 100
        names(weight) <- x$slab
        if (na.act == "na.omit") 
            weight <- weight[x$not.na]
        if (na.act == "na.fail" && any(!x$not.na)) 
            stop("Missing values in weights.")
        return(weight)
    }
    if (type == "matrix") {
        Wfull <- matrix(NA, nrow = x$k.f, ncol = x$k.f)
        Wfull[x$not.na, x$not.na] <- diag(wi)
        rownames(Wfull) <- x$slab
        colnames(Wfull) <- x$slab
        if (na.act == "na.omit") 
            Wfull <- Wfull[x$not.na, x$not.na, drop = FALSE]
        if (na.act == "na.fail" && any(!x$not.na)) 
            stop("Missing values in results.")
        return(Wfull)
    }
}
weights.rma.mv <-
function (object, type = "diagonal", ...) 
{
    if (!is.element("rma.mv", class(object))) 
        stop("Argument 'object' must be an object of class \"rma.mv\".")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    type <- match.arg(type, c("diagonal", "matrix"))
    x <- object
    if (is.null(x$W)) {
        W <- chol2inv(chol(x$M))
    }
    else {
        W <- x$W
    }
    if (type == "diagonal") {
        wi <- as.vector(diag(W))
        weight <- rep(NA, x$k.f)
        weight[x$not.na] <- wi/sum(wi) * 100
        names(weight) <- x$slab
        if (na.act == "na.omit") 
            weight <- weight[x$not.na]
        if (na.act == "na.fail" && any(!x$not.na)) 
            stop("Missing values in weights.")
        return(weight)
    }
    if (type == "matrix") {
        Wfull <- matrix(NA, nrow = x$k.f, ncol = x$k.f)
        Wfull[x$not.na, x$not.na] <- W
        rownames(Wfull) <- x$slab
        colnames(Wfull) <- x$slab
        if (na.act == "na.omit") 
            Wfull <- Wfull[x$not.na, x$not.na, drop = FALSE]
        if (na.act == "na.fail" && any(!x$not.na)) 
            stop("Missing values in results.")
        return(Wfull)
    }
}
weights.rma.peto <-
function (object, type = "diagonal", ...) 
{
    if (!is.element("rma.peto", class(object))) 
        stop("Argument 'object' must be an object of class \"rma.peto\".")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    type <- match.arg(type, c("diagonal", "matrix"))
    x <- object
    n1i <- x$ai + x$bi
    n2i <- x$ci + x$di
    Ni <- x$ai + x$bi + x$ci + x$di
    xt <- x$ai + x$ci
    yt <- x$bi + x$di
    wi <- xt * yt * (n1i/Ni) * (n2i/Ni)/(Ni - 1)
    if (type == "diagonal") {
        weight <- rep(NA, x$k.f)
        weight[x$not.na] <- wi/sum(wi) * 100
        names(weight) <- x$slab
        if (na.act == "na.omit") 
            weight <- weight[x$not.na]
        if (na.act == "na.fail" && any(!x$not.na)) 
            stop("Missing values in weights.")
        return(weight)
    }
    if (type == "matrix") {
        Wfull <- matrix(NA, nrow = x$k.f, ncol = x$k.f)
        Wfull[x$not.na, x$not.na] <- diag(wi)
        rownames(Wfull) <- x$slab
        colnames(Wfull) <- x$slab
        if (na.act == "na.omit") 
            Wfull <- Wfull[x$not.na, x$not.na, drop = FALSE]
        if (na.act == "na.fail" && any(!x$not.na)) 
            stop("Missing values in results.")
        return(Wfull)
    }
}
weights.rma.uni <-
function (object, type = "diagonal", ...) 
{
    if (!is.element("rma.uni", class(object))) 
        stop("Argument 'object' must be an object of class \"rma.uni\".")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    type <- match.arg(type, c("diagonal", "matrix"))
    x <- object
    if (is.null(x$weights)) {
        W <- diag(1/(x$vi + x$tau2), nrow = x$k, ncol = x$k)
    }
    else {
        W <- diag(x$weights, nrow = x$k, ncol = x$k)
    }
    if (type == "diagonal") {
        wi <- as.vector(diag(W))
        weight <- rep(NA, x$k.f)
        weight[x$not.na] <- wi/sum(wi) * 100
        names(weight) <- x$slab
        if (na.act == "na.omit") 
            weight <- weight[x$not.na]
        if (na.act == "na.fail" && any(!x$not.na)) 
            stop("Missing values in weights.")
        return(weight)
    }
    if (type == "matrix") {
        Wfull <- matrix(NA, nrow = x$k.f, ncol = x$k.f)
        Wfull[x$not.na, x$not.na] <- W
        rownames(Wfull) <- x$slab
        colnames(Wfull) <- x$slab
        if (na.act == "na.omit") 
            Wfull <- Wfull[x$not.na, x$not.na, drop = FALSE]
        if (na.act == "na.fail" && any(!x$not.na)) 
            stop("Missing values in results.")
        return(Wfull)
    }
}
