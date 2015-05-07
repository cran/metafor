confint.rma.uni <-
function (object, parm, level, fixed = FALSE, random = TRUE, 
    digits, transf, targs, verbose = FALSE, control, ...) 
{
    if (!is.element("rma.uni", class(object))) 
        stop("Argument 'object' must be an object of class \"rma.uni\".")
    x <- object
    k <- x$k
    p <- x$p
    yi <- x$yi
    vi <- x$vi
    X <- x$X
    Y <- cbind(yi)
    weights <- x$weights
    if (missing(level)) 
        level <- x$level
    if (missing(digits)) 
        digits <- x$digits
    if (missing(transf)) 
        transf <- FALSE
    if (missing(targs)) 
        targs <- NULL
    if (missing(control)) 
        control <- list()
    if (!fixed && !random) 
        stop("At least one of the arguments 'fixed' and 'random' must be TRUE.")
    alpha <- ifelse(level > 1, (100 - level)/100, 1 - level)
    con <- list(tol = .Machine$double.eps^0.25, maxiter = 1000, 
        tau2.min = ifelse(is.null(x$control$tau2.min), 0, x$control$tau2.min), 
        tau2.max = ifelse(is.null(x$control$tau2.max), 100, x$control$tau2.max), 
        verbose = FALSE)
    con[pmatch(names(control), names(con))] <- control
    if (verbose) 
        con$verbose <- verbose
    verbose <- con$verbose
    if (random) {
        if (k == 1) 
            stop("Stopped because k = 1.")
        lb.conv <- TRUE
        ub.conv <- TRUE
        if (x$method != "GENQ") {
            if (!x$allvipos) 
                stop("Cannot compute confidence interval for the amount of (residual)\n  heterogeneity with non-positive sampling variances in the data.")
            crit.u <- qchisq(alpha/2, k - p, lower.tail = FALSE)
            crit.l <- qchisq(alpha/2, k - p, lower.tail = TRUE)
            if (.QE.func(con$tau2.min, Y = Y, vi = vi, X = X, 
                k = k, objective = 0, verbose = FALSE) < crit.l) {
                tau2.lb <- NA
                tau2.ub <- NA
            }
            else {
                if (.QE.func(con$tau2.min, Y = Y, vi = vi, X = X, 
                  k = k, objective = 0, verbose = FALSE) > crit.u) {
                  tau2.lb <- try(uniroot(.QE.func, interval = c(con$tau2.min, 
                    con$tau2.max), tol = con$tol, maxiter = con$maxiter, 
                    Y = Y, vi = vi, X = X, k = k, objective = crit.u, 
                    verbose = verbose, digits = digits)$root, 
                    silent = TRUE)
                  if (!is.numeric(tau2.lb)) {
                    tau2.lb <- NA
                    lb.conv <- FALSE
                  }
                }
                else {
                  tau2.lb <- con$tau2.min
                }
                tau2.ub <- try(uniroot(.QE.func, interval = c(ifelse(is.na(tau2.lb), 
                  con$tau2.min, tau2.lb), con$tau2.max), tol = con$tol, 
                  maxiter = con$maxiter, Y = Y, vi = vi, X = X, 
                  k = k, objective = crit.l, verbose = verbose, 
                  digits = digits)$root, silent = TRUE)
                if (!is.numeric(tau2.ub)) {
                  tau2.ub <- NA
                  ub.conv <- FALSE
                }
            }
        }
        if (x$method == "GENQ") {
            if (!requireNamespace("CompQuadForm", quietly = TRUE)) 
                stop("Please install the 'CompQuadForm' package when method='QGEN'.")
            A <- diag(weights, nrow = k, ncol = k)
            stXAX <- .invcalc(X = X, W = A, k = k)
            P <- A - A %*% X %*% stXAX %*% t(X) %*% A
            Q <- crossprod(Y, P) %*% Y
            S <- diag(sqrt(vi), nrow = k, ncol = k)
            lambda <- Re(eigen(S %*% P %*% S, symmetric = TRUE, 
                only.values = TRUE)$values)
            test <- CompQuadForm::farebrother(Q, lambda[1:(k - 
                p)])$res
            if (test > 1 - alpha/2) {
                tau2.lb <- NA
                tau2.ub <- NA
            }
            else {
                if (test < alpha/2) {
                  tau2.lb <- try(uniroot(.GENQ.func, c(0, con$tau2.max), 
                    P = P, vi = vi, Q = Q, alpha = alpha/2, k = k, 
                    p = p, getlower = TRUE, verbose = verbose, 
                    digits = digits)$root, silent = TRUE)
                  if (!is.numeric(tau2.lb)) {
                    tau2.lb <- NA
                    lb.conv <- FALSE
                  }
                }
                else {
                  tau2.lb <- 0
                }
                tau2.ub <- try(uniroot(.GENQ.func, c(ifelse(is.na(tau2.lb), 
                  0, tau2.lb), con$tau2.max), P = P, vi = vi, 
                  Q = Q, alpha = alpha/2, k = k, p = p, getlower = FALSE, 
                  verbose = verbose, digits = digits)$root, silent = TRUE)
                if (!is.numeric(tau2.ub)) {
                  tau2.ub <- NA
                  ub.conv <- FALSE
                }
            }
        }
        if (!lb.conv) 
            warning("Error in iterative search for the lower bound.")
        if (!ub.conv) 
            warning("Error in iterative search for the upper bound.")
        if (!lb.conv || !ub.conv) 
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
        if (x$knha || x$robust) {
            crit <- qt(alpha/2, df = k - p, lower.tail = FALSE)
        }
        else {
            crit <- qnorm(alpha/2, lower.tail = FALSE)
        }
        b <- x$b
        ci.lb <- c(x$b - crit * x$se)
        ci.ub <- c(x$b + crit * x$se)
        if (is.function(transf)) {
            if (is.null(targs)) {
                b <- sapply(b, transf)
                ci.lb <- sapply(ci.lb, transf)
                ci.ub <- sapply(ci.ub, transf)
            }
            else {
                b <- sapply(b, transf, targs)
                ci.lb <- sapply(ci.lb, transf, targs)
                ci.ub <- sapply(ci.ub, transf, targs)
            }
        }
        res.fixed <- cbind(estimate = b, ci.lb = ci.lb, ci.ub = ci.ub)
        rownames(res.fixed) <- rownames(x$b)
        colnames(res.fixed) <- c("estimate", "ci.lb", "ci.ub")
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
