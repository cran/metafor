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
