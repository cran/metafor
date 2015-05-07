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
