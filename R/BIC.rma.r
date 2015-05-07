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
