coef.summary.rma <-
function (object, ...) 
{
    if (!is.element("summary.rma", class(object))) 
        stop("Argument 'object' must be an object of class \"summary.rma\".")
    x <- object
    res.table <- cbind(estimate = x$b, se = x$se, zval = x$zval, 
        pval = x$pval, ci.lb = x$ci.lb, ci.ub = x$ci.ub)
    colnames(res.table) <- c("estimate", "se", "zval", "pval", 
        "ci.lb", "ci.ub")
    if (x$knha || x$robust) 
        colnames(res.table)[3] <- "tval"
    res.table <- data.frame(res.table)
    return(res.table)
}
