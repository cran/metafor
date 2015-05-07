confint.rma.mv <-
function (object, parm, level, digits, transf, targs, ...) 
{
    if (!is.element("rma.mv", class(object))) 
        stop("Argument 'object' must be an object of class \"rma.mv\".")
    stop("Method not yet implemented for objects of class \"rma.mv\". Sorry!")
}
