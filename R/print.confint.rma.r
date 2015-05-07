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
