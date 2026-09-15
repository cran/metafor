# c(m) calculation function for the bias correction of SMDs and related measures

cmicalc <- function(mi, method="exact") {

   # the exact formula can overflow if mi is 'large' (if mi >= 344)
   # cmi <- gamma(mi/2)/(sqrt(mi/2)*gamma((mi-1)/2))
   # could catch those cases and apply the approximate formula (which is accurate then)
   # is.na <- is.na(cmi)
   # cmi[is.na] <- 1 - 3/(4*mi[is.na] - 1)
   # instead, use log() of the exact formula (i.e., using lgamma()) and then exponentiate, which avoids this issue

   if (is.logical(method)) {
      method <- ifelse(isTRUE(method), "exact", "none")
   } else {
      opts <- c("exact", "approx")
      method <- opts[pmatch(method, opts)]
      if (is.na(method)) {
         mstyle <- .get.mstyle()
         stop(mstyle$stop("Correction method must be either 'exact' or 'approx'."), call.=FALSE)
      }
   }

   if (method == "exact")
      cmi <- ifelse(mi <= 1, NA_real_, exp(lgamma(mi/2) - log(sqrt(mi/2)) - lgamma((mi-1)/2)))
   if (method == "approx")
      cmi <- ifelse(mi <= 1, NA_real_, 1 - 3 / (4*mi - 1))
   if (method == "none")
      cmi <- rep(1, length(mi))

   return(cmi)

}
