\name{print.escalc}
\alias{print.escalc}
\alias{summary.escalc}
\alias{print.summary.escalc}
\title{Print and Summary Methods for escalc Objects}
\description{
   Print and summary methods for objects of class \code{"escalc"}.
}
\usage{
\method{print}{escalc}(x, digits, \dots)
\method{summary}{escalc}(object, out.names=c("sei","zi","ci.lb","ci.ub"), var.names,
        append=TRUE, replace=TRUE, level=95, digits, transf=FALSE, \dots)
\method{print}{summary.escalc}(x, digits, \dots)
}
\arguments{
   \item{x}{an object of class \code{"escalc"} or \code{"summary.escalc"}.}
   \item{object}{an object of class \code{"escalc"}.}
   \item{digits}{integer specifying the number of decimal places to which the printed results should be rounded (the default is to take the value from the object if possible).}
   \item{out.names}{character string with four elements, specifying the variable names for the standard errors, test statistics, and lower/upper confidence interval bounds.}
   \item{var.names}{character string with two elements, specifying the variable names for the observed outcomes and the sampling variances (the default is to take the value from the object if possible).}
   \item{append}{logical indicating whether the data frame specified via the \code{object} argument should be returned together with the additional variables that are calculated by the \code{summary} function (default is \code{TRUE}).}
   \item{replace}{logical indicating whether existing values for \code{sei}, \code{zi}, \code{ci.lb}, and \code{ci.ub} in the data frame should be replaced or not. Only relevant when the data frame already contains these variables. If \code{replace=TRUE} (the default), all of the existing values will be overwritten. If \code{replace=FALSE}, only \code{NA} values will be replaced.}
   \item{level}{numerical value between 0 and 100 specifying the confidence interval level (the default is 95).}
   \item{transf}{optional argument specifying the name of a function that should be used to transform the observed outcomes and interval bounds (e.g., \code{transf=exp}). Defaults to \code{FALSE}, which means that no transformation is used (any additional arguments needed for the function specified here can be passed via \code{\ldots}).}
   \item{\dots}{other arguments.}
}
\value{
   The \code{print.escalc} function formats and prints the data frame, so that the observed outcomes and sampling variances are rounded (to the number of digits specified).

   The \code{summary.escalc} function creates an object of class \code{c("summary.escalc","data.frame")}. The object is a data frame containing the original data (if \code{append=TRUE}) and the following components:
   \item{yi}{observed outcomes or effect size estimates (transformed if \code{transf} is not \code{FALSE}).}
   \item{vi}{corresponding (estimated) sampling variances.}
   \item{sei}{standard errors of the observed outcomes or effect size estimates.}
   \item{zi}{test statistics for testing \eqn{H_0: \theta_i = 0} (i.e., \code{yi/sei}).}
   \item{ci.lb}{lower confidence interval bounds (transformed if \code{transf} is not \code{FALSE}).}
   \item{ci.ub}{upper confidence interval bounds (transformed if \code{transf} is not \code{FALSE}).}

   Note that the actual variable names above depend on the \code{out.names} (and \code{var.names}) arguments. If the data frame already contains variables with names as specified by the \code{out.names} argument, the values for these variables will be overwritten when \code{replace=TRUE} (which is the default). By setting \code{replace=FALSE}, only values that are \code{NA} will be replaced.

   The \code{print.summary.escalc} function formats and prints the data frame.
}
\note{
   If some transformation function has been specified for the \code{transf} argument, then \code{yi}, \code{ci.lb}, and \code{ci.ub} will be transformed accordingly. However, \code{vi} and \code{sei} then still reflect the sampling variances and standard errors of the untransformed values.

   The \code{summary.escalc} function computes \code{level} \% Wald-type confidence intervals, which may or may not be the most accurate method for computing confidence intervals for the chosen outcome or effect size measure.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \cr
   package homepage: \url{http://www.metafor-project.org/} \cr
   author homepage: \url{http://www.wvbauer.com/}
}
\references{
   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \url{http://www.jstatsoft.org/v36/i03/}.
}
\seealso{
   \code{\link{escalc}}
}
\examples{
### load BCG vaccine data
data(dat.bcg)

### calculate log relative risks and corresponding sampling variances
dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
dat

### apply summary function
summary(dat)
summary(dat, transf=exp)
}
\keyword{print}
