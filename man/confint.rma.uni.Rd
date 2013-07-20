\name{confint.rma.uni}
\alias{confint.rma.uni}
\alias{confint.rma.mh}
\alias{confint.rma.peto}
\alias{confint.rma.glmm}
\title{Compute Confidence Intervals for 'rma' Objects}
\description{
   The function calculates a confidence interval for the amount of (residual) heterogeneity in random- and mixed-effects models that were fitted with the \code{\link{rma.uni}} function. Confidence intervals for the model coefficients can also be obtained.
}
\usage{
\method{confint}{rma.uni}(object, parm, level=object$level, fixed=FALSE, random=TRUE,
        digits=object$digits, verbose=FALSE, control, \dots)
\method{confint}{rma.mh}(object, parm, level=object$level, digits=object$digits, \dots)
\method{confint}{rma.peto}(object, parm, level=object$level, digits=object$digits, \dots)
\method{confint}{rma.glmm}(object, parm, level=object$level, digits=object$digits, \dots)
}
\arguments{
   \item{object}{an object of class \code{"rma.uni"}, \code{"rma.mh"}, or \code{"rma.peto"}. The method is not yet implemented for objects of class \code{"rma.glmm"}.}
   \item{parm}{this argument is here for compatability with the generic function \code{\link{confint}}, but is ignored.}
   \item{level}{numerical value between 0 and 100 specifying the confidence interval level (the default is to take the value from the object).}
   \item{fixed}{logical indicating whether confidence intervals for the model coefficients should be returned (default is \code{FALSE}).}
   \item{random}{logical indicating whether a confidence interval for the amount of (residual) heterogeneity should be returned (default is \code{TRUE}).}
   \item{digits}{integer specifying the number of decimal places to which the results should be rounded (the default is to take the value from the object).}
   \item{verbose}{logical indicating whether output should be generated for the iterative algorithm used as part of the Q-profile method (default is \code{FALSE}). See \sQuote{Note}.}
   \item{control}{list of control values for the iterative algorithm. If unspecified, default values are defined inside the function. See \sQuote{Note}.}
   \item{\dots}{other arguments.}
}
\details{
   The confidence intervals for the model coefficients are simply the usual Wald-type intervals which are also shown when printing the fitted object.

   The confidence interval for the amount of (residual) heterogeneity is obtained iteratively via the Q-profile method as described by Hartung and Knapp (2005) and Viechtbauer (2007). The method provides an exact confidence interval for \eqn{\tau^2} in random- and mixed-effects models. The square root of the interval bounds is also returned for easier interpretation. Confidence intervals for \eqn{I^2} and \eqn{H^2} are also provided (Higgins & Thompson, 2002). Since \eqn{I^2} and \eqn{H^2} are just monotonic transformation of \eqn{\tau^2}, the confidence intervals for \eqn{I^2} and \eqn{H^2} are also exact.
}
\value{
   An object of class \code{"confint.rma"}. The object is a list with either one or two elements (named \code{fixed} and \code{random}) with the following elements:
   \item{estimate}{estimate of the model coefficient or variance component.}
   \item{ci.lb}{lower bound of the confidence interval.}
   \item{ci.ub}{upper bound of the confidence interval.}

   For fixed-effects models, the variance component estimates are \code{NA}, but the confidence interval bounds are still provided.

   The results are formated and printed with the \code{\link{print.confint.rma}} function.
}
\note{
   The iterative algorithm used as part of the Q-profile method makes use of the \code{\link{uniroot}} function. By default, the desired accuracy is set equal to \code{.Machine$double.eps^0.25} and the maximum number of iterations to \code{1000}. The upper bound of the interval searched is set to 50 (which should be large enough for most cases). The desired accuracy (\code{tol}), maximum number of iterations (\code{maxiter}), and upper bound (\code{tau2.max}) can be adjusted with \code{control=list(tol=value, maxiter=value, tau2.max=value)}. One can also adjust the lower bound of the interval searched with \code{control=list(tau2.min=value)} (the default is to take that value from the object, which is 0 by default). You should only play around with this latter value if you know what you are doing. If \code{verbose=TRUE}, output is generated for the iterative algorithm.

   It is possible that the lower and upper confidence interval bounds both fall below zero (or whatever value was chosen for \code{tau2.min}). Since both values then fall outside of the parameter space, the confidence interval then just consists of the null set.

   Usually, the estimate of \eqn{\tau^2} from the random/mixed-effects model will fall within the confidence interval provided by the Q-profile method. However, this is not guaranteed. Depending on the method used to estimate \eqn{\tau^2} and the width of the confidence interval, it can happen in rare cases that the confidence interval does not actually contain the estimate (trying to explain this to reviewers can be tricky). However, using the empirical Bayes estimator of \eqn{\tau^2} when fitting the model (i.e., using \code{method="EB"}) guarantees that the estimate of \eqn{\tau^2} falls within the confidence interval.

   The Q-profile method is exact under the assumptions of the random- and mixed-effects models (i.e., normally distributed observed and true outcomes and known sampling variances). These assumptions are usually only approximately true, turning the confidence interval for \eqn{\tau^2} also into an approximation.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \cr
   package homepage: \url{http://www.metafor-project.org/} \cr
   author homepage: \url{http://www.wvbauer.com/}
}
\references{
   Hartung, J., & Knapp, G. (2005). On confidence intervals for the among-group variance in the one-way random effects model with unequal error variances. \emph{Journal of Statistical Planning and Inference}, \bold{127}, 157--177.

   Higgins, J. P. T., & Thompson, S. G. (2002). Quantifying heterogeneity in a meta-analysis. \emph{Statistics in Medicine}, \bold{21}, 1539--1558.

   Mantel, N., & Haenszel, W. (1959). Statistical aspects of the analysis of data from retrospective studies of disease. \emph{Journal of the National Cancer Institute}, \bold{22}, 719--748.

   Viechtbauer, W. (2007). Confidence intervals for the amount of heterogeneity in meta-analysis. \emph{Statistics in Medicine}, \bold{26}, 37--52.

   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \url{http://www.jstatsoft.org/v36/i03/}.

   Yusuf, S., Peto, R., Lewis, J., Collins, R., & Sleight, P. (1985). Beta blockade during and after myocardial infarction: An overview of the randomized trials. \emph{Progress in Cardiovascular Disease}, \bold{27}, 335--371.
}
\seealso{
   \code{\link{rma.uni}}, \code{\link{rma.mh}}, \code{\link{rma.peto}}
}
\examples{
### load BCG vaccine data
data(dat.bcg)

### meta-analysis of the log relative risks using a random-effects model
res <- rma(ai=tpos, bi=tneg, ci=cpos, di=cneg,
           data=dat.bcg, measure="RR", method="REML")

### confidence interval for the total amount of heterogeneity
confint(res)

### mixed-effects model with absolute latitude in the model
res <- rma(ai=tpos, bi=tneg, ci=cpos, di=cneg, mods=ablat,
           data=dat.bcg, measure="RR", method="REML")

### confidence interval for the residual amount of heterogeneity
confint(res)
}
\keyword{models}
