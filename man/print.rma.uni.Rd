\name{print.rma.uni}
\alias{print.rma.uni}
\alias{print.rma.mh}
\alias{print.rma.peto}
\alias{print.rma.glmm}
\alias{summary.rma}
\alias{print.summary.rma}
\title{Print and Summary Methods for rma Objects}
\description{
   Print and summary methods for objects of class \code{"rma.uni"}, \code{"rma.mh"}, \code{"rma.peto"}, and \code{"rma.glmm"}.
}
\usage{
\method{print}{rma.uni}(x, digits=x$digits, showfit=FALSE, signif.stars=getOption("show.signif.stars"),
      signif.legend=signif.stars, \dots)
\method{print}{rma.mh}(x, digits=x$digits, showfit=FALSE, \dots)
\method{print}{rma.peto}(x, digits=x$digits, showfit=FALSE, \dots)
\method{print}{rma.glmm}(x, digits=x$digits, showfit=FALSE, signif.stars=getOption("show.signif.stars"),
      signif.legend=signif.stars, \dots)
\method{summary}{rma}(object, digits=object$digits, showfit=TRUE, \dots)
\method{print}{summary.rma}(x, digits=x$digits, showfit=TRUE, signif.stars=getOption("show.signif.stars"),
      signif.legend=signif.stars, \dots)
}
\arguments{
   \item{x}{an object of class \code{"rma.uni"}, \code{"rma.mh"}, \code{"rma.peto"}, \code{"rma.glmm"}, or \code{"summary.rma"} (for \code{print}).}
   \item{object}{an object of class \code{"rma"} (for \code{summary}).}
   \item{digits}{integer specifying the number of decimal places to which the printed results should be rounded (the default is to take the value from the object).}
   \item{showfit}{logical indicating whether the fit statistics and information criteria should be printed (default is \code{FALSE} for \code{print} and \code{TRUE} for \code{summary}).}
   \item{signif.stars}{logical indicating whether p-values should be encoded visually with \sQuote{significance stars}. Defaults to the \code{show.signif.stars} slot of \code{\link{options}}.}
   \item{signif.legend}{logical indicating whether the legend for the \sQuote{significance stars} should be printed. Defaults to the value for \code{signif.stars}.}
   \item{\dots}{other arguments.}
}
\details{
   The output includes:
   \itemize{
   \item the log likelihood, deviance, AIC, and BIC value (when setting \code{showfit=TRUE} or by default for \code{summary}).
   \item the amount of (residual) heterogeneity in the random/mixed-effects model (i.e., the estimate of \eqn{\tau^2} and its square root). Suppressed for fixed-effects models. The (asymptotic) standard error of the estimate of \eqn{\tau^2} is also provided (where possible).
   \item the \eqn{I^2} statistic. For a random-effects model, \eqn{I^2} estimates (in percent) how much of the total variability in the effect size estimates (which is composed of heterogeneity plus sampling variability) can be attributed to heterogeneity among the true effects. For a mixed-effects model, \eqn{I^2} estimates how much of the unaccounted variability (which is composed of residual heterogeneity plus sampling variability) can be attributed to residual heterogeneity.
   \item the \eqn{H^2} statistic. For a random-effects model, \eqn{H^2} estimates the ratio of the total amount of variability in the effect size estimates to the amount of sampling variability. For a mixed-effects model, \eqn{H^2} estimates the ratio of the unaccounted variability in the effect size estimates to the amount of sampling variability.
   \item the \eqn{R^2} statistic, which estimates the amount of heterogeneity accounted for by the moderators included in the model and can be regarded as a pseudo \eqn{R^2} statistic (Raudenbush, 2009). Only for objects of class \code{"rma.uni"} and mixed-effects models (i.e., for models including moderators). This is suppressed (and set to \code{NULL}) for models without moderators, fixed-effects models, or if the model does not contain an intercept. \code{NA} if the amount of heterogeneity is equal to zero to begin with.
   \item the amount of study level variability (only for objects of class \code{"rma.glmm"} and when using a model that models study level differences as a random effect).
   \item the results of the test for (residual) heterogeneity. This is the usual Q-test for heterogeneity when not including moderators in the model and the QE-test for residual heterogeneity when moderators are included. For objects of class \code{"rma.glmm"}, the results from a Wald-type test and a likelihood ratio test are provided (see \code{\link{rma.glmm}} for more details).
   \item the results of the omnibus (Wald-type) test of the coefficients in the model (the indices of the coefficients tested are also indicated). Suppressed if the model includes only one coefficient (e.g., only an intercept, like in the fixed- and random-effects model).
   \item a table with the estimated coefficients, corresponding standard errors, test statistics, p-values, and confidence interval bounds.
   \item the Cochran-Mantel-Haenszel test and Tarone's test for heterogeneity (only when analyzing odds ratios using the Mantel-Haenszel method, i.e., \code{"rma.mh"}). 
   }
}
\value{
   The \code{print} functions do not return an object. The \code{summary} function returns the object passed to it (with additional class \code{"summary.rma"}).
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \cr
   package homepage: \url{http://www.metafor-project.org/} \cr
   author homepage: \url{http://www.wvbauer.com/}
}
\references{
   Raudenbush, S. W. (2009). Analyzing effect sizes: Random effects models. In H. Cooper, L. V. Hedges, & J. C. Valentine (Eds.), \emph{The handbook of research synthesis and meta-analysis} (2nd ed., pp. 295--315). New York: Russell Sage Foundation.

   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \url{http://www.jstatsoft.org/v36/i03/}.
}
\seealso{
   \code{\link{rma.uni}}, \code{\link{rma.mh}}, \code{\link{rma.peto}}, \code{\link{rma.glmm}}
}
\keyword{print}
