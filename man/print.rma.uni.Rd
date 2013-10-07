\name{print.rma.uni}
\alias{print.rma.uni}
\alias{print.rma.mh}
\alias{print.rma.peto}
\alias{print.rma.glmm}
\alias{print.rma.mv}
\alias{summary.rma}
\alias{print.summary.rma}
\title{Print and Summary Methods for 'rma' Objects}
\description{
   Print and summary methods for objects of class \code{"rma.uni"}, \code{"rma.mh"}, \code{"rma.peto"}, \code{"rma.glmm"}, and \code{"rma.glmm"}.
}
\usage{
\method{print}{rma.uni}(x, digits=x$digits, showfit=FALSE, signif.stars=getOption("show.signif.stars"),
      signif.legend=signif.stars, \dots)
\method{print}{rma.mh}(x, digits=x$digits, showfit=FALSE, \dots)
\method{print}{rma.peto}(x, digits=x$digits, showfit=FALSE, \dots)
\method{print}{rma.glmm}(x, digits=x$digits, showfit=FALSE, signif.stars=getOption("show.signif.stars"),
      signif.legend=signif.stars, \dots)
\method{print}{rma.mv}(x, digits=x$digits, showfit=FALSE, signif.stars=getOption("show.signif.stars"),
      signif.legend=signif.stars, \dots)
\method{summary}{rma}(object, digits=object$digits, showfit=TRUE, \dots)
\method{print}{summary.rma}(x, digits=x$digits, showfit=TRUE, signif.stars=getOption("show.signif.stars"),
      signif.legend=signif.stars, \dots)
}
\arguments{
   \item{x}{an object of class \code{"rma.uni"}, \code{"rma.mh"}, \code{"rma.peto"}, \code{"rma.glmm"}, \code{"rma.mv"}, or \code{"summary.rma"} (for \code{print}).}
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
   \item the log-likelihood, deviance, AIC, BIC, and AICc value (when setting \code{showfit=TRUE} or by default for \code{summary}).
   \item for objects of class \code{"rma.uni"} and \code{"rma.glmm"}, the amount of (residual) heterogeneity in the random/mixed-effects model (i.e., the estimate of \eqn{\tau^2} and its square root). Suppressed for fixed-effects models. The (asymptotic) standard error of the estimate of \eqn{\tau^2} is also provided (where possible).
   \item for objects of \code{"rma.mv"}, a table providing information about the variance components and correlations in the model. For \eqn{\sigma^2} components, the estimate and its square root are provided, in addition to the number of values/levels, whether the component was fixed or estimated, and the name of the grouping variable/factor. If the \code{R} argument was used to specify known correlation matrices, this is also indicated. For models with an \code{~ inner | outer} formula term, the name of the inner and outer grouping variable/factor are given and the number of values/levels of these variables/factors. In addition, for each \eqn{\tau^2} component, the estimate and its square root are provided, the number of effects or outcomes observed at each level of the inner grouping variable/factor (only for \code{struct="HCS"} and \code{struct="UN"}), and whether the component was fixed or estimated. Finally, either the estimate of \eqn{\rho} (for \code{struct="CS"} or \code{struct="HCS"}) or the entire estimated correlation matrix (for \code{struct="UN"}) between the levels of the inner grouping variable/factor is provided, again with information whether a particular correlation was fixed or estimated, and how often each combination of levels of the inner grouping variable/factor was observed across the levels of the outer grouping variable/factor.
   \item for objects of class \code{"rma.uni"} and \code{"rma.glmm"}, the \eqn{I^2} statistic. For a random-effects model, \eqn{I^2} estimates (in percent) how much of the total variability in the effect size estimates (which is composed of heterogeneity plus sampling variability) can be attributed to heterogeneity among the true effects. For a mixed-effects model, \eqn{I^2} estimates how much of the unaccounted variability (which is composed of residual heterogeneity plus sampling variability) can be attributed to residual heterogeneity. See \sQuote{Note} for how \eqn{I^2} is computed.
   \item for objects of class \code{"rma.uni"} and \code{"rma.glmm"}, the \eqn{H^2} statistic. For a random-effects model, \eqn{H^2} estimates the ratio of the total amount of variability in the effect size estimates to the amount of sampling variability. For a mixed-effects model, \eqn{H^2} estimates the ratio of the unaccounted variability in the effect size estimates to the amount of sampling variability. See \sQuote{Note} for how \eqn{H^2} is computed.
   \item for objects of class \code{"rma.uni"}, the \eqn{R^2} statistic, which estimates the amount of heterogeneity accounted for by the moderators included in the model and can be regarded as a pseudo \eqn{R^2} statistic (Raudenbush, 2009). Only provided when fitting a mixed-effects models (i.e., for models including moderators). This is suppressed (and set to \code{NULL}) for models without moderators, fixed-effects models, or if the model does not contain an intercept. \code{NA} if the amount of heterogeneity is equal to zero to begin with. See \sQuote{Note} for how \eqn{R^2} is computed.
   \item for objects of class \code{"rma.glmm"}, the amount of study level variability (only when using a model that models study level differences as a random effect).
   \item the results of the test for (residual) heterogeneity. This is the usual \eqn{Q}-test for heterogeneity when not including moderators in the model and the \eqn{Q_E}-test for residual heterogeneity when moderators are included. For objects of class \code{"rma.glmm"}, the results from a Wald-type test and a likelihood ratio test are provided (see \code{\link{rma.glmm}} for more details).
   \item the results of the omnibus (Wald-type) test of the coefficients in the model (the indices of the coefficients tested are also indicated). Suppressed if the model includes only one coefficient (e.g., only an intercept, like in the fixed- and random-effects model).
   \item a table with the estimated coefficients, corresponding standard errors, test statistics, p-values, and confidence interval bounds.
   \item the Cochran-Mantel-Haenszel test and Tarone's test for heterogeneity (only when analyzing odds ratios using the Mantel-Haenszel method, i.e., \code{"rma.mh"}).
   }
}
\value{
   The \code{print} functions do not return an object. The \code{summary} function returns the object passed to it (with additional class \code{"summary.rma"}).
}
\note{
   The \eqn{I^2} statistic is computed with \deqn{I^2 = 100\% \times \frac{\hat{\tau}^2}{\hat{\tau}^2 + s^2},}{I^2 = 100\% hat(\tau)^2 / (hat(\tau)^2 + s^2),} where \eqn{\hat{\tau}^2}{hat(\tau)^2} is the estimated value of \eqn{\tau^2} (either the total or the residual amount of heterogeneity) and \deqn{s^2 = \frac{\sum w_i (k-1)}{(\sum w_i)^2 - \sum w_i^2},}{s^2 = (\sum w_i (k-1))((\sum w_i)^2 - \sum w_i^2),} where \eqn{w_i} is the inverse of the sampling variance of the \eqn{i}th study (\eqn{s^2} is equation 9 in Higgins & Thompson, 2002, and can be regarded as the \sQuote{typical} within-study variance of the observed effect sizes or outcomes). The \eqn{H^2} statistic is computed with \deqn{H^2 = \frac{\hat{\tau}^2 + s^2}{s^2}.}{H^2 = (hat(\tau)^2 + s^2) / s^2.}

   Therefore, depending on the estimator of \eqn{\tau^2} used, the values of \eqn{I^2} and \eqn{H^2} will change. For random-effects models, \eqn{I^2} and \eqn{H^2} are typically computed in practice with \eqn{I^2 = (Q-(k-1))/Q} and \eqn{H^2 = Q/(k-1)}, where \eqn{Q} denotes the statistic for the test of heterogeneity and \eqn{k} the number of studies (i.e., observed effects or outcomes) included in the meta-analysis (analogous equations also exist for mixed-effects models). The equations used in the metafor package to compute these statistics are more general and have the advantage that the values of \eqn{I^2} and \eqn{H^2} will be consistent with the estimated value of \eqn{\tau^2} (i.e., if \eqn{\hat{\tau}^2 = 0}{hat(\tau)^2 = 0}, then \eqn{I^2 = 0} and \eqn{H^2 = 1} and if \eqn{\hat{\tau}^2 > 0}{hat(\tau)^2 > 0}, then \eqn{I^2 > 0} and \eqn{H^2 > 1}).

   The two definitions of \eqn{I^2} and \eqn{H^2} actually coincide when using the DerSimonian-Laird estimator of \eqn{\tau^2} (i.e., the commonly used equations are actually special cases of the more general definitions given above). Therefore, if you prefer the more conventional definitions of these statistics, use \code{method="DL"} when fitting the random/mixed-effects model with the \code{\link{rma.uni}} function.

   The pseudo \eqn{R^2} statistic (Raudenbush, 2009) is computed with \deqn{R^2 = \frac{\hat{\tau}_{RE}^2 - \hat{\tau}_{ME}^2}{\hat{\tau}_{RE}^2},}{R^2 = (hat(\tau)_RE^2 - hat(\tau)_ME^2) / hat(\tau)_RE^2,} where \eqn{\hat{\tau}_{RE}^2}{hat(\tau)_RE^2} denotes the estimated value of \eqn{\tau^2} based on the random-effects model (i.e., the total amount of heterogeneity) and \eqn{\hat{\tau}_{ME}^2}{hat(\tau)_ME^2} denotes the estimated value of \eqn{\tau^2} based on the mixed-effects model (i.e., the residual amount of heterogeneity). It can happen that \eqn{\hat{\tau}_{RE}^2 < \hat{\tau}_{ME}^2}{hat(\tau)_RE^2 < hat(\tau)_ME^2}, in which case \eqn{R^2} is set to zero. Again, the value of \eqn{R^2} will change depending on the estimator of \eqn{\tau^2} used. Also note that this statistic is only computed when the mixed-effects model includes an intercept (so that the random-effects model is clearly nested within the mixed-effects model). You can also use the \code{\link{anova.rma.uni}} function to compute \eqn{R^2} for any two models that are known to be nested.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \cr
   package homepage: \url{http://www.metafor-project.org/} \cr
   author homepage: \url{http://www.wvbauer.com/}
}
\references{
   Higgins, J. P. T., & Thompson, S. G. (2002). Quantifying heterogeneity in a meta-analysis. \emph{Statistics in Medicine}, \bold{21}, 1539--1558.

   Raudenbush, S. W. (2009). Analyzing effect sizes: Random effects models. In H. Cooper, L. V. Hedges, & J. C. Valentine (Eds.), \emph{The handbook of research synthesis and meta-analysis} (2nd ed., pp. 295--315). New York: Russell Sage Foundation.

   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \url{http://www.jstatsoft.org/v36/i03/}.
}
\seealso{
   \code{\link{rma.uni}}, \code{\link{rma.mh}}, \code{\link{rma.peto}}, \code{\link{rma.glmm}}, \code{\link{rma.mv}}
}
\keyword{print}
