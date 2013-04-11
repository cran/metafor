\name{anova.rma.uni}
\alias{anova.rma.uni}
\title{Compare Fit Statistics and Likelihoods of rma.uni Objects}
\description{
   The function provides a full versus reduced model comparison of two objects of class \code{"rma.uni"}. Model fit statistics for the two models are provided. A likelihood ratio test comparing the two models is also performed.
}
\usage{
\method{anova}{rma.uni}(object, object2, digits=object$digits, \dots)
}
\arguments{
   \item{object}{an object of class \code{"rma.uni"}.}
   \item{object2}{an object of class \code{"rma.uni"}.}
   \item{digits}{integer specifying the number of decimal places to which the printed results should be rounded (the default is to take the value from the object).}
   \item{\dots}{other arguments.}
}
\details{
   The two models must be based on the same set of data and should be nested for the likelihood ratio test to make sense.
}
\value{
   An object of class \code{"anova.rma.uni"}. The object is a list containing the following components:
   \item{fit.stats.f}{log likelihood, deviance, AIC, and BIC for the full model.}
   \item{fit.stats.r}{log likelihood, deviance, AIC, and BIC for the reduced model.}
   \item{p.f}{number of parameters in the full model.}
   \item{p.r}{number of parameters in the reduced model.}
   \item{LRT}{likelihood ratio test statistic.}
   \item{pval}{p-value for the likelihood ratio test.}
   \item{QE.f}{test statistic for the test of (residual) heterogeneity from the full model.}
   \item{QE.r}{test statistic for the test of (residual) heterogeneity from the reduced model.}
   \item{tau2.f}{tau2 value from the full model.}
   \item{tau2.r}{tau2 value from the reduced model.}
   \item{VAF}{amount of (residual) heterogeneity in the reduced model that is accounted for in the full model. \code{NA} for fixed-effects models or if the amount of heterogeneity in the reduced model is equal to zero. This can be regarded as a pseudo-R^2 statistic (Raudenbush, 2009).}

   The results are formated and printed with the \code{\link{print.anova.rma.uni}} function.
}
\note{
   Note that likelihood ratio tests are not meaningful when using REML estimation and the two models have different fixed effects.

   In principle, one can also consider likelihood ratio tests for (residual) heterogeneity in random- and mixed-effects models. The full model should then be fitted with either \code{method="ML"} or \code{method="REML"} and the reduced model with \code{method="FE"}. The p-value from that test is based on a chi-square distribution with 1 degree of freedom, but actually needs to be adjusted for the fact that the parameter (i.e., \eqn{\tau^2}) falls on the boundary of the parameter space under the null hypothesis. Moreover, the Q-test usually keeps better control of the Type I error rate and therefore should be preferred (see Viechtbauer, 2007, for more details).
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \cr
   package homepage: \url{http://www.metafor-project.org/} \cr
   author homepage: \url{http://www.wvbauer.com/}
}
\references{
   Hardy, R. J., & Thompson, S. G. (1996). A likelihood approach to meta-analysis with random effects. \emph{Statistics in Medicine}, \bold{15}, 619--629.

   Huizenga, H. M., Visser, I., & Dolan, C. V. (2011). Testing overall and moderator effects in random effects meta-regression. \emph{British Journal of Mathematical and Statistical Psychology}, \bold{64}, 1--19.

   Raudenbush, S. W. (2009). Analyzing effect sizes: Random effects models. In H. Cooper, L. V. Hedges, & J. C. Valentine (Eds.), \emph{The handbook of research synthesis and meta-analysis} (2nd ed., pp. 295--315). New York: Russell Sage Foundation.

   Viechtbauer, W. (2007). Hypothesis tests for population heterogeneity in meta-analysis. \emph{British Journal of Mathematical and Statistical Psychology}, \bold{60}, 29--60.

   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \url{http://www.jstatsoft.org/v36/i03/}.
}
\seealso{
   \code{\link{rma.uni}}, \code{\link{print.anova.rma.uni}}
}
\examples{
### load BCG vaccine data
data(dat.bcg)

### calculate log relative risks and corresponding sampling variances
dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

### random-effects model
res2 <- rma(yi, vi, data=dat, method="ML")

### mixed-effects model with two moderators (absolute latitude and publication year)
res1 <- rma(yi, vi, mods = ~ ablat + year, data=dat, method="ML")

anova(res1, res2)
}
\keyword{models}
