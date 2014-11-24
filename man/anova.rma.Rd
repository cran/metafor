\name{anova.rma}
\alias{anova.rma}
\title{Likelihood Ratio and Wald-Type Tests for 'rma' Objects}
\description{
   For two (nested) models of class \code{"rma.uni"} or \code{"rma.mv"}, the function provides a full versus reduced model comparison in terms of model fit statistics and a likelihood ratio test. When a single model is specified, a Wald-type test of one or more model coefficients is carried out.
}
\usage{
\method{anova}{rma}(object, object2, btt, digits, \dots)
}
\arguments{
   \item{object}{an object of class \code{"rma.uni"} or \code{"rma.mv"}.}
   \item{object2}{an (optional) object of class \code{"rma.uni"} or \code{"rma.mv"}. See \sQuote{Details}.}
   \item{btt}{optional vector of indices specifying which coefficients to include in the Wald-type test. See \sQuote{Details}.}
   \item{digits}{integer specifying the number of decimal places to which the printed results should be rounded (if unspecified, the default is to take the value from the object).}
   \item{\dots}{other arguments.}
}
\details{
   When a single model is specified, the function provides a Wald-type test of one or more model coefficients. For a fixed- or random-effects model (i.e., a model without moderators), this is just the test of the single coefficient of the model. For models including moderators, an omnibus test of all the model coefficients is conducted that excludes the intercept (the first coefficient) if it is included in the model. If no intercept is included in the model, then the omnibus test includes all of the coefficients in the model including the first. Alternatively, one can manually specify the indices of the coefficients to test via the \code{btt} argument. For example, with \code{btt=c(3,4)}, only the third and fourth coefficient from the model would be included in the test (if an intercept is included in the model, then it corresponds to the first coefficient in the model).

   When specifying two models for comparison, the function provides a likelihood ratio test comparing the two models. The two models must be based on the same set of data, must be of the same class, and should be nested for the likelihood ratio test to make sense.
}
\value{
   An object of class \code{"anova.rma"}. When a single model is specified, the object is a list containing the following components:
   \item{QM}{test statistic for the Wald-type test of model coefficients.}
   \item{QMp}{p-value for the Wald-type test of model coefficients.}
   \item{btt}{indices of the coefficients tested by the Wald-type test.}
   \item{k}{number of outcomes included in the model fitting.}
   \item{p}{number of coefficients in the model (including the intercept).}
   \item{m}{number of coefficients included in the Wald-type test of coefficients.}
   \item{\dots}{some additional elements/values.}
   When two models are specified, the object is a list containing the following components:
   \item{fit.stats.f}{log-likelihood, deviance, AIC, BIC, and AICc for the full model.}
   \item{fit.stats.r}{log-likelihood, deviance, AIC, BIC, and AICc for the reduced model.}
   \item{p.f}{number of parameters in the full model.}
   \item{p.r}{number of parameters in the reduced model.}
   \item{LRT}{likelihood ratio test statistic.}
   \item{pval}{p-value for the likelihood ratio test.}
   \item{QE.f}{test statistic for the test of (residual) heterogeneity from the full model.}
   \item{QE.r}{test statistic for the test of (residual) heterogeneity from the reduced model.}
   \item{tau2.f}{estimated \if{latex}{\eqn{\tau^2}}\if{html}{\eqn{\tau²}}\if{text}{\eqn{\tau^2}} value from the full model. \code{NA} for \code{"rma.mv"} objects.}
   \item{tau2.r}{estimated \if{latex}{\eqn{\tau^2}}\if{html}{\eqn{\tau²}}\if{text}{\eqn{\tau^2}} value from the reduced model. \code{NA} for \code{"rma.mv"} objects.}
   \item{R2}{amount of (residual) heterogeneity in the reduced model that is accounted for in the full model. \code{NA} for fixed-effects models, if the amount of heterogeneity in the reduced model is equal to zero, or for \code{"rma.mv"} objects. This can be regarded as a pseudo \if{latex}{\eqn{R^2}}\if{html}{\eqn{R²}}\if{text}{\eqn{R^2}} statistic (Raudenbush, 2009).}
   \item{\dots}{some additional elements/values.}

   The results are formated and printed with the \code{\link{print.anova.rma}} function.
}
\note{
   Note that likelihood ratio tests are not meaningful when using REML estimation and the two models have different fixed effects.

   In principle, one can also consider likelihood ratio tests for amount of (residual) heterogeneity in random- and mixed-effects models. The full model should then be fitted with either \code{method="ML"} or \code{method="REML"} and the reduced model with \code{method="FE"}. The p-value from that test is based on a chi-square distribution with 1 degree of freedom, but actually needs to be adjusted for the fact that the parameter (i.e., \if{latex}{\eqn{\tau^2}}\if{html}{\eqn{\tau²}}\if{text}{\eqn{\tau^2}}) falls on the boundary of the parameter space under the null hypothesis. Moreover, the Q-test usually keeps better control of the Type I error rate and therefore should be preferred (see Viechtbauer, 2007, for more details).
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \cr
   package website: \url{http://www.metafor-project.org/} \cr
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
   \code{\link{rma.uni}}, \code{\link{print.anova.rma}}
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

### Wald-type test of the two moderators
anova(res1)

### corresponding likelihood ratio test
anova(res1, res2)
}
\keyword{models}
