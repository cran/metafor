\name{regtest}
\alias{regtest}
\alias{regtest.rma}
\alias{regtest.default}
\title{Regression Test for Funnel Plot Asymmetry}
\description{The function can be used to carry out (various versions of) Egger's regression test for funnel plot asymmetry.}
\usage{
regtest(x, \dots)

\method{regtest}{rma}(x, model="rma", predictor="sei", ret.fit=FALSE, digits, \dots)

\method{regtest}{default}(x, vi, sei, ni, subset, model="rma", predictor="sei", ret.fit=FALSE, digits, \dots)
}
\arguments{
   \item{x}{an object of class \code{"rma"} or a vector with the observed effect sizes or outcomes.}
   \item{vi}{vector with the corresponding sampling variances.}
   \item{sei}{vector with the corresponding standard errors.}
   \item{ni}{vector with the corresponding sample sizes.}
   \item{subset}{optional (logical or numeric) vector indicating the subset of studies that should be included in the test. Only relevant when passing a vector via \code{x}.}
   \item{model}{either \code{"rma"} or \code{"lm"} to indicate the type of model to use for the regression test. See \sQuote{Details}.}
   \item{predictor}{either \code{"sei"} \code{"vi"}, \code{"ni"}, \code{"ninv"}, \code{"sqrtni"}, or \code{"sqrtninv"} to indicate the predictor to use for the regression test. See \sQuote{Details}.}
   \item{ret.fit}{logical indicating whether the full results from the fitted model should also be returned.}
   \item{digits}{integer specifying the number of decimal places to which the printed results should be rounded (the default is 4).}
   \item{\dots}{other arguments.}
}
\details{
   Various tests for funnel plot asymmetry have been suggested in the literature, including the rank correlation test by Begg and Mazumdar (1994) and the regression test by Egger et al. (1997). Extensions, modifications, and further developments of the regression test are described (among others) by Macaskill, Walter, and Irwig (2001), Sterne and Egger (2005), Harbord, Egger, and Sterne (2006), Peters et al. (2006), \enc{Rücker}{Ruecker} et al. (2008), and Moreno et al. (2009). The various versions of the regression test differ in terms of the model (either a weighted regression model with a multiplicative dispersion term or a fixed/random/mixed-effects meta-regression model is used), in terms of the predictor variable that the observed outcomes are hypothesized to be related to when publication bias is present (suggested predictors include the standard error, the sampling variance, and the sample size -- or transformations thereof), and in terms of the outcome measure used (e.g., for \eqn{2 \times 2}{2x2} table data, one has the choice between various outcome measures). The idea behind the various tests is the same though: If there is a relationship between the observed outcomes and the chosen predictor, then this usually implies asymmetry in the funnel plot, which in turn may be an indication of publication bias.

   The \code{regtest.rma} function can be used to carry out various versions of the regression test. The model is chosen via the \code{model} argument, with \code{model="lm"} for weighted regression with a multiplicative dispersion term or \code{model="rma"} for the meta-analytic models. In the latter case, arguments such as \code{method}, \code{weighted}, and \code{test} used during the initial model fitting are also used for the regression test. Therefore, if one wants to conduct the regression test with a random/mixed-effects model, one should first fit a random-effects model with the \code{rma} function and then use the \code{regtest.rma} function on the fitted model object.

   The predictor is chosen via the \code{predictor} argument:
   \itemize{
   \item \code{predictor="sei"} for the standard error,
   \item \code{predictor="vi"} for the sampling variance,
   \item \code{predictor="ni"} for the sample size,
   \item \code{predictor="ninv"} for the inverse of the sample size,
   \item \code{predictor="sqrtni"} for the square root transformed sample size, and
   \item \code{predictor="sqrtninv"} for the inverse of the square root transformed sample size.
   } For predictors based on the sample size, the object \code{x} obviously must contain the information about the sample sizes. This will automatically be the case when \code{measure} was \emph{not} equal to \code{"GEN"} or the \code{ni} values were explicitly specified during the initial model fitting.

   If the model passed to the \code{regtest.rma} function already included one or more moderators, then \code{regtest.rma} will add the chosen \code{predictor} to the moderator(s) already included in the model. This way, one can test for funnel plot asymmetry after accounting first for the influence of the moderator(s) already included.

   One can also pass a vector of observed effect sizes or outcomes (via \code{x}) and the corresponding sampling variances via \code{vi} (or the standard errors via \code{sei}) directly to the function (in this case, the \code{regtest.default} function is used). When the \code{predictor} is the sample size or a transformation thereof, then \code{ni} needs to be specified here as well.

   The outcome measure used for the regression test is simply determined by what measure was used in fitting the original model (or what values are passed to \code{regtest.default}).
}
\value{
   An object of class \code{"regtest.rma"}. The object is a list containing the following components:
   \item{model}{the model used for the regression test.}
   \item{predictor}{the predictor used for the regression test.}
   \item{zval}{the value of the test statistic.}
   \item{pval}{the corresponding p-value}
   \item{dfs}{the degrees of freedom of the test stastistic (if the test is based on a t-distribution).}
   \item{fit}{the full results from the fitted model.}

   The results are formatted and printed with the \code{\link{print.regtest.rma}} function.
}
\note{
   The classical ``Egger test'' is obtained by setting \code{model="lm"} and \code{predictor="sei"}. For the random/mixed-effects version of the Egger test, one should first fit a random-effects model to the data and then set \code{model="rma"} and \code{predictor="sei"} when using the \code{regtest.rma} function. See Sterne and Egger (2005) for details on these two types of models/tests.

   All of the tests do not directly test for publication bias, but for a relationship between the observed outcomes and the chosen predictor. If such a relationship is present, then this usually implies asymmetry in the funnel plot, which in turn may be an indication of publication bias. However, it is important to keep in mind that there can be other reasons besides publication bias that could lead to asymmetry in the funnel plot.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{http://www.metafor-project.org/}
}
\references{
   Begg, C. B., & Mazumdar, M. (1994). Operating characteristics of a rank correlation test for publication bias. \emph{Biometrics}, \bold{50}, 1088--1101.

   Egger, M., Davey Smith, G., Schneider, M., & Minder, C. (1997). Bias in meta-analysis detected by a simple, graphical test. \emph{British Medical Journal}, \bold{315}, 629--634.

   Harbord, R. M., Egger, M., & Sterne, J. A. C. (2006). A modified test for small-study effects in meta-analyses of controlled trials with binary endpoints. \emph{Statistics in Medicine}, \bold{25}, 3443--3457.

   Macaskill, P., Walter, S. D., & Irwig, L. (2001). A comparison of methods to detect publication bias in meta-analysis. \emph{Statistics in Medicine}, \bold{20}, 641--654.

   Moreno, S. G., Sutton, A. J., Ades, A. E., Stanley, T. D., Abrams, K. R., Peters, J. L., & Cooper, N. J. (2009). Assessment of regression-based methods to adjust for publication bias through a comprehensive simulation study. \emph{BMC Medical Research Methodology}, \bold{9}:2.

   Peters, J. L., Sutton, A. J., Jones, D. R., Abrams, K. R., & Rushton, L. (2006). Comparison of two methods to detect publication bias in meta-analysis. \emph{Journal of the American Medical Association}, \bold{295}, 676--680.

   \enc{Rücker}{Ruecker}, G., Schwarzer, G., & Carpenter, J. (2008). Arcsine test for publication bias in meta-analyses with binary outcomes. \emph{Statistics in Medicine}, \bold{27}, 746--763.

   Sterne, J. A. C., & Egger, M. (2005). Regression methods to detect publication and other bias in meta-analysis. In H. R. Rothstein, A. J. Sutton, & M. Borenstein (Eds.) \emph{Publication bias in meta-analysis: Prevention, assessment, and adjustments} (pp. 99--110). Chichester, England: Wiley.

   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \url{https://www.jstatsoft.org/v036/i03}.
}
\seealso{
   \code{\link{ranktest}}
}
\examples{
### calculate log risk ratios and corresponding sampling variances
dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

### fit random-effects model
res <- rma(yi, vi, data=dat)

### classical Egger test
regtest(res, model="lm")

### random/mixed-effects version of the Egger test
regtest(res)

### same tests, but passing outcomes directly
regtest(dat$yi, dat$vi, model="lm")
regtest(dat$yi, dat$vi)

### examples using the sample size (or a transformation thereof) as predictor
regtest(res, model="lm", predictor="ni")
regtest(res, model="lm", predictor="ninv")
regtest(res, model="rma", predictor="ni")
regtest(res, model="rma", predictor="ninv")

### if dat$yi is computed with escalc(), sample size information is stored in attributes
dat$yi

### then this will work
regtest(dat$yi, dat$vi, predictor="ni")

### otherwise have to supply sample sizes manually
dat$ni <- with(dat, tpos + tneg + cpos + cneg)
regtest(dat$yi, dat$vi, dat$ni, predictor="ni")

### testing for asymmetry after accounting for the influence of a moderator
res <- rma(yi, vi, mods = ~ ablat, data=dat)
regtest(res, model="lm")
regtest(res)
}
\keyword{htest}
