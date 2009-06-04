\name{ranktest}
\alias{ranktest}
\title{Rank Correlation Test for Funnel Plot Asymmetry}
\description{Rank correlation test for funnel plot asymmetry by Begg and Mazumdar (1994).}
\usage{
ranktest(yi, vi, sei, data=NULL, subset=NULL)
}
\arguments{
	\item{yi}{a vector with the observed effect sizes or outcomes.}
	\item{vi}{a vector with the corresponding sampling variances.}
	\item{sei}{a vector with the corresponding standard errors.\cr(note: only one of the two, \code{vi} or \code{sei}, needs to be specified)}
	\item{data}{an optional data frame containing the variables given to the arguments above.}
	\item{subset}{an optional vector indicating the subset of studies that should be used for the calculation. This can be a logical vector of length \eqn{k} or a numeric vector indicating the indices of the observations to include.}
}
\details{
	The function carries out the rank correlation test as described in Begg and Mazumdar (1994). The test can be used to examine whether the observed outcomes and the corresponding sampling variances are correlated. An asymmetric funnel plot would give rise to such a correlation and may be indicative of publication bias.

	The function is actually a wrapper for the \code{\link{cor.test}} function with \code{method="kendall"}. If possible, an exact p-value is provided; otherwise, a large-sample approximation is used.
}
\value{
	See the documentation for the \code{\link{cor.test}} function.
}
\note{
	Another test that is often used to detect funnel plot asymmetry was suggested by Egger et al. (1997) and is therefore usually called Egger's regression test. Extensions and further developments of this test are described by Sterne and Egger (2005) and Moreno et al. (2009). Other tests for funnel plot asymmetry were suggested, among others, by Macaskill, Walter, and Irwig (2001) and Peters et al. (2006). The various tests can be easily carried out with the \code{\link{lm}} and \code{\link{rma.uni}} functions. Examples are shown below.
}
\author{Wolfgang Viechtbauer; \email{wvb@www.wvbauer.com}; \url{http://www.wvbauer.com/}}
\references{
	Begg, C. B. & Mazumdar, M. (1994) Operating characteristics of a rank correlation test for publication bias. \emph{Biometrics}, \bold{50}, 1088--1101.

	Egger, M., Davey Smith, G., Schneider, M. & Minder, C. (1997) Bias in meta-analysis detected by a simple, graphical test. \emph{British Medical Journal}, \bold{315}, 629--634.

	Macaskill, P., Walter, S. D. & Irwing, L. (2001) A comparison of methods to detect publication bias in meta-analysis. \emph{Statistics in Medicine}, \bold{20}, 641--654.

	Moreno, S. G., Sutton, A. J., Ades, A. E., Stanley, T. D., Abrams, K. R., Peters, J. L. & Cooper, N. J. (2009) Assessment of regression-based methods to adjust for publication bias through a comprehensive simulation study. \emph{BMC Medical Research Methodology}, \bold{9}:2.

	Peters, J. L., Sutton, A. J., Jones, D. R., Abrams, K. R. & Rushton, L. (2006) Comparison of two methods to detect publication bias in meta-analysis. \emph{Journal of the American Medical Association}, \bold{295}, 676--680.

	Sterne, J. A. C. & Egger, M. (2005). Regression methods to detect publication and other bias in meta-analysis. In H. R. Rothstein, A. J. Sutton & M. Borenstein (Eds.) (2005) \emph{Publication bias in meta-analysis: Prevention, assessment, and adjustments} (pp. 99--110). Chichester, England: Wiley.
}
\seealso{
	\code{\link{funnel}}, \code{\link{fsn}}, \code{\link{trimfill}}
}
\examples{
### load BCG vaccine data
data(dat.bcg)

### calculate log risk rates and corresponding sampling variances
dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
dat.bcg$yi <- dat$yi
dat.bcg$vi <- dat$vi

ranktest(yi, vi, data=dat.bcg)

### Egger's regression test (multiplicative model)
summary(lm(yi ~ sqrt(vi), weights=1/vi, data=dat.bcg))

### Egger's regression test (fixed-effects model)
rma(yi, vi, mods=sqrt(vi), data=dat.bcg, method="FE")

### Egger's regression test (random-effects model)
rma(yi, vi, mods=sqrt(vi), data=dat.bcg, method="DL")

### Egger's test using the sampling variance as moderator
summary(lm(yi ~ vi, weights=1/vi, data=dat.bcg))
rma(yi, vi, mods=vi, data=dat.bcg, method="FE")
rma(yi, vi, mods=vi, data=dat.bcg, method="DL")

### add total sample size to data frame
dat.bcg$n.tot <- dat.bcg$tpos + dat.bcg$tneg + dat.bcg$cpos + dat.bcg$cneg

### Macaskill's test (using total sample size)
summary(lm(yi ~ n.tot, weights=1/vi, data=dat.bcg))

### Peters' test (using inverse of total sample size)
summary(lm(yi ~ I(1/n.tot), weights=1/vi, data=dat.bcg))

### note: one may consider using a pooled estimate of the sampling variance 
### to reduce the correlation between yi and vi; see Macaskill et al. (2005) 
### and Peters et al. (2006) for more details when using odds ratios as the
### outcome measure

### Egger's regression test using a mixed-effects model that also
### includes absolute latitude as a moderator
rma(yi, vi, mods=cbind(ablat, sei=sqrt(vi)), data=dat.bcg, method="DL")
}
\keyword{htest}
