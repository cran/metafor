\name{fitstats.rma}
\alias{fitstats.rma}
\alias{logLik.rma}
\alias{deviance.rma}
\alias{AIC.rma}
\alias{BIC.rma}
\title{Fit Statistics and Information Criteria for rma Objects}
\description{
	Functions to extract the log likelihood, deviance, AIC, and BIC values from objects of class \code{"rma"}.
}
\usage{
\method{fitstats}{rma}(object, REML=NULL, \dots)
\method{logLik}{rma}(object, REML, \dots)
\method{deviance}{rma}(object, REML, \dots)
\method{AIC}{rma}(object, \dots, k=2)
\method{BIC}{rma}(object, \dots)
}
\arguments{
	\item{object}{an object of class \code{"rma"}.}
	\item{REML}{logical indicating whether the regular or restricted likelihood function should be used to obtain the fit statistics and information criteria. Defaults to the method of estimation used, that is \code{TRUE} if \code{object} was fitted with \code{method = "REML"} and \code{FALSE} otherwise.}
	\item{k}{numeric specifying the penalty per parameter to be used. The default \code{k=2} is the classical AIC. See \code{\link{AIC}} for more details.}
	\item{\dots}{other arguments.}
}
\value{
	For \code{fitstats.rma}, a column vector with the (restricted) log likelihood, deviance, AIC, and BIC values.

	For \code{logLik.rma}, an object of class \code{"logLik"}, providing the (restricted) log likelihood of the model evaluated at the estimated coefficient(s).

	For \code{deviance.rma}, a numeric value with the corresponding deviance.

	For \code{AIC.rma}, a numeric value with the corresponding AIC.

	For \code{BIC.rma}, a numeric value with the corresponding BIC.
}
\note{
	When a random/mixed-effects model is fitted, \eqn{\tau^2}{tau^2} is counted as an additional parameter in the calculation of the AIC and BIC.
}
\author{Wolfgang Viechtbauer; \email{wvb@www.wvbauer.com}; \url{http://www.wvbauer.com/}}
\references{
	Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \url{http://www.jstatsoft.org/v36/i03/}.
}
\seealso{
	\code{\link{rma.uni}}, \code{\link{rma.mh}}, \code{\link{rma.peto}}, \code{\link{anova.rma.uni}}
}
\examples{
### load BCG vaccine data
data(dat.bcg)

### meta-analysis of the log relative risks using a random-effects model
res1 <- rma(ai=tpos, bi=tneg, ci=cpos, di=cneg, 
            data=dat.bcg, measure="RR", method="ML")

### mixed-effects model with two moderators (latitude and publication year)
res2 <- rma(ai=tpos, bi=tneg, ci=cpos, di=cneg, mods=cbind(ablat, year),
            data=dat.bcg, measure="RR", method="ML")

fitstats(res1)
fitstats(res2)

logLik(res1)
logLik(res2)

deviance(res1)
deviance(res2)

AIC(res1)
AIC(res2)

BIC(res1)
BIC(res2)
}
\keyword{models}
