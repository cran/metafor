\name{logLik.rma}
\alias{logLik.rma}
\title{Log-Likelihood of rma Objects}
\description{
   The function extracts the (restricted) log likelihood for objects of class \code{"rma"}.
}
\usage{
\method{logLik}{rma}(object, REML, \dots)
}
\arguments{
	\item{object}{an object of class \code{"rma"}.}
	\item{REML}{an optional logical value. If \code{TRUE}, the restricted log-likelihood is returned, else, if \code{FALSE}, the log-likelihood is returned. Defaults to the method of estimation used, that is \code{TRUE} if \code{object} was fitted with \code{method = "REML"} and \code{FALSE} otherwise.}
	\item{\dots}{other arguments.}
}
\value{
	An object of class \code{"logLik"}. The object provides the (restricted) log likelihood of the model evaluated at the estimated coefficient(s).
}
\note{
	For a more complete set of fit statistics and information criteria, see \code{\link{fitstats.rma}}.
}
\author{Wolfgang Viechtbauer; \email{wvb@www.wvbauer.com}; \url{http://www.wvbauer.com/}}
\seealso{
	\code{\link{fitstats.rma}}, \code{\link{rma.uni}}, \code{\link{rma.mh}}, \code{\link{rma.peto}}
}
\examples{
### load BCG vaccine data
data(dat.bcg)

### meta-analysis of the log relative risks using a mixed-effects model
### with two moderators (absolute latitude and publication year)
res <- rma(ai=tpos, bi=tneg, ci=cpos, di=cneg, mods=cbind(ablat, year),
           data=dat.bcg, measure="RR", method="REML")
logLik(res)
}
\keyword{models}
