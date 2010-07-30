\name{weights.rma.uni}
\alias{weights.rma.uni}
\alias{weights.rma.mh}
\alias{weights.rma.peto}
\title{Weights for rma Objects}
\description{
   The function extracts the weights (in \%) given to the observed effects or outcomes during the model fitting from objects of class \code{"rma.uni"}, \code{"rma.mh"}, and \code{"rma.peto"}.
}
\usage{
\method{weights}{rma.uni}(object, \dots)
\method{weights}{rma.mh}(object, \dots)
\method{weights}{rma.peto}(object, \dots)
}
\arguments{
	\item{object}{an object of class \code{"rma.uni"}, \code{"rma.mh"}, or \code{"rma.peto"}.}
	\item{\dots}{other arguments.}
}
\value{
	A vector with the weights.
}
\author{Wolfgang Viechtbauer; \email{wvb@www.wvbauer.com}; \url{http://www.wvbauer.com/}}
\references{
	Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \url{http://www.jstatsoft.org/v36/i03/}.
}
\seealso{
	\code{\link{rma.uni}}, \code{\link{rma.mh}}, \code{\link{rma.peto}}, \code{\link{influence.rma.uni}}
}
\examples{
### load BCG vaccine data
data(dat.bcg)

### meta-analysis of the log relative risks using a mixed-effects model
### with two moderators (absolute latitude and publication year)
res <- rma(ai=tpos, bi=tneg, ci=cpos, di=cneg, mods=cbind(ablat, year),
           data=dat.bcg, measure="RR", method="REML")
weights(res)
}
\keyword{models}
