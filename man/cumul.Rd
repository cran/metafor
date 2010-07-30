\name{cumul}
\alias{cumul}
\title{Cumulative Meta-Analysis}
\description{
   The function \code{cumul} is generic. For suitable model objects, it repeatedly fits the model, adding one observation at a time to the model.
}
\usage{
cumul(x, \dots)
}
\arguments{
	\item{x}{an object of class \code{"rma.uni"}, \code{"rma.mh"}, or \code{"rma.peto"}.}
	\item{\dots}{other arguments.}
}
\details{
	Currently, there are methods for handling objects of class \code{"rma.uni"}, \code{"rma.mh"}, and \code{"rma.peto"} with the \code{cumul} function. Accordingly, the corresponding methods are called \code{\link{cumul.rma.uni}}, \code{\link{cumul.rma.mh}}, and \code{\link{cumul.rma.peto}}. See the documentation for those functions for more details.
}
\author{Wolfgang Viechtbauer; \email{wvb@www.wvbauer.com}; \url{http://www.wvbauer.com/}}
\references{
	Chalmers, T. C. & Lau, J. (1993). Meta-analytic stimulus for changes in clinical trials. \emph{Statistical Methods in Medical Research}, \bold{2}, 161--172.

	Lau, J., Schmid, C. H. & Chalmers, T. C. (1995). Cumulative meta-analysis of clinical trials builds evidence for exemplary medical care. \emph{Journal of Clinical Epidemiology}, \bold{48}, 45--57.

	Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \url{http://www.jstatsoft.org/v36/i03/}.
}
\seealso{
	\code{\link{cumul.rma.uni}}, \code{\link{cumul.rma.mh}}, \code{\link{cumul.rma.peto}}
}
\keyword{methods}
