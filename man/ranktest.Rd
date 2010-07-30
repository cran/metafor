\name{ranktest}
\alias{ranktest}
\title{Rank Correlation Test for Funnel Plot Asymmetry}
\description{
   The function \code{ranktest} is generic. It can be used to carry out the rank correlation test for funnel plot asymmetry as described by Begg and Mazumdar (1994).
}
\usage{
ranktest(x, \dots)
}
\arguments{
	\item{x}{an object of class \code{"rma"}.}
	\item{\dots}{other arguments.}
}
\details{
	Currently, there is only a method for handling objects of class \code{"rma"} with the \code{ranktest} function. Accordingly, the corresponding method is called \code{\link{ranktest.rma}}. See the documentation for that function for more details.
}
\author{Wolfgang Viechtbauer; \email{wvb@www.wvbauer.com}; \url{http://www.wvbauer.com/}}
\references{
	Begg, C. B. & Mazumdar, M. (1994). Operating characteristics of a rank correlation test for publication bias. \emph{Biometrics}, \bold{50}, 1088--1101.

	Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \url{http://www.jstatsoft.org/v36/i03/}.
}
\seealso{
	\code{\link{ranktest.rma}}
}
\keyword{htest}
