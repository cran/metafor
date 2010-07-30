\name{print.infl.rma.uni}
\alias{print.infl.rma.uni}
\title{Print Method for infl.rma.uni Objects}
\description{
   Print method for objects of class \code{"infl.rma.uni"}.
}
\usage{
\method{print}{infl.rma.uni}(x, digits=x$digits, \dots)
}
\arguments{
	\item{x}{an object of class \code{"infl.rma.uni"}.}
	\item{digits}{an integer specifying the number of decimal places to which the printed results should be rounded (the default is to take the value from the object).}
	\item{\dots}{other arguments.}
}
\value{
	See the documentation of the \code{\link{influence.rma.uni}} function for details on what is printed.
}
\author{Wolfgang Viechtbauer; \email{wvb@www.wvbauer.com}; \url{http://www.wvbauer.com/}}
\references{
	Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \url{http://www.jstatsoft.org/v36/i03/}.

	Viechtbauer, W. & Cheung, M. W.-L. (in press). Outlier and influence diagnostics for meta-analysis. \emph{Research Synthesis Methods}.
}
\seealso{
	\code{\link{influence.rma.uni}}
}
\keyword{print}
