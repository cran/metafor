\name{print.predict.rma}
\alias{print.predict.rma}
\title{Print Method for predict.rma Objects}
\description{
   Print method for objects of class \code{"predict.rma"}.
}
\usage{
\method{print}{predict.rma}(x, digits=x$digits, \dots)
}
\arguments{
	\item{x}{an object of class \code{"predict.rma"}.}
	\item{digits}{an integer specifying the number of decimal places to which the printed results should be rounded (the default is to take the value from the object).}
	\item{\dots}{other arguments.}
}
\value{
	See the documentation of the \code{\link{predict.rma}} function for details on what is printed.
}
\author{Wolfgang Viechtbauer; \email{wvb@www.wvbauer.com}; \url{http://www.wvbauer.com/}}
\references{
	Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \url{http://www.jstatsoft.org/v36/i03/}.
}
\keyword{print}
