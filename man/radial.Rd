\name{radial}
\alias{radial}
\alias{galbraith}
\title{Radial Plots}
\description{
   The function \code{radial} is generic. It can be used to create radial (also called Galbraith) plots.
}
\usage{
radial(x, \dots)
galbraith(x, \dots)
}
\arguments{
	\item{x}{an object of class \code{"rma"}.}
	\item{\dots}{other arguments.}
}
\details{
	Currently, there is only a method for handling objects of class \code{"rma"} with the \code{radial} function. Accordingly, the corresponding method is called \code{\link{radial.rma}}. See the documentation for that function for more details.
}
\author{Wolfgang Viechtbauer; \email{wvb@www.wvbauer.com}; \url{http://www.wvbauer.com/}}
\references{
	Galbraith, R. F. (1988). Graphical display of estimates having differing standard errors. \emph{Technometrics}, \bold{30}, 271--281.

	Galbraith, R. F. (1988). A note on graphical presentation of estimated odds ratios from several clinical trials. \emph{Statistics in Medicine}, \bold{7}, 889--894.

	Galbraith, R. F (1994). Some applications of radial plots. \emph{Journal of the American Statistical Association}, \bold{89}, 1232--1242.

	Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \url{http://www.jstatsoft.org/v36/i03/}.
}
\seealso{
	\code{\link{radial.rma}}
}
\keyword{hplot}
