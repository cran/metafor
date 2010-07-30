\name{funnel}
\alias{funnel}
\title{Funnel Plots}
\description{
   The function \code{funnel} is generic. It can be used to create funnel plots.
}
\usage{
funnel(x, \dots)
}
\arguments{
	\item{x}{an object of class \code{"rma"}.}
	\item{\dots}{other arguments.}
}
\details{
	Currently, there is only a method for handling objects of class \code{"rma"} with the \code{funnel} function. Accordingly, the corresponding method is called \code{\link{funnel.rma}}. See the documentation for that function for more details.
}
\author{Wolfgang Viechtbauer; \email{wvb@www.wvbauer.com}; \url{http://www.wvbauer.com/}}
\references{
	Light, R. J. & Pillemer, D. B. (1984). \emph{Summing up: The science of reviewing research}. Cambridge, MA: Harvard University Press.

	Sterne, J. A. C. & Egger, M. (2001). Funnel plots for detecting bias in meta-analysis: Guidelines on choice of axis. \emph{Journal of Clinical Epidemiology}, \bold{54}, 1046--1055.

	Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \url{http://www.jstatsoft.org/v36/i03/}.
}
\seealso{
	\code{\link{funnel.rma}}
}
\keyword{hplot}
