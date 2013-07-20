\name{baujat}
\alias{baujat}
\title{Baujat Plots}
\description{
   The function \code{baujat} is generic. It can be used to create Baujat plots.
}
\usage{
baujat(x, \dots)
}
\arguments{
   \item{x}{an object of class \code{"rma.uni"}, \code{"rma.mh"}, or \code{"rma.peto"}.}
   \item{\dots}{other arguments.}
}
\details{
   Currently, there are methods for handling objects of class \code{"rma.uni"}, \code{"rma.mh"}, and \code{"rma.peto"} with the \code{baujat} function. Accordingly, the corresponding methods are called \code{\link{baujat.rma.uni}}, \code{\link{baujat.rma.mh}}, and \code{\link{baujat.rma.peto}}. See the documentation for those functions for more details.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \cr
   package homepage: \url{http://www.metafor-project.org/} \cr
   author homepage: \url{http://www.wvbauer.com/}
}
\references{
   Baujat, B., Mahe, C., Pignon, J.-P., & Hill, C. (2002). A graphical method for exploring heterogeneity in meta-analyses: Application to a meta-analysis of 65 trials. \emph{Statistics in Medicine}, \bold{21}(18), 2641--2652.

   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \url{http://www.jstatsoft.org/v36/i03/}.
}
\seealso{
   \code{\link{baujat.rma.uni}}, \code{\link{baujat.rma.mh}}, \code{\link{baujat.rma.peto}}
}
\keyword{hplot}
