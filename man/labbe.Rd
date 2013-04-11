\name{labbe}
\alias{labbe}
\title{L'Abbe Plots}
\description{
   The function \code{labbe} is generic. It can be used to create L'Abbe plots.
}
\usage{
labbe(x, \dots)
}
\arguments{
   \item{x}{an object of class \code{"rma"}.}
   \item{\dots}{other arguments.}
}
\details{
   Currently, there is only a method for handling objects of class \code{"rma"} with the \code{labbe} function. Accordingly, the corresponding method is called \code{\link{labbe.rma}}. See the documentation for that function for more details.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \cr
   package homepage: \url{http://www.metafor-project.org/} \cr
   author homepage: \url{http://www.wvbauer.com/}
}
\references{
   L'Abbe, K. A., Detsky, A. S., & O'Rourke, K. (1987). Meta-analysis in clinical research. \emph{Annals of Internal Medicine}, \bold{107}, 224--233.

   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \url{http://www.jstatsoft.org/v36/i03/}.
}
\seealso{
   \code{\link{labbe.rma}}
}
\keyword{hplot}
