\name{fitstats}
\alias{fitstats}
\title{Extract Fit Statistics and Information Criteria}
\description{
   The function \code{fitstats} is generic. It can be used to extract the log likelihood, deviance, and information criteria (e.g., AIC and BIC) for suitable objects.
}
\usage{
fitstats(object, \dots)
}
\arguments{
   \item{object}{an object for which fit statistics and information criteria can be calculated. See \sQuote{Details}.}
   \item{\dots}{other arguments.}
}
\details{
   Currently, there is only a method for handling objects of class \code{"rma"} with the \code{fitstats} function. Accordingly, the corresponding method is called \code{\link{fitstats.rma}}. See the documentation for that function for more details.
}
\value{
   Log likelihood, deviance, AIC, and BIC values either under the regular or restricted likelihood.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \cr
   package homepage: \url{http://www.metafor-project.org/} \cr
   author homepage: \url{http://www.wvbauer.com/}
}
\references{
   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \url{http://www.jstatsoft.org/v36/i03/}.
}
\seealso{
   \code{\link{fitstats.rma}}
}
\keyword{methods}
