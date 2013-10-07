\name{print.hc.rma.uni}
\alias{print.hc.rma.uni}
\title{Print Method for 'hc.rma.uni' Objects}
\description{
   Print method for objects of class \code{"hc.rma.uni"}.
}
\usage{
\method{print}{hc.rma.uni}(x, digits=x$digits, \dots)
}
\arguments{
   \item{x}{an object of class \code{"hc.rma.uni"}.}
   \item{digits}{integer specifying the number of decimal places to which the printed results should be rounded (the default is to take the value from the object).}
   \item{\dots}{other arguments.}
}
\details{
   The output is a data frame with two rows, the first (labeled \code{rma}) corresponding to the results based on the usual estimation method, the second (labeled \code{hc}) corresponding to the results based on the method by Henmi and Copas (2010). The data frame includes the following variables:
   \itemize{
   \item the method used to estimate \eqn{\tau^2} (always \code{DL} for \code{hc})
   \item the estimated amount of heterogeneity
   \item the estimated average true effect
   \item the corresponding standard errror (\code{NA} when \code{transf} argument has been used)
   \item the lower and upper confidence interval bounds
   }
}
\value{
   The function returns the data frame invisibly.
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
   \code{\link{hc.rma.uni}}
}
\keyword{print}