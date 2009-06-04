\name{fsn}
\alias{fsn}
\title{Fail-Safe N}
\description{Function to calculate Rosenthal's fail-safe N.}
\usage{
fsn(yi, vi, sei, data=NULL, subset=NULL, alpha=.05)
}
\arguments{
	\item{yi}{a vector with the observed effect sizes or outcomes.}
	\item{vi}{a vector with the corresponding sampling variances.}
	\item{sei}{a vector with the corresponding standard errors.\cr(note: only one of the two, \code{vi} or \code{sei}, needs to be specified)}
	\item{data}{an optional data frame containing the variables given to the arguments above.}
	\item{subset}{an optional vector indicating the subset of studies that should be used for the calculation. This can be a logical vector of length \eqn{k} or a numeric vector indicating the indices of the observations to include.}
	\item{alpha}{one-tailed alpha level (.05 by default).}
}
\value{
	The fail-safe N value.
}
\details{
	The function calculates the number of studies averaging null results that would have to be added to the given set of observed outcomes to reduce the combined (one-sided) significance level to a desired alpha level (usually .05). The calculation is based on Stouffer's method to combine p-values and is described in Rosenthal (1979).
}
\author{Wolfgang Viechtbauer; \email{wvb@www.wvbauer.com}; \url{http://www.wvbauer.com/}}
\references{
	Rosenthal, R. (1979) The "file drawer problem" and tolerance for null results. \emph{Psychological Bulletin}, \bold{86}, 638--641.
}
\seealso{
	\code{\link{ranktest}}, \code{\link{trimfill}}
}
\examples{
### load BCG vaccine data
data(dat.bcg)

### calculate log risk rates and corresponding sampling variances
dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
dat.bcg$yi <- dat$yi
dat.bcg$vi <- dat$vi

fsn(yi, vi, data=dat.bcg)
}
\keyword{htest}
