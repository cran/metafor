\name{forest}
\alias{forest}
\title{Forest Plots}
\description{
   The function \code{forest} is generic. It can be used to create forest plots.
}
\usage{
forest(x, ...)
}
\arguments{
	\item{x}{either an object of class \code{"rma"} or a vector with the observed effect size or outcomes. See \sQuote{Details}.}
	\item{...}{other arguments.}
}
\details{
	Currently, methods exist for two types of situations.

	In the first case, object \code{x} is a fitted model object coming from the \code{\link{rma.uni}}, \code{\link{rma.mh}}, or \code{\link{rma.peto}} functions. The corresponding method is then \code{\link{forest.rma}}.

	Alternatively, object \code{x} can be a vector with observed effect size or outcomes. The corresponding method is then \code{\link{forest.default}}.
}
\author{Wolfgang Viechtbauer; \email{wvb@www.wvbauer.com}; \url{http://www.wvbauer.com/}}
\seealso{
	\code{\link{forest.rma}}, \code{\link{forest.default}}
}
\keyword{hplot}
