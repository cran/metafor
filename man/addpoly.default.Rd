\name{addpoly.default}
\alias{addpoly.default}
\title{Add Polygons to Forest Plot}
\description{Function to add one or more polygons to a forest plot.}
\usage{
\method{addpoly}{default}(x, vi, sei, rows=-1, level=95, digits=2, annotate=TRUE, 
        mlab=NULL, transf=FALSE, atransf=FALSE, targs=NULL, 
        col="black", efac=1, cex=NULL, \dots)
}
\arguments{
	\item{x}{a vector with the values at which the polygons should be drawn.}
	\item{vi}{a vector with the corresponding variances.}
	\item{sei}{a vector with the corresponding standard errors.\cr(note: only one of the two, \code{vi} or \code{sei}, needs to be specified)}
	\item{rows}{a vector specifying the rows (or more generally, the horizontal positions) for plotting the polygons (defaults is \code{-1}). Can also be a single value specifying the row (horizontal position) of the first polygon (the remaining polygons are then plotted below this starting row).}
	\item{level}{a numerical value between 0 and 100 specifying the confidence interval level (default is 95).}
	\item{digits}{integer value specifying the number of decimal places to which the annotations should be rounded (default is 2).}
	\item{annotate}{logical specifying whether annotations should be added to the plot for the polygons that are drawn (default is \code{TRUE}).}
	\item{mlab}{optional character vector with the same length as \code{x} giving labels for the polygons that are drawn. Defaults to \code{NULL}, which suppresses labels.}
	\item{transf}{an optional argument specifying the name of a function that should be used to transform the \code{x} values and confidence interval bounds. Defaults to \code{FALSE}, which means that no transformation is used.}
	\item{atransf}{an optional argument specifying the name of a function that should be used to transform the annotations. Defaults to \code{FALSE}, which means that no transformation is used.}
	\item{targs}{optional arguments needed by the function specified via \code{transf} or \code{atransf}.}
	\item{col}{color of the polygons that are drawn (default is \code{"black"}).}
	\item{efac}{vertical expansion factor for the polygons. The default value of 1 should usually work okay.}
	\item{cex}{an optional symbol expansion factor. If \code{NULL} (default), the function tries to set this to a sensible value.}
	\item{\dots}{other arguments.}
}
\details{
	The function can be used to add one or more polygons to an existing forest plot created with the \code{\link{forest}} function. For example, summary estimates based on a subgrouping of the studies or from models involving moderators can be added to the plot this way. See examples below.

	The arguments \code{transf}, \code{atransf}, \code{efac}, and \code{cex} should always be set equal to the same values used to create the forest plot.
}
\author{Wolfgang Viechtbauer; \email{wvb@www.wvbauer.com}; \url{http://www.wvbauer.com/}}
\references{
	Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \url{http://www.jstatsoft.org/v36/i03/}.
}
\seealso{
	\code{\link{forest.rma}}, \code{\link{forest.default}}
}
\examples{
### load BCG vaccine data
data(dat.bcg)

### meta-analysis of the log relative risks using a mixed-
### effects model with absolute latitude as a moderator
res <- rma(ai=tpos, bi=tneg, ci=cpos, di=cneg, mods=cbind(ablat),
           slab=paste(author, year, sep=", "), data=dat.bcg, 
           measure="RR", method="REML")

### forest plot of the observed relative risks
forest(res, addfit=FALSE, atransf=exp, ylim=c(-3.5,16))

### predicted log average relative risks for 10, 30, and 50 degrees absolute latitude
x <- predict(res, newmods=c(10, 30, 50))

### add predicted average relative risks to forest plot
addpoly(x$pred, sei=x$se, atransf=exp, 
        mlab=c("10 Degrees", "30 Degrees", "50 Degrees"))
abline(h=0)

### forest plot with subgrouping of studies and summaries per subgroup
res <- rma(ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg, measure="RR", 
           slab=paste(author, year, sep=", "), method="REML")
forest(res, slab=paste(dat.bcg$author, dat.bcg$year, sep=", "), 
       xlim=c(-16, 6), at=log(c(.05, .25, 1, 4)), atransf=exp,
       ilab=cbind(dat.bcg$tpos, dat.bcg$tneg, dat.bcg$cpos, dat.bcg$cneg), 
       ilab.xpos=c(-9.5,-8,-6,-4.5), cex=.75, ylim=c(-1, 27), 
       order=order(dat.bcg$alloc), rows=c(3:4,9:15,20:23),
       mlab="RE Model for All Studies")
op <- par(cex=.75, font=4)
text(-16, c(24,16,5), c("Systematic Allcoation", "Random Allocation", 
                        "Alternate Allocation"), pos=4)
par(font=2)
text(c(-9.5,-8,-6,-4.5), 26, c("TB+", "TB-", "TB+", "TB-"))
text(c(-8.75,-5.25),     27, c("Vaccinated", "Control"))
text(-16,                26, "Author(s) and Year",     pos=4)
text(6,                  26, "Relative Risk [95\% CI]", pos=2)
par(op)
res <- rma(ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg, measure="RR", 
           subset=(alloc=="systematic"), method="REML")
addpoly(res, row=18.5, cex=.75, atransf=exp, mlab="RE Model for Subgroup")
res <- rma(ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg, measure="RR", 
           subset=(alloc=="random"), method="REML")
addpoly(res, row=7.5, cex=.75, atransf=exp, mlab="RE Model for Subgroup")
res <- rma(ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg, measure="RR", 
           subset=(alloc=="alternate"), method="REML")
addpoly(res, row=1.5, cex=.75, atransf=exp, mlab="RE Model for Subgroup")
}
\keyword{aplot}
