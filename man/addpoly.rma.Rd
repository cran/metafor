\name{addpoly.rma}
\alias{addpoly.rma}
\title{Add Polygons to Forest Plots Based on 'rma' Objects}
\description{Function to add a polygon to a forest plot showing the summary estimate with correspondong confidence interval based on an object of class \code{"rma"}.}
\usage{
\method{addpoly}{rma}(x, row=-2, level=x$level, digits=2,
        annotate=TRUE, mlab, transf=FALSE, atransf=FALSE, targs,
        col="black", efac=1, cex, \dots)
}
\arguments{
   \item{x}{an object of class \code{"rma"}.}
   \item{row}{value specifying the row (or more generally, the horizontal position) for plotting the polygon (default is \code{-2}).}
   \item{level}{numerical value between 0 and 100 specifying the confidence interval level (the default is to take the value from the object).}
   \item{digits}{integer specifying the number of decimal places to which the annotations should be rounded (default is 2).}
   \item{annotate}{logical specifying whether annotations for the summary estimate should be added to the plot (default is \code{TRUE}).}
   \item{mlab}{optional character string giving a label for the summary estimate polygon. If unspecified, the function sets a default label.}
   \item{transf}{optional argument specifying the name of a function that should be used to transform the summary estimate and confidence interval bound. Defaults to \code{FALSE}, which means that no transformation is used.}
   \item{atransf}{optional argument specifying the name of a function that should be used to transform the annotations. Defaults to \code{FALSE}, which means that no transformation is used.}
   \item{targs}{optional arguments needed by the function specified via \code{transf} or \code{atransf}.}
   \item{col}{color of the polygon that is drawn (default is \code{"black"}).}
   \item{efac}{vertical expansion factor for the polygon. The default value of 1 should usually work okay.}
   \item{cex}{optional symbol expansion factor. If unspecified, the function tries to set this to a sensible value.}
   \item{\dots}{other arguments.}
}
\details{
   The function can be used to add a polygon to an existing forest plot created with the \code{\link{forest}} function. The polygon shows the summary estimate based on a fixed- or random-effects model. Using this function, summary estimates based on different types of models can be shown in the same plot. Also, summary estimates based on a subgrouping of the studies can be added to the plot this way. See examples below.

   The arguments \code{transf}, \code{atransf}, \code{efac}, and \code{cex} should always be set equal to the same values used to create the forest plot.
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
   \code{\link{forest.rma}}, \code{\link{forest.default}}
}
\examples{
### load BCG vaccine data
data(dat.bcg)

### meta-analysis of the log relative risks using the Mantel-Haenszel method
res <- rma.mh(ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg,
              slab=paste(author, year, sep=", "), measure="RR")

### forest plot of the observed relative risks with summary estimate
forest(res, atransf=exp, ylim=c(-2.5,16))

### meta-analysis of the log relative risks using a random-effects model
res <- rma(ai=tpos, bi=tneg, ci=cpos, di=cneg,
           data=dat.bcg, measure="RR", method="REML")

### add summary estimate from the random-effects model to forest plot
addpoly(res, atransf=exp)

### forest plot with subgrouping of studies and summaries per subgroup
res <- rma(ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg, measure="RR",
           slab=paste(author, year, sep=", "), method="REML")
forest(res, xlim=c(-16, 6), at=log(c(.05, .25, 1, 4)), atransf=exp,
       ilab=cbind(dat.bcg$tpos, dat.bcg$tneg, dat.bcg$cpos, dat.bcg$cneg),
       ilab.xpos=c(-9.5,-8,-6,-4.5), cex=.75, ylim=c(-1, 27),
       order=order(dat.bcg$alloc), rows=c(3:4,9:15,20:23),
       mlab="RE Model for All Studies")
op <- par(cex=.75, font=4)
text(-16, c(24,16,5), c("Systematic Allocation", "Random Allocation",
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
