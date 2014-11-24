\name{baujat}
\alias{baujat}
\alias{baujat.rma.uni}
\alias{baujat.rma.mh}
\alias{baujat.rma.peto}
\title{Baujat Plots for 'rma' Objects}
\description{Function to create Baujat plots for objects of class \code{"rma.uni"}, \code{"rma.mh"}, and \code{"rma.peto"}.}
\usage{
baujat(x, \dots)

\method{baujat}{rma.uni}(x, xlim, ylim, xlab, ylab, cex, grid=TRUE, \ldots)

\method{baujat}{rma.mh}(x, xlim, ylim, xlab, ylab, cex, grid=TRUE, \ldots)

\method{baujat}{rma.peto}(x, xlim, ylim, xlab, ylab, cex, grid=TRUE, \ldots)
}
\arguments{
   \item{x}{an object of class \code{"rma.uni"}, \code{"rma.mh"}, or \code{"rma.peto"}.}
   \item{xlim}{x-axis limits. If unspecified, the function tries to set the x-axis limits to some sensible values.}
   \item{ylim}{y-axis limits. If unspecified, the function tries to set the y-axis limits to some sensible values.}
   \item{xlab}{title for the x-axis. If unspecified, the function tries to set an appropriate axis title.}
   \item{ylab}{title for the y-axis. If unspecified, the function tries to set an appropriate axis title.}
   \item{cex}{optional character expansion factor. If unspecified, the function tries to set this to a sensible value.}
   \item{grid}{logical indicating whether a grid should be added to the plot.}
   \item{\dots}{other arguments.}
}
\details{
   Baujat et al. (2002) proposed a diagnostic plot to detect sources of heterogeneity in meta-analytic data. The plot shows the contribution of each study to the overall Q-test statistic for heterogeneity on the x-axis versus the influence of each study (defined as the standardized squared difference between the overall estimate based on a fixed-effects model with and without the study included in the model fitting) on the y-axis. The same type of plot can be produced by first fitting a fixed-effects model with either the \code{\link{rma.uni}} (using \code{method="FE"}), \code{\link{rma.mh}}, or \code{\link{rma.peto}} functions and then passing the fitted model object to the \code{baujat} function.

   For models fitted with the \code{\link{rma.uni}} function (which may involve moderators and/or may be random/mixed-effects models), the idea underlying this type of plot can be generalized as follows: The x-axis then corresponds to the squared Pearson residual of a study, while the y-axis corresponds to the standardized squared difference between the predicted/fitted value for the study with and without the study included in the model fitting. Therefore, for a fixed-effect with moderators model, the x-axis corresponds to the contribution of the study to the QE-test statistic for residual heterogeneity.

   The points are labeled according to the \code{slab} argument (see model fitting functions). If \code{slab} was unspecified, then the numbers 1 through \eqn{k} are used by default.
}
\value{
   A data frame with components:
   \item{x}{the x coordinates of the points that were plotted.}
   \item{y}{the y coordinates of the points that were plotted.}
   \item{slab}{study labels of the points that were plotted.}
   Note that the data frame is returned invisibly.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \cr
   package website: \url{http://www.metafor-project.org/} \cr
   author homepage: \url{http://www.wvbauer.com/}
}
\references{
   Baujat, B., Mahe, C., Pignon, J.-P., & Hill, C. (2002). A graphical method for exploring heterogeneity in meta-analyses: Application to a meta-analysis of 65 trials. \emph{Statistics in Medicine}, \bold{21}(18), 2641--2652.

   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \url{http://www.jstatsoft.org/v36/i03/}.
}
\seealso{
   \code{\link{rma.uni}}, \code{\link{rma.mh}}, \code{\link{rma.peto}}, \code{\link{influence.rma.uni}}, \code{\link{radial}}
}
\examples{
### load data from Pignon et al. (2000)
dat <- get(data(dat.pignon2000))

### compute estimated log hazard ratios and sampling variances
dat$yi <- with(dat, OmE/V)
dat$vi <- with(dat, 1/V)

### meta-analysis based on all 65 trials
res <- rma(yi, vi, data=dat, method="FE", slab=id)

### create Baujat plot
baujat(res)
}
\keyword{hplot}
