\name{dat.collins1985a}
\docType{data}
\alias{dat.collins1985a}
\title{Studies on the Treatment of Upper Gastrointestinal Bleeding by a Histamine H2 Antagonist}
\description{Results from studies examining the presence of persistent or recurrent bleedings in patients receiving either a histamine H2 antagonist or placebo.}
\usage{dat.collins1985a}
\format{The data frame contains the following columns:
\tabular{lll}{
\bold{id}    \tab \code{numeric}   \tab study number \cr
\bold{ref}   \tab \code{numeric}   \tab reference number \cr
\bold{year}  \tab \code{numeric}   \tab year of publication \cr
\bold{nti}   \tab \code{numeric}   \tab number of patients in treatment group \cr
\bold{xti}   \tab \code{numeric}   \tab number of patients in treatment group with persistent or recurrent bleedings \cr
\bold{nci}   \tab \code{numeric}   \tab number of patients in placebo group \cr
\bold{xci}   \tab \code{numeric}   \tab number of patients in placebo group with persistent or recurrent bleedings
}
}
\details{
   The data contained in this dataset were obtained from Table I in van Houwelingen, Zwinderman, and Stijnen (1993).
}
\source{
   Collins, R., & Langman, M. (1985). Treatment with histamine H2 antagonists in acute upper gastrointestinal hemorrhage. \emph{New England Journal of Medicine}, \bold{313}, 660--666.

   van Houwelingen, H. C., Zwinderman, K. H., & Stijnen, T. (1993). A bivariate approach to meta-analysis. \emph{Statistics in Medicine}, \bold{12}, 2273--2284.
}
\references{
   Curtis, P. S., & Wang, X. (1998). A meta-analysis of elevated CO2 effects on woody plant mass, form, and physiology. \emph{Oecologia}, \bold{113}, 299--313.
}
\examples{
### load data
data(dat.collins1985a)

### calculate (log) odds ratio and sampling variance
dat <- escalc(measure="OR", ai=xci, n1i=nci, ci=xti, n2i=nti,
              data=dat.collins1985a, to="all")
summary(dat, digits=2, transf=exp)

### meta-analysis of log odds ratios using Peto's method
res <- rma.peto(ai=xci, n1i=nci, ci=xti, n2i=nti, data=dat)
summary(res)

\dontrun{
### meta-analysis of log odds ratios using conditional logistic regression model
res <- rma.glmm(measure="OR", ai=xci, n1i=nci, ci=xti, n2i=nti,
                data=dat, model="CM.EL", method="FE")
summary(res)

### plot the log likelihoods of the odds ratios
llplot(measure="OR", ai=xci, n1i=nci, ci=xti, n2i=nti, data=dat,
       lwd=1, refline=NA, xlim=c(-4,4), drop00=FALSE)

### meta-analysis of log odds ratios using conditional logistic regression model
res <- rma.glmm(measure="OR", ai=xci, n1i=nci, ci=xti, n2i=nti,
                data=dat, model="CM.EL", method="ML")
summary(res)
}
}
\keyword{datasets}
