\name{dat.warfarin}
\docType{data}
\alias{dat.warfarin}
\title{Studies on the Effectiveness of Warfarin for Preventing Strokes}
\description{Results from 6 clinical trials examining the effectiveness of adjusted-dose warfarin for preventing strokes in patients with atrial fibrillation.}
\usage{dat.warfarin}
\format{The data frame contains the following columns:
\tabular{lll}{
\bold{trial}    \tab \code{numeric}   \tab trial number \cr
\bold{study}    \tab \code{character} \tab study name (abbreviated) \cr
\bold{year}     \tab \code{numeric}   \tab publication year \cr
\bold{x1i}      \tab \code{numeric}   \tab number of strokes in the warfarin group \cr
\bold{n1i}      \tab \code{numeric}   \tab number of patients in the warfarin group \cr
\bold{t1i}      \tab \code{numeric}   \tab total person-time (years) in the warfarin group \cr
\bold{x2i}      \tab \code{numeric}   \tab number of strokes in the placebo/control group \cr
\bold{n2i}      \tab \code{numeric}   \tab number of patients in the placebo/control group \cr
\bold{t2i}      \tab \code{numeric}   \tab total person-time (years) in the placebo/control group \cr
\bold{compgrp}  \tab \code{character} \tab type of comparison group (placebo or control) \cr
\bold{prevtype} \tab \code{character} \tab type of prevention (primary or secondary) \cr
\bold{trinr}    \tab \code{character} \tab target range for the international normalized ratio (INR) 
}
}
\details{
	The 6 studies provide data with respect to the number of strokes in the warfarin and the comparison (placebo or control) group. In addition, the number of patients and the total person-time (in years) is provided for the two groups. The goal of the meta-analysis was to examine the effectiveness of adjusted-dose warfarin for preventing strokes in patients with atrial fibrillation.
}
\source{
	Hart, R. G., Benavente, O., McBride, R., & Pearce, L. A. (1999). Antithrombotic therapy to prevent stroke in patients with atrial fibrillation: A meta-analysis. \emph{Annals of Internal Medicine}, \bold{131}, 492--501.
}
\examples{
### load warfarin data
data(dat.warfarin)

### calculate log incidence rate ratios and corresponding sampling variances
dat <- escalc(x1i=x1i, x2i=x2i, t1i=t1i, t2i=t2i, measure="IRR", 
              data=dat.warfarin, append=TRUE)
dat

### meta-analysis of log incidence rate ratios using a random-effects model
res <- rma(yi, vi, data=dat)
res

### average incidence rate ratio with 95\% CI
predict(res, transf=exp)

### forest plot with extra annotations
forest(res, xlim=c(-14, 6), at=log(c(.05, .25, 1, 4)), atransf=exp, 
       slab=paste(dat$study, " (", dat$year, ")", sep=""), 
       ilab=cbind(paste(dat$x1i, "/", dat$t1i, sep=" "), 
       paste(dat$x2i, "/", dat$t2i, sep=" ")), 
       ilab.xpos=c(-8,-5), cex=.85)
op <- par(cex=.85, font=2)
text(-14, 7.5, "Study (Year)", pos=4)
text(6,   7.5, "IRR [95\% CI]", pos=2)
text(c(-8,-5), 8.0, c("Strokes /", "Strokes /"))
text(c(-8,-5), 7.5, c("Person-Time", "Person-Time"))
text(c(-8,-5), 8.5, c("Warfarin", "Control"))
segments(x0=-9, y0=8.25, x1=-4, y1=8.25)
par(op)

### meta-analysis of incidence rate differences using a random-effects model
res <- rma(x1i=x1i, x2i=x2i, t1i=t1i, t2i=t2i, measure="IRD", data=dat.warfarin)
res
}
\keyword{datasets}
