\name{dat.berkey1998}
\docType{data}
\alias{dat.berkey1998}
\title{Studies on Treatments for Periodontal Disease}
\description{Results from 5 trials comparing surgical and non-surgical treatments for medium-severity periodontal disease one year after treatment.}
\usage{dat.berkey1998}
\format{The data frame contains the following columns:
\tabular{lll}{
\bold{trial}      \tab \code{numeric}    \tab trial number \cr
\bold{author}     \tab \code{character}  \tab study author(s) \cr
\bold{year}       \tab \code{numeric}    \tab publication year \cr
\bold{ni}         \tab \code{numeric}    \tab number of patients \cr
\bold{outcome}    \tab \code{character}  \tab patient outcome (PD = probing depth; AL = attachment level) \cr
\bold{yi}         \tab \code{numeric}    \tab observed mean difference in outcome (surgical versus non-surgical) \cr
\bold{v1i}        \tab \code{numeric}    \tab variance or covariance between observed effects \cr
\bold{v2i}        \tab \code{numeric}    \tab variance or covariance between observed effects
}
}
\details{
   The dataset includes the results from 5 trials that compared surgical and non-surgical methods for the treatment of medium-severity periodontal disease. Reported outcomes include the probing depth (PD) and attachment level (AL) one year after the treatment. The effect size measure used for this meta-analysis was the (raw) mean difference, calculated in such a way that positive values indicate that surgery was more effective than non-surgical treatment (so, the results from the various trials indicate that surgery is preferable for reducing the probing depth, while non-surgical treatment is preferable for increasing the attachment level). Since each trial provides effect size estimates for both outcomes, the estimates are correlated. A multivariate model can be used to meta-analyze the two outcomes simultaneously.
}
\source{
   Berkey, C. S., Antczak-Bouckoms, A., Hoaglin, D. C., Mosteller, F., & Pihlstrom, B. L. (1995). Multiple-outcomes meta-analysis of treatments for periodontal disease. \emph{Journal of Dental Research}, \bold{74}, 1030--1039.

   Berkey, C. S., Hoaglin, D. C., Antczak-Bouckoms, A., Mosteller, F., & Colditz, G. A. (1998). Meta-analysis of multiple outcomes by regression with random effects. \emph{Statistics in Medicine}, \bold{17}, 2537--2550.
}
\examples{
### load data
dat <- get(data(dat.berkey1998))

### for bdiag() function
require(Matrix)

### construct variance-covariance matrix of the observed outcomes
V <- as.matrix(bdiag(lapply(split(dat[,c("v1i", "v2i")], dat$trial), as.matrix)))

### fit multiple outcomes (meta-regression) model (with ML estimation)
res <- rma.mv(yi, V, mods = ~ outcome + outcome:I(year - 1983) - 1, 
              random = ~ outcome | trial, struct="UN", data=dat, method="ML")
print(res, digits=3)
}
\keyword{datasets}
