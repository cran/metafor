\name{regtest}
\alias{regtest}
\title{Regression Tests for Funnel Plot Asymmetry}
\description{
   The function \code{regtest} is generic. It can be used to carry out various tests for funnel plot asymmetry, including Egger's regression test and variations thereof.
}
\usage{
regtest(x, \dots)
}
\arguments{
	\item{x}{an object of class \code{"rma"}.}
	\item{\dots}{other arguments.}
}
\details{
	Currently, there is only a method for handling objects of class \code{"rma"} with the \code{regtest} function. Accordingly, the corresponding method is called \code{\link{regtest.rma}}. See the documentation for that function for more details.
}
\author{Wolfgang Viechtbauer; \email{wvb@www.wvbauer.com}; \url{http://www.wvbauer.com/}}
\references{
	Begg, C. B. & Mazumdar, M. (1994). Operating characteristics of a rank correlation test for publication bias. \emph{Biometrics}, \bold{50}, 1088--1101.

	Egger, M., Davey Smith, G., Schneider, M. & Minder, C. (1997). Bias in meta-analysis detected by a simple, graphical test. \emph{British Medical Journal}, \bold{315}, 629--634.

	Harbord, R. M., Egger, M. & Sterne, J. A. C. (2006). A modified test for small-study effects in meta-analyses of controlled trials with binary endpoints. \emph{Statistics in Medicine}, \bold{25}, 3443--3457.

	Macaskill, P., Walter, S. D. & Irwig, L. (2001). A comparison of methods to detect publication bias in meta-analysis. \emph{Statistics in Medicine}, \bold{20}, 641--654.

	Moreno, S. G., Sutton, A. J., Ades, A. E., Stanley, T. D., Abrams, K. R., Peters, J. L. & Cooper, N. J. (2009). Assessment of regression-based methods to adjust for publication bias through a comprehensive simulation study. \emph{BMC Medical Research Methodology}, \bold{9}:2.

	Peters, J. L., Sutton, A. J., Jones, D. R., Abrams, K. R. & Rushton, L. (2006). Comparison of two methods to detect publication bias in meta-analysis. \emph{Journal of the American Medical Association}, \bold{295}, 676--680.

	Ruecker, G., Schwarzer, G. & Carpenter, J. (2008). Arcsine test for publication bias in meta-analyses with binary outcomes. \emph{Statistics in Medicine}, \bold{27}, 746--763.

	Sterne, J. A. C. & Egger, M. (2005). Regression methods to detect publication and other bias in meta-analysis. In H. R. Rothstein, A. J. Sutton & M. Borenstein (Eds.) (2005) \emph{Publication bias in meta-analysis: Prevention, assessment, and adjustments} (pp. 99--110). Chichester, England: Wiley.

	Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \url{http://www.jstatsoft.org/v36/i03/}.
}
\seealso{
	\code{\link{regtest.rma}}
}
\keyword{htest}
