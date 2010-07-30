\name{escalc}
\alias{escalc}
\alias{escalc.default}
\alias{escalc.formula}
\title{Calculate Effect Size and Outcome Measures}
\description{
   The function can be used to calculate various effect size or outcome measures (and the corresponding sampling variances) that are commonly used in meta-analyses.
}
\usage{
escalc(measure, formula, \dots)

\method{escalc}{default}(measure, formula, ai, bi, ci, di, n1i, n2i, 
       m1i, m2i, sd1i, sd2i, xi, mi, ri, ni, data, 
       add=1/2, to="only0", vtype="LS", append=FALSE, ...)

\method{escalc}{formula}(measure, formula, weights, data, 
       add=1/2, to="only0", vtype="LS", ...)
}
\arguments{
	\item{measure}{a character string indicating which effect size or outcome measure should be calculated. See \sQuote{Details} for possible options and how the data should then be specified.}
	\item{formula}{when using the formula interface of the function (see \sQuote{Details} below), a model formula specifying the data structure should be specified via this argument. When not using the formula interface, this argument can be ignored and the data requierd to calculate the effect sizes or outcomes are then passed to the function via the following set of argument. See \sQuote{Details} below.}
	\item{weights}{vector of weights to specify the group sizes or cell frequencies (only needed when using the formula interface). See \sQuote{Details} below.}
	\item{ai}{vector to specify 2x2 table frequencies (upper left cell).}
	\item{bi}{vector to specify 2x2 table frequencies (upper right cell).}
	\item{ci}{vector to specify 2x2 table frequencies (lower left cell).}
	\item{di}{vector to specify 2x2 table frequencies (lower right cell).}
	\item{n1i}{vector to specify group sizes or row total (first group/row).}
	\item{n2i}{vector to specify group sizes or row total (second group/row).}
	\item{m1i}{vector to specify means (first group).}
	\item{m2i}{vector to specify means (second group).}
	\item{sd1i}{vector to specify standard deviations (first group).}
	\item{sd2i}{vector to specify standard deviations (second group).}
	\item{xi}{vector to specify the frequencies of the event of interest.}
	\item{mi}{vector to specify the frequencies of the complement of the event of interest.}
	\item{ri}{vector to specify the raw correlation coefficients.}
	\item{ni}{vector to specify the sample sizes.}
	\item{data}{an optional data frame containing the variables given to the arguments above.}
	\item{add}{See \sQuote{Details}.}
	\item{to}{See \sQuote{Details}.}
	\item{vtype}{See \sQuote{Details}.}
	\item{append}{logical indicating whether the data frame specified via the \code{data} argument (if one has been specified) should be returned together with the effect sizes and sampling variances (default is \code{FALSE}).}
	\item{\dots}{other arguments.}
}
\details{
	There are two interfaces to using the \code{escalc} function, the default and a formula interface. The two interfaces are described below.
	
	\subsection{Default Interface}{

		The default interface works as follows. The argument \code{measure} is a character string specifying which outcome measure should be calculated (see below for the various options), arguments \code{ai} through \code{ni} are then used to supply the needed information to calculate the various measures (depending on the outcome measure, different arguments need to be supplied), and \code{data} can be used to specify a data frame containing the variables given to the previous arguments. The \code{add} and \code{to} arguments may be needed when dealing with 2x2 table data that contain cells with zeros. Finally, the \code{vtype} argument is used to specify how to calculate the sampling variance estimate (see below).

		\subsection{Effect Size and Outcome Measures for 2x2 Table Data}{

			Meta-analyses in the health/medical sciences are often based on studies providing data in terms of 2x2 tables. In particular, assume that we have \eqn{k} tables of the form:
			\tabular{lccc}{
					  \tab outcome 1 \tab outcome 2 \tab total \cr
			group 1 \tab \code{ai} \tab \code{bi} \tab \code{n1i} \cr
			group 2 \tab \code{ci} \tab \code{di} \tab \code{n2i} \cr
			} where \code{ai}, \code{bi}, \code{ci}, and \code{di} denote the cell frequencies and \code{n1i} and \code{n2i} the row totals. For example, in a set of randomized clinical trials, group 1 and group 2 may refer to the treatment and placebo/control group, with outcome 1 denoting some event of interest (e.g., remission) and outcome 2 its complement. In a set of case-control studies, group 1 and group 2 may refer to the group of cases and the group of controls, with outcome 1 denoting, for example, exposure to some risk factor and outcome 2 non-exposure. The 2x2 table may also be the result of cross-sectional (i.e., multinomial) sampling, so that none of the table margins (except the total sample size \code{n1i+n2i}) are fixed through the study design.

			Depending on the type of design (sampling method), a meta-analysis of 2x2 table data can be based on one of several different outcome measures, including the odds ratio, the relative risk (also called risk ratio), the risk difference, and the arcsine transformed risk difference. The phi coefficient, Yule's Q, and Yule's Y are additional measures of association for 2x2 table data (but they may not be the most ideal choices for meta-analyses of such data). For these measures, one needs to supply either \code{ai}, \code{bi}, \code{ci}, and \code{di} or alternatively \code{ai}, \code{ci}, \code{n1i}, and \code{n2i}. The options for the \code{measure} argument are then:
			\itemize{
			\item \code{"RR"}: The \emph{log relative risk} is equal to the log of \code{(ai/n1i)/(ci/n2i)}.
			\item \code{"OR"}: The \emph{log odds ratio} is equal to the log of \code{(ai*di)/(bi*ci)}.
			\item \code{"RD"}: The \emph{risk difference} is equal to \code{(ai/n1i)-(ci/n2i)}.
			\item \code{"AS"}: The \emph{arcsine transformed risk difference} is equal to \code{ asin(sqrt(ai/n1i)) -} \code{asin(sqrt(ci/n2i))}. See Ruecker et al. (2009) for a discussion of this and other outcome measures for 2x2 table data.
			\item \code{"PETO"}: The \emph{log odds ratio estimated with Peto's method} (see Yusuf et al., 1985) is equal to \code{(ai-si*n1i/ni)/((si*ti*n1i*n2i)/(ni^2*(ni-1)))}, where \code{si=ai+ci}, \code{ti=bi+di}, and \code{ni=n1i+n2i}. Note that this measure technically assumes that the true odds ratio is equal to 1 in all tables.
			\item \code{"PHI"}: The \emph{phi coefficient} is equal to \code{(ai*di-bi*ci)/sqrt(n1i*n2i*si*ti)}, where \code{si=ai+ci} and \code{ti=bi+di}.
			\item \code{"YUQ"}: \emph{Yule's Q} is equal to \code{(oi-1)/(oi+1)}, where \code{oi} is the odds ratio.
			\item \code{"YUY"}: \emph{Yule's Y} is equal to \code{(sqrt(oi)-1)/(sqrt(oi)+1)}, where \code{oi} is the odds ratio.
			} Note that the log is taken of the relative risk and the odds ratio, which makes these outcome measures symmetric around 0 and helps to make the distribution of these outcome measure closer to normal.
			
			Cell entries with a zero can be problematic, especially for the relative risk and the odds ratio. Adding a small constant to the cells of the 2x2 tables is a common solution to this problem. When \code{to="all"}, the value of \code{add} is added to each cell of the 2x2 tables in all \eqn{k} tables. When \code{to="only0"}, the value of \code{add} is added to each cell of the 2x2 tables only in those tables with at least one cell equal to 0. When \code{to="if0all"}, the value of \code{add} is added to each cell of the 2x2 tables in all \eqn{k} tables, but only when there is at least one 2x2 table with a zero entry. Setting \code{to="none"} or \code{add=0} has the same effect: No adjustment to the observed table frequencies is made. Depending on the outcome measure and the data, this may lead to division by zero inside of the function (when this occurs, the resulting \code{Inf} value is recoded to \code{NA}).

		}

		\subsection{Raw and Standardized Mean Difference}{

			The raw mean difference and standardized mean difference are useful effect size measures when meta-analyzing a set of studies comparing two experimental groups (e.g., treatment and control groups) or two naturally occurring groups (e.g., men and women) with respect to some quantitative (and ideally normally distributed) dependent variable. For these outcome measures, \code{m1i} and \code{m2i} specify the means of the two groups, \code{sd1i} and \code{sd2i} the standard deviations of the scores in the two groups, and \code{n1i} and \code{n2i} the sample sizes of the two groups.
			\itemize{
			\item \code{"MD"}: The \emph{raw mean difference} is equal to \code{m1i-m2i}.
			\item \code{"SMD"}: The \emph{standardized mean difference} is equal to \code{(m1i-m2i)/spi}, where \code{spi} is the pooled standard deviation of the two groups (which is calculated inside of the function). The standardized mean difference is automatically corrected for its slight positive bias within the function (see Hedges & Olkin, 1985). When \code{vtype="LS"}, the sampling variances are calculated based on the large sample approximation. Alternatively, the unbiased estimates of the sampling variances can be obtained with \code{vtype="UB"}.
			}
		
		}

		\subsection{Raw and Transformed Correlation Coefficient}{

			Another frequently used outcome measure in meta-analyses is the correlation coefficient, which is used to measure the strength of the (linear) relationship between two quantitative variables. Here, one needs to specify \code{ri}, the vector with the raw correlation coefficients, and \code{ni}, the corresponding sample sizes.
			\itemize{
			\item \code{"COR"}: The \emph{raw correlation coefficient} is simply equal to \code{ri} as supplied to the function. When \code{vtype="LS"}, the sampling variances are calculated based on the large sample approximation. Alternatively, an approximation to the unbiased estimates of the sampling variances can be obtained with \code{vtype="UB"} (see Hedges, 1989).	
			\item \code{"UCOR"}: The \emph{unbiased estimate of the correlation coefficient} is obtained by correcting the raw correlation coefficient for its slight negative bias (based on equation 2.7 in Olkin & Pratt, 1958). Again, \code{vtype="LS"} and \code{vtype="UB"} can be used to choose between the large sample approximation or approximately unbiased estimates of the sampling variances.
			\item \code{"ZCOR"}: Fisher's r-to-z transformation is a variance stabilizing transformation for correlation coefficients with the added benefit of also being a rather effective normalizing transformation (Fisher, 1921). The \emph{Fisher's r-to-z transformed correlation coefficient} is equal to \code{1/2*log((1+ri)/(1-ri))}.
			}

		}
		
		\subsection{Proportions and Transformations Thereof}{

			When the studies provide data for single groups with respect to a dichotomous dependent variable, then the raw proportion, the logit transformed proportion, the arcsine transformed proportion, and the Freeman-Tukey double arcsine transformed proportion are useful outcome measures (the log transformed proportion is also a possibility, but not frequently used in meta-analyses). Here, one needs to specify \code{xi} and \code{ni}, denoting the number of individuals experiencing the event of interest and the total number of individuals, respectively. Instead of specifying \code{ni}, one can use \code{mi} to specify the number of individuals that do not experience the event of interest.
			\itemize{
			\item \code{"PR"}: The \emph{raw proportion} is equal to \code{xi/ni}.
			\item \code{"PLN"}: The \emph{log transformed proportion} is equal to the log of \code{xi/ni}.
			\item \code{"PLO"}: The \emph{logit transformed proportion} is equal to the log of \code{xi/(ni-xi)}.
			\item \code{"PAS"}: The arcsine transformation is a variance stabilizing transformation for proportions and is equal to \code{asin(sqrt(xi/ni))}.
			\item \code{"PFT"}: Yet another variance stabilizing transformation for proportions was suggested by Freeman & Tukey (1950). The \emph{Freeman-Tukey double arcsine transformed proportion} is equal to \code{1/2*(asin(sqrt(xi/(ni+1))) + asin(sqrt((xi+1)/(ni+1))))}.
			} Again, zero cell entries can be problematic. When \code{to="all"}, the value of \code{add} is added to \code{xi} and \code{mi} in all \eqn{k} studies. When \code{to="only0"}, the value of \code{add} is added only for studies where the \code{xi} or \code{mi} is equal to 0. When \code{to="if0all"}, the value of \code{add} is added in all \eqn{k} studies, but only when there is at least one study with a zero value for \code{xi} or \code{mi}. Setting \code{to="none"} or \code{add=0} again means that no adjustment to the observed values is made.
		}

	}

	\subsection{Formula Interface}{

		The formula interface works as follows. As above, the argument \code{measure} is a character string specifying which outcome measure should be calculated. The \code{formula} argument is then used to specify the data structure as a multipart formula together with the \code{weights} argument for the group sizes or cell frequencies. The \code{data} argument can be used to specify a data frame containing the variables in the formula and the \code{weights} variable. The \code{add}, \code{to}, and \code{vtype} arguments work as described above.

		The \code{formula} argument takes the form \code{outcome ~ group | study}, where \code{group} is a factor specifying the group variable and \code{study} is a factor specifying the study factor.

		\subsection{Effect Size and Outcome Measures for 2x2 Table Data}{

			For 2x2 table data, \code{group} is a two-level factor specifying the rows of the tables and \code{outcome} is a two-level factor specifying the columns of the tables (the two possible outcomes). The \code{weights} argument is then used to specify the frequencies in the various cells.

		}

		\subsection{Raw and Standardized Mean Difference}{

			For these outcome measures, \code{group} is a two-level factor specifying the group factor and the left-hand side of the formula is composed of two parts, with the first variable for the means and the second variable for the standard deviations (i.e., \code{means + sds ~ group | study}). The \code{weights} argument is then used to specify the sample sizes in the groups.

		}

		\subsection{Raw and Transformed Correlation Coefficient}{

			For these outcome measures, \code{group} is a one-level factor and \code{outcome} is used to specify the observed correlations. The \code{weights} argument is again used to specify the sample sizes.

		}

		\subsection{Proportions and Transformations Thereof}{

			For these outcome measures, \code{group} is a one-level factor and \code{outcome} is a two-level factor specifying the columns of the tables (the two possible outcomes). The \code{weights} argument is then used to specify the frequencies in the various cells.

		}
	
	}
}
\value{
	A data frame with the following elements:
	\item{yi}{value of the effect size or outcome measure.}
	\item{vi}{corresponding (estimated) sampling variance.}
	If \code{append=TRUE} and a data frame was specified via the \code{data} argument, then \code{yi} and \code{vi} are append to this data frame.
}
\note{
	For standard meta-analyses using the typical (wide-format) data layout (i.e., one row in the dataset per study), the default interface is typically easier to use. The advantage of the formula interface is that it can, in principle, handle more complicated data structures (e.g., studies with more than two treatment groups or more than two outcomes). While such functionality is currently not implemented, this may be the case in the future.
}	
\author{Wolfgang Viechtbauer; \email{wvb@www.wvbauer.com}; \url{http://www.wvbauer.com/}}
\references{
	Cooper, H. C., Hedges, L. V. & Valentine, J. C. (Eds.) (2009). \emph{The handbook of research synthesis and meta-analysis} (2nd ed.). New York: Russell Sage Foundation.

	Fisher, R. A. (1921). On the \dQuote{probable error} of a coefficient of correlation deduced from a small sample. \emph{Metron}, \bold{1}, 1--32.

	Freeman, M. F. & Tukey, J. W. (1950). Transformations related to the angular and the square root. \emph{Annals of Mathematical Statistics}, \bold{21}, 607--611.

	Hedges, L. V. (1989). An unbiased correction for sampling error in validity generalization studies. \emph{Journal of Applied Psychology}, \bold{74}, 469--477.

	Hedges, L. V. & Olkin, I. (1985). \emph{Statistical methods for meta-analysis}. San Diego, CA: Academic Press.

	Ruecker, G., Schwarzer, G., Carpenter, J. & Olkin, I. (2009). Why add anything to nothing? The arcsine difference as a measure of treatment effect in meta-analysis with zero cells. \emph{Statistics in Medicine}, \bold{28}, 721--738.

	Olkin, I. & Pratt, J. W. (1958). Unbiased estimation of certain correlation coefficients. \emph{Annals of Mathematical Statistics}, \bold{29}, 201--211.

	Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \url{http://www.jstatsoft.org/v36/i03/}.

	Yusuf, S., Peto, R., Lewis, J., Collins, R. & Sleight, P. (1985). Beta blockade during and after myocardial infarction: An overview of the randomized trials. \emph{Progress in Cardiovascular Disease}, \bold{27}, 335--371.
}
\seealso{
	\code{\link{rma.uni}}, \code{\link{rma.mh}}, \code{\link{rma.peto}}
}
\examples{
### load BCG vaccine data
data(dat.bcg)

### calculate log relative risks and corresponding sampling variances
dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, 
              data=dat.bcg, append=TRUE)
dat

### using formula interface (first rearrange data into required format)
k <- length(dat.bcg$trial)
dat.fm <- data.frame(study=factor(rep(1:k, each=4)))
dat.fm$grp  <- factor(rep(c("T","T","C","C"), k), levels=c("T", "C"))
dat.fm$out  <- factor(rep(c("+","-","+","-"), k), levels=c("+", "-"))
dat.fm$freq <- with(dat.bcg, c(rbind(tpos, tneg, cpos, cneg)))
dat.fm

escalc(out ~ grp | study, weights=freq, data=dat.fm, measure="RR")
}
\keyword{datagen}
