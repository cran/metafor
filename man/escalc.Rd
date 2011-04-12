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
       x1i, x2i, t1i, t2i, m1i, m2i, sd1i, sd2i, 
       xi, mi, ri, ni, ti, data, 
       add=1/2, to="only0", vtype="LS", append=FALSE, ...)

\method{escalc}{formula}(measure, formula, weights, data, 
       add=1/2, to="only0", vtype="LS", ...)
}
\arguments{
   \item{measure}{a character string indicating which effect size or outcome measure should be calculated. See \sQuote{Details} for possible options and how the data should then be specified.}
   \item{formula}{when using the formula interface of the function (see \sQuote{Details} below), a model formula specifying the data structure should be specified via this argument. When not using the formula interface, this argument can be ignored and the data required to calculate the effect sizes or outcomes are then passed to the function via the following set of arguments. See \sQuote{Details}.}
   \item{weights}{vector of weights to specify the group sizes or cell frequencies (only needed when using the formula interface). See \sQuote{Details}.}
   \item{ai}{vector to specify the 2x2 table frequencies (upper left cell). See \sQuote{Details}.}
   \item{bi}{vector to specify the 2x2 table frequencies (upper right cell). See \sQuote{Details}.}
   \item{ci}{vector to specify the 2x2 table frequencies (lower left cell). See \sQuote{Details}.}
   \item{di}{vector to specify the 2x2 table frequencies (lower right cell). See \sQuote{Details}.}
   \item{n1i}{vector to specify the group sizes or row totals (first group/row). See \sQuote{Details}.}
   \item{n2i}{vector to specify the group sizes or row totals (second group/row). See \sQuote{Details}.}
   \item{x1i}{vector to specify the number of cases (first group). See \sQuote{Details}.}
   \item{x2i}{vector to specify the number of cases (second group). See \sQuote{Details}.}
   \item{t1i}{vector to specify the total person-times (first group). See \sQuote{Details}.}
   \item{t2i}{vector to specify the total person-times (second group). See \sQuote{Details}.}
   \item{m1i}{vector to specify the means (first group). See \sQuote{Details}.}
   \item{m2i}{vector to specify the means (second group). See \sQuote{Details}.}
   \item{sd1i}{vector to specify the standard deviations (first group). See \sQuote{Details}.}
   \item{sd2i}{vector to specify the standard deviations (second group). See \sQuote{Details}.}
   \item{xi}{vector to specify the frequencies of the event of interest. See \sQuote{Details}.}
   \item{mi}{vector to specify the frequencies of the complement of the event. See \sQuote{Details}.}
   \item{ri}{vector to specify the raw correlation coefficients. See \sQuote{Details}.}
   \item{ni}{vector to specify the sample sizes. See \sQuote{Details}.}
   \item{ti}{vector to specify the total person-times. See \sQuote{Details}.}
   \item{data}{an optional data frame containing the variables given to the arguments above.}
   \item{add}{a non-negative number indicating the amount to add to zero cells, counts, or frequencies. See \sQuote{Details}.}
   \item{to}{a string indicating when the values under \code{add} should be added (either \code{"all"}, \code{"only0"}, \code{"if0all"}, or \code{"none"}). See \sQuote{Details}.}
   \item{vtype}{a string indicating the type of sampling variances to calculate (either \code{"LS"} or \code{"UB"}). See \sQuote{Details}.}
   \item{append}{logical indicating whether the data frame specified via the \code{data} argument (if one has been specified) should be returned together with the effect sizes and sampling variances (default is \code{FALSE}).}
   \item{\dots}{other arguments.}
}
\details{
   There are two interfaces to using the \code{escalc} function, the default and a formula interface. The two interfaces are described below.

   \subsection{Default Interface}{
   
      The default interface works as follows. The argument \code{measure} is a character string specifying which outcome measure should be calculated (see below for the various options), arguments \code{ai} through \code{ni} are then used to specify the needed information to calculate the various measures (depending on the outcome measure, different arguments need to be supplied), and \code{data} can be used to specify a data frame containing the variables given to the previous arguments. The \code{add} and \code{to} arguments may be needed when dealing with 2x2 table data that contain cells with zeros. Finally, the \code{vtype} argument is used to specify how to calculate the sampling variance estimate (see below).

      \subsection{Effect Size and Outcome Measures for 2x2 Table Data}{

         Meta-analyses in the health/medical sciences are often based on studies providing data in terms of 2x2 tables. In particular, assume that we have \eqn{k} tables of the form:
         \tabular{lccc}{
                 \tab outcome 1 \tab outcome 2 \tab total      \cr
         group 1 \tab \code{ai} \tab \code{bi} \tab \code{n1i} \cr
         group 2 \tab \code{ci} \tab \code{di} \tab \code{n2i}
         } where \code{ai}, \code{bi}, \code{ci}, and \code{di} denote the cell frequencies and \code{n1i} and \code{n2i} the row totals. For example, in a set of randomized clinical trials (RCTs) or cohort studies, group 1 and group 2 may refer to the treatment (exposed) and placebo/control (not exposed) group, with outcome 1 denoting some event of interest (e.g., death) and outcome 2 its complement. In a set of case-control studies, group 1 and group 2 may refer to the group of cases and the group of controls, with outcome 1 denoting, for example, exposure to some risk factor and outcome 2 non-exposure. The 2x2 table may also be the result of cross-sectional (i.e., multinomial) sampling, so that none of the table margins (except the total sample size \code{n1i+n2i}) are fixed through the study design.

         Depending on the type of design (sampling method), a meta-analysis of 2x2 table data can be based on one of several different outcome measures, including the odds ratio, the relative risk (also called risk ratio), the risk difference, and the arcsine transformed risk difference (for example, for case-control, the odds ratio is the measure of choice, while for RCTs and cohort studies, all of these measures may be applicable). The phi coefficient, Yule's Q, and Yule's Y are additional measures of association for 2x2 table data (although they are not frequently used in meta-analyses). For these outcome measures, one needs to specify either \code{ai}, \code{bi}, \code{ci}, and \code{di} or alternatively \code{ai}, \code{ci}, \code{n1i}, and \code{n2i}. The options for the \code{measure} argument are then:
         \itemize{
         \item \code{"RR"}: The \emph{log relative risk} is equal to the log of \code{(ai/n1i)/(ci/n2i)}.
         \item \code{"OR"}: The \emph{log odds ratio} is equal to the log of \code{(ai*di)/(bi*ci)}.
         \item \code{"RD"}: The \emph{risk difference} is equal to \code{(ai/n1i)-(ci/n2i)}.
         \item \code{"AS"}: The arcsine transformation is a variance stabilizing transformation for proportions. The \emph{arcsine transformed risk difference} is equal to \code{asin(sqrt(ai/n1i)) - asin(sqrt(ci/n2i))}. See Ruecker et al. (2009) for a discussion of this and other outcome measures for 2x2 table data.
         \item \code{"PETO"}: The \emph{log odds ratio estimated with Peto's method} (see Yusuf et al., 1985) is equal to \code{(ai-si*n1i/ni)/((si*ti*n1i*n2i)/(ni^2*(ni-1)))}, where \code{si=ai+ci}, \code{ti=bi+di}, and \code{ni=n1i+n2i}.
         \item \code{"PHI"}: The \emph{phi coefficient} is equal to \code{(ai*di-bi*ci)/sqrt(n1i*n2i*si*ti)}, where \code{si=ai+ci} and \code{ti=bi+di}.
         \item \code{"YUQ"}: \emph{Yule's Q} is equal to \code{(oi-1)/(oi+1)}, where \code{oi} is the odds ratio.
         \item \code{"YUY"}: \emph{Yule's Y} is equal to \code{(sqrt(oi)-1)/(sqrt(oi)+1)}, where \code{oi} is the odds ratio.
         } Note that the log is taken of the relative risk and the odds ratio, which makes these outcome measures symmetric around 0 and helps to make the distribution of these outcome measure closer to normal.
         
         Cell entries with a zero can be problematic, especially for the relative risk and the odds ratio. Adding a small constant to the cells of the 2x2 tables is a common solution to this problem. When \code{to="all"}, the value of \code{add} is added to each cell of all 2x2 tables. When \code{to="only0"} (the default), the value of \code{add} (the default is 1/2) is added to each cell of the 2x2 tables with at least one cell equal to 0. When \code{to="if0all"}, the value of \code{add} is added to each cell of all 2x2 tables, but only when there is at least one 2x2 table with a zero cell. Setting \code{to="none"} or \code{add=0} has the same effect: No adjustment to the observed table frequencies is made. Depending on the outcome measure and the data, this may lead to division by zero inside of the function (when this occurs, the resulting \code{Inf} value is recoded to \code{NA}).

         An example dataset corresponding to data of this type is provided in \code{\link{dat.bcg}}.

      }

      \subsection{Incidence Rate Ratios and Differences}{

         Epidemiological studies often compare the incidence rates (i.e., the rate of occurrence of a particular outcome, e.g., a certain disease, over a particular period of time) of two different groups (e.g., exposed, not exposed). In particular, assume that we have \eqn{k} tables of the form:
         \tabular{lcc}{
                 \tab cases      \tab person-time \cr
         group 1 \tab \code{x1i} \tab \code{t1i} \cr
         group 2 \tab \code{x2i} \tab \code{t2i}
         } where \code{x1i} and \code{x2i} denote the number of cases in the first and the second group, respectively, and \code{t1i} and \code{t2i} the corresponding total person-times at risk. Commonly used effect size or outcome measures in this context are the ratio or the difference between the two incidence rates. The options for the \code{measure} argument are then:
         \itemize{
         \item \code{"IRR"}: The \emph{log incidence rate ratio} is equal to the log of \code{(x1i/t1i)/(x2i/t2i)}.
         \item \code{"IRD"}: The \emph{incidence rate difference} is equal to \code{(x1i/t1i)-(x2i/t2i)}.
         \item \code{"IRSD"}: The square-root transformation is a variance stabilizing transformation for incidence rates. The \emph{square-root transformed incidence rate difference} is equal to \code{sqrt(x1i/t1i)-sqrt(x2i/t2i)}.
         } Note that the log is taken of the incidence rate ratio, which makes this outcome measure symmetric around 0 and helps to make its distribution closer to normal.

         Studies with zero cases in one or both groups can be problematic, especially for the incidence rate ratio. Adding a small constant to the number of cases is a common solution to this problem. When \code{to="all"}, the value of \code{add} is added to \code{x1i} and \code{x2i} in all \code{k} studies. When \code{to="only0"} (the default), the value of \code{add} (the default is 1/2) is added to \code{x1i} and \code{x2i} only in the studies that have zero cases in one or both groups. When \code{to="if0all"}, the value of \code{add} is added to \code{x1i} and \code{x2i} in all \code{k} studies, but only when there is at least one study with zero cases in one or both groups. Setting \code{to="none"} or \code{add=0} has the same effect: No adjustment to the observed number of cases is made. Depending on the outcome measure and the data, this may lead to division by zero inside of the function (when this occurs, the resulting \code{Inf} value is recoded to \code{NA}).

         An example dataset corresponding to data of this type is provided in \code{\link{dat.warfarin}}.
      }

      \subsection{Mean Differences, Standardized Mean Differences, and Ratio of Means}{

         The raw mean difference, the standardized mean difference, and the ratio of means (also called response ratio) are useful effect size measures when meta-analyzing a set of studies comparing two experimental groups (e.g., treatment and control groups) or two naturally occurring groups (e.g., men and women) with respect to some quantitative (and ideally normally distributed) dependent variable. For these outcome measures, \code{m1i} and \code{m2i} are used to specify the means of the two groups, \code{sd1i} and \code{sd2i} the standard deviations of the scores in the two groups, and \code{n1i} and \code{n2i} the sample sizes of the two groups. The options for the \code{measure} argument are then:
         \itemize{
         \item \code{"MD"}: The \emph{raw mean difference} is equal to \code{m1i-m2i}.
         \item \code{"SMD"}: The \emph{standardized mean difference} is equal to \code{(m1i-m2i)/spi}, where \code{spi} is the pooled standard deviation of the two groups (which is calculated inside of the function based on \code{sd1i} and \code{sd2i}). The standardized mean difference is automatically corrected for its slight positive bias within the function (see Hedges & Olkin, 1985). When \code{vtype="LS"}, the sampling variances are calculated based on the large sample approximation. Alternatively, the unbiased estimates of the sampling variances can be obtained with \code{vtype="UB"}.
         \item \code{"ROM"}: The \emph{log transformed ratio of means} is equal to \code{log(m1i/m2i)}. See Hedges et al. (1999) for more details on this outcome measure.
         } Note that the log is taken of the ratio of means, which makes this outcome measures symmetric around 0 and helps to make the distribution of this measure closer to normal (however, if \code{m1i} and \code{m2i} have opposite signs, this outcome measure cannot be computed).

         An example dataset corresponding to data of this type is provided in \code{\link{dat.los}} (for mean differences and standardized mean differences). An example dataset showing the use of the ratio of means measure is provided in \code{\link{dat.co2}}.
      
      }

      \subsection{Raw and Transformed Correlation Coefficients}{

         Another frequently used outcome measure in meta-analyses is the correlation coefficient, which is used to measure the strength of the (linear) relationship between two quantitative variables. Here, one needs to specify \code{ri}, the vector with the raw correlation coefficients, and \code{ni}, the corresponding sample sizes. The options for the \code{measure} argument are then:
         \itemize{
         \item \code{"COR"}: The \emph{raw correlation coefficient} is simply equal to \code{ri} as supplied to the function. When \code{vtype="LS"}, the sampling variances are calculated based on the large sample approximation. Alternatively, an approximation to the unbiased estimates of the sampling variances can be obtained with \code{vtype="UB"} (see Hedges, 1989).   
         \item \code{"UCOR"}: The \emph{unbiased estimate of the correlation coefficient} is obtained by correcting the raw correlation coefficient for its slight negative bias (based on equation 2.7 in Olkin & Pratt, 1958). Again, \code{vtype="LS"} and \code{vtype="UB"} can be used to choose between the large sample approximation or approximately unbiased estimates of the sampling variances.
         \item \code{"ZCOR"}: Fisher's r-to-z transformation is a variance stabilizing transformation for correlation coefficients with the added benefit of also being a rather effective normalizing transformation (Fisher, 1921). The \emph{Fisher's r-to-z transformed correlation coefficient} is equal to \code{1/2*log((1+ri)/(1-ri))}.
         }

         An example dataset corresponding to data of this type is provided in \code{\link{dat.empint}}.

      }
      
      \subsection{Proportions and Transformations Thereof}{

         When the studies provide data for single groups with respect to a dichotomous dependent variable, then the raw proportion, the logit transformed proportion, the arcsine transformed proportion, and the Freeman-Tukey (double arcsine) transformed proportion are useful outcome measures (the log transformed proportion is also a possibility, but not frequently used in meta-analyses). Here, one needs to specify \code{xi} and \code{ni}, denoting the number of individuals experiencing the event of interest and the total number of individuals, respectively. Instead of specifying \code{ni}, one can use \code{mi} to specify the number of individuals that do not experience the event of interest. The options for the \code{measure} argument are then:
         \itemize{
         \item \code{"PR"}: The \emph{raw proportion} is equal to \code{xi/ni}.
         \item \code{"PLN"}: The \emph{log transformed proportion} is equal to the log of \code{xi/ni}.
         \item \code{"PLO"}: The \emph{logit transformed proportion} is equal to the log of \code{xi/(ni-xi)} (i.e., the log of the odds).
         \item \code{"PAS"}: The arcsine transformation is a variance stabilizing transformation for proportions. The \emph{arcsine transformed proportion} is equal to \code{asin(sqrt(xi/ni))}.
         \item \code{"PFT"}: Another variance stabilizing transformation for proportions was suggested by Freeman & Tukey (1950). The \emph{Freeman-Tukey double arcsine transformed proportion} is equal to \code{1/2*(asin(sqrt(xi/(ni+1))) + asin(sqrt((xi+1)/(ni+1))))}.
         } Zero cell entries can be problematic for certain outcome measures. When \code{to="all"}, the value of \code{add} is added to \code{xi} and \code{mi} in all \eqn{k} studies. When \code{to="only0"} (the default), the value of \code{add} (the default is 1/2) is added to \code{xi} and \code{mi} only for studies where \code{xi} or \code{mi} is equal to 0. When \code{to="if0all"}, the value of \code{add} is added in all \eqn{k} studies, but only when there is at least one study with a zero value for \code{xi} or \code{mi}. Setting \code{to="none"} or \code{add=0} again means that no adjustment to the observed values is made.
      }

      \subsection{Incidence Rates and Transformations Thereof}{

         Instead of proportions, we may also be interested in aggregating individual incidence rates. Here, one needs to specify \code{xi} and \code{ti}, denoting the number of individuals experiencing the event of interest and the total person-time at risk, respectively. The options for the \code{measure} argument are then:
         \itemize{
         \item \code{"IR"}: The \emph{raw incidence rate} is equal to \code{xi/ti}.
         \item \code{"IRLN"}: The \emph{log transformed incidence rate} is equal to the log of \code{xi/ti}.
         \item \code{"IRS"}: The square-root transformation is a variance stabilizing transformation for incidence rates. The \emph{square-root transformed incidence rate} is equal to \code{sqrt(xi/ti)}.
         \item \code{"IRFT"}: Another variance stabilizing transformation for incidence rates can be based on Freeman & Tukey (1950). The \emph{Freeman-Tukey transformed incidence rate} is equal to \code{sqrt(xi/ti) + sqrt(xi/ti+1/ti)}.
         } Studies with zero cases can be problematic, especially for the log transformed incidence rate. Adding a small constant to the number of cases is a common solution to this problem. When \code{to="all"}, the value of \code{add} is added to \code{xi} in all \code{k} studies. When \code{to="only0"} (the default), the value of \code{add} (the default is 1/2) is added to \code{xi} only in the studies that have zero cases. When \code{to="if0all"}, the value of \code{add} is added to \code{xi} in all \code{k} studies, but only when there is at least one study with zero cases. Setting \code{to="none"} or \code{add=0} has the same effect: No adjustment to the observed number of cases is made. Depending on the outcome measure and the data, this may lead to division by zero inside of the function (when this occurs, the resulting \code{Inf} value is recoded to \code{NA}).
      }

   }
   
   \subsection{Formula Interface}{

      The formula interface works as follows. As above, the argument \code{measure} is a character string specifying which outcome measure should be calculated. The \code{formula} argument is then used to specify the data structure as a multipart formula. The \code{data} argument can be used to specify a data frame containing the variables in the formula. The \code{add}, \code{to}, and \code{vtype} arguments work as described above.

      \subsection{Effect Size and Outcome Measures for 2x2 Table Data}{

         For 2x2 table data, the \code{formula} argument takes the form \code{outcome ~ group | study}, where \code{group} is a two-level factor specifying the rows of the tables, \code{outcome} is a two-level factor specifying the columns of the tables (the two possible outcomes), and \code{study} is a factor specifying the study factor. The \code{weights} argument is used to specify the frequencies in the various cells.

      }

      \subsection{Incidence Rate Ratios and Differences}{

         For these outcome measures, the \code{formula} argument takes the form \code{cases/times ~ group | study}, where \code{group} is a two-level factor specifying the group factor and \code{study} is a factor specifying the study factor. The left-hand side of the formula is composed of two parts, with the first variable for the number of cases and the second variable for the person-time at risk.

      }

      \subsection{Mean Differences, Standardized Mean Differences, and Ratio of Means}{

         For these outcome measures, the \code{formula} argument takes the form \code{means/sds ~ group | study}, where \code{group} is a two-level factor specifying the group factor and \code{study} is a factor specifying the study factor. The left-hand side of the formula is composed of two parts, with the first variable for the means and the second variable for the standard deviations. The \code{weights} argument is used to specify the sample sizes in the groups.

      }

      \subsection{Raw and Transformed Correlation Coefficients}{

         For these outcome measures, the \code{formula} argument takes the form \code{outcome ~ 1 | study}, where \code{outcome} is used to specify the observed correlations and \code{study} is a factor specifying the study factor. The \code{weights} argument is used to specify the sample sizes.

      }

      \subsection{Proportions and Transformations Thereof}{

         For these outcome measures, the \code{formula} argument takes the form \code{outcome ~ 1 | study}, where \code{outcome} is a two-level factor specifying the columns of the tables (the two possible outcomes) and \code{study} is a factor specifying the study factor. The \code{weights} argument is used to specify the frequencies in the various cells.

      }

         \subsection{Incidence Rates and Transformations Thereof}{

         For these outcome measures, the \code{formula} argument takes the form \code{cases/times ~ 1 | study}, where \code{study} is a factor specifying the study factor. The left-hand side of the formula is composed of two parts, with the first variable for the number of cases and the second variable for the person-time at risk.

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
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \cr
   project homepage: \url{http://www.metafor-project.org/} \cr
   author homepage: \url{http://www.wvbauer.com/}
}
\references{
   Cooper, H. C., Hedges, L. V., & Valentine, J. C. (Eds.) (2009). \emph{The handbook of research synthesis and meta-analysis} (2nd ed.). New York: Russell Sage Foundation.

   Fisher, R. A. (1921). On the \dQuote{probable error} of a coefficient of correlation deduced from a small sample. \emph{Metron}, \bold{1}, 1--32.

   Freeman, M. F. & Tukey, J. W. (1950). Transformations related to the angular and the square root. \emph{Annals of Mathematical Statistics}, \bold{21}, 607--611.

   Hedges, L. V. (1989). An unbiased correction for sampling error in validity generalization studies. \emph{Journal of Applied Psychology}, \bold{74}, 469--477.

   Hedges, L. V. & Olkin, I. (1985). \emph{Statistical methods for meta-analysis}. San Diego, CA: Academic Press.

   Hedges, L. V., Gurevitch, J., & Curtis, P. S. (1999). The meta-analysis of response ratios in experimental ecology. \emph{Ecology}, \bold{80}, 1150--1156.
   
   Ruecker, G., Schwarzer, G., Carpenter, J., & Olkin, I. (2009). Why add anything to nothing? The arcsine difference as a measure of treatment effect in meta-analysis with zero cells. \emph{Statistics in Medicine}, \bold{28}, 721--738.

   Olkin, I. & Pratt, J. W. (1958). Unbiased estimation of certain correlation coefficients. \emph{Annals of Mathematical Statistics}, \bold{29}, 201--211.

   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \url{http://www.jstatsoft.org/v36/i03/}.

   Yusuf, S., Peto, R., Lewis, J., Collins, R., & Sleight, P. (1985). Beta blockade during and after myocardial infarction: An overview of the randomized trials. \emph{Progress in Cardiovascular Disease}, \bold{27}, 335--371.
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
dat.fm      <- data.frame(study=factor(rep(1:k, each=4)))
dat.fm$grp  <- factor(rep(c("T","T","C","C"), k), levels=c("T","C"))
dat.fm$out  <- factor(rep(c("+","-","+","-"), k), levels=c("+","-"))
dat.fm$freq <- with(dat.bcg, c(rbind(tpos, tneg, cpos, cneg)))
dat.fm
escalc(out ~ grp | study, weights=freq, data=dat.fm, measure="RR")
}
\keyword{datagen}
