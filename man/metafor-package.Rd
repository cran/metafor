\name{metafor-package}
\alias{metafor-package}
\alias{metafor}
\docType{package}
\title{Metafor: A Meta-Analysis Package for R}
\description{
	The \pkg{metafor} package consists of a collection of functions for conducting meta-analyses in \R. The package includes functions for calculating various effect size or outcome measures frequently used in meta-analyses (e.g., risk differences, risk ratios, odds ratios, standardized mean differences, Fisher's z-transformed correlation coefficients) and then allows the user to fit fixed-, random-, and mixed-effects models to these data. By including study-level covariates (\sQuote{moderators}) in these models, so-called \sQuote{meta-regression} analyses can be carried out. For meta-analyses of 2x2 table data, the Mantel-Haenszel and Peto's method are also implemented.

	Various methods are available to assess model fit, to identify outliers and influential studies, and for conducting sensitivity analyses (e.g., standardized residuals, Cook's distances, leave-one-out analyses). Advanced techniques for conducting hypothesis tests and obtaining confidence intervals (e.g., for the average effect size or for the model coefficients in meta-regression models) have also been implemented (e.g., permutation tests, the Knapp and Hartung method).

	The package also provides functions for creating forest, funnel, radial, normal quantile-quantile, and L'Abbe plots. The presence of publication bias (or more accurately, for funnel plot asymmetry) and its potential impact on the results can be examined via the rank correlation and Egger's regression test and by applying the trim and fill method.
}
\section{The rma.uni Function}{
	The various meta-analytic models that are usually used in practice are special cases of the linear (mixed-effects) model. The \code{\link{rma.uni}} function (with alias \code{\link{rma}}) provides a general framework for fitting such models. The function can be used in conjunction with any of the usual effect size or outcome measures used in meta-analyses (e.g., log odds ratios, log relative risks, risk differences, mean differences, standardized mean differences, raw correlation coefficients, correlation coefficients transformed with Fisher's r-to-z transformation, and so on). For details on these effect size or outcome measures, see the documentation of the \code{\link{escalc}} function. The notation and models underlying the \code{\link{rma.uni}} function are explained below.

	For a set of \eqn{i = 1, \ldots, k}{i = 1, ..., k} independent studies, let \eqn{y_i} denote the observed value of the effect size or outcome measure in the \eqn{i^{th}}{ith} study. Let \eqn{\theta_i}{theta_i} denote the corresponding (unknown) true effect or outcome in the \eqn{i^{th}}{ith} study, such that \eqn{y_i | \theta_i \sim N(\theta_i, v_i)}{y_i | theta_i ~ N(theta_i, v_i)}. In other words, the observed effects or outcomes are assumed to be unbiased and normally distributed estimates of the corresponding true effects or outcomes with sampling variances equal to \eqn{v_i}. The \eqn{v_i} values are assumed to be known. Depending on the outcome measure used, a bias correction, normalizing, and/or variance stabilizing transformation may be necessary to ensure that these assumptions are (approximately) true (e.g., the log transformation for odds ratios, Fisher's r-to-z transformation for correlations, the bias correction for standardized mean differences described by Hedges & Olkin, 1985; see section \code{\link{escalc}} for more details).

	The \bold{fixed-effects model} conditions on the true effects or outcomes and therefore provides a \emph{conditional inference} about the set of \eqn{k} studies included in the meta-analysis. This implies that the fitted model provides an estimate of \eqn{\bar{\theta}_w = \sum_{i=1}^k w_i \theta_i / \sum_{i=1}^k w_i}{theta_w = sum w_i theta_i / sum w_i}, that is, the \emph{weighted average} of the true effects in the set of \eqn{k} studies, with weights equal to \eqn{w_i = 1/v_i} (this is what is often described as the \sQuote{inverse-variance} method in the meta-analytic literature). One can also employ an unweighted estimation method, which provides an estimate of the \emph{unweighted average} of the true effects in the set of \eqn{k} studies (i.e., an estimate of \eqn{\bar{\theta}_u = \sum_{i=1}^k \theta_i / k}{theta_u = sum theta_i / k}).

	Moderators can be included in the fixed-effects model, yielding a \bold{fixed-effects with moderators model}. Again, since the model conditions on the set of \eqn{k} studies included in the meta-analysis, the regression coefficients from the fitted model estimate the weighted least squares relationship between the true effects and the moderator variables within the set of \eqn{k} studies included in the meta-analysis (again using weights equal to \eqn{w_i = 1/v_i}). The (unweighted) least squares relationship between the true effects and the moderator variables can be obtained when using the unweighted estimation method.

	The \bold{random-effects model} does not condition on the true effects. Instead, the \eqn{k} studies included in the meta-analysis are assumed to be a random selection from a hypothetical population of studies. One can envision this hypothetical population as an essentially infinite set of studies comprising all of the studies that have been conducted, that could have been conducted, or that may be conducted in the future. We assume \eqn{theta_i \sim N(\mu, \tau^2)}{theta_i ~ N(mu, tau^2)}, that is, the true effects or outcomes in the population of studies are assumed to be normally distributed with \eqn{\mu}{mu} denoting the average effect and \eqn{\tau^2}{tau^2} denoting the variance of the true effects in the population (\eqn{\tau^2}{tau^2} is therefore often referred to as the \sQuote{amount of heterogeneity} in the true effects). The fitted model provides an estimate of \eqn{\mu}{mu} and \eqn{\tau^2}{tau^2}. Consequently, the random-effects model provides an \emph{unconditional inference} about the average effect in the population of studies from which the \eqn{k} studies included in the meta-analysis are assumed to be a random selection.

	When including moderator variables in the random-effects model, we obtain what is typically called a \bold{mixed-effects model} in the meta-analytic literature. The coefficients from the fitted model then estimate the relationship between the average true effect or outcome in the population of studies and the moderator variables included in the model. The value of \eqn{\tau^2}{tau^2} in the mixed-effects model denotes the \sQuote{amount of residual heterogeneity} in the true effects or outcomes (i.e., the amount of variability among the true effects or outcomes that is not accounted for by the moderators included in the model).

	One can also choose between weighted and unweighted estimation in the context of the random- and mixed-effects model, although the parameters that are estimated remain the same regardless of the estimation method used (as opposed to the fixed-effects model case, where the parameter estimated is different for weighted and unweighted estimation).

	Contrary to what is often stated in the literature, it is important to realize that the fixed-effects model does \emph{not} assume that the true effects or outcomes are homogeneous (i.e., that \eqn{\theta_i}{theta_i} is equal to some common value \eqn{\theta}{theta} for all \eqn{k} studies). In other words, fixed-effects models provide perfectly valid inferences under heterogeneity, as long as one is restricting these inferences to the set of studies included in the meta-analysis (more specifically, to sets of \eqn{k} studies with true effects equal to the true effects of the \eqn{k} studies included in the meta-analysis). On the other hand, the random-effects model provides an inference about the average effect in the entire population of studies from which the included studies are assumed to be a random selection.

	In the special case that the true effects are actually homogeneous, the distinction between fixed- and random-effects models disappears, since homogeneity implies that \eqn{\mu = \bar{\theta}_w = \bar{\theta}_u \equiv \theta}{mu = theta_w = theta_u = theta}. However, since there is no infallible method to test whether the true effects are really homogeneous or not, a researcher should decide on the type of inference desired before examining the data and choose the model accordingly. For more details on the distinction between fixed- and random-effects models, see Hedges and Vevea (1998) and Laird and Mosteller (1990).
}
\section{The rma.mh Function}{
	The Mantel-Haenszel method provides an alternative approach to fitting the fixed-effects model when dealing with studies providing data in the form of 2x2 tables or in the form of incidence rates (Mantel & Haenszel, 1959). The method is particularly advantageous when aggregating a large number of studies/tables with small sample sizes (the so-called sparse data or increasing strata case). The Mantel-Haenszel method is implemented in the \code{\link{rma.mh}} function. It can be used in combination with odds ratios, relative risks, risk differences, and incidence rate ratios. The Mantel-Haenszel method is always based on a weighted estimation approach.
}
\section{The rma.peto Function}{
	Yet another method that can be used in the context of a meta-analysis of 2x2 tables is Peto's method (see Yusuf et al., 1985), implemented in the \code{\link{rma.peto}} function. It is a weighted estimation approach for the combination of odds ratios.
}
\section{Future Plans and Updates}{
	The \pkg{metafor} package is a work in progress and is updated on a regular basis with new functions and options. With \code{\link{metafor.news}}, you can read the \file{NEWS} file of the package after installation. Comments, feedback, and suggestions for improvements are very welcome.

	And since this is a frequently-asked-question: Functions for more complex syntheses (e.g., network meta-analyses, multivariate meta-analyses, meta-analyses with correlated outcomes) are currently under development and will be incorporated into the package in the future.
}
\section{Citing the Package}{
	To cite the package, please use the following reference:

	Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \url{http://www.jstatsoft.org/v36/i03/}.
}
\section{Getting Started with the Package}{
	The paper mentioned above is a good starting place for those interesting in using the \pkg{metafor} package. The purpose of the article is to provide a general overview of the package and its capabilities (as of version 1.4-0). Not all of the functions and options are described in the paper, but it should provide a pretty thorough introduction. The paper can be freely downloaded from the URL given above.

	Alternatively (or in addition to reading the paper), carefully read this page and then the help pages for the \code{\link{escalc}} and the \code{\link{rma.uni}} functions (or the \code{\link{rma.mh}} and \code{\link{rma.peto}} functions if you intend to use the Mantel-Haenszel or Peto's method). The model fitting functions (i.e., \code{\link{rma.uni}}, \code{\link{rma.mh}}, and \code{\link{rma.peto}}) provide links to the help pages of many additional functions, which can be used after fitting a model.
}
\author{
	Wolfgang Viechtbauer \email{wvb@metafor-project.org} \cr
	project homepage: \url{http://www.metafor-project.org/} \cr
	author homepage: \url{http://www.wvbauer.com/} \cr

	Please post any questions about the package to the R-help mailing list (\url{https://stat.ethz.ch/mailman/listinfo/r-help}). Make sure you first read and follow the posting guide (\url{http://www.r-project.org/posting-guide.html}). To make sure that my mail filter catches your post, please include the package name (i.e., \sQuote{metafor}) in your post.
}
\references{
	Cooper, H. C., Hedges, L. V., & Valentine, J. C. (Eds.) (2009). \emph{The handbook of research synthesis and meta-analysis} (2nd ed.). New York: Russell Sage Foundation.

	Hedges, L. V. & Olkin, I. (1985). \emph{Statistical methods for meta-analysis}. San Diego, CA: Academic Press.

	Hedges, L. V. & Vevea, J. L. (1998). Fixed- and random-effects models in meta-analysis. \emph{Psychological Methods}, \bold{3}, 486--504.

	Laird, N. M. & Mosteller, F. (1990). Some statistical methods for combining experimental results. \emph{International Journal of Technology Assessment in Health Care}, \bold{6}, 5--30.

	Mantel, N. & Haenszel, W. (1959). Statistical aspects of the analysis of data from retrospective studies of disease. \emph{Journal of the National Cancer Institute}, \bold{22}, 719--748.

	Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \url{http://www.jstatsoft.org/v36/i03/}.

	Yusuf, S., Peto, R., Lewis, J., Collins, R., & Sleight, P. (1985). Beta blockade during and after myocardial infarction: An overview of the randomized trials. \emph{Progress in Cardiovascular Disease}, \bold{27}, 335--371.
}
\keyword{package}
