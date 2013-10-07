\name{rma.mv}
\alias{rma.mv}
\title{Meta-Analysis via Multivariate/Multilevel Linear (Mixed-Effects) Models}
\description{Function to fit meta-analytic multivariate/multilevel fixed- and random/mixed-effects models with or without moderators via linear (mixed-effects) models. See below and the documentation of the \pkg{\link{metafor-package}} for more details on these models.}
\usage{
rma.mv(yi, V, mods, random, struct="CS",
       intercept=TRUE, data, slab, subset,
       method="REML", tdist=FALSE, level=95, digits=4, btt,
       R, sigma2, tau2, rho, verbose=FALSE, control)
}
\arguments{
   \item{yi}{vector of length \eqn{k} with the observed effect sizes or outcomes. See \sQuote{Details}.}
   \item{V}{vector of length \eqn{k} with the corresponding sampling variances or a \eqn{k \times k}{k x k} variance-covariance matrix of the sampling errors. See \sQuote{Details}.}
   \item{mods}{optional argument to include one or more moderators in the model. A single moderator can be given as a vector of length \eqn{k} specifying the values of the moderator. Multiple moderators are specified by giving a matrix with \eqn{k} rows and \eqn{p'} columns. Alternatively, a model \code{\link{formula}} can be used to specify the model. See \sQuote{Details}.}
   \item{random}{either a single one-sided formula or list of one-sided formulas to specify the random-effects structure of the model. See \sQuote{Details}.}
   \item{struct}{character string to specify the variance structure of an \code{~ inner | outer} formula in the \code{random} argument. Either \code{"CS"} for compound symmetry, \code{"HCS"} for heteroscedastic compound symmetry, or \code{"UN"} for an unstructured variance-covariance matrix. See \sQuote{Details}.}
   \item{intercept}{logical indicating whether an intercept term should be added to the model (default is \code{TRUE}).}
   \item{data}{optional data frame containing the data supplied to the function.}
   \item{slab}{optional vector with unique labels for the \eqn{k} studies.}
   \item{subset}{optional vector indicating the subset of studies that should be used for the analysis. This can be a logical vector of length \eqn{k} or a numeric vector indicating the indices of the observations to include.}
   \item{method}{character string specifying whether a fixed- or a random/mixed-effects model should be fitted. A fixed-effects model (with or without moderators) is fitted when using \code{method="FE"}. Random/mixed-effects models are fitted by setting \code{method} equal to \code{"ML"} (for maximum-likelihood estimation) or \code{"REML"} (for restricted maximum-likelihood estimation). Default is \code{"REML"}. See \sQuote{Details}.}
   \item{tdist}{logical specifying whether test statistics and confidence intervals should be based on the normal (when \code{FALSE}, the default) or the t-distribution (when \code{TRUE}). See \sQuote{Details}.}
   \item{level}{numerical value between 0 and 100 specifying the confidence interval level (default is 95).}
   \item{digits}{integer specifying the number of decimal places to which the printed results should be rounded (default is 4).}
   \item{btt}{optional vector of indices specifying which coefficients to include in the omnibus test of moderators. See \sQuote{Details}.}
   \item{R}{an optional named list of known correlation matrices corresponding to (some of) the components specified via the \code{random} argument. See \sQuote{Details}.}
   \item{sigma2}{optional numerical vector (of the same length as the number of random intercept components specified via the \code{random} argument) to fix the corresponding \eqn{\sigma^2} value(s). A specific \eqn{\sigma^2} value can be fixed by setting the corresponding element of this argument to the desired value. A specific \eqn{\sigma^2} value will be estimated if the corresponding element is set equal to \code{NA}. See \sQuote{Details}.}
   \item{tau2}{optional numerical value (for \code{struct="CS"}) or vector (for \code{struct="HCS"} or \code{struct="UN"}) to fix the amount of (residual) heterogeneity in the levels of the inner factor corresponding to an \code{~ inner | outer} formula specified in the \code{random} argument. A numeric value fixes a particular \eqn{\tau^2} value, while \code{NA} means that the value should be estimated. See \sQuote{Details}.}
   \item{rho}{optional numerical value (for \code{struct="CS"} or \code{struct="HCS"}) or vector (for \code{struct="UN"}) to fix the correlation between levels of the inner factor corresponding to an \code{~ inner | outer} formula specified in the \code{random} argument. A numeric value fixes a particular \eqn{\rho} value, while \code{NA} means that the value should be estimated. See \sQuote{Details}.}
   \item{verbose}{logical indicating whether output should be generated on the progress of the model fitting (default is \code{FALSE}). Can also be an integer. Values > 1 generate more verbose output. See \sQuote{Note}.}
   \item{control}{optional list of control values for the estimation algorithms. If unspecified, default values are defined inside the function. See \sQuote{Note}.}
}
\details{

   \bold{Specifying the Data}

   The function can be used in conjunction with any of the usual effect size or outcome measures used in meta-analyses (e.g., log odds ratios, log relative risks, risk differences, mean differences, standardized mean differences, raw correlation coefficients, correlation coefficients transformed with Fisher's r-to-z transformation, and so on). Simply specify the observed outcomes via the \code{yi} argument and the corresponding sampling variances via the \code{V} argument. In case the sampling errors are correlated, then one can specify the entire variance-covariance matrix of the sampling errors via the \code{V} argument. The \code{\link{escalc}} function can be used to compute a wide variety of effect size and outcome measures (and the corresponding sampling variances) based on summary statistics. Equations for computing the covariance between sampling errors for a variety of different effect size or outcome measures can be found in Gleser and Olkin (2009). For raw and Fisher's r-to-z transformed correlations, one can find suitable equations, for example, in Steiger (1980).

   \bold{Specifying Fixed Effects}

   When \code{method="FE"}, a fixed-effects model is fitted to the data (note: arguments \code{random}, \code{struct}, \code{sigma2}, \code{tau2}, \code{rho}, and \code{R} are not relevant then and are ignored). The model is then simply given by \eqn{\mathbf{y} \sim N(\mathbf{1} \beta_0, \mathbf{V})}{y ~ N(1 \beta_0, V)}, where \eqn{\mathbf{y}}{y} is the (column) vector with the observed effect sizes or outcomes, \eqn{\mathbf{1}}{1} is a column vector of 1's, \eqn{\beta_0} is the (average) true effect size or outcome, and \eqn{\mathbf{V}}{V} is the variance-covariance matrix of the sampling errors (if a vector of sampling variances is provided via the \code{V} argument, then \eqn{\mathbf{V}}{V} is assumed to be diagonal).

   One or more moderators can be included in the model via the \code{mods} argument. A single moderator can be given as a (row or column) vector of length \eqn{k} specifying the values of the moderator. Multiple moderators are specified by giving an appropriate model matrix (i.e., \eqn{\mathbf{X}}{X}) with \eqn{k} rows and \eqn{p'} columns (e.g., using \code{mods = cbind(mod1, mod2, mod3)}, where \code{mod1}, \code{mod2}, and \code{mod3} correspond to the names of the variables for the three moderator variables). The intercept is added to the model matrix by default unless \code{intercept=FALSE}.

   Alternatively, one can use the standard \code{\link{formula}} syntax to specify the model. In this case, the \code{mods} argument should be set equal to a one-sided formula of the form \code{mods = ~ model} (e.g., \code{mods = ~ mod1 + mod2 + mod3}). Interactions, polynomial terms, and factors can be easily added to the model in this manner. When specifying a model formula via the \code{mods} argument, the \code{intercept} argument is ignored. Instead, the inclusion/exclusion of the intercept term is controlled by the specified formula (e.g., \code{mods = ~ mod1 + mod2 + mod3 - 1} would lead to the removal of the intercept term).

   With moderators included, the model is then given by \eqn{\mathbf{y} \sim N(\mathbf{X} \mathbf{\beta}, \mathbf{V})}{y ~ N(X \beta, V)}, where \eqn{\mathbf{X}} denotes the model matrix containing the moderator values (and possibly the intercept) and \eqn{\mathbf{\beta}}{\beta} is a column vector containing the corresponding model coefficients. Fixed-effects models with or without moderators are fitted via generalized/weighted least squares estimation.

   \bold{Specifying Random Effects}

   When \code{method="ML"} or \code{method="REML"}, one can fit random/mixed-effects models to the data by specifying the random effects structure via the \code{random} argument. The \code{random} argument is either a single one-sided formula or a list of one-sided formulas. A formula specified via this argument can be of the form \code{~ 1 | id}. Such a formula adds random effects corresponding to the grouping variable/factor \code{id} to the model. Effects or outcomes with the same value/level of the \code{id} variable/factor receive the same random effect, while effects or outcomes with different values/levels of the \code{id} variable/factor are assumed to be independent. The variance component corresponding to such a formula is denoted by \eqn{\sigma^2}. An arbitrary number of such formulas can be specified as a list of formulas (with variance components \eqn{\sigma^2_1}, \eqn{\sigma^2_2}, and so on).

   Such random effects components are useful to model clustering (i.e., correlation) induced by a multilevel structure in the data (e.g., effects derived from the same paper, lab, research group, or species may be more similar to each other than effects derived from different papers, labs, research groups, or species).

   In addition or alternatively to specifying one or multiple \code{~ 1 | id} terms, the \code{random} argument can also contain one (and only one!) formula of the form \code{~ inner | outer}. Effects or outcomes with different values/levels of the \code{outer} grouping variable/factor are assumed to be independent, while effects or outcomes with the same value/level of the \code{outer} grouping variable/factor receive correlated random effects corresponding to the levels of the \code{inner} grouping variable/factor. The \code{struct} argument is used to specify the variance structure corresponding to the \code{inner} variable/factor. With \code{struct="CS"}, a compound symmetric structure is assumed (i.e., a single variance component \eqn{\tau^2} corresponding to all values/levels of the \code{inner} variable/factor and a single correlation coefficient \eqn{\rho} for the correlation between different values/levels). With \code{struct="HCS"}, a heteroscedastic compound symmetric structure is assumed (with variance components \eqn{\tau^2_1}, \eqn{\tau^2_2}, and so on, corresponding to the values/levels of the \code{inner} variable/factor and a single correlation coefficient \eqn{\rho} for the correlation between different values/levels). Finally, with \code{struct="UN"}, an unstructured variance-covariance matrix is assumed (with variance components \eqn{\tau^2_1}, \eqn{\tau^2_2}, and so on, corresponding to the values/levels of the \code{inner} variable/factor and correlation coefficients \eqn{\rho_{12}}, \eqn{\rho_{13}}, \eqn{\rho_{23}}, and so on, for the various combinations of the values/levels of the \code{inner} variable/factor). \if{html,latex}{For example, for an \code{inner} grouping variable/factor with three levels, the three structures correspond to:}

   \if{html}{\figure{structs.png}{options: width=700}}
   \if{latex}{\figure{structs.png}{options: width=5.4in}}

   With the \code{outer} factor corresponding to a study id variable and the \code{inner} factor corresponding to a variable indicating the treatment type or study arm, such a random effects component could be used to estimate how strongly different true effects or outcomes within the same study are correlated and/or whether the amount of heterogeneity differs across different treatment types/arms. The meta-analytic bivariate model (van Houwelingen, Arends, & Stijnen, 2002) can also be fitted in this manner (see the examples below).

   When the \code{random} argument contains a formula of the form \code{~ 1 | id}, one can use the (optional) argument \code{R} to specify a corresponding known correlation matrix of the random effects (i.e., \code{R = list(id = Cor)}, where \code{Cor} is the correlation matrix). In that case, effects or outcomes with the same value/level of the \code{id} variable/factor receive the same random effect, while effects or outcomes with different values/levels of the \code{id} variable/factor receive random effects that are correlated as specified in the corresponding correlation matrix given via the \code{R} argument. The column/row names of the correlation matrix given via the \code{R} argument must therefore contain all of the values/levels of the \code{id} variable/factor. When the \code{random} argument contains multiple formulas of the form \code{~ 1 | id}, one can specify known correlation matrices for none, some, or all of those terms (e.g., \code{random = list(~ 1 | id1, ~ 1 | id2), R = list(id1 = Cor1)} or \code{random = list(~ 1 | id1, ~ 1 | id2), R = list(id1 = Cor1, id2 = Cor2)}, where \code{Cor1} and \code{Cor2} are the correlation matrices corresponding to the grouping variables/factors \code{id1} and \code{id2}, respectively).

   Random effects with a known (or at least approximately known) correlation structure are useful in a variety of contexts. For example, such components can be used to account for the correlations induced by a shared phylogenetic history among organisms (e.g., plants, fungi, animals). In that case, \code{~ 1 | id} is used to specify the organisms and argument \code{R} is used to specify the phylogenetic correlation matrix of the organisms studied in the meta-analysis. The corresponding variance component then indicates how much variance/heterogeneity is attributable to the specified phylogeny. As another example, in a genetic meta-analysis studying disease association for several single nucleotide polymorphisms (SNPs), linkage disequilibrium (LD) among the SNPs can induce an approximately known degree of correlation among the effects. In that case, \code{~ 1 | id} is used to specify the SNP and \code{R} the corresponding LD correlation map.

   \bold{Fixing Variance Components and/or Correlations}

   Arguments \code{sigma2}, \code{tau2}, and \code{rho} can be used to fix particular variance components and/or correlations at a given value. This is useful for sensitivity analyses (e.g., for plotting the regular/restricted log-likelihood as a function of a particular variance component or correlation) or for imposing a desired variance-covariance structure on the data.

   For example, if \code{random = list(~ 1 | id1, ~ 1 | id2)}, then \code{sigma2} must be of length 2 (corresponding to \eqn{\sigma^2_1} and \eqn{\sigma^2_2}) and a fixed value can be assigned to either or both variance components. Setting a particular component to \code{NA} means that the component will be estimated by the function.

   Argument \code{tau2} is only relevant when the \code{random} argument contains an \code{~ inner | outer} formula. In that case, if the \code{tau2} argument is used, it must be either of length 1 (for \code{struct="CS"}) or of the same length as the number of levels of the inner factor (for \code{struct="HCS"} or \code{struct="UN"}). A numeric value in the \code{tau2} argument then fixes the corresponding variance component to that value, while \code{NA} means that the component will be estimated. Similarly, if argument \code{rho} is used, it must be either of length 1 (for \code{struct="CS"} or \code{struct="HCS"}) or of length \eqn{l(l-1)/2} (for \code{struct="UN"}), where \eqn{l} denotes the number of levels of the inner factor. Again, a numeric value fixes the corresponding correlation, while \code{NA} means that the correlation will be estimated. For example, with \code{struct="CS"} and \code{rho=0}, the variance-covariance matrix of the inner factor will be diagonal with \eqn{\tau^2} along the diagonal.

   \bold{Omnibus Test of Parameters}

   For models including moderators, an omnibus test of all the model coefficients is conducted that excludes the intercept (the first coefficient) if it is included in the model. If no intercept is included in the model, then the omnibus test includes all of the coefficients in the model including the first. Alternatively, one can manually specify the indices of the coefficients to test via the \code{btt} argument. For example, use \code{btt=c(3,4)} to only include the third and fourth coefficient from the model in the test (if an intercept is included in the model, then it corresponds to the first coefficient in the model).

   \bold{Categorical Moderators}

   Categorical moderator variables can be included in the model via the \code{mods} argument in the same way that appropriately (dummy) coded categorical independent variables can be included in linear models. One can either do the dummy coding manually or use a model formula together with the \code{\link{factor}} function to let \R handle the coding automatically.

   \bold{Tests and Confidence Intervals}

   By default, the test statistics of the individual coefficients in the model (and the corresponding confidence intervals) are based on the normal distribution, while the omnibus test is based on a chi-square distribution with \eqn{m} degrees of freedom (\eqn{m} being the number of coefficients tested). As an alternative, one can set \code{tdist=TRUE}, which slightly mimics the Knapp and Hartung (2003) method by using a t-distribution with \eqn{k-p} degrees of freedom for tests of individual coefficients and confidence intervals and an F-distribution with \eqn{m} and \eqn{k-p} degrees of freedom (\eqn{p} being the total number of model coefficients including the intercept if it is present) for the omnibus test statistic.

   \bold{Test for (Residual) Heterogeneity}

   A test for (residual) heterogeneity is automatically carried out by the function. Without moderators in the model, this test is the generalized/weighted least squares extension of Cochran's \eqn{Q}-test, which tests whether the variability in the observed effect sizes or outcomes is larger than one would expect based on sampling variability (and the given covariances among the sampling errors) alone. A significant test suggests that the true effects or outcomes are heterogeneous. When moderators are included in the model, this is the \eqn{Q_E}-test for residual heterogeneity, which tests whether the variability in the observed effect sizes or outcomes that is not accounted for by the moderators included in the model is larger than one would expect based on sampling variability (and the given covariances among the sampling errors) alone.
}
\value{
   An object of class \code{c("rma.mv","rma")}. The object is a list containing the following components:
   \item{b}{estimated coefficients of the model.}
   \item{se}{standard errors of the coefficients.}
   \item{zval}{test statistics of the coefficients.}
   \item{pval}{p-values for the test statistics.}
   \item{ci.lb}{lower bound of the confidence intervals for the coefficients.}
   \item{ci.ub}{upper bound of the confidence intervals for the coefficients.}
   \item{vb}{variance-covariance matrix of the estimated coefficients.}
   \item{sigma2}{estimated \eqn{\sigma^2} value(s).}
   \item{tau2}{estimated \eqn{\tau^2} value(s).}
   \item{rho}{estimated \eqn{\rho} value(s).}
   \item{k}{number of studies included in the model.}
   \item{p}{number of coefficients in the model (including the intercept).}
   \item{m}{number of coefficients included in the omnibus test of coefficients.}
   \item{QE}{test statistic for the test of (residual) heterogeneity.}
   \item{QEp}{p-value for the test of (residual) heterogeneity.}
   \item{QM}{test statistic for the omnibus test of coefficients.}
   \item{QMp}{p-value for the omnibus test of coefficients.}
   \item{int.only}{logical that indicates whether the model is an intercept-only model.}
   \item{yi, V, X}{the vector of outcomes, the corresponding variance-covariance matrix of the sampling errors, and the model matrix of the model.}
   \item{fit.stats}{a list with the log-likelihood, deviance, AIC, BIC, and AICc values.}
   \item{\dots}{some additional elements/values.}

   The results of the fitted model are neatly formated and printed with the \code{\link{print.rma.mv}} function. If fit statistics should also be given, use \code{\link{summary.rma}} (or use the \code{\link{fitstats.rma}} function to extract them).

   For random/mixed-effects models, the \code{\link{profile.rma.mv}} function can be used to obtain a plot of the (restricted) log-likelihood as a function of a specific variance component or correlation parameter of the model.
}
\note{
   Model fitting is done via numerical optimization over the model parameters. By default, \code{\link{optim}} is used with \code{method="Nelder-Mead"}. A different method can be chosen via the \code{control} argument (e.g., \code{control=list(optmethod="BFGS")}). One can also chose a different optimizer via the \code{optimizer} argument (the default is \code{control=list(optimizer="optim")}). Besides \code{\link{optim}}, one can also choose \code{\link{nlminb}} or one of the optimizers in the \code{minqa} package (i.e., \code{\link[minqa]{uobyqa}}, \code{\link[minqa]{newuoa}}, or \code{\link[minqa]{bobyqa}}). Note that the optimizer name should be given as a character string (i.e., in quotes). Any additional control parameters for the optimizer can also be specified via the \code{control} argument (e.g., \code{control=list(optimizer="nlminb", iter.max=500)}).

   At the moment, the starting values are not chosen in a terribly clever way and may be far off. As a result, the optimizer may be slow to converge or may even get stuck in a logical maximum. One can set the starting values manually for the various variance components and correlations in the model via the \code{control} argument by specifying the vectors \code{sigma2.init}, \code{tau2.init}, and/or \code{rho.init} as needed. Especially for complex models, it is a good idea to try out different starting values to make sure that the same estimates are obtained.

   Information on the progress of the optimization algorithm can be obtained by setting \code{verbose=TRUE} or with \code{control=list(verbose=TRUE)}.

   Whether particular variance components and/or correlations are actually identifiable needs to be carefully examined when fitting complex models. The function does some limited checking internally to fix variances and/or correlations at zero when it is clear that insufficient information is available to estimate a particular parameter (e.g., if a particular factor has only a single level, the corresponding variance component cannot be estimated). However, it is strongly advised in general to do post model fitting checks to make sure that the likelihood surface around the ML/REML estimates is not flat for some combination of the parameter estimates (which would imply that the estimates are essentially arbitrary). For example, one can plot the (restricted) log-likelihood as a function of each variance component and correlation in the model to make sure that each profile plot shows a clear peak at the corresponding ML/REML estimate. The \code{\link{profile.rma.mv}} function can be used for this purpose.

   Finally, note that the model fitting is not done in a very efficient manner at the moment, which is partly a result of allowing for correlations among the effect sizes or outcomes across the entire dataset (when using the \code{R} argument). As a result, the function works directly with the entire \eqn{k \times k}{k x k} (marginal) variance-covariance matrix of the effect sizes or outcomes (instead of working with smaller blocks in a block-diagonal structure). As a result, model fitting can be slow for large \eqn{k}. Using \code{verbose=TRUE} will provide information on how the model fitting is progressing.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \cr
   package homepage: \url{http://www.metafor-project.org/} \cr
   author homepage: \url{http://www.wvbauer.com/}
}
\references{
   Gleser, L. J., & Olkin, I. (2009). Stochastically dependent effect sizes. In H. Cooper, L. V. Hedges, & J. C. Valentine (Eds.), \emph{The handbook of research synthesis and meta-analysis} (2nd ed., pp. 357--376). New York: Russell Sage Foundation.

   van Houwelingen, H. C., Arends, L. R., & Stijnen, T. (2002). Advanced methods in meta-analysis: Multivariate approach and meta-regression. \emph{Statistics in Medicine}, \bold{21}, 589--624.

   Steiger, J. H. (1980). Tests for comparing elements of a correlation matrix. \emph{Psychological Bulletin}, \bold{87}, 245--251.

   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \url{http://www.jstatsoft.org/v36/i03/}.
}
\seealso{
   \code{\link{rma.uni}}, \code{\link{rma.mh}}, \code{\link{rma.peto}}, \code{\link{rma.glmm}}
}
\examples{
### load BCG vaccine data
data(dat.bcg)

### calculate log odds ratios and corresponding sampling variances
dat <- escalc(measure="OR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

### random-effects model using rma.uni()
rma(yi, vi, data=dat)

### change data into long format
dat.long <- to.long(measure="OR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

### set levels of group variable ("exp" = experimental/vaccinated; "con" = control/non-vaccinated)
levels(dat.long$group) <- c("exp", "con")

### set "con" to reference level
dat.long$group <- relevel(dat.long$group, ref="con")

### calculate log odds and corresponding sampling variances
dat.long <- escalc(measure="PLO", xi=out1, mi=out2, data=dat.long)

### bivariate random-effects model using rma.mv()
res <- rma.mv(yi, vi, mods = ~ group, random = ~ group | study, struct="UN", data=dat.long)
res

### see help(dat.berkey1998) for another example of a multivariate model
}
\keyword{models}
