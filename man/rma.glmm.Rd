\name{rma.glmm}
\alias{rma.glmm}
\title{Meta-Analysis via the Generalized Linear (Mixed-Effects) Model}
\description{Function to fit the meta-analytic fixed- and random-effects models with or without moderators via the generalized linear (mixed-effects) model. See the documentation of the \pkg{\link{metafor-package}} for more details on these models.}
\usage{
rma.glmm(ai, bi, ci, di, n1i, n2i, x1i, x2i, t1i, t2i, xi, mi, ti, ni,
         mods, measure, intercept=TRUE, data, slab, subset,
         add=1/2, to="only0", drop00=TRUE, vtype="LS",
         model="UM.FS", method="ML", tdist=FALSE,
         level=95, digits=4, btt, nAGQ=100, verbose=FALSE, control)
}
\arguments{
   \item{ai}{see below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{bi}{see below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{ci}{see below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{di}{see below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{n1i}{see below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{n2i}{see below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{x1i}{see below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{x2i}{see below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{t1i}{see below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{t2i}{see below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{xi}{see below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{mi}{see below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{ti}{see below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{ni}{see below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{mods}{optional argument to include one or more moderators in the model. A single moderator can be given as a vector of length \eqn{k} specifying the values of the moderator. Multiple moderators are specified by giving a matrix with \eqn{k} rows and \eqn{p'} columns. Alternatively, a model \code{\link{formula}} can be used to specify the model. See \sQuote{Details}.}
   \item{measure}{character string indicating the outcome measure to use for the meta-analysis. Possible options are the odds ratio (\code{"OR"}), the incidence rate ratio (\code{"IRR"}), the logit transformed proportion (\code{"PLO"}), or the log transformed incidence rate (\code{"IRLN"}).}
   \item{intercept}{logical, indicating whether an intercept term should be added to the model (default is \code{TRUE}).}
   \item{data}{optional data frame containing the data supplied to the function.}
   \item{slab}{optional vector with unique labels for the \eqn{k} studies.}
   \item{subset}{optional vector indicating the subset of studies that should be used for the analysis. This can be a logical vector of length \eqn{k} or a numeric vector indicating the indices of the observations to include.}
   \item{add}{non-negative number indicating the amount to add to zero cells, counts, or frequencies when calculating the individual outcomes. See below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{to}{character string indicating when the values under \code{add} should be added (either \code{"only0"}, \code{"all"}, \code{"if0all"}, or \code{"none"}). See below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{drop00}{logical indicating whether studies with no cases/events (or only cases) in both groups should be dropped. See the documentation of the \code{\link{escalc}} function for more details.}
   \item{vtype}{character string indicating the type of sampling variances to calculate when calculating the individual outcomes. See below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{model}{character string specifying the general model type to use for the analysis (either \code{"UM.FS"} (the default), \code{"UM.RS"}, \code{"CM.EL"}, or \code{"CM.AL"}). See \sQuote{Details}.}
   \item{method}{character string specifying whether a fixed- or a random/mixed-effects model should be fitted. A fixed-effects model (with or without moderators) is fitted when using \code{method="FE"}. Random/mixed-effects models are fitted by setting \code{method="ML"} (the default). See \sQuote{Details}.}
   \item{tdist}{logical specifying whether test statistics and confidence intervals should be based on the normal (when \code{FALSE}, the default) or the t-distribution (when \code{TRUE}). See \sQuote{Details}.}
   \item{level}{numerical value between 0 and 100 specifying the confidence interval level (default is 95).}
   \item{digits}{integer specifying the number of decimal places to which the printed results should be rounded (default is 4).}
   \item{btt}{optional vector of indices specifying which coefficients to include in the omnibus test of moderators. See \sQuote{Details}.}
   \item{nAGQ}{positive integer specifying the number of points per axis for evaluating the adaptive Gauss-Hermite approximation to the log likelihood. This defaults to 100, which should provide sufficient accuracy in the evaluation of the log likelihood in most cases, but at the expense of speed. Setting this to 1 corresponds to the Laplacian approximation. See \sQuote{Note}.}
   \item{verbose}{logical indicating whether output should be generated for the fitting algorithms used (default is \code{FALSE}). See \sQuote{Note}.}
   \item{control}{optional list of control values for the estimation algorithms. If unspecified, default values are defined inside the function. See \sQuote{Note}.}
}
\details{

   \bold{Specifying the Data}

   The function can be used in conjunction with the following effect size or outcome measures:
   \itemize{
   \item \code{measure="OR"} for odds ratios (analyzed in log units)
   \item \code{measure="IRR"} for incidence rate ratios (analyzed in log units)
   \item \code{measure="PLO"} for logit transformed proportions (i.e., log odds)
   \item \code{measure="IRLN"} for log transformed incidence rates.
   }
   The \code{\link{escalc}} function describes the data/arguments that should be specified/used for these measures.

   \bold{Specifying the Model}

   A variety of general model types are available when analyzing 2x2 table data (i.e., when \code{measure="OR"}) or two-group event count data (i.e., when \code{measure="IRR"}).
   \itemize{
   \item \code{model="UM.FS"} for an unconditional generalized linear mixed-effects model with fixed study effects
   \item \code{model="UM.RS"} for an unconditional generalized linear mixed-effects model with random study effects
   \item \code{model="CM.AL"} for a conditional generalized linear mixed-effects model (approximate likelihood)
   \item \code{model="CM.EL"} for a conditional generalized linear mixed-effects model (exact likelihood).
   }
   For \code{measure="OR"}, models \code{"UM.FS"} and \code{"UM.RS"} are essentially (mixed-effects) logistic regression models, while for \code{measure="IRR"}, these models are (mixed-effects) Poisson regression models. A choice must be made on how to model study level variability (i.e., differences in outcomes across studies irrespective of group membership). One can choose between using fixed study effects (which means that \eqn{k} dummy variables are added to the model) or random study effects (which means adding a random effect at the study level).

   The conditional model (\code{model="CM.EL"}) avoids having to model study level variability by conditioning on the total numbers of cases/events in each study. For \code{measure="OR"}, this leads to a non-central hypergeometric distribution for the data within each study and the corresponding model is then a (mixed-effects) conditional logistic model. Fitting this model can be difficult and computationally expensive. When the number of cases in each study is small relative to the group sizes, one can approximate the exact likelihood by a binomial distribution, which leads to a regular (mixed-effects) logistic regression model (\code{model="CM.AL"}). For \code{measure="IRR"}, the conditional model leads directly to a binomial distribution for the data within each study and the resulting model is again a (mixed-effects) logistic regression model (no approximate likelihood model is needed here).

   When analyzing proportions (i.e., when \code{measure="PLO"}) or incidence rates (i.e., when \code{measure="IRLN"}) of individual groups, the model type is always a (mixed-effects) logistic or Poisson regression model, respectively (i.e., the \code{model} argument is not relevant here).

   Aside from choosing the general model type, one has to decide whether to fit a fixed- or random-effects model to the data. A \emph{fixed-effects model} is fitted by setting \code{method="FE"}. A \emph{random-effects model} is fitted by setting \code{method="ML"} (the default). Note that random-effects models with dichotomous data are often referred to as \sQuote{binomial-normal} models in the meta-analytic literature. Analogously, for event count data, such models could be referred to as \sQuote{Poisson-normal} models.

   One or more moderators can be included in all of these models via the \code{mods} argument. A single moderator can be given as a (row or column) vector of length \eqn{k} specifying the values of the moderator. Multiple moderators are specified by giving an appropriate design matrix with \eqn{k} rows and \eqn{p'} columns (e.g., using \code{mods = cbind(mod1, mod2, mod3)}, where \code{mod1}, \code{mod2}, and \code{mod3} correspond to the names of the variables for the three moderator variables). The intercept is included in the model by default unless \code{intercept=FALSE}.

   Alternatively, one can use standard \code{\link{formula}} syntax to specify the model. In this case, the \code{mods} argument should be set equal to a one-sided formula of the form \code{mods = ~ model} (e.g., \code{mods = ~ mod1 + mod2 + mod3}). Interactions, polynomial terms, and factors can be easily added to the model in this manner. When specifying a model formula via the \code{mods} argument, the \code{intercept} argument is ignored. Instead, the inclusion/exclusion of the intercept term is controlled by the specified formula (e.g., \code{mods = ~ mod1 + mod2 + mod3 - 1} would lead to the removal of the intercept term). With moderators, a \emph{fixed-effects with moderators model} is then fitted by setting \code{method="FE"}, while a \emph{mixed-effects model} is fitted by setting \code{method="ML"}.

   \bold{Fixed-, Saturated-, and Random/Mixed-Effects Models}

   When fitting a particular model, actually up to three different models are fitted within the function:
   \itemize{
   \item the fixed-effects model (i.e., where \eqn{\tau^2} is set to 0),
   \item the saturated model (i.e., the model with a deviance of 0), and
   \item the random/mixed-effects model (i.e., where \eqn{\tau^2} is estimated) (only if \code{method="ML"}).
   } The saturated model is obtained by adding as many dummy variables to the model as needed so that the model deviance is equal to zero. Even when \code{method="ML"}, the fixed-effects and saturated models are fitted, as they are used to compute the test statistics for the Wald-type and likelihood ratio tests for (residual) heterogeneity (see below).

   \bold{Omnibus Test of Parameters}

   For models including moderators, an omnibus test of all the model coefficients is conducted that excludes the intercept (the first coefficient) if it is included in the model. If no intercept is included in the model, then the omnibus test includes all of the coefficients in the model including the first. Alternatively, one can manually specify the indices of the coefficients to test via the \code{btt} argument. For example, use \code{btt=c(3,4)} to only include the third and fourth coefficient from the model in the test (if an intercept is included in the model, then it corresponds to the first coefficient in the model).

   \bold{Categorical Moderators}

   Categorical moderator variables can be included in the model in the same way that appropriately (dummy) coded categorical independent variables can be included in linear models. One can either do the dummy coding manually or use a model formula together with the \code{\link{factor}} function to let \R handle the coding automatically.

   \bold{Tests and Confidence Intervals}

   By default, the test statistics of the individual coefficients in the model (and the corresponding confidence intervals) are based on the normal distribution, while the omnibus test is based on a chi-square distribution with \eqn{m} degrees of freedom (\eqn{m} being the number of coefficients tested). As an alternative, one can set \code{tdist=TRUE}, which slightly mimics the Knapp and Hartung (2003) method by using a t-distribution with \eqn{k-p} degrees of freedom for tests of individual coefficients and confidence intervals and an F-distribution with \eqn{m} and \eqn{k-p} degrees of freedom (\eqn{p} being the total number of model coefficients including the intercept if it is present) for the omnibus test statistic.

   \bold{Tests for (Residual) Heterogeneity}

   Two different tests for (residual) heterogeneity are automatically carried out by the function. The first is a Wald-type test, which tests the coefficients corresponding to the dummy variables added in the saturated model for significance. The second is a likelihood ratio test, which tests the same set of coefficients, but does so by comparing the deviance of the fixed-effects and the saturated model. These two tests are not identical for the types of models fitted by the \code{rma.glmm} function and may even lead to conflicting conclusions.

   \bold{Individual Outcomes}

   The various models do not require the calculation of the individual outcome values and directly make use of the table/event counts. Zero cells/events are not a problem (except in extreme cases, such as when one of the two outcomes never occurs or when there are no events in any of the studies). Therefore, it is unnecessary to add some constant to the cell counts (or the number of events) when there are zero cells/events. However, for plotting and various other functions, it is necessary to calculate the individual outcome values for the \eqn{k} studies. Here, zero cells/events can be problematic, so adding a constant value to the cell counts (or the number of events) ensures that all \eqn{k} values can be calculated. The \code{add} and \code{to} arguments are used to specify what value should be added to the cell frequencies (or the number of events) and under what circumstances when calculating the individual outcome values. The documentation of the \code{\link{escalc}} function explains how the \code{add} and \code{to} arguments work. Note that \code{drop00} is set to \code{TRUE} by default, since studies where \code{ai=ci=0} or \code{bi=di=0} or studies where \code{x1i=x2i=0} are uninformative about the size of the effect (the counts for such studies are set to \code{NA}).
}
\value{
   An object of class \code{c("rma.glmm","rma")}. The object is a list containing the following components:
   \item{b}{estimated coefficients of the model.}
   \item{se}{standard errors of the coefficients.}
   \item{zval}{test statistics of the coefficients.}
   \item{pval}{p-values for the test statistics.}
   \item{ci.lb}{lower bound of the confidence intervals for the coefficients.}
   \item{ci.ub}{upper bound of the confidence intervals for the coefficients.}
   \item{vb}{variance-covariance matrix of the estimated coefficients.}
   \item{tau2}{estimated amount of (residual) heterogeneity. Always \code{0} when \code{method="FE"}.}
   \item{sigma2}{estimated amount of study level variability (only for \code{model="UM.RS"}).}
   \item{k}{number of studies included in the model.}
   \item{p}{number of coefficients in the model (including the intercept).}
   \item{m}{number of coefficients included in the omnibus test of coefficients.}
   \item{QE.Wld}{Wald-type test statistic for the test of (residual) heterogeneity.}
   \item{QEp.Wld}{p-value for the Wald-type test of (residual) heterogeneity.}
   \item{QE.LRT}{likelihood ratio test statistic for the test of (residual) heterogeneity.}
   \item{QEp.LRT}{p-value for the likelihood ratio test of (residual) heterogeneity.}
   \item{QM}{test statistic for the omnibus test of coefficients.}
   \item{QMp}{p-value for the omnibus test of coefficients.}
   \item{I2}{value of \eqn{I^2}.}
   \item{H2}{value of \eqn{H^2}.}
   \item{int.only}{logical that indicates whether the model is an intercept-only model.}
   \item{yi, vi, X}{the vector of outcomes, the corresponding sampling variances, and the design matrix of the model.}
   \item{fit.stats}{a list with the log likelihood, deviance, AIC, and BIC values.}
   \item{\dots}{some additional elements/values.}

   The results of the fitted model are neatly formated and printed with the \code{\link{print.rma.uni}} function. If fit statistics should also be given, use \code{\link{summary.rma}} (or use the \code{\link{fitstats.rma}} function to extract them).
}
\note{
   Fitting the various types of models requires several different iterative algorithms:
   \itemize{
   \item For \code{model="UM.FS"} and \code{model="CM.AL"}, iteratively reweighted least squares (IWLS) as implemented in the \code{\link{glm}} function is used for fitting the fixed-effects and the saturated models. For \code{method="ML"}, adaptive Gauss-Hermite quadrature as implemented in the \code{\link[lme4]{lmer}} function is used. The same applies when \code{model="CM.EL"} is used in combination with \code{measure="IRR"} or when \code{measure="PLO"} or \code{measure="IRLN"} (regardless of which general model type is then specified).
   \item For \code{model="UM.RS"}, adaptive Gauss-Hermite quadrature as implemented in the \code{\link[lme4]{lmer}} function is used to fit all of the models.
   \item For \code{model="CM.EL"} and \code{measure="OR"}, the quasi-Newton method (\code{"BFGS"}) as implemented in the \code{\link{optim}} function is used for fitting the fixed-effects and the saturated models. For \code{method="ML"}, the same algorithm is used, together with adaptive quadrature as implemented in the \code{\link{integrate}} function (for the integration over the density of the non-central hypergeometric distribution).
   } When \code{model="CM.EL"} and \code{measure="OR"}, actually \code{model="CM.AL"} is first used to obtain starting values for \code{\link{optim}}, so either 4 (if \code{method="FE"}) or 6 (if \code{method="ML"}) models need to be fitted in total.

   Some important control parameters for the various algorithms can be adjusted via the \code{control} argument. In particular,
   \itemize{
   \item for \code{\link{glm}}, the \code{epsilon} (default is \code{1e-8}) and \code{maxit} (default is \code{25}) arguments,
   \item for \code{\link[lme4]{lmer}}, the \code{maxIter} (default is \code{300}) and \code{maxFN} (default is \code{900}) arguments,
   \item for \code{\link{optim}}, the \code{method} (default is \code{"BFGS"}) and \code{reltol} (default is \code{1e-8}) arguments,
   \item for \code{\link{integrate}}, the \code{rel.tol} (default is \code{1e-8}) and \code{subdivisions} (default is \code{100}) arguments.
   } Also, for \code{\link[lme4]{lmer}}, the \code{nAGQ} argument is used to specify the number of quadrature points. The default value is 100, which should be more than enough for most cases, but this can be quite slow. Setting this to 1 corresponds to the Laplacian approximation.

   Information on the evolution of the various algorithms is obtained by setting \code{verbose=TRUE} or with \code{control=list(verbose=TRUE)}. Since fitting the various models can be computationally expensive, this option is quite useful to determine how the fitting is progressing.

   For \code{model="CM.EL"} and \code{measure="OR"}, optimization involves repeated calculation of the density of the non-central hypergeometric distribution. When \code{method="ML"}, this also requires integration over the same density. This is currently implemented in a rather brute-force manner and may not be numerically stable, especially when models with moderators are fitted. Stability can be improved by scaling the moderators in a similar manner (i.e., don't use a moderator that is coded 0 and 1, while another uses values in the 1000s). Sensitivity analyses are highly recommended here, to ensure that the results do not depend on the scaling of the moderators.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \cr
   package homepage: \url{http://www.metafor-project.org/} \cr
   author homepage: \url{http://www.wvbauer.com/}
}
\references{
   Agresti, A. (2002). \emph{Categorical data analysis} (2nd. ed). Hoboken, NJ: Wiley.

   Bagos, P. G., & Nikolopoulos, G. K. (2009). Mixed-effects Poisson regression models for meta-analysis of follow-up studies with constant or varying durations. \emph{The International Journal of Biostatistics}, \bold{5}(1), article 21.

   van Houwelingen, H. C., Zwinderman, K. H., & Stijnen, T. (1993). A bivariate approach to meta-analysis. \emph{Statistics in Medicine}, \bold{12}, 2273--2284.

   Stijnen, T., Hamza, T. H., & Ozdemir, P. (2010). Random effects meta-analysis of event outcome in the framework of the generalized linear mixed model with applications in sparse data. \emph{Statistics in Medicine}, \bold{29}, 3046--3067.

   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \url{http://www.jstatsoft.org/v36/i03/}.
}
\seealso{
   \code{\link{rma.uni}}, \code{\link{rma.mh}}, \code{\link{rma.peto}}
}
\examples{
### load BCG vaccine data
data(dat.bcg)

### random-effects model using rma.uni() (standard RE model analysis)
rma(measure="OR", ai=tpos, bi=tneg, ci=cpos, di=cneg,
    data=dat.bcg, method="REML", verbose=TRUE)

### random-effects models using rma.glmm() (requires 'lme4' package)

### unconditional model with fixed study effects
### note: fitting this model is slow with the default of nAGQ=100
\dontrun{
rma.glmm(measure="OR", ai=tpos, bi=tneg, ci=cpos, di=cneg, 
         data=dat.bcg, model="UM.FS", verbose=TRUE)
}

### unconditional model with random study effects
### note: fitting this model is slow with the default of nAGQ=100
\dontrun{
rma.glmm(measure="OR", ai=tpos, bi=tneg, ci=cpos, di=cneg,
         data=dat.bcg, model="UM.RS", verbose=TRUE)
}

### conditional model with approximate likelihood
rma.glmm(measure="OR", ai=tpos, bi=tneg, ci=cpos, di=cneg,
         data=dat.bcg, model="CM.AL", verbose=TRUE)

### conditional model with exact likelihood
### note: fitting this model is slow
\dontrun{
rma.glmm(measure="OR", ai=tpos, bi=tneg, ci=cpos, di=cneg,
         data=dat.bcg, model="CM.EL", verbose=TRUE)
}
}
\keyword{models}
