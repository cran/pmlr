\name{pmlr}
\alias{pmlr}
\alias{summary.pmlr}
\title{Penalized maximum likelihood estimation for multinomial logistic regression using the Jeffreys prior}
\description{
Extends the approach proposed by Firth (1993) for bias reduction of MLEs in exponential family models to the multinomial logistic regression model with general covariate types.  Modification of the logistic regression score function to remove first-order bias is equivalent to penalizing the likelihood by the Jeffreys prior, and yields penalized maximum likelihood estimates (PLEs) that always exist.  Hypothesis testing is conducted via likelihood ratio statistics.  Profile confidence intervals (CI) are constructed for the PLEs.
}
\usage{
pmlr(formula, data, weights = NULL, alpha = 0.05, penalized = TRUE, 
  method = c("likelihood", "wald")[1], joint = FALSE) 
}
\arguments{                             
  \item{formula}{an object of class \code{\link{formula}} (or one that can be coerced to that class): a symbolic description of the model to be fitted.  Typically, \code{formula} has the form \code{response ~ terms} where \code{response} is either a factor with \eqn{J+1} levels (the first is used as the baseline category) or a \eqn{J}-column indicator matrix, and \code{terms} is a series of terms which specifies a linear predictor for the response.}
  \item{data}{a data frame containing the variables in the model}
  \item{weights}{an optional vector of weights to be used in the fitting process. Should be a numeric vector of counts in the case of frequency data.}
  \item{alpha}{the significance level (default is \eqn{\alpha = 0.05}).}
  \item{penalized}{a logical variable indicating whether penalized maximum likelihood should be used (default)}
  \item{method}{a character string specifying whether p-values and confidence intervals should be based on the profile likelihood or the Wald method. Must be one of \dQuote{likelihood} (default) or \dQuote{wald}.}
  \item{joint}{a logical variable indicating whether joint hypothesis tests should be performed in addition to individual parameter tests. If \code{TRUE}, \eqn{H_0: \beta_{1i} = \cdots = \beta_{Ji} = 0} and \eqn{H_0: \beta_{1i} = \cdots = \beta_{Ji}} are also tested for covariate \eqn{i, i = 1, \ldots, P}.}

}
\details{
Logistic regression is one of the most widely used regression models in practice, but alternatives to conventional maximum likelihood estimation methods may be more appropriate for small or sparse samples.  It is well known that the usual maximum likelihood estimates (MLEs) of the log odds ratios are biased in finite samples, and there is a non-zero probability that an MLE is infinite (i.e., does not exist).  This corresponds to the problem of separation (Lesaffre and Albert, 1989).   

In this package, we extend the approach proposed by Firth (1993) for bias reduction of MLEs in exponential family models to the multinomial logistic regression model with general covariate types.  Modification of the logistic regression score function to remove first order bias is equivalent to penalizing the likelihood by the Jeffreys prior, and yields penalized likelihood estimates (PLEs) that always exist, even in samples in which MLEs are infinite.  PLEs are an attractive alternative in small to moderate-sized samples, and are preferred to exact conditional MLEs when there are continuous covariates.     
 
We consider a multicategory outcome \eqn{y} that is a multinomial variable with \eqn{J + 1} categories.  For each category \eqn{j (j = 1, \ldots, J)} there is a regression function in which the log odds of response in category \eqn{j}, relative to category 0, is a linear function of regression parameters and a vector \eqn{\bold{x}} of \eqn{P} covariates (including a constant):  \eqn{\log\{\mathrm{prob}(y = j| \bold{x})/\mathrm{prob}(y = 0 | \bold{x})\} = \bold{\beta}_{j}^T \bold{x}}.  Let \eqn{\bold{y}_i} be a \eqn{J \times 1} vector of indicators for the observed response category for observation \eqn{i}, with the corresponding \eqn{J \times 1} vector of probabilities \eqn{\bold{\Theta}_{i} = (\Theta_{i1}, \ldots, \Theta_{iJ})^T}.  The vector of MLEs, \eqn{\hat{B} = vec[(\bold{\hat{\beta}}_1, \ldots, \bold{\hat{\beta}}_J)^T]}, is estimated from observations \eqn{(\bold{y}_i,\bold{x}_i), i = 1, \ldots, n}, by solving the score equations of the log-likelihood \eqn{l(B)}.  We denote the score function by \eqn{U(B)} and the \eqn{PJ \times PJ} Fisher information matrix by \eqn{A = A(B)}.

The order \eqn{n^{-1}} bias of estimates based on the usual likelihood \eqn{L(B)} is removed by applying the penalty \eqn{|A|^{1/2}}, and basing estimation on the penalized likelihood \eqn{L^*(B) = L(B) |A|^{1/2}}.  The vector \eqn{\hat{B}^*} of penalized estimates (PLEs) is the solution to the score equations of the penalized log-likelihood \eqn{l^*(B) = l(B) + \frac{1}{2} \log{|A|}}.  The introduction of bias into the score function through the penalty removes the leading term in the asymptotis bias of the MLEs.  The modified score function proposed by Firth for the binomial logistic model extends directly to the multinomial model as \eqn{U^*(B) = U(B) - A \, b(B)}.  The bias term \eqn{b(B)} is the leading term in the asymptotic bias of the multinomial MLEs, obtained from the Taylor series expansion of the log-likelihood (Cox and Snell, 1968) and is a function of the matrix of third derivatives of \eqn{l(B)} with respect to \eqn{B} (Bull et al., 2002).  

The PLEs are obtained by a modified Fisher scoring algorithm.  Using \eqn{t} to denote the iteration number, the modified iterative updating equations are:
\deqn{B^*_{(t+1)} = B^*_{(t)} + A^{-1}_{(t)}U^*(B^*_{(t)}) = B^*_{(t)} + b(B^*_{(t)}) + A^{-1}_{(t)}U(B^*_{(t)}).}
Thus, in comparison to the usual updating equation used to obtain the MLEs, the score function modification operates by applying the asymptotic bias correction at each step in the iterative process.  This prevents estimates from going off to infinity and failing to converge when there is separation in the data.  

Symmetric Wald-type CIs for \eqn{\beta_{jp}} can be constructed using \eqn{Var^{* 1/2}(\hat{\beta}^*_{jp})}, obtained from the inverse of \eqn{A^*}, however, performance is expected to be poor in situations where separation is likely to occur.  Asymmetric CIs for the PLEs can be constructed from the profile log-likelihood for \eqn{\beta_{jp}}, which is the function \eqn{l^*_0(B(s))}, where \eqn{B(s)} is the argument that maximizes \eqn{l^*} under the single-parameter constraint \eqn{H_0: \beta_{jp} = s}.  The \eqn{100(1 - \alpha)\%} CI for \eqn{\beta_{jp}} is given by all parameter values that are compatible with the data (i.e., all \eqn{s} such that the likelihood ratio statistic \eqn{LR_P(s) \le q}, where \eqn{q} is the \eqn{1 - \alpha} percentile of the \eqn{\chi^2} distribution).  This is equivalent to \eqn{l^*_0(B(s)) \ge l^*(\hat{B}^*) - \frac{1}{2} q}.  The endpoints of the interval are then found by numerically solving the equality for values of \eqn{s}.  Based on the algorithm employed in SAS PROC LOGISTIC for MLEs, our method for finding these roots does not require computing \eqn{l^*_0(B(s))}, which in itself would involve maximizing \eqn{l^*(B)}, but proceeds directly by solving the constrained optimization problem: maximize \eqn{l^*(B)} such that \eqn{l^*_(B) = l^*(\hat{B}^*) - \frac{1}{2} q} and \eqn{\beta_{jp} = s}.  We, however, modify the starting values in the iterative scheme used by SAS to obtain a new algorithm that is slower, but simple and more robust (see Bull et al. 2007 for details).
}
\value{\code{pmlr} returns an object of class \dQuote{\code{pmlr}} with the following components:
  \item{coefficients}{an array containing the coefficients of the \eqn{P} parameters for the \eqn{J} categories}
  \item{var}{an array containing the variance-covariance matrices for the \eqn{J} categories}
  \item{stat}{an array containing the test statistics from the individual parameter tests}
  \item{pval}{an array containing the p-values from the individual parameter tests}
  \item{CI.lower}{an array containing the lower confidence limits from the individual parameter tests}
  \item{CI.upper}{an array containing the upper confidence limits from the individual parameter tests}
  \item{separation}{an array indicating coefficients for which separation occured} 
  \item{stat.joint}{an array containing the test statistics from the joint tests}
  \item{pval.joint}{an array containing the p-values from the joint tests}
  \item{call}{the matched call}
  \item{method}{a character string indicating the method used to compute p-values and confidence limits}
  \item{joint}{a logical variable indicating whether the joint hypothesis tests were performed}
}
\note{This implementation is not a true scoring or Newton-type algorithm because it updates with the inverse of \eqn{A}, the Fisher information matrix for the MLEs, rather than the information for the PLEs, \eqn{A^*}, which includes an additional term corresponding to the second derivatives of the penalty: \eqn{\frac{1}{2} \log |A|}.  As a result, convergence using the modified scoring algorithm for the PLEs is slower than a scoring algorithm based on \eqn{A^*}, which converges at a quadratic rate.  For well behaved and larger datasets this usually means no more than 2 or 3 steps beyond that required for the MLEs.  Starting values of \eqn{\beta{jp} = 0} are used and are generally satisfactory.  For smaller datasets (i.e., less than 50 observations) and especially datasets in which there are infinite MLEs, convergence is slower and could take up to 35 or 40 iterations.  In datasets with separations, starting values other than 0 can lead to divergence problems.}

\note{The first argument to the \code{pmlr} function is a formula of the form \code{response ~ terms} where \code{response} can be either a \eqn{J}-column indicator matrix or a factor with \eqn{J + 1} levels.  In the case of frequency data (e.g., the \code{hepatitis} dataset) the baseline category is determined by the \eqn{J}-column indicator matrix with baseline category 0.  In the case of data with individual records (e.g., the \code{enzymes} dataset) the baseline category is determined by the outcome, which is coded as a factor.  The first level (i.e., \code{levels(x)[1]} is taken to be the baseline.  Note that the default levels attribute of a factor is the unique set of values taken \code{as.character(x)}, sorted into increasing order of \code{x}.  For example, the enzymes dataset has categories 1, 2, and 3 with baseline category 1.  The baseline category can be changed to category 2 via \code{enzymes$Group <- factor(enzymes$Group, levels-c("2","1","3"))}.}

\references{
Bull, S. B., Greenwood, C. M. T., and Hauck, W. W. (1997) Jacknife bias reduction for polychotomous logistic regression. \emph{Statistics in Medicine}, \bold{16}, 545--560; Bull, S. B. (1997) Correction. \emph{Statistics in Medicine}, \bold{16}, 2928.  

Bull, S. B., Mak, C. and Greenwood, C. M. T. (2002) A modified score function estimator for multinomial logistic regression in small samples. \emph{Computational Statistics & Data Analysis}, \bold{39}, 57--74.

Bull, S. B., Lewinger, J. P., Lee, S. S. F. (2005) Penalized maximum likelihood estimation for multinomial logistic regression using the Jeffreys prior.  \emph{Technical Report No. 0505}, Department of Statistics, University of Toronto.

Bull, S. B., Lewinger, J. P. and Lee, S. S. F. (2007) Confidence intervals for multinomial logistic regression in sparse data. \emph{Statistics in Medicine}, \bold{26}, 903--918. 

Cox, D. R. and Snell, E. J. (1968) A general definition of residuals. \emph{Journal of the Royal Statistical Society, Series B}, \bold{30}, 248--275.

Firth, D. (1993) Bias reduction of maximum likelihood estimates. \emph{Biometrika}, \bold{80}, 27--38.

Lesaffre, E. and Albert, A. (1989) Partial separation in logistic discrimination. \emph{Journal of the Royal Statistical Society, Series B}         , \bold{51}, 109--116.

SAS Institute Inc. (1999) \emph{SAS OnlineDoc, Version 8, The LOGISTIC Procedure}, Confidence Intervals for Parameters, Chapter 39, Section 26, Cary, NC.
}
\examples{
# As reported in Bull et al. (2007)

data(hepatitis)
fit <- pmlr(cbind(HCV, nonABC) ~ group + time + group:time, 
  data = hepatitis, weights = counts)
summary(fit)

data(enzymes)
# Exclude patients in Group 4 (post-necrotic cirrhosis)
enzymes <- enzymes[enzymes$Group != 4,]
# Center and scale covariates
AST <- scale(log(enzymes$AST))
ALT <- scale(log(enzymes$ALT))
GLDH <- scale(log(enzymes$GLDH))
OCT <- scale(log(enzymes$OCT))
enzymes <- data.frame(Patient = enzymes$Patient,
  Group = enzymes$Group, AST, ALT, GLDH, OCT)
# Remove 10 observations to create separation
enzymes <- enzymes[-c(9, 18, 33, 58, 61, 77, 94, 97, 99, 100),]

# Multinomial: acute viral hepatitis and aggressive chronic hepatits
# vs. persistent chronic hepatitis
# Assign Group 2 (persistent chronic hepatitis) as baseline category
enzymes$Group <- factor(enzymes$Group, levels=c("2","1","3"))
fit <- pmlr(Group ~ AST + GLDH, data = enzymes)
summary(fit)

# Binomial: Acute viral hepatitis vs. persistent chronic hepatitis
# Exclude patients in Group 3 (agressive chronic hepatitis)
enzymes.1vs2 <- enzymes[enzymes$Group != 3,]
# Assign Group 2 (persistent chronic hepatitis) as baseline category
enzymes.1vs2$Group <- factor(enzymes.1vs2$Group, levels=c("2","1"))
fit <- pmlr(Group ~ AST + GLDH, data = enzymes.1vs2)
summary(fit)

# Binomial: Aggressive chronic hepatitis vs. persistent chronic hepatitis
# Exclude patients in Group 1 (acute viral hepatitis)
enzymes.3vs2 <- enzymes[enzymes$Group != 1,]
# Assign Group 2 (persistent chronic hepatitis) as baseline category
enzymes.3vs2$Group <- factor(enzymes.3vs2$Group, levels=c("2","3"))
fit <- pmlr(Group ~ AST + GLDH, data = enzymes.3vs2)
summary(fit)
}
\keyword{models}
\keyword{regression}
\keyword{htest}
