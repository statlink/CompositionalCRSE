\name{jack.boot}
\alias{jack.boot}
\title{
Cluster jackknife covariance matrix of the regression coefficients in compositional regression models
}

\description{
Delete-one-cluster (or ordinary delete-one) jackknife estimate of the covariance matrix of the
regression coefficients of several regression models for compositional data, with optional
parallel computing.
}

\usage{
jack.boot(y, x, type = "kl", a = 1, id = NULL, center = "mean", ncores = 1)
}

\arguments{
\item{y}{
A matrix with the compositional data (the response variable). The rows must sum to 1.
}
\item{x}{
A vector or a numerical matrix with the predictor variables.
Factors or character variables are not allowed.
}
\item{type}{
The regression model to be fitted. Possible values are: "lm" for multivariate
linear regression on the additive log-ratio transformed compositional data,
"kl" for the Kullback-Leibler divergence based regression (\code{\link[Compositional]{kl.compreg}}),
"diri" for Dirichlet regression (\code{\link[Compositional]{diri.reg.irls}}), "zadr"
for zero-adjusted Dirichlet regression (\code{zadr.irls} of the package \bold{CompositionalZADR}),
"areg" for the \eqn{\alpha}-regression (\code{areg} of the package \bold{CompositionalSR}) and
"cln" for the conditional logistic normal regression (\code{cln.reg} of the package
\bold{Compositionalcln}).
}
\item{a}{
The value of the power transformation \eqn{\alpha}, used only when type = "areg".
This is treated as fixed and is not re-selected within each jackknife sample.
}
\item{id}{
A vector with the cluster indicators, one value for every observation. If NULL (default)
each observation is treated as its own cluster and the ordinary delete-one jackknife
is performed.
}
\item{center}{
The point around which the jackknife estimates are centred. If "mean" (default), the mean of the
leave-one-cluster-out estimates is used (the classical jackknife). If "full", the estimate
from the full sample is used, which gives a slightly larger (more conservative) covariance matrix.
}
\item{ncores}{
The number of cores to use. If greater than 1, the jackknife samples are run in parallel
using \code{\link[parallel]{parSapply}}, and the package required by the selected model is
loaded on every worker.
}
}

\details{
Each cluster (the whole group of observations sharing the same value of \code{id}) is removed
in turn and the model is refitted to the remaining data, giving \eqn{G} leave-one-cluster-out
estimates \eqn{\hat{\bm{\beta}}_{(g)}}, where \eqn{G} is the number of clusters. The covariance
matrix is then
\deqn{\frac{G-1}{G}\sum_{g=1}^G\left(\hat{\bm{\beta}}_{(g)}-\bm{c}\right)\left(\hat{\bm{\beta}}_{(g)}-\bm{c}\right)^\top,}
where \eqn{\bm{c}} is either the mean of the \eqn{\hat{\bm{\beta}}_{(g)}} or the full sample estimate
(argument \code{center}). Since the clusters are kept intact, the within-cluster correlation is
accounted for and no working correlation structure needs to be specified. The procedure is
deterministic, so no number of resamples is required, and exactly \eqn{G} refits are performed.
When \code{id} is NULL this reduces to the ordinary delete-one jackknife with \eqn{n} refits,
which may be slow for the iterative models when the sample size is large.

The regression coefficients are vectorised column-wise, so the first \eqn{p+1} elements
correspond to the coefficients (constant and \eqn{p} predictor variables) of the second
component, the next \eqn{p+1} to the third component and so on. The row and column names
of the covariance matrix reflect this ordering. For type = "cln" the last component is the
reference, so the blocks correspond to the first up to the penultimate component.

Leave-one-cluster-out samples for which the model cannot be fitted are discarded with a warning,
and \eqn{G} becomes the number of successful fits. For type = "lm" the compositional data must
contain no zero values.

The cluster jackknife is more conservative than the cluster bootstrap and the sandwich estimator,
and is recommended when the number of clusters is small or the cluster sizes are unbalanced.
}

\value{
A list including:
\item{be}{
The matrix of the estimated regression coefficients, fitted to the full data.
}
\item{covb}{
The jackknife covariance matrix of the vectorised regression coefficients.
}
\item{G}{
The number of clusters (successful leave-one-cluster-out fits) used.
}
}

\references{
Efron B. and Tibshirani R.J. (1993). An Introduction to the Bootstrap. Chapman and Hall/CRC.

Bell R.M. and McCaffrey D.F. (2002). Bias reduction in standard errors for linear regression
with multi-stage samples. Survey Methodology, 28(2): 169--181.

MacKinnon J.G., Nielsen M.O. and Webb M.D. (2023). Cluster-robust inference: A guide to
empirical practice. Journal of Econometrics, 232(2): 272--299.

Tsagris M. and Stewart C. (2018). A Dirichlet regression model for compositional data with zeros.
Lobachevskii Journal of Mathematics, 39(3): 398--412.

Alzeley O. and Tsagris M. (2026). Modelling compositional data with structural zero values
https://arxiv.org/pdf/2608.29954

Aitchison J. (1986). The statistical analysis of compositional data, pages 272--274.

Murteira Jose MR, and Joaquim JS Ramalho (2016). Regression analysis of multivariate fractional data.
Econometric Reviews 35(4): 515-552.

Alenazi A. A. (2022). f-divergence regression models for compositional data.
Pakistan Journal of Statistics and Operation Research, 18(4): 867--882.

Maier Marco J. (2014) DirichletReg: Dirichlet Regression for Compositional Data in R.
Research Report Series/Department of Statistics and Mathematics, 125. WU Vienna University of
Economics and Business, Vienna.

Mardia K.V., Kent J.T., and Bibby J.M. (1979). Multivariate analysis. Academic press.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

\seealso{
\code{ \link{clust.boot} }
}

\examples{
set.seed(1)
G <- 30  ;  m <- 5  ;  n <- G * m
id <- rep(1:G, each = m)
x <- matrix( rnorm(n * 2), ncol = 2 )
u <- matrix( rnorm(G * 2, sd = 0.5), ncol = 2 )[id, ]
be <- rbind( c(0.5, -0.3), c(1, 0), c(0, -0.8) )
mu <- cbind( 1, exp( cbind(1, x) \%*\% be + u ) )
mu <- mu / rowSums(mu)
y <- t( apply(mu, 1, function(p) rgamma(3, 20 * p) ) )
y <- y / rowSums(y)

## cluster jackknife
mod <- jack.boot(y, x, type = "kl", id = id)
mod$be
sqrt( diag(mod$covb) )

## centred at the full sample estimate
mod2 <- jack.boot(y, x, type = "kl", id = id, center = "full")
sqrt( diag(mod2$covb) )
}
