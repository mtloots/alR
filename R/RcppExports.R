# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' @rdname alEfit
#' @return alE: A list with the following components (see \code{\link{optim}}):
#' \itemize{
#' \item par: The estimated parameters.
#' \item abstol: The absolute tolerance level (default 1e-15).
#' \item fail: An integer code indicating convergence.
#' \item fncount: Number of function evaluations.
#' }
#'
#' @examples
#' x <- rnorm(1000)
#' alE(x,0.025, 0.975, TRUE, -1)
#' alE(x,c(0.025, 0.5), c(0.5, 0.975), TRUE, -1)
#' alE(x,0.025, 0.975, FALSE, -1)
#' alE(x,c(0.025, 0.5), c(0.5, 0.975), FALSE, -1)
#'
#' @export
alE <- function(x, q1, q2, dc, type, distribution = 1L) {
    .Call(`_alR_alE`, x, q1, q2, dc, type, distribution)
}

#' @rdname alEfit
#' @return alEfitdist: A matrix of parameter estimates resulting from the estimated arc lengths over the specified interval(s), i.e. the bootstrap distribution for the estimated parameters resulting from the chosen sample arc length statistic.
#' @examples
#' \dontrun{
#' alEfitdist(x, 0.025, 0.975, TRUE, -1, 100)
#' alEfitdist(x, 0.025, 0.975, FALSE, -1, 100)
#' }
#' @export
alEfitdist <- function(x, q1, q2, dc, type, bootstraps, distribution = 1L) {
    .Call(`_alR_alEfitdist`, x, q1, q2, dc, type, bootstraps, distribution)
}

#' @rdname alEtest
#' @return alEdist: A vector (matrix) of arc lengths over the specified interval(s), i.e. the simulated distribution for the chosen sample arc length statistic.
#' @examples
#' \dontrun{
#' alEdist(50, 100, 2, 3.5, 0.025, 0.975, TRUE, TRUE, -1)
#' alEdist(50, 100, 2, 3.5, c(0.025,0.5), c(0.5,0.975), TRUE, TRUE, -1)
#' alEdist(50, 100, 2, 3.5, 0.025, 0.975, TRUE, FALSE, -1)
#' alEdist(50, 100, 2, 3.5, c(0.025,0.5), c(0.5,0.975), TRUE, FALSE, -1)
#' alEdist(50, 100, 2, 3.5, qnorm(0.025,2,3.5),
#' qnorm(0.975, 2, 3.5), FALSE, FALSE, -1)
#' alEdist(50, 100, 2, 3.5, c(qnorm(0.025, 2, 3.5),2),
#' c(2,qnorm(0.975, 2, 3.5)), FALSE, FALSE, -1)
#' }
#' @export
alEdist <- function(n, bootstraps, mu, sigma, q1, q2, quantile, dc, type, distribution = 1L) {
    .Call(`_alR_alEdist`, n, bootstraps, mu, sigma, q1, q2, quantile, dc, type, distribution)
}

#' Objective function for KDE arc length matching.
#'
#' The arc lengths over specified intervals, in the domain of kernel density estimates, are matched.
#'
#' @param beta Vector of regression parameters.
#' @param gamma Design matrix.
#' @param aly A numerical vector of arc length segments to be matched to (same length as \code{beta}).
#' @param q1 A vector of points (not quantiles) specifying the lower limit of the arc length segments.
#' @param q2 A vector of points (not quantiles) specifying the upper limit of the arc length segments.
#' @param type An integer specifying the bandwidth selection method, see \code{\link{bw}}.
#'
#' @return Square root of the sum of squared differences between \code{gamma}*\code{beta} and \code{aly} (Eucledian distance).
#'
#' @export
alrKDE <- function(beta, gamma, aly, q1, q2, type) {
    .Call(`_alR_alrKDE`, beta, gamma, aly, q1, q2, type)
}

#' @rdname bw
#' @return Silverman: Bandwidth estimator based on Silverman's rule of thumb.
#'
#' @examples
#' set.seed(1)
#' x <- rnorm(100)
#' Silverman(x)
#'
#' @export
Silverman <- function(x) {
    .Call(`_alR_Silverman`, x)
}

#' @rdname bw
#' @return Silverman2: Bandwidth estimator based on Silverman's adapted rule of thumb.
#'
#' @examples
#' Silverman2(x)
#'
#' @export
Silverman2 <- function(x) {
    .Call(`_alR_Silverman2`, x)
}

#' Kernel density bandwidth estimators.
#'
#' Calculate the bandwidth estimator using various methods.
#'
#' @param x A vector of data points.
#' @param type One of the following options:
#' \itemize{
#' \item -1: Silverman's rule of thumb.
#' \item -2: Silverman's adapted rule of thumb.
#' \item >0: The real number is returned without any calculations.
#' }
#'
#' @return bw: Bandwidth estimator based on the selected method.
#'
#' @examples
#' bw(x, -1)
#' bw(x, -2)
#' bw(x, 0.5)
#'
#' @export
bw <- function(x, type) {
    .Call(`_alR_bw`, x, type)
}

#' Arc length of Gaussian PDF.
#'
#' Calculate the arc length for a univariate Gaussian probability density function over a specified interval.
#'
#' The arc length of a univariate Gaussian probability density function is approximated using the numerical integration C code implimented for R's integrate functions, i.e. using Rdqags.  For this approximation, subdiv = 100 (100 subdivisions), and eps_abs = eps_rel = 1e-10, i.e. the absolute and relative errors respectively.
#'
#' @param mu A real number specifying the location parameter.
#' @param sigma A positive real number specifying the scale parameter.
#' @param q1 The point (or vector for \code{GaussInt2}) specifying the lower limit of the arc length integral.
#' @param q2 The point (or vector for \code{GaussInt2}) specifying the upper limit of the arc length integral.
#' @param quantile Logical, TRUE/FALSE, whether \code{q1} and \code{q2} are quantiles, or actual points in the domain.
#'
#' @return GaussInt: A list with the following components:
#' \itemize{
#' \item value: The resultant arc length.
#' \item abs.err: The absolute error between iterations.
#' subdivisions: Number of subdivisions used in the numerical approximation.
#' \item neval: Number of function evaluations used by the numerical approximation.
#' }
#'
#' @examples
#' library(alR)
#' mu <- 2
#' sigma <- 3.5
#' GaussInt(mu, sigma, 0.025, 0.975, TRUE)
#' GaussInt(mu, sigma, -1.96, 1.96, FALSE)
#'
#' @export
GaussInt <- function(mu, sigma, q1, q2, quantile) {
    .Call(`_alR_GaussInt`, mu, sigma, q1, q2, quantile)
}

#' @rdname GaussInt
#' @return GaussInt2: A vector having length equal to that of the vector of lower quantile bounds, containing the arc lengths requested for a Gaussian probability density function.
#'
#' @examples
#' GaussInt2(mu, sigma, c(0.025, 0.5), c(0.5, 0.975), TRUE)
#' GaussInt2(mu, sigma, c(-1.96, 0), c(0, 1.96), FALSE)
#'
#' @export
GaussInt2 <- function(mu, sigma, q1, q2, quantile) {
    .Call(`_alR_GaussInt2`, mu, sigma, q1, q2, quantile)
}

#' The Generalised Pareto Distribution (GPD).
#'
#' Density, distribution function, quantile function and random generation for the generalised Pareto distribution (GPD) with location parameter \code{mu}, scale parameter \code{sigma}, and shape parameter \code{alpha}. .
#'
#' @param x,q A data point, or quantile, at which the GPD should be evaluated (may be vectors).
#' @param p A probability, at which the GPD should be evaluated.
#' @param mu The location parameter of the GPD.
#' @param sigma The scale parameter of the GPD (should be strictly positive).
#' @param alpha The shape parameter of the GPD (valid for all real numbers).
#' @param n The number of random numbers which should be generated for the GPD.
#' @rdname GPD
#' @return dGPD: The estimated value of the density function of the GPD at the point \code{x}.
#' @examples
#' library(alR)
#' dGPD(0.5, 0, 1, 2)
#' @export
dGPD <- function(x, mu, sigma, alpha) {
    .Call(`_alR_dGPD`, x, mu, sigma, alpha)
}

#' @rdname GPD
#' @examples
#' pGPD(0.5, 0, 1, 2)
#' @return pGPD: The value of the cumulative distribution function of the GPD at the point \code{x}.
#' @export
pGPD <- function(q, mu, sigma, alpha) {
    .Call(`_alR_pGPD`, q, mu, sigma, alpha)
}

#' @rdname GPD
#' @examples
#' qGPD(0.5, 0, 1, 2)
#' @return qGPD: The \code{x}th quantile of the GPD.
#' @export
qGPD <- function(p, mu, sigma, alpha) {
    .Call(`_alR_qGPD`, p, mu, sigma, alpha)
}

#' @rdname GPD
#' @examples
#' rGPD(10, 0, 1, 2)
#' @return rGPD: \code{n} random numbers from the GPD.
#' @export
rGPD <- function(n, mu, sigma, alpha) {
    .Call(`_alR_rGPD`, n, mu, sigma, alpha)
}

#' Arc length of GPD PDF.
#'
#' Calculate the arc length for a univariate generalised Pareto probability density function over a specified interval.
#'
#' The arc length of a univariate generalised Pareto probability density function is approximated using the numerical integration C code implimented for R's integrate functions, i.e. using Rdqags.  For this approximation, subdiv = 100 (100 subdivisions), and eps_abs = eps_rel = 1e-10, i.e. the absolute and relative errors respectively.
#'
#' @param mu A real number specifying the location parameter.
#' @param sigma A strictly positive real number specifying the scale parameter.
#' @param alpha A real number specifying the shape parameter.
#' @param q1 The point (or vector for \code{GPDInt2}) specifying the lower limit of the arc length integral.
#' @param q2 The point (or vector for \code{GPDInt2}) specifying the upper limit of the arc length integral.
#' @param quantile Logical, TRUE/FALSE, whether \code{q1} and \code{q2} are quantiles, or actual points in the domain.
#'
#' @return GPDInt: A list with the following components:
#' \itemize{
#' \item value: The resultant arc length.
#' \item abs.err: The absolute error between iterations.
#' subdivisions: Number of subdivisions used in the numerical approximation.
#' \item neval: Number of function evaluations used by the numerical approximation.
#' }
#'
#' @examples
#' library(alR)
#' mu <- 0
#' sigma <- 1
#' alpha <- 2
#' GPDInt(mu, sigma, alpha, 0.025, 0.975, TRUE)
#' GPDInt(mu, sigma, alpha, 0.001, 0.5, FALSE)
#'
#' @export
GPDInt <- function(mu, sigma, alpha, q1, q2, quantile) {
    .Call(`_alR_GPDInt`, mu, sigma, alpha, q1, q2, quantile)
}

#' @rdname GPDInt
#' @return GPDInt2: A vector having length equal to that of the vector of lower quantile bounds, containing the arc lengths requested for a generalised Pareto probability density function.
#'
#' @examples
#' GPDInt2(mu, sigma, alpha, c(0.025, 0.5), c(0.5, 0.975), TRUE)
#' GPDInt2(mu, sigma, alpha, c(-1.96, 0), c(0, 1.96), FALSE)
#'
#' @export
GPDInt2 <- function(mu, sigma, alpha, q1, q2, quantile) {
    .Call(`_alR_GPDInt2`, mu, sigma, alpha, q1, q2, quantile)
}

#' The qML Method for the Generalised Pareto Distribution (GPD).
#'
#' Estimation of the parameters of the generalised Pareto distribution (GPD) with location parameter \code{mu}, scale parameter \code{sigma}, and shape parameter \code{alpha}.  using the quasi-maximum-likelihood method.
#'
#' @param x A vector of data points.
#' @rdname GPDqML
#' @return GPDMLE:  A named vector containing the maximum likelihood estimates for the GPD fitted to the data \code{x}.
#' @examples
#' x <- rGPD(100, 0, 1, 2)
#' GPDMLE(x)
#' @export
GPDMLE <- function(x) {
    .Call(`_alR_GPDMLE`, x)
}

#' @rdname GPDqML
#' @return GPDqML:  A named vector containing the estimated parameters for the GPD fitted to the data \code{x}, using quasi-maximum likelihood.
#' @examples
#' GPDqML(x)
#' @export
GPDqML <- function(x) {
    .Call(`_alR_GPDqML`, x)
}

#' @rdname GPDqML
#' @return GPDMSE:  A named vector containing the estimated parameters for the GPD fitted to the data \code{x}, using maximum spacings estimation (MSE).
#' @examples
#' \dontrun{
#' GPDMSE(x)
#' }
#' @export
GPDMSE <- function(x) {
    .Call(`_alR_GPDMSE`, x)
}

#' Four-parameter kappa distribution.
#'
#' Functions for the four-parameter kappa distribution.
#'
#' @param x A data point, or quantile, at which the four-parameter kappa distribution should be evaluated.
#' @param mu A real value representing the location of the distribution.
#' @param sigma A positive real number representing the scale parameter of the distribution.
#' @param h,k Real numbers representing shape parameters of the distribution.
#' @param n Number of random variates to generate.
#' @rdname kappa4
#' @examples
#' pkappa4(1, 1, 2, 0.5, 2)
#' @return pkappa4: The cumulative distribution function at the point \code{x}.
#' @export
pkappa4 <- function(x, mu, sigma, h, k) {
    .Call(`_alR_pkappa4`, x, mu, sigma, h, k)
}

#' @rdname kappa4
#' @return dkappa4: The density function at the point x.
#' @examples
#' dkappa4(1, 1, 2, 0.5, 2)
#' @export
dkappa4 <- function(x, mu, sigma, h, k) {
    .Call(`_alR_dkappa4`, x, mu, sigma, h, k)
}

#' @rdname kappa4
#' @examples
#' qkappa4(0.25, 1, 2, 0.5, 2)
#' @return qkappa4: The \code{x}th quantile of the distribution.
#' @export
qkappa4 <- function(x, mu, sigma, h, k) {
    .Call(`_alR_qkappa4`, x, mu, sigma, h, k)
}

#' @rdname kappa4
#' @examples
#' rkappa4(10, 1, 2, 0.5, 2)
#' @return rkappa4: Randomly generated numbers from the distribution.
#' @export
rkappa4 <- function(n, mu, sigma, h, k) {
    .Call(`_alR_rkappa4`, n, mu, sigma, h, k)
}

#' @rdname kappa4
#' @return dddkappa4: The second derivative of dkappa4.
#' @examples
#' dddkappa4(1, 1, 2, 0.5, 2)
#' @export
dddkappa4 <- function(x, mu, sigma, h, k) {
    .Call(`_alR_dddkappa4`, x, mu, sigma, h, k)
}

#' @rdname kappa4
#' @return kappa4cond: The resultant induction period (IP).
#' @examples
#' kappa4cond(1, 2, 0.5, 2)
#' @export
kappa4cond <- function(mu, sigma, h, k) {
    .Call(`_alR_kappa4cond`, mu, sigma, h, k)
}

#' @rdname kappa4
#' @return kappa4tc: A list with the following components:
#' \itemize{
#' \item $par: The k shape parameter corresponding to a given h parameter for the time-conductivity problem.
#' \item $abstol: The absolute tolerance for the numerical optimisation.
#' \item $fail: A code relating to the optimisation routine.
#' \item $fncount: Number of function calls.
#' }
#' @examples
#' kappa4tc(-4, 0, 1)
#' @export
kappa4tc <- function(h, mu, sigma) {
    .Call(`_alR_kappa4tc`, h, mu, sigma)
}

#' Arc length of four-parameter kappa CDF.
#'
#' Calculate the arc length for a univariate four-parameter kappa cumulative distribution function over a specified interval.
#'
#' The arc length of a univariate four-parameter kappa cumulative distribution function is approximated using the numerical integration C code implimented for R's integrate functions, i.e. using Rdqags.  For this approximation, subdiv = 100 (100 subdivisions), and eps_abs = eps_rel = 1e-10, i.e. the absolute and relative errors respectively.
#'
#' @param mu A real number specifying the location parameter.
#' @param sigma A positive real number specifying the scale parameter.
#' @param h,k Real numbers specifying the two shape parameters.
#' @param tau A real number between 0 and 1, corresponding to the CDF value at the point of truncation.
#' @param q1 The point (or vector for \code{kappa4Int2}) specifying the lower limit of the arc length integral.
#' @param q2 The point (or vector for \code{kappa4Int2}) specifying the upper limit of the arc length integral.
#' @param quantile Logical, TRUE/FALSE, whether \code{q1} and \code{q2} are quantiles, or actual points in the domain.
#'
#' @return kappa4Int: A list with the following components:
#' \itemize{
#' \item value: The resultant arc length.
#' \item abs.err: The absolute error between iterations.
#' subdivisions: Number of subdivisions used in the numerical approximation.
#' \item neval: Number of function evaluations used by the numerical approximation.
#' }
#'
#' @examples
#' library(alR)
#' mu <- 4
#' sigma <- 0.4
#' h <- -4
#' tau <- 1 ## no truncation
#' k <- kappa4tc(-4, 0, 1)$par
#' kappa4Int(mu, sigma, h, k, tau, 0.025, 0.975, TRUE)
#' p1 <- qkappa4(0.025, mu, sigma, h, k)
#' p2 <- qkappa4(0.975, mu, sigma, h, k)
#' kappa4Int(mu, sigma, h, k, tau, p1, p2, FALSE)
#'
#' @export
kappa4Int <- function(mu, sigma, h, k, tau, q1, q2, quantile) {
    .Call(`_alR_kappa4Int`, mu, sigma, h, k, tau, q1, q2, quantile)
}

#' @rdname kappa4Int
#' @return kappa4Int2: A vector having length equal to that of the vector of lower quantile bounds, containing the arc lengths requested for a four-parameter kappa cumulative distribution function.
#'
#' @examples
#' kappa4Int2(mu, sigma, h, k, tau, c(0.025, 0.5), c(0.5, 0.975), TRUE)
#' p12 <- qkappa4(0.5, mu, sigma, h, k)
#' kappa4Int2(mu, sigma, h, k, tau, c(p1, p12), c(p12, p2), FALSE)
#'
#' @export
kappa4Int2 <- function(mu, sigma, h, k, tau, q1, q2, quantile) {
    .Call(`_alR_kappa4Int2`, mu, sigma, h, k, tau, q1, q2, quantile)
}

#' Sample arc length statistic.
#'
#' The arc length over a specified interval is calculated for use in non-linear estimation.
#'
#' @param x Numeric vector of independent outcomes.
#' @param y Numeric vector of dependent outcomes \eqn{y=F(x)}.
#' @param q1,q2 Quantiles (between 0 and 1) over which the arc length segment is to be computed.
#' @return kappa4IntApprox: The resultant arc length.
#' @param quantile Logical, TRUE/FALSE, whether \code{q1} and \code{q2} are quantiles, or actual points in the domain.
#' @examples
#' x <- rnorm(100)
#' y <- pnorm(x)
#' kappa4IntApprox(x, y, 0.025, 0.975, TRUE)
#' kappa4IntApprox(x, y, -1.96, 1.96, FALSE)
#'
#' @export
kappa4IntApprox <- function(x, y, q1, q2, quantile) {
    .Call(`_alR_kappa4IntApprox`, x, y, q1, q2, quantile)
}

#' @rdname kappa4IntApprox
#' @return kappa4IntApprox2: A vector having length equal to that of the vector of lower quantile bounds, containing the discrete arc length segments over the specified intervals.
#'
#' @examples
#' kappa4IntApprox2(x, y, c(0.025, 0.5), c(0.5, 0.975), TRUE)
#' kappa4IntApprox2(x, y, c(-1.96, 0), c(0, 1.96), FALSE)
#'
#' @export
kappa4IntApprox2 <- function(x, y, q1, q2, quantile) {
    .Call(`_alR_kappa4IntApprox2`, x, y, q1, q2, quantile)
}

#' Sigmoidal curve fitting.
#'
#' Support functions for fitting four-parameter kappa sigmoidal curves.
#'
#' @param parms A numeric vector of parameters to be estimated.
#' @param  xvec A numeric vector of independent observations.
#' @param y A numeric vector of dependent observations.
#' @param x_min The minimum xvec value.
#' @param x_max The maximum xvec value.
#' @param q1,q2 Numeric vectors, for the lower and upper bounds of the intervals over which arc lengths are to be computed.
#' @param al_samp The sample arc length statistic.
#'
#' @return kappa4NLSobj: The nonlinear least squares objective function.
#' @export
kappa4NLSobj <- function(parms, xvec, y, x_min, x_max) {
    .Call(`_alR_kappa4NLSobj`, parms, xvec, y, x_min, x_max)
}

#' @rdname kappa4NLSobj
#' @return kappa4NLScon: A vector with three conditions evaluated.
#' @export
kappa4NLScon <- function(parms, xvec, y, x_min, x_max) {
    .Call(`_alR_kappa4NLScon`, parms, xvec, y, x_min, x_max)
}

#' @rdname kappa4NLSobj
#' @return kappa4NLShin: A vector specifying a single nonlinear inequality constraint.
#' @export
kappa4NLShin <- function(parms, xvec, y, x_min, x_max) {
    .Call(`_alR_kappa4NLShin`, parms, xvec, y, x_min, x_max)
}

#' @rdname kappa4NLSobj
#' @return kappa4NLSheq: A vector specifying two nonlinear equality constraints.
#' @export
kappa4NLSheq <- function(parms, xvec, y, x_min, x_max) {
    .Call(`_alR_kappa4NLSheq`, parms, xvec, y, x_min, x_max)
}

#' @rdname kappa4NLSobj
#' @return kappa4ALobj: The arc length objective function.
#' @export
kappa4ALobj <- function(parms, al_samp, x_min, x_max, q1, q2) {
    .Call(`_alR_kappa4ALobj`, parms, al_samp, x_min, x_max, q1, q2)
}

#' @rdname kappa4NLSobj
#' @return kappa4ALcon: A vector with three conditions evaluated.
#' @export
kappa4ALcon <- function(parms, al_samp, x_min, x_max, q1, q2) {
    .Call(`_alR_kappa4ALcon`, parms, al_samp, x_min, x_max, q1, q2)
}

#' @rdname kappa4NLSobj
#' @return kappa4ALhin: A vector specifying a single nonlinear inequality constraint.
#' @export
kappa4ALhin <- function(parms, al_samp, x_min, x_max, q1, q2) {
    .Call(`_alR_kappa4ALhin`, parms, al_samp, x_min, x_max, q1, q2)
}

#' @rdname kappa4NLSobj
#' @return kappa4ALheq: A vector specifying two nonlinear equality constraints.
#' @export
kappa4ALheq <- function(parms, al_samp, x_min, x_max, q1, q2) {
    .Call(`_alR_kappa4ALheq`, parms, al_samp, x_min, x_max, q1, q2)
}

#' Gaussian kernel density estimator.
#'
#' Estimate a density function using a kernel density estimator with a Gaussian kernel.
#'
#' The cumulative distribution function is calculated using the numerical integration C code implimented for R's integrate functions, i.e. using Rdqagi.  For this approximation, subdiv = 100 (100 subdivisions), and eps_abs = eps_rel = 1e-10, i.e. the absolute and relative errors respectively.
#'
#' The quantiles of the Gaussian kernel density estimator are calculated using Brent's method.  This method requires an interval in which a solution is saught.  The objective funcion for which a zero is saught is \code{\link{dkdeGauss}}-\code{x}, where \code{x} is the quantile saught.  The first interval in which a solution is searched for, corresponds to the range of \code{mu}, and is expanded in multiples thereof in consequtive steps.  The maximum number of iterations is set at 1000, and the accuracy saught between iterations, is set at 1e-10.
#'
#' @param x A data point, or quantile, at which the kernel density estimator should be evaluated.
#' @param mu A vector of data points on which the kernel density estimator is based.
#' @param h The kernel density estimator bandwidth.
#' @rdname kdeGauss
#' @return dkdeGauss: The estimated value of the density function at the point x.
#' @examples
#' library(alR)
#' x <- rnorm(100)
#' h_x <- bw(x, type=1)
#' dkdeGauss(0, x, h_x)
#' @export
dkdeGauss <- function(x, mu, h) {
    .Call(`_alR_dkdeGauss`, x, mu, h)
}

#' @rdname kdeGauss
#' @examples
#' pkdeGauss(0, x, h_x)
#' @return pkdeGauss: The estimated value of the cumulative distribution function at the point \code{x}.
#' @export
pkdeGauss <- function(x, mu, h) {
    .Call(`_alR_pkdeGauss`, x, mu, h)
}

#' @rdname kdeGauss
#' @examples
#' qkdeGauss(0.5, x, h_x)
#' @return qkdeGauss: A list with the following components:
#' \itemize{
#' \item result: The \code{x}th quantile of the Gaussian kernel density estimator.
#' \item value: The value of the cumulative distribution function of the Gaussian kernel density estimator at the \code{x}th quantile.
#' \item obj.fun: The value of the objective function resulting from Brent's method; should be less than 1e-10.
#' \item iterations: Number of iterations for Brent's method in order to achieve the desired accuracy.
#' \item steps: Number of range expansions of the search boundaries for Brent's method.
#' }
#' @export
qkdeGauss <- function(x, mu, h) {
    .Call(`_alR_qkdeGauss`, x, mu, h)
}

#' Arc length of Gaussian KDE.
#'
#' Calculate the arc length for a univariate Gaussian kernel density estimator over a specified interval.
#'
#' For \code{kdeGaussInt} and \code{kdeGaussInt2}, the arc length of a univariate Gaussian kernel density estimator is approximated using the numerical integration C code implimented for R's integrate functions, i.e. using Rdqags.  For this approximation, subdiv = 100 (100 subdivisions), and eps_abs = eps_rel = 1e-10, i.e. the absolute and relative errors respectively.
#'
#' For \code{kdeGaussIntApprox}, the arc length is approximated by constructing the KDE, and then calculated as the sum of a finite collection of straight lines, based on the Pythagorean theorem.
#'
#' @param mu A vector of data points on which the kernel density estimator is based.
#' @param h The kernel density estimator bandwidth.
#' @param q1 The point (or vector for \code{kdeGaussInt2}) specifying the lower limit of the arc length integral.
#' @param q2 The point (or vector for \code{kdeGaussInt2}) specifying the upper limit of the arc length integral.
#' @param quantile Logical, TRUE/FALSE, whether \code{q1} and \code{q2} are quantiles, or actual points in the domain.
#'
#' @return kdeGaussInt: A list with the following components:
#' \itemize{
#' \item value: The resultant arc length.
#' \item abs.err: The absolute error between iterations.
#' subdivisions: Number of subdivisions used in the numerical approximation.
#' \item neval: Number of function evaluations used by the numerical approximation.
#' }
#'
#' @examples
#' library(alR)
#' mu <- rnorm(100)
#' h <- bw(mu, type=1)
#' kdeGaussInt(mu, h, 0.025, 0.975, TRUE)
#' kdeGaussInt(mu, h, -1.96, 1.96, FALSE)
#'
#' @export
kdeGaussInt <- function(mu, h, q1, q2, quantile) {
    .Call(`_alR_kdeGaussInt`, mu, h, q1, q2, quantile)
}

#' @rdname kdeGaussInt
#' @return kdeGaussInt2: A vector having length equal to that of the vector of lower quantile bounds, containing the arc lengths requested for a Gaussian kernel density estimator.
#'
#' @examples
#' kdeGaussInt2(mu, h, c(0.025, 0.5), c(0.5, 0.975), TRUE)
#' kdeGaussInt2(mu, h, c(-1.96, 0), c(0, 1.96), FALSE)
#'
#' @export
kdeGaussInt2 <- function(mu, h, q1, q2, quantile) {
    .Call(`_alR_kdeGaussInt2`, mu, h, q1, q2, quantile)
}

#' @rdname kdeGaussInt
#' @return kdeGaussIntApprox: The resultant arc length.
#' @examples
#' kdeGaussIntApprox(mu, h, 0.025, 0.975, TRUE)
#' kdeGaussIntApprox(mu, h, -1.96, 1.96, FALSE)
#'
#' @export
kdeGaussIntApprox <- function(mu, h, q1, q2, quantile) {
    .Call(`_alR_kdeGaussIntApprox`, mu, h, q1, q2, quantile)
}

#' @rdname kdeGaussInt
#' @return kdeGaussIntApprox2: A vector having length equal to that of the vector of lower quantile bounds, containing the discrete arc lengths requested for a Gaussian kernel density estimator.
#'
#' @examples
#' kdeGaussIntApprox2(mu, h, c(0.025, 0.5), c(0.5, 0.975), TRUE)
#' kdeGaussIntApprox2(mu, h, c(-1.96, 0), c(0, 1.96), FALSE)
#'
#' @export
kdeGaussIntApprox2 <- function(mu, h, q1, q2, quantile) {
    .Call(`_alR_kdeGaussIntApprox2`, mu, h, q1, q2, quantile)
}

#' Objective function for KDE moment matching.
#'
#' The moments of kernel density estimates are matched.
#'
#' @param beta Vector of regression parameters.
#' @param gamma Design matrix.
#' @param momy A numerical vector of moments to be matched to (same length as \code{beta}).
#' @param kdeGaussMom A R function object to be passed that calculates exact moments.
#' @param type An integer specifying the bandwidth selection method used, see \code{\link{bw}}.
#'
#' @return Square root of sum of squared differences between \code{gamma}*\code{beta} and \code{momy} (Eucledian distance).
#'
#' @export
momKDE <- function(beta, gamma, momy, kdeGaussMom, type) {
    .Call(`_alR_momKDE`, beta, gamma, momy, kdeGaussMom, type)
}

#' Point corresponding to quantile.
#' Calculate the \eqn{x} point corresponding to a quantile \eqn{F(x)} using linear interpolation.
#'
#' @param x A numeric vector, specifying the \eqn{x} values.
#' @param y A numeric vector, specifying the \eqn{F(x)} values.
#' @param q A real number between 0 and 1 inclusive, specifying the desired quantile.
#'
#' @return The interpolated quantile, \eqn{x}, corresponding to \eqn{q=F(x)}.
#' @examples
#' x <- rnorm(100)
#' y <- pnorm(x)
#' qlin(x, y, 0.5)
#' @export
qlin <- function(x, y, q) {
    .Call(`_alR_qlin`, x, y, q)
}

#' Empirical sample quantile.
#' Calculate empirical sample quantile.
#'
#' @param x A numeric vector, specifying the sample for which the quantile is to be calculated.
#' @param q A real number between 0 and 1 inclusive, specifying the desired quantile.
#'
#' @return The empirical quantile of the provided sample.
#' @examples
#' x<-rnorm(100)
#' qsamp(x, 0.5)
#' @export
qsamp <- function(x, q) {
    .Call(`_alR_qsamp`, x, q)
}

#' Sorted vector index.
#'
#' The sorted vector is returned along with the original index of the vector it belonged to.
#'
#' @param x A real-valued vector to be sorted.
#'
#' @return A list with two components:
#' \itemize{
#' \item sorted: The sorted version of \code{x}.
#' \item index: The index of the \eqn{i^{th}} element in \code{x}.
#' }
#'
#' @examples
#' pairSort(c(5, 2, 6))
#' @export
pairSort <- function(x) {
    .Call(`_alR_pairSort`, x)
}

#' LULU smoother.
#'
#' Performs LULU smoothing of the provided vector.
#'
#' @param x A real-valued vector.
#'
#' @return The LULU-smoothed version of \code{x}.
#'
#' @examples
#' x <- rnorm(10)
#' lulu(x)
#' @export
lulu <- function(x) {
    .Call(`_alR_lulu`, x)
}

