% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{GPDInt}
\alias{GPDInt}
\alias{GPDInt2}
\title{Arc length of GPD PDF.}
\usage{
GPDInt(mu, sigma, alpha, q1, q2, quantile)

GPDInt2(mu, sigma, alpha, q1, q2, quantile)
}
\arguments{
\item{mu}{A real number specifying the location parameter.}

\item{sigma}{A strictly positive real number specifying the scale parameter.}

\item{alpha}{A real number specifying the shape parameter.}

\item{q1}{The point (or vector for \code{GPDInt2}) specifying the lower limit of the arc length integral.}

\item{q2}{The point (or vector for \code{GPDInt2}) specifying the upper limit of the arc length integral.}

\item{quantile}{Logical, TRUE/FALSE, whether \code{q1} and \code{q2} are quantiles, or actual points in the domain.}
}
\value{
GPDInt: A list with the following components:
\itemize{
\item value: The resultant arc length.
\item abs.err: The absolute error between iterations.
subdivisions: Number of subdivisions used in the numerical approximation.
\item neval: Number of function evaluations used by the numerical approximation.
}

GPDInt2: A vector having length equal to that of the vector of lower quantile bounds, containing the arc lengths requested for a generalised Pareto probability density function.
}
\description{
Calculate the arc length for a univariate generalised Pareto probability density function over a specified interval.
}
\details{
The arc length of a univariate generalised Pareto probability density function is approximated using the numerical integration C code implimented for R's integrate functions, i.e. using Rdqags.  For this approximation, subdiv = 100 (100 subdivisions), and eps_abs = eps_rel = 1e-10, i.e. the absolute and relative errors respectively.
}
\examples{
library(alR)
mu <- 0
sigma <- 1
alpha <- 2
GPDInt(mu, sigma, alpha, 0.025, 0.975, TRUE)
GPDInt(mu, sigma, alpha, 0.001, 0.5, FALSE)

GPDInt2(mu, sigma, alpha, c(0.025, 0.5), c(0.5, 0.975), TRUE)
GPDInt2(mu, sigma, alpha, c(-1.96, 0), c(0, 1.96), FALSE)

}
