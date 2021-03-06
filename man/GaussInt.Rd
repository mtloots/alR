% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{GaussInt}
\alias{GaussInt}
\alias{GaussInt2}
\title{Arc length of Gaussian PDF.}
\usage{
GaussInt(mu, sigma, q1, q2, quantile)

GaussInt2(mu, sigma, q1, q2, quantile)
}
\arguments{
\item{mu}{A real number specifying the location parameter.}

\item{sigma}{A positive real number specifying the scale parameter.}

\item{q1}{The point (or vector for \code{GaussInt2}) specifying the lower limit of the arc length integral.}

\item{q2}{The point (or vector for \code{GaussInt2}) specifying the upper limit of the arc length integral.}

\item{quantile}{Logical, TRUE/FALSE, whether \code{q1} and \code{q2} are quantiles, or actual points in the domain.}
}
\value{
GaussInt: A list with the following components:
\itemize{
\item value: The resultant arc length.
\item abs.err: The absolute error between iterations.
subdivisions: Number of subdivisions used in the numerical approximation.
\item neval: Number of function evaluations used by the numerical approximation.
}

GaussInt2: A vector having length equal to that of the vector of lower quantile bounds, containing the arc lengths requested for a Gaussian probability density function.
}
\description{
Calculate the arc length for a univariate Gaussian probability density function over a specified interval.
}
\details{
The arc length of a univariate Gaussian probability density function is approximated using the numerical integration C code implimented for R's integrate functions, i.e. using Rdqags.  For this approximation, subdiv = 100 (100 subdivisions), and eps_abs = eps_rel = 1e-10, i.e. the absolute and relative errors respectively.
}
\examples{
library(alR)
mu <- 2
sigma <- 3.5
GaussInt(mu, sigma, 0.025, 0.975, TRUE)
GaussInt(mu, sigma, -1.96, 1.96, FALSE)

GaussInt2(mu, sigma, c(0.025, 0.5), c(0.5, 0.975), TRUE)
GaussInt2(mu, sigma, c(-1.96, 0), c(0, 1.96), FALSE)

}
