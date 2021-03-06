% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mmKDEshort.R
\name{mmKDEshort}
\alias{mmKDEshort}
\title{Moment matching for kernel density estimators.}
\usage{
mmKDEshort(formula, data = list(), xin, type, ...)
}
\arguments{
\item{formula}{An LHS ~ RHS formula, specifying the linear model to be estimated.}

\item{data}{A data.frame which contains the variables in \code{formula}.}

\item{xin}{Numeric vector of length equal to the number of independent variables, of initial values, for the parameters to be estimated.}

\item{type}{An integer specifying the bandwidth selection method used, see \code{\link{bw}}.}

\item{...}{Arguments to be passed on to the control argument of the \code{\link{optim}} function.}
}
\value{
mmKDEshort: A list with the following components:
\itemize{
\item coefficients: A vector of estimated coefficients.
\item error: The value of the objective function.
}
}
\description{
Estimate parameters of a linear model by matching the moments of kernel density estimators.
}
\details{
A shortened version of \code{\link{mmKDE}}.
}
\examples{
x <- 1:10
y <- x+rnorm(10)
XIn <- lm(y~x)
mmKDEshort(y~x, xin=coef(XIn), type=-1)

}
