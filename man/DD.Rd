% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kdeGaussMom.R
\name{DD}
\alias{DD}
\title{Calculate derivatives.}
\usage{
DD(expr, name, order)
}
\arguments{
\item{expr}{An R expression of the function to be differentiated.}

\item{name}{A qoted string containing the name of the variable with respect to which the derivative is calculated.}

\item{order}{The order of the derivative saught.}
}
\value{
An R expression of the calculated derivative.
}
\description{
Calculate higher order derivatives for use in exact moment calculations.
}
\examples{
DD(expression(exp(m*t+(1/2)*(s^2)*(t^2))), "t", 1) ## m
DD(expression(exp(m*t+(1/2)*(s^2)*(t^2))), "t", 2) ## m^2+s^2

}
