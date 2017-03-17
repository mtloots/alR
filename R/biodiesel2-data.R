#' Data on the Measurement of the oxidative stability of biodiesel.
#'
#' Data from an experiment at 0.33\%, 0\%, and 0.67\% concentrations of antioxidants Orox PK, Naugard P and Anox
#' 20 respectively, for measuring the oxidative stability of biodiesel.  4080 measurements were recorded.
#' The experiment was conducted by Isbe van der Westhuizen en Yong Fah Mian using equipment from the CSIR, and was
#'  sponsored by Walter Focke from the Institute of Applied Materials, at the University of Pretoria.
#'
#' @docType data
#'
#' @usage data(biodiesel2)
#'
#' @format A \code{data.frame} with three variables:
#' \describe{
#' \item{x}{The time measured in hours.}
#' \item{y}{The conductivity value measured in S/cm.}
#' \item{S60}{A logical indicator giving a 60\% split for validation purposes.}
#' }
#'
#' @keywords datasets
#'
#' @examples
#' library(alR)
#' str(biodiesel2)
#' summary(biodiesel2)
#' plot(y~x, data=biodiesel2)
'biodiesel2'