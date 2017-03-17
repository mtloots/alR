#' Trimmed jackknife distributions.
#'
#' Remove jackknife estimates associated with extreme jackknife distribution values, resulting from arc length methods.
#'
#' @param object An \code{alKDEjack}, \code{mmKDEjack}, \code{kappa4alJack}, or \code{kappa4nlsJack}.
#' @param trim #' @param trim The fraction (0 to 0.5) of observations to be trimmed from each end of the jackknife distribution before the mean and standard deviation are computed. Values of \code{trim} outside that range are taken as the nearest endpoint.
#' @param na.rm A logical value indicating whether NA values should be stripped before the computation proceeds.
#'
#' @return The \code{object} is returned with the estimated coefficients corresponding to the trimmed bootstrap distribution.
#'
#' @export
fixJack <- function(object, trim, na.rm=TRUE)
{
object$coefDist <- apply(object$coefDist, 2, trim, trim=trim, na.rm=na.rm)
coefMean <- colMeans(object$coefDist)
n <- nrow(object$coefDist)
coefVar <- (1/(n))*apply(object$coefDist, 2, sd)
object$se <- sqrt(coefVar)
object$bias <- (n-1)*(coefMean-object$coefficients)
object$jcoefficients <- object$coefficients-object$bias

object
}