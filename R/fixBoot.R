#' Trimmed bootstrap distributions.
#'
#' Remove bootstrap estimates associated with extreme bootstrap distribution values, resulting from arc length methods.
#'
#' @param object An \code{alKDEboot}, \code{mmKDEboot}, \code{kappa4alBoot}, or \code{kappa4nlsBoot}.
#' @param trim #' @param trim The fraction (0 to 0.5) of observations to be trimmed from each end of the jackknife distribution before the mean and standard deviation are computed. Values of \code{trim} outside that range are taken as the nearest endpoint.
#' @param na.rm A logical value indicating whether NA values should be stripped before the computation proceeds.
#'
#' @return The \code{object} is returned with the estimated coefficients corresponding to the trimmed bootstrap distribution.
#'
#' @export
fixBoot <- function(object, trim=0, na.rm=TRUE)
{
object$coefDist <- apply(object$coefDist, 2, trim, trim=trim, na.rm=na.rm)
object$bcoefficients <- colMeans(object$coefDist)
object$se <- apply(object$coefDist, 2, sd)
object
}