#' Calculate derivatives.
#'
#' Calculate higher order derivatives for use in exact moment calculations.
#'
#' @param expr An R expression of the function to be differentiated.
#' @param name A qoted string containing the name of the variable with respect to which the derivative is calculated.
#' @param order The order of the derivative saught.
#'
#' @return An R expression of the calculated derivative.
#' @importFrom stats D
#'
#' @examples
#' DD(expression(exp(m*t+(1/2)*(s^2)*(t^2))), "t", 1) ## m
#' DD(expression(exp(m*t+(1/2)*(s^2)*(t^2))), "t", 2) ## m^2+s^2
#'
DD <- function(expr, name, order)
{
if(order==1) D(expr, name)
else DD(D(expr, name), name, order-1)
}

#' Non-central Gaussian KDE moments.
#'
#' Calculate the non-central moments for a univariate Gaussian kernel density estimator.
#'
#' @param n The highest order of non-central moments saught.
#' @param mu A vector of data points on which the kernel density estimator is based.
#' @param h The kernel density estimator bandwidth.
#'
#' @return A vector of length \code{n} of non-central moments requested for a Gaussian kernel density estimator.
#'
#' @examples
#' library(alR)
#' x <- rnorm(100)
#' h_x <- bw(x, type=-1)
#' kdeGaussMom(3, x, h_x)
#'
#' @export
kdeGaussMom <- function(n, mu, h)
{
(1/length(mu))*do.call(rbind, lapply(1:n, function(i) sum(eval(DD(expression(exp(m*t+(1/2)*(s^2)*(t^2))), "t", i), list(t=0, m=mu, s=h)))))
}