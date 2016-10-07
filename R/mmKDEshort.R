#' Moment matching for kernel density estimators.
#'
#' Estimate parameters of a linear model by matching the moments of kernel density estimators.
#'
#' A shortened version of \code{mmKDE}.
#'
#' @param formula An LHS ~ RHS formula, specifying the linear model to be estimated.
#' @param data A data.frame which contains the variables in \code{formula}.
#' @param lower,upper Numeric vectors, of length equal to the number of independent variables, for the lower and upper bounds for the parameters to be estimated.
#' @param itermax Number of iterations for the Differential Evolution algorithm.
#' @param type An integer specifying the bandwidth selection method used, see \code{\link{bw}}.
#' @param ... Arguments to be passed on to \code{DEoptim.control()} of the Differential Evolution algorithm.
#'
#' @importFrom DEoptim DEoptim DEoptim.control
#' @importFrom stats model.frame model.matrix model.response
#' @return mmKDEshort: A list with the following components:
#' \itemize{
#' \item coefficients: A vector of estimated coefficients.
#' \item error: The value of the objective function.
#' }
#' 
#' @examples
#' x <- 1:10
#' y <- x+rnorm(10)
#' mmKDEshort(y~x, lower=c(-2,2), upper=c(2,2), itermax=50, type=-1)
#'
#' @export
mmKDEshort <- function(formula, data=list(), lower, upper, itermax, type, ...)
{
mf <- model.frame(formula=formula, data=data)
X <- model.matrix(attr(mf, "terms"), data=mf)
y <- model.response(mf)

h_y <- bw(y, type)
MOMy <- kdeGaussMom(ncol(X), y, h_y)

mom <- DEoptim(momKDE, lower=lower, upper=upper, control=DEoptim.control(trace=FALSE, itermax=itermax, strategy=6, ...), gamma=X, momy=MOMy, kdeGaussMom, type=type)

coefficients <- mom$optim$bestmem

labels <- if(attr(attr(mf, "terms"), "intercept") == 1) c("Intercept", attr(attr(mf, "terms"), "term.labels")) else attr(attr(mf, "terms"), "term.labels")

names(coefficients) <- labels

error <- mom$optim$bestval

momout <- list(coefficients=coefficients,
error=error)

momout
}