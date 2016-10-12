#' Moment matching for kernel density estimators.
#'
#' Estimate parameters of a linear model by matching the moments of kernel density estimators.
#'
#' A shortened version of \code{\link{mmKDE}}.
#'
#' @param formula An LHS ~ RHS formula, specifying the linear model to be estimated.
#' @param data A data.frame which contains the variables in \code{formula}.
#' @param xin Numeric vector of length equal to the number of independent variables, of initial values, for the parameters to be estimated.
#' @param type An integer specifying the bandwidth selection method used, see \code{\link{bw}}.
#' @param ... Arguments to be passed on to the control argument of the \code{\link{optim}} function.
#'
#' @importFrom stats model.frame model.matrix model.response optim
#' @return mmKDEshort: A list with the following components:
#' \itemize{
#' \item coefficients: A vector of estimated coefficients.
#' \item error: The value of the objective function.
#' }
#' 
#' @examples
#' x <- 1:10
#' y <- x+rnorm(10)
#' XIn <- lm(y~x)
#' mmKDEshort(y~x, xin=coef(XIn), type=-1)
#'
#' @export
mmKDEshort <- function(formula, data=list(), xin, type, ...)
{
mf <- model.frame(formula=formula, data=data)
X <- model.matrix(attr(mf, "terms"), data=mf)
y <- model.response(mf)

h_y <- bw(y, type)
MOMy <- kdeGaussMom(ncol(X), y, h_y)

mom <- optim(xin, momKDE, gamma=X, momy=MOMy, kdeGaussMom=kdeGaussMom, type=type, control=list(...))

coefficients <- mom$par

labels <- if(attr(attr(mf, "terms"), "intercept") == 1) c("Intercept", attr(attr(mf, "terms"), "term.labels")) else attr(attr(mf, "terms"), "term.labels")

names(coefficients) <- labels

error <- mom$value

momout <- list(coefficients=coefficients,
error=error)

momout
}