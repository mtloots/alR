#' Arc length matching for kernel density estimators.
#'
#' Estimate parameters of a linear model by matching the arc lengths of kernel density estimators.
#'
#' A shortened version of \code{\link{alKDE}}.
#'
#' @param formula An LHS ~ RHS formula, specifying the linear model to be estimated.
#' @param data A data.frame which contains the variables in \code{formula}.
#' @param xin Numeric vector of length equal to the number of independent variables, of initial values, for the parameters to be estimated.
#' @param q1,q2 Numeric vectors, for the lower and upper bounds of the intervals over which arc lengths are to be computed.
#' @param type An integer specifying the bandwidth selection method used, see \code{\link{bw}}.
#' @param ... Arguments to be passed on to the control argument of the \code{\link{optim}} function.
#'
#' @importFrom stats model.frame model.matrix model.response optim
#' @return alKDEshort: A list with the following components:
#' \itemize{
#' \item coefficients: A vector of estimated coefficients.
#' \item error: The value of the objective function.
#' }
#' 
#' @examples
#' x <- 1:10
#' y <- x+rnorm(10)
#' XIn <- lm(y~x)
#' alKDEshort(y~x, xin=coef(XIn), q1=c(0.1,0.5), q2=c(0.5,0.9), type=-1)
#'
#' @export
alKDEshort <- function(formula, data=list(), xin, q1, q2, type, ...)
{
mf <- model.frame(formula=formula, data=data)
X <- model.matrix(attr(mf, "terms"), data=mf)
y <- model.response(mf)

h_y <- bw(y, type)
p1 <- sapply(1:length(q1), function(i) qkdeGauss(q1[i], y, h_y)$result)
p2 <- sapply(1:length(q2), function(i) qkdeGauss(q2[i], y, h_y)$result)

ALy <- kdeGaussInt2(y, h_y, p1, p2, FALSE)

al <- optim(xin, alrKDE, gamma=X, aly=ALy, q1=p1, q2=p2, type=type, control=list(...))

coefficients <- al$par

labels <- if(attr(attr(mf, "terms"), "intercept") == 1) c("Intercept", attr(attr(mf, "terms"), "term.labels")) else attr(attr(mf, "terms"), "term.labels")

names(coefficients) <- labels

error <- al$value

alout <- list(coefficients=coefficients,
error=error)

alout
}