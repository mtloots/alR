#' Arc length matching for kernel density estimators.
#'
#' Estimate parameters of a linear model by matching the arc lengths of kernel density estimators.
#'
#' A shortened version of \code{alKDE}.
#'
#' @param formula An LHS ~ RHS formula, specifying the linear model to be estimated.
#' @param data A data.frame which contains the variables in \code{formula}.
#' @param lower,upper Numeric vectors, of length equal to the number of independent variables, for the lower and upper bounds for the parameters to be estimated.
#' @param q1,q2 Numeric vectors, of length equal to the number of independent variables, for the lower and upper bounds of the intervals over which arc lengths are to be computed.
#' @param itermax Number of iterations for the Differential Evolution algorithm.
#' @param ... Arguments to be passed on to \code{DEoptim.control()} of the Differential Evolution algorithm.
#'
#' @importFrom RcppDE DEoptim DEoptim.control
#' @importFrom stats model.frame model.matrix model.response
#' @return alKDEshort: A list with the following components:
#' \itemize{
#' \item coefficients: A vector of estimated coefficients.
#' \item error: The value of the objective function.
#' }
#' 
#' @examples
#' x <- 1:10
#' y <- x+rnorm(10)
#' alKDEshort(y~x, lower=c(-2,2), upper=c(2,2), q1=c(0.1,0.5), q2=c(0.5,0.9), itermax=50)
#'
#' @export
alKDEshort <- function(formula, data=list(), lower, upper, q1, q2, itermax, ...)
{
mf <- model.frame(formula=formula, data=data)
X <- model.matrix(attr(mf, "terms"), data=mf)
y <- model.response(mf)

p1 <- sapply(1:length(q1), function(i) qkdeGauss(q1[i], y, h_y)$result)
p2 <- sapply(1:length(q2), function(i) qkdeGauss(q2[i], y, h_y)$result)

h_y <- Silverman(y)
ALy <- kdeGaussInt2(y, h_y, p1, p2, FALSE)

al <- DEoptim(alrKDE, lower=lower, upper=upper, control=DEoptim.control(trace=FALSE, itermax=itermax, strategy=2, ...), gamma=X, aly=ALy, q1=p1, q2=p2)

coefficients <- al$optim$bestmem

labels <- if(attr(attr(mf, "terms"), "intercept") == 1) c("Intercept", attr(attr(mf, "terms"), "term.labels")) else attr(attr(mf, "terms"), "term.labels")

names(coefficients) <- labels

error <- al$optim$bestval

alout <- list(coefficients=coefficients,
error=error)

alout
}