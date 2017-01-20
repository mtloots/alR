#' Sigmoidal curve fitting.
#'
#' A framework for nonlinear least squares fitting of the four-parameter kappa sigmoidal function.
#'
#' A shortened version of \code{\link{kappa4nls}}.
#'
#' @param formula An LHS ~ RHS formula, specifying the linear model to be estimated.
#' @param data A data.frame which contains the variables in \code{formula}.
#' @param xin Numeric vector of length 3 containing initial values, for \eqn{\sigma}, \eqn{h}, and \eqn{k}.
#' @param ... Arguments to be passed on to the outer control list of \code{\link{constrOptim.nl}}.
#'
#' @importFrom stats coef fitted model.frame model.matrix model.response printCoefmat
#' @importFrom alabama constrOptim.nl
#'
#' @return kappa4nlsShort: A list with the following components:
#' \itemize{
#' \item coefficients: A vector of estimated coefficients.
#' \item error: The value of the objective function.
#' }
#' 
#' @examples
#' k <- kappa4tc(-4)$par
#' x <- seq(qkappa4(0, 4, 0.4, -4, k), qkappa4(0.7, 4, 0.4, -4, k), length.out=100)
#' y <- sapply(x, function(i) pkappa4(i, 4, 0.4, -4, k))
#' kappa4nlsShort(y~x, xin=c(0.1, -3, -0.1))
#'
#' @export
kappa4nlsShort <- function(formula, data=list(), xin, ...)
{
mf <- model.frame(formula=formula, data=data)
X <- model.matrix(attr(mf, "terms"), data=mf)
y <- model.response(mf)

intercept <- if(attr(attr(mf, "terms"), "intercept") == 1) TRUE else FALSE

if(intercept)
{
x <- X[,2]
}
else
{
x <- X[,1]
}

nls <- alabama::constrOptim.nl(par=xin, fn=kappa4NLSobj, hin=kappa4NLShin, heq=kappa4NLSheq, control.outer=list(eps=1e-25, itmax=50000, trace=FALSE), control.optim=list(maxit=50000, abstol=1e-25, reltol=1e-25), xvec=x, y=y/max(y), x_min=min(x), x_max=max(x))

if(nls$par[2] <= 0)
{
mu <- min(x)-(nls$par[1]/nls$par[3])
}
else
{
mu <- min(x)-nls$par[1]*(1-(nls$par[2])^(-nls$par[3]))/nls$par[3]
}

coefficients <- c(mu, nls$par)
names(coefficients) <- c("mu", "sigma", "h", "k")

error <- nls$value

nlsout <- list(coefficients=coefficients,
error=error)
nlsout
}