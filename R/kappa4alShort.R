#' Sigmoidal curve fitting.
#'
#' A framework for arc length fitting of the four-parameter kappa sigmoidal function.
#'
#' A shortened version of \code{\link{kappa4al}}.
#'
#' @param formula An LHS ~ RHS formula, specifying the linear model to be estimated.
#' @param data A data.frame which contains the variables in \code{formula}.
#' @param xin Numeric vector of length 3 containing initial values, for \eqn{\sigma}, \eqn{h}, and \eqn{k}.
#' @param q1,q2 Numeric vectors, for the lower and upper bounds of the intervals over which arc lengths are to be computed.
#' @param ... Arguments to be passed on to the outer control list of \code{\link{constrOptim.nl}}.
#'
#' @importFrom stats coef fitted model.frame model.matrix model.response printCoefmat
#' @importFrom alabama constrOptim.nl
#'
#' @return kappa4alShort: A list with the following components:
#' \itemize{
#' \item coefficients: A vector of estimated coefficients.
#' \item error: The value of the objective function.
#' }
#' 
#' @examples
#' k <- kappa4tc(-4, 0, 1)$par
#' x <- seq(qkappa4(0, 4, 0.4, -4, k), qkappa4(0.7, 4, 0.4, -4, k), length.out=100)
#' y <- sapply(x, function(i) pkappa4(i, 4, 0.4, -4, k))
#' kappa4alShort(y~x, xin=c(0.1, -3, -0.1), q1=c(0.1, 0.5), q2=c(0.5, 0.9))
#'
#' @export
kappa4alShort <- function(formula, data=list(), xin, q1, q2, ...)
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

p1 <- sapply(1:length(q1), function(i) x[length(which((y/max(y)) <= q1[i]))+1])
p2 <- sapply(1:length(q2), function(i) x[length(which((y/max(y)) <= q2[i]))+1])

al_samp <- kappa4IntApprox2(x, y/max(y), p1, p2, FALSE)

al <- alabama::constrOptim.nl(par=xin, fn=kappa4ALobj, hin=kappa4ALhin, heq=kappa4ALheq, control.outer=list(eps=1e-25, itmax=50000, trace=FALSE), control.optim=list(maxit=50000, abstol=1e-25, reltol=1e-25), al_samp=al_samp, x_min=min(x), x_max=max(x), q1=p1, q2=p2)

k <- al$par[3]

if(al$par[2] <= 0)
{
mu <- min(x)-(al$par[1]/k)
}
else
{
mu <- min(x)-al$par[1]*(1-(al$par[2])^(-k))/k
}

coefficients <- c(mu, al$par)
names(coefficients) <- c("mu", "sigma", "h", "k")

error <- al$value

alout <- list(coefficients=coefficients,
error=error)
alout
}