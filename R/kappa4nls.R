#' Sigmoidal curve fitting.
#'
#' A framework for nonlinear least squares fitting of the four-parameter kappa sigmoidal function.
#'
#' @param formula An LHS ~ RHS formula, specifying the linear model to be estimated.
#' @param data A data.frame which contains the variables in \code{formula}.
#' @param lower A vector of lower constraints for the parameters to be estimated; defaults to c(0, -5, -5).
#' @param upper A vector of upper constraints for the parameters to be estimated; defaults to c(10, 1, 1).
#' @param tol Error tolerance level; defaults to 1e-15.
#' @param maxiter The maximum number of iterations allowed; defaults to 50000.
#' @param ... Arguments to be passed on to the differential evolution function \code{\link{JDEoptim}}.
#'
#' @return A generic S3 object with class kappa4nls.
#' @importFrom stats coef fitted model.frame model.matrix model.response printCoefmat var
#' @importFrom DEoptimR JDEoptim
#'
#' @export
kappa4nls <- function(formula, data=list(), lower, upper, tol, maxiter, ...) UseMethod("kappa4nls")

#' @describeIn kappa4nls default method for kappa4nls.
#'
#' @return kappa4nls.default: A list with all components from \code{\link{JDEoptim}}, as well as:
#' \itemize{
#' \item intercept: Did the model contain an intercept TRUE/FALSE?
#' \item coefficients: A vector of estimated coefficients.
#' \item error: The value of the objective function.
#' \item fitted.values: A vector of estimated values.
#' \item residuals: The residuals resulting from the fitted model.
#' \item call: The call to the function.
#' }
#' 
#' @examples
#' k <- kappa4tc(-4, 0, 1)$par
#' x <- seq(qkappa4(0, 4, 0.4, -4, k), qkappa4(0.7, 4, 0.4, -4, k), length.out=100)
#' y <- sapply(x, function(i) pkappa4(i, 4, 0.4, -4, k))
#' kappa4nls.default(y~x, tol=1e-5)
#'
#' @export
kappa4nls.default <- function(formula, data=list(), lower=c(0, -5, -5), upper=c(10, 1, 1), tol=1e-15, maxiter=50000, ...)
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

nls <- DEoptimR::JDEoptim(lower=lower, upper=upper, fn=kappa4NLSobj, constr=kappa4NLScon, meq=2, tol=tol, maxiter=maxiter, x=x, y=y/max(y), x_min=min(x), x_max=max(x))

nls$intercept <- intercept

if(nls$par[2] <= 0)
{
mu <- min(x)-(nls$par[1]/nls$par[3])
}
else
{
mu <- min(x)-nls$par[1]*(1-(nls$par[2])^(-nls$par[3]))/nls$par[3]
}

nls$coefficients <- c(mu, nls$par)
names(nls$coefficients) <- c("mu", "sigma", "h", "k")

nls$error <- nls$value

nls$fitted.values <- sapply(x, function(i) pkappa4(i, mu, nls$par[1], nls$par[2], nls$par[3]))/pkappa4(max(x), mu, nls$par[1], nls$par[2], nls$par[3])
nls$residuals <- y/max(y)-nls$fitted.values
nls$call <- match.call()

class(nls) <- "kappa4nls"
nls
}

#' @describeIn kappa4nls print method for kappa4nls.
#'
#' @param x A kappa4nls object.
#'
#' @export
print.kappa4nls <- function(x, ...)
{
cat("\nCall:\n")
print(x$call)
cat("\nCoefficients:\n")
print(x$coefficients, digits=5)
}

#' @describeIn kappa4nls summary method for kappa4nls.
#'
#' @param object A kappa4nls object.
#'
#' @return summary.kappa4nls: A list of class summary.kappa4nls with the following components:
#' \itemize{
#' \item call: Original call to \code{kappa4nls} function.
#' \item coefficients: A vector with parameter estimates.
#' \item r.squared: The \eqn{r^{2}} coefficient.
#' \item sigma: The residual standard error.
#' \item error: Value of the objective function.
#' \item residSum: Summary statistics for the distribution of the residuals.
#' }
#'
#' @export
summary.kappa4nls <- function(object, ...)
{
TAB <- cbind(Estimate = coef(object))

rownames(TAB) <- names(object$coefficients)
    colnames(TAB) <- c("Estimate")

y <- object$residuals+object$fitted.values

r.squared <- 1-var(object$residuals)/var(y)

nls <- list(call=object$call,
coefficients=TAB,
r.squared=r.squared,
sigma=sqrt(sum((object$residuals)^2)),
error=object$error,
residSum=summary(object$residuals, digits=5)[-4])

class(nls) <- "summary.kappa4nls"
nls
}

#' @describeIn kappa4nls print method for summary.kappa4nls.
#'
#' @return print.summary.kappa4nls: The object passed to the function is returned invisibly.
#'
#' @export
print.summary.kappa4nls <- function(x, ...)
{
cat("\nCall:\n")
print(x$call)
cat("\nResiduals:\n")
print(x$residSum)

cat("\n")

printCoefmat(x$coefficients, P.values=FALSE, has.Pvalue=FALSE)
digits <- max(3, getOption("digits") - 3)
cat("\nResidual standard error: ", formatC(x$sigma, digits=digits), sep="")
cat("\nMultiple R-squared: ", formatC(x$r.squared, digits=digits), sep="")
cat("\tValue of objective function: ",formatC(x$error, digits=digits, format="f"), "\n", sep="")
invisible(x)
}

#' @describeIn kappa4nls formula method for kappa4nls.
#' @export
kappa4nls.formula <- function(formula, data=list(), lower=c(0, -5, -5), upper=c(10, 1, 1), tol=1e-15, maxiter=50000, ...)
{
mf <- model.frame(formula=formula, data=data)
x <- model.matrix(attr(mf, "terms"), data=mf)
y <- model.response(mf)

nls <- kappa4nls.default(formula, data=data, lower=lower, upper=upper, tol=tol, maxiter=maxiter, ...)
nls$call <- match.call()
nls$formula <- formula
nls$intercept <- attr(attr(mf, "terms"), "intercept")
nls
}

#' @describeIn kappa4nls predict method for kappa4nls.
#'
#' @param newdata The data on which the estimated model is to be fitted.
#'
#' @return predict.kappa4nls: A vector of predicted values resulting from the estimated model.
#'
#' @examples
#' u <- seq(qkappa4(0.1, 4, 0.4, -4, k), qkappa4(0.8, 4, 0.4, -4, k), length.out=100)
#' v <- sapply(u, function(i) pkappa4(i, 4, 0.4, -4, k))
#' nls <- kappa4nls(y~x, tol=1e-5)
#' predict(nls, newdata=data.frame(y=v, x=u))
#'
#' @export
predict.kappa4nls <- function(object, newdata=NULL, ...)
{
if(is.null(newdata))
{
y <- fitted(object)
}
else
{
if(!is.null(object$formula))
{
x <- model.matrix(object$formula, newdata)
}
else
{
x <- newdata
}

if(object$intercept)
{
X <- x[,2]
}
else
{
X <- x[,1]
}

y <- sapply(X, function(i) pkappa4(i, coef(object)[1], coef(object)[2], coef(object)[3], coef(object)[4]))/pkappa4(max(X), coef(object)[1], coef(object)[2], coef(object)[3], coef(object)[4])
}
y
}