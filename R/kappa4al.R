#' Sigmoidal curve fitting.
#'
#' A framework for arc length fitting of the four-parameter kappa sigmoidal function.
#'
#' @param formula An LHS ~ RHS formula, specifying the linear model to be estimated.
#' @param data A data.frame which contains the variables in \code{formula}.
#' @param lower A vector of lower constraints for the parameters to be estimated; defaults to c(0, -5, -5).
#' @param upper A vector of upper constraints for the parameters to be estimated; defaults to c(10, 1, 1).
#' @param q1,q2 Numeric vectors, for the lower and upper bounds of the intervals over which arc lengths are to be computed.
#' @param tol Error tolerance level; defaults to 1e-15.
#' @param maxiter The maximum number of iterations allowed; defaults to 50000.
#' @param ... Arguments to be passed on to the differential evolution function \code{\link{JDEoptim}}.
#'
#' @return A generic S3 object with class kappa4al.
#' @importFrom stats coef fitted model.frame model.matrix model.response printCoefmat var
#' @importFrom DEoptimR JDEoptim
#'
#' @export
kappa4al <- function(formula, data=list(), lower, upper, q1, q2, tol, maxiter, ...) UseMethod("kappa4al")

#' @describeIn kappa4al default method for kappa4al.
#'
#' @return kappa4al.default: A list with all components from \code{\link{JDEoptim}}, as well as:
#' \itemize{
#' \item intercept: Did the model contain an intercept TRUE/FALSE?
#' \item coefficients: A vector of estimated coefficients.
#' \item error: The value of the objective function.
#' \item fitted.values: A vector of estimated values.
#' \item residuals: The residuals resulting from the fitted model.
#' \item call: The call to the function.
#' \item ALFHat: Arc length segments of the empirical CDF (calculated from data).
#' \item ALF: Arc length segments of the CDF of the four-parameter kappa distribution (theoretical).
#' p1: The vector of sample quantiles in the data corresponding to \code{q1}.
#' p2: The vector of sample quantiles in the data corresponding to \code{q2}.
#' }
#' 
#' @examples
#' k <- kappa4tc(-4, 0, 1)$par
#' x <- seq(qkappa4(0, 4, 0.4, -4, k), qkappa4(0.7, 4, 0.4, -4, k), length.out=100)
#' y <- sapply(x, function(i) pkappa4(i, 4, 0.4, -4, k))
#' kappa4nls.default(y~x, q1=c(0.025, 0.5), q2=c(0.5, 0.975), tol=1e-5)
#'
#' @export
kappa4al.default <- function(formula, data=list(), lower=c(0, -5, -5), upper=c(10, 1, 1), q1, q2, tol=1e-15, maxiter=50000, ...)
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

al <- DEoptimR::JDEoptim(lower=lower, upper=upper, fn=kappa4ALobj, constr=kappa4ALcon, meq=2, tol=tol, maxiter=maxiter, al_samp=al_samp, x_min=min(x), x_max=max(x), q1=p1, q2=p2)

al$intercept <- intercept

k <- al$par[3]

if(al$par[2] <= 0)
{
mu <- min(x)-(al$par[1]/k)
}
else
{
mu <- min(x)-al$par[1]*(1-(al$par[2])^(-k))/k
}

al$coefficients <- c(mu, al$par)
names(al$coefficients) <- c("mu", "sigma", "h", "k")

al$error <- al$value

al$fitted.values <- sapply(x, function(i) pkappa4(i, mu, al$par[1], al$par[2], k))/pkappa4(max(x), mu, al$par[1], al$par[2], k)
al$residuals <- (y/max(y))-al$fitted.values
al$call <- match.call()

al$ALFHat <- al_samp
al$ALF <- kappa4Int2(mu, al$par[1], al$par[2], k, pkappa4(max(x), mu, al$par[1], al$par[2], k), p1, p2, FALSE)
al$p1 <- p1
al$p2 <- p2

class(al) <- "kappa4al"
al
}

#' @describeIn kappa4al print method for kappa4al.
#'
#' @param x A kappa4al object.
#'
#' @export
print.kappa4al <- function(x, ...)
{
cat("\nCall:\n")
print(x$call)
cat("\nCoefficients:\n")
print(x$coefficients, digits=5)
}

#' @describeIn kappa4al summary method for kappa4al.
#'
#' @param object A kappa4al object.
#'
#' @return summary.kappa4al: A list of class summary.kappa4al with the following components:
#' \itemize{
#' \item call: Original call to the \code{kappa4al} function.
#' \item coefficients: A vector with parameter estimates.
#' \item arclengths: A matrix of the arc length segments of the dependent and independent variables that were matched.
#' \item r.squared: The \eqn{r^{2}} coefficient.
#' \item sigma: The residual standard error.
#' \item error: Value of the objective function.
#' \item residSum: Summary statistics for the distribution of the residuals.
#' }
#'
#' @export
summary.kappa4al <- function(object, ...)
{
TAB <- cbind(Estimate = coef(object))

rownames(TAB) <- names(object$coefficients)
    colnames(TAB) <- c("Estimate")

alTAB <- cbind(LHS = object$ALF, RHS = object$ALFHat)

rownames(alTAB) <- paste("[", round(object$p1, 5), ", ", round(object$p2, 5), "]", sep="")
    colnames(alTAB) <- c("Theoretical", "Sample")

y <- object$residuals+object$fitted.values

r.squared <- 1-var(object$residuals)/var(y)

al <- list(call=object$call,
coefficients=TAB,
arclengths = alTAB,
r.squared=r.squared,
sigma=sqrt(sum((object$residuals)^2)),
error=object$error,
residSum=summary(object$residuals, digits=5)[-4])

class(al) <- "summary.kappa4al"
al
}

#' @describeIn kappa4al print method for summary.kappa4al.
#'
#' @return print.summary.kappa4al: The object passed to the function is returned invisibly.
#'
#' @export
print.summary.kappa4al <- function(x, ...)
{
cat("\nCall:\n")
print(x$call)
cat("\nResiduals:\n")
print(x$residSum)

cat("\nKappa4 CDF Arc Lengths:\n")
print(x$arclengths)
cat("\n")

printCoefmat(x$coefficients, P.values=FALSE, has.Pvalue=FALSE)
digits <- max(3, getOption("digits") - 3)
cat("\nResidual standard error: ", formatC(x$sigma, digits=digits), sep="")
cat("\nMultiple R-squared: ", formatC(x$r.squared, digits=digits), sep="")
cat("\tValue of objective function: ",formatC(x$error, digits=digits, format="f"), "\n", sep="")

invisible(x)
}

#' @describeIn kappa4al formula method for kappa4al.
#' @export
kappa4al.formula <- function(formula, data=list(), lower=c(0, -5, -5), upper=c(10, 1, 1), q1, q2, tol=1e-15, maxiter=50000, ...)
{
mf <- model.frame(formula=formula, data=data)
x <- model.matrix(attr(mf, "terms"), data=mf)
y <- model.response(mf)

al <- kappa4al.default(formula, data=data, lower=lower, upper=upper, q1=q1, q2=q2, tol=tol, maxiter=maxiter, ...)
al$call <- match.call()
al$formula <- formula
al$intercept <- attr(attr(mf, "terms"), "intercept")
al
}

#' @describeIn kappa4al predict method for kappa4al.
#'
#' @param newdata The data on which the estimated model is to be fitted.
#'
#' @return predict.kappa4al: A vector of predicted values resulting from the estimated model.
#'
#' @examples
#' u <- seq(qkappa4(0.1, 4, 0.4, -4, k), qkappa4(0.8, 4, 0.4, -4, k), length.out=100)
#' v <- sapply(u, function(i) pkappa4(i, 4, 0.4, -4, k))
#' al <- kappa4al(y~x, q1=c(0.025, 0.5), q2=c(0.5, 0.975), tol=1e-5)
#' predict(al, newdata=data.frame(y=v, x=u))
#'
#' @export
predict.kappa4al <- function(object, newdata=NULL, ...)
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