#' Sigmoidal curve fitting.
#'
#' Bias corrected jackknife estimates, along with standard  errors and confidence intervals, of a nonlinear model, resulting from  nonlinear least squares fitting of the four-parameter kappa sigmoidal function.
#'
#' @param formula An LHS ~ RHS formula, specifying the linear model to be estimated.
#' @param data A data.frame which contains the variables in \code{formula}.
#' @param xin Numeric vector of length 3 containing initial values, for \eqn{\sigma}, \eqn{h}, and \eqn{k}.
#' @param lower A vector of lower constraints for the parameters to be estimated; defaults to c(0, -5, -5).
#' @param upper A vector of upper constraints for the parameters to be estimated; defaults to c(10, 1, 1).
#' @param tol Error tolerance level; defaults to 1e-15.
#' @param maxiter The maximum number of iterations allowed; defaults to 50000.
#' @param jackName The name of the .rds file to store the kappa4nlsJack object.  May include a path.
#' @param ... Arguments to be passed on to the differential evolution function \code{\link{JDEoptim}}.
#'
#' @return A generic S3 object with class kappa4nlsJack.
#'
#' @import pbdMPI
#' @importFrom stats coef fitted model.frame model.matrix model.response printCoefmat var
#'
#' @export
kappa4nlsJack <- function(formula, data=list(), xin, lower, upper, tol, maxiter, jackName, ...) UseMethod("kappa4nlsJack")

#' @describeIn kappa4nlsJack default method for kappa4nlsJack.
#'
#' @return kappa4nlsJack.default: A list object (saved using \code{saveRDS} in the specified location) with the following components:
#' \itemize{
#' \item intercept: Did the model contain an intercept TRUE/FALSE?
#' \item coefficients: A vector of estimated coefficients.
#' \item coefDist The jackknife parameter distribution.
#' \item jcoefficients: A vector of bias-corrected coefficients, resulting from jackknife estimation.
#' \item bias: The corrections applied in obtaining the bias-corrected estimates.
#' \item se: The standard errors for the estimates resulting from jackknife estimation.
#' \item error: The value of the objective function.
#' \item errorList: A vector of values of the objective function at jackknife points.
#' \item fitted.values: A vector of estimated values.
#' \item residuals: The residuals resulting from the fitted model.
#' \item call: The call to the function.
#' \item time: Min, mean and max time incurred by the computation, as obtained from \code{\link{comm.timer}}.
#' }
#' 
#' @export
kappa4nlsJack.default <- function(formula, data=list(), xin, lower=c(0, -5, -5), upper=c(10, 1, 1), tol=1e-15, maxiter=50000, jackName, ...)
{
ret.time <- comm.timer({
ret <- task.pull(1:nrow(data), function(jid) kappa4nlsShort(formula, data[-jid], xin))
})

if(comm.rank() == 0)
{
jack <- kappa4nls(formula, data, lower, upper, tol, maxiter, ...)

jack$call <- match.call()

jack$time <- ret.time
ret.coef <- do.call(rbind, lapply(1:nrow(data), function(x) ret[[x]]$coefficients))

jack$coefDist <- ret.coef

coefMean <- colMeans(ret.coef)
coefVar <- ((nrow(data)-1)/nrow(data))*colSums((ret.coef-coefMean)^2)
jack$se <- sqrt(coefVar)
jack$bias <- (nrow(data)-1)*(coefMean-jack$coefficients)
jack$jcoefficients <- jack$coefficients-jack$bias

jack$errorList <- as.vector(do.call(rbind, lapply(1:nrow(data), function(x) ret[[x]]$error)))

class(jack) <- "kappa4nlsJack"
saveRDS(jack, paste(jackName, ".rds", sep=""))
}

finalize()
}

#' @describeIn kappa4nlsJack print method for kappa4nlsJack.
#'
#' @param x A kappa4nlsJack object.
#'
#' @export
print.kappa4nlsJack <- function(x, ...)
{
cat("Call:\n")
print(x$call)
cat("\nCoefficients:\n")
print(x$coefficients, digits=5)
}

#' @describeIn kappa4nlsJack summary method for kappa4nlsJack.
#'
#' @param object A kappa4nlsJack object.
#'
#' @return summary.kappa4nlsJack: A list of class summary.kappa4nlsJack with the following components:
#' \itemize{
#' \item call: Original call to the \code{kappa4nlsJack} function.
#' \item coefficients: A matrix with estimates, estimated errors, and 95\% parameter confidence intervals (based on the inverse empirical distribution function).
#' \item r.squared: The \eqn{r^{2}} coefficient.
#' \item sigma: The residual standard error.
#' \item error: Value of the objective function.
#' \item time: Min, mean and max time incurred by the computation, as obtained from \code{\link{comm.timer}}.
#' \item residSum: Summary statistics for the distribution of the residuals.
#' \item errorSum: Summary statistics for the distribution of the value of the objective function.
#' }
#'
#' @export
summary.kappa4nlsJack <- function(object, ...)
{
ci <- do.call(rbind, lapply(1:ncol(object$coefDist), function(j) c(qsamp(object$coefDist[,j], 0.025), qsamp(object$coefDist[,j], 0.975))))

pval <- do.call(rbind, lapply(1:length(object$coefficients), function(i) {
lpv <- mean(object$coefDist[,i]-object$jcoefficients[i] <= object$coefficients[i])
2*min(lpv, 1-lpv)
}))

TAB <- cbind(Estimate = coef(object),
StdErr=object$se,
ci=ci,
j.value = object$jcoefficients,
p.value = pval)

rownames(TAB) <- names(object$coefficients)
    colnames(TAB) <- c("Estimate", "StdErr", "LCI", "UCI", "j.value", "p.value")

y <- object$residuals+object$fitted.values
r.squared <- 1-var(object$residuals)/var(y)

k4 <- list(call=object$call,
coefficients=TAB,
r.squared=r.squared,
sigma=sqrt(sum((object$residuals)^2)),
error=object$error,
time=object$time,
residSum=summary(object$residuals, digits=5)[-4],
errorSum=summary(object$errorList, digits=5)[-4])

class(k4) <- "summary.kappa4nlsJack"
k4
}

#' @describeIn kappa4nlsJack print method for summary.kappa4nlsJack.
#'
#' @return print.summary.kappa4nlsJack: The object passed to the function is returned invisibly.
#'
#' @export
print.summary.kappa4nlsJack <- function(x, ...)
{
cat("Call:\n")
print(x$call)
cat("\nResiduals:\n")
print(x$residSum)

cat("\nValue of objective function:\n")
print(x$errorSum)
cat("\n")

printCoefmat(x$coefficients, P.values=TRUE, has.Pvalue=TRUE)
digits <- max(3, getOption("digits")-3)

cat("\nResidual standard error: ", formatC(x$sigma, digits=digits), "\n", sep="")
cat("Multiple R-squared: ", formatC(x$r.squared, digits=digits),
"\nValue of objective function: ",formatC(x$error,digits=digits, format="f"), "\n", sep="")

cat("\nElapsed time:\n", sep="")
print(x$time)
invisible(x)
}

#' @describeIn kappa4nlsJack formula method for kappa4nlsJack.
#' @export
kappa4nlsJack.formula <- function(formula, data=list(), xin, lower, upper, tol, maxiter, jackName, ...)
{
mf <- model.frame(formula=formula, data=data)
x <- model.matrix(attr(mf, "terms"), data=mf)
y <- model.response(mf)

k4 <- kappa4nlsJack.default(formula, data=data, xin=xin, lower=c(0, -5, -5), upper=c(10, 1, 1), tol=1e-15, maxiter=50000, jackName, ...)
k4$call <- match.call()
k4$formula <- formula
k4$intercept <- attr(attr(mf, "terms"), "intercept")
k4
}

#' @describeIn kappa4nlsJack predict method for kappa4nlsJack.
#'
#' @param newdata The data on which the estimated model is to be fitted.
#'
#' @return predict.kappa4nlsJack: A vector of predicted values resulting from the estimated model.
#'
#' @export
predict.kappa4nlsJack <- function(object, newdata=NULL, ...)
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