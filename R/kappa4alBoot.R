#' Sigmoidal curve fitting.
#'
#' Bootstrap estimates, along with standard  errors and confidence intervals, of a nonlinear model, resulting from  arc length fitting of the four-parameter kappa sigmoidal function.
#'
#' On systems where the pbMPI package is available, this code will run in parallel.
#'
#' @param formula An LHS ~ RHS formula, specifying the linear model to be estimated.
#' @param data A data.frame which contains the variables in \code{formula}.
#' @param xin Numeric vector of length 3 containing initial values, for \eqn{\sigma}, \eqn{h}, and \eqn{k}.
#' @param lower A vector of lower constraints for the parameters to be estimated; defaults to c(0, -5, -5).
#' @param upper A vector of upper constraints for the parameters to be estimated; defaults to c(10, 1, 1).
#' @param q1,q2 Numeric vectors, for the lower and upper bounds of the intervals over which arc lengths are to be computed.
#' @param tol Error tolerance level; defaults to 1e-15.
#' @param maxiter The maximum number of iterations allowed; defaults to 50000.
#' @param bootstraps An integer giving the number of bootstrap samples.
#' @param bootName The name of the .rds file to store the kappa4alBoot object.  May include a path.
#' @param ... Arguments to be passed on to the differential evolution function \code{\link{JDEoptim}}.
#'
#' @return A generic S3 object with class kappa4alBoot.
#'
#' @importFrom stats coef fitted model.frame model.matrix model.response printCoefmat sd var
#'
#' @export
kappa4alBoot <- function(formula, data=list(), xin, lower, upper, q1, q2, tol, maxiter, bootstraps, bootName, ...) UseMethod("kappa4alBoot")

#' @describeIn kappa4alBoot default method for kappa4alBoot.
#'
#' @return kappa4alBoot.default: A list object (saved using \code{saveRDS} in the specified location) with the following components:
#' \itemize{
#' \item intercept: Did the model contain an intercept TRUE/FALSE?
#' \item coefficients: A vector of estimated coefficients.
#' \item bcoefficients: A vector of bootstrap coefficients, resulting from bootstrap estimation.
#' \item se: The standard errors for the estimates resulting from bootstrap estimation.
#' \item error: The value of the objective function.
#' \item errorList: A vector of values of the objective function for each bootstrap sample.
#' \item fitted.values: A vector of estimated values.
#' \item residuals: The residuals resulting from the fitted model.
#' \item call: The call to the function.
#' \item time: Min, mean and max time incurred by the computation, as obtained from \code{comm.timer}, or that obtained from \code{\link{system.time}}.
#' }
#' 
#' @export
kappa4alBoot.default <- function(formula, data=list(), xin, lower=c(0, -5, -5), upper=c(10, 1, 1), q1, q2, tol=1e-15, maxiter=50000, bootstraps, bootName, ...)
{
if(requireNamespace("pbdMPI", quietly=TRUE))
{
pbdMPI::comm.set.seed(123, diff=TRUE)
N <- nrow(data)

ret.time <- pbdMPI::comm.timer({
ret <- pbdMPI::task.pull(1:bootstraps, function(jid)
{
id <- sample(1:N, N, replace=TRUE)
id <- id[order(id)]
kappa4alShort(formula, data[id,], xin, q1, q2)
})
})

if(pbdMPI::comm.rank() == 0)
{
boot <- kappa4al(formula, data, lower, upper, q1, q2, tol, maxiter, ...)

boot$call <- match.call()

boot$time <- ret.time
ret.coef <- do.call(rbind, lapply(1:bootstraps, function(x) ret[[x]]$coefficients))

boot$coefDist <- ret.coef

boot$bcoefficients <- colMeans(ret.coef)
boot$se <- apply(ret.coef, 2, sd)

boot$errorList <- as.vector(do.call(rbind, lapply(1:bootstraps, function(x) ret[[x]]$error)))

class(boot) <- "kappa4alBoot"
saveRDS(boot, paste(bootName, ".rds", sep=""))
}

pbdMPI::finalize()
}
else
{
set.seed(123)
N <- nrow(data)

ret.time <- system.time({
ret <- lapply(1:bootstraps, function(jid)
{
id <- sample(1:N, N, replace=TRUE)
id <- id[order(id)]
kappa4alShort(formula, data[id,], xin, q1, q2)
})
})

boot <- kappa4al(formula, data, lower, upper, q1, q2, tol, maxiter, ...)

boot$call <- match.call()

boot$time <- ret.time
ret.coef <- do.call(rbind, lapply(1:bootstraps, function(x) ret[[x]]$coefficients))

boot$coefDist <- ret.coef

boot$bcoefficients <- colMeans(ret.coef)
boot$se <- apply(ret.coef, 2, sd)

boot$errorList <- as.vector(do.call(rbind, lapply(1:bootstraps, function(x) ret[[x]]$error)))

class(boot) <- "kappa4alBoot"
saveRDS(boot, paste(bootName, ".rds", sep=""))
}
}

#' @describeIn kappa4alBoot print method for kappa4alBoot.
#'
#' @param x A kappa4alBoot object.
#'
#' @export
print.kappa4alBoot <- function(x, ...)
{
cat("Call:\n")
print(x$call)
cat("\nCoefficients:\n")
print(x$coefficients, digits=5)
}

#' @describeIn kappa4alBoot summary method for kappa4alBoot.
#'
#' @param object A kappa4alBoot object.
#'
#' @return summary.kappa4alBoot: A list of class summary.kappa4alBoot with the following components:
#' \itemize{
#' \item call: Original call to the \code{kappa4alBoot} function.
#' \item coefficients: A matrix with estimates, estimated errors, and 95\% parameter confidence intervals (based on the inverse empirical distribution function).
#' \item arclengths: A matrix of the arc length segments that were matched, for the dependent and independent variables.
#' \item r.squared: The \eqn{r^{2}} coefficient.
#' \item sigma: The residual standard error.
#' \item error: Value of the objective function.
#' \item time: Min, mean and max time incurred by the computation, as obtained from \code{comm.timer}, or that obtained from \code{\link{system.time}}.
#' \item residSum: Summary statistics for the distribution of the residuals.
#' \item errorSum: Summary statistics for the distribution of the value of the objective function.
#' }
#'
#' @export
summary.kappa4alBoot <- function(object, ...)
{
ci <- do.call(rbind, lapply(1:ncol(object$coefDist), function(j) c(qsamp(object$coefDist[,j], 0.025), qsamp(object$coefDist[,j], 0.975))))

pval <- do.call(rbind, lapply(1:length(object$coefficients), function(i) {
lpv <- mean(object$coefDist[,i]-object$bcoefficients[i] <= object$coefficients[i])
2*min(lpv, 1-lpv)
}))

TAB <- cbind(Estimate = coef(object),
StdErr=object$se,
ci=ci,
b.value = object$bcoefficients,
p.value = pval)

rownames(TAB) <- names(object$coefficients)
    colnames(TAB) <- c("Estimate", "StdErr", "LCI", "UCI", "b.value", "p.value")

alTAB <- cbind(LHS = object$ALF, RHS = object$ALFHat)

rownames(alTAB) <- paste("[", round(object$p1, 5), ", ", round(object$p2, 5), "]", sep="")
colnames(alTAB) <- c("Theoretical", "Sample")

y <- object$residuals+object$fitted.values
r.squared <- 1-var(object$residuals)/var(y)

k4 <- list(call=object$call,
coefficients=TAB,
arclengths = alTAB,
r.squared=r.squared,
sigma=sqrt(sum((object$residuals)^2)),
error=object$error,
time=object$time,
residSum=summary(object$residuals, digits=5)[-4],
errorSum=summary(object$errorList, digits=5)[-4])

class(k4) <- "summary.kappa4alBoot"
k4
}

#' @describeIn kappa4alBoot print method for summary.kappa4alBoot.
#'
#' @return print.summary.kappa4alBoot: The object passed to the function is returned invisibly.
#'
#' @export
print.summary.kappa4alBoot <- function(x, ...)
{
cat("Call:\n")
print(x$call)
cat("\nResiduals:\n")
print(x$residSum)

cat("\nKappa4 CDF Arc Lengths:\n")
print(x$arclengths)

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

#' @describeIn kappa4alBoot formula method for kappa4alBoot.
#' @export
kappa4alBoot.formula <- function(formula, data=list(), xin, lower, upper, q1, q2, tol, maxiter, bootstraps, bootName, ...)
{
mf <- model.frame(formula=formula, data=data)
x <- model.matrix(attr(mf, "terms"), data=mf)
y <- model.response(mf)

k4 <- kappa4alBoot.default(formula, data=data, xin=xin, lower=c(0, -5, -5), upper=c(10, 1, 1), q1=q1, q2=q2, tol=1e-15, maxiter=50000, bootstraps, bootName, ...)
k4$call <- match.call()
k4$formula <- formula
k4$intercept <- attr(attr(mf, "terms"), "intercept")
k4
}

#' @describeIn kappa4alBoot predict method for kappa4alBoot.
#'
#' @param newdata The data on which the estimated model is to be fitted.
#'
#' @return predict.kappa4alBoot: A vector of predicted values resulting from the estimated model.
#'
#' @export
predict.kappa4alBoot <- function(object, newdata=NULL, ...)
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