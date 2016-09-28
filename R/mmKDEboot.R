#' Moment matching for kernel density estimators.
#'
#' Bootstrap estimates, along with standard  errors and confidence intervals, of a linear model, resulting from moment matching of kernel density estimates.
#'
#' @param formula An LHS ~ RHS formula, specifying the linear model to be estimated.
#' @param data A data.frame which contains the variables in \code{formula}.
#' @param lower,upper Numeric vectors, of length equal to the number of independent variables, for the lower and upper bounds for the parameters to be estimated.
#' @param itermax Number of iterations for the Differential Evolution algorithm.
#' @param type An integer specifying the bandwidth selection method used, see \code{\link{bw}}.
#' @param bootstraps An integer giving the number of bootstrap samples.
#' @param bootName The name of the .rds file to store the mmKDEboot object.  May include a path.
#' @param ... Arguments to be passed on to \code{DEoptim.control()} of the Differential Evolution algorithm.
#'
#' @return A generic S3 object with class mmKDEboot.
#'
#' @import pbdMPI
#' @importFrom DEoptim DEoptim DEoptim.control
#' @importFrom stats coef ecdf fitted model.frame model.matrix model.response printCoefmat quantile
#'
#' @export
mmKDEboot <- function(formula, data=list(), lower, upper, itermax, type, bootstraps, bootName, ...) UseMethod("mmKDEboot")

#' @describeIn mmKDEboot default method for mmKDEboot.
#'
#' @return mmKDEboot.default: A list object (saved using \code{saveRDS} in the specified location) with the following components:
#' \itemize{
#' \item intercept: Did the model contain an intercept TRUE/FALSE?
#' \item coefficients: A vector of estimated coefficients.
#' \item coefDist The bootstrap parameter distribution.
#' \item bcoefficients: A vector of bootstrap coefficients, resulting from bootstrap estimation.
#' \item df: Degrees of freedom of the model.
#' \item se: The standard errors for the estimates resulting from bootstrap estimation.
#' \item error: The value of the objective function.
#' \item errorList: A vector of values of the objective function for each bootstrap sample.
#' \item fitted.values: A vector of estimated values.
#' \item residuals: The residuals resulting from the fitted model.
#' \item call: The call to the function.
#' \item h_y: The KDE bandwidth estimator for the dependent variable.
#' \item h_X: The KDE bandwidth estimator for the independent variables, i.e. \eqn{\mathbf{X}\underline{\hat{\beta}}}.
#' \item MOMy: The first \eqn{n} non central moments of the dependent variable, where \eqn{n} is the number of columns in the design matrix.
#' \item MOMX: The first \eqn{n} non central moments of the independent variables \eqn{\mathbf{X}\underline{\hat{\beta}}}, where \eqn{n} is the number of columns in the design matrix.
#' \item time: Min, mean and max time incurred by the computation, as obtained from \code{\link{comm.timer}}.
#' }
#' 
#' @export
mmKDEboot.default <- function(formula, data=list(), lower, upper, itermax, type, bootstraps, bootName, ...)
{
comm.set.seed(123, diff=TRUE)
N <- nrow(data)

ret.time <- comm.timer({
ret <- task.pull(1:bootstraps, function(jid)
{
id <- sample(1:N, N, replace=TRUE)
mmKDEshort(formula, data[id,], lower, upper, itermax, type, ...)
})
})

if(comm.rank() == 0)
{
mf <- model.frame(formula=formula, data=data)
X <- model.matrix(attr(mf, "terms"), data=mf)
y <- model.response(mf)

boot <- mmKDE(formula, data, lower, upper, itermax, type, ...)

boot$call <- match.call()

boot$time <- ret.time
ret.coef <- do.call(rbind, lapply(1:bootstraps, function(x) ret[[x]]$coefficients))

boot$coefDist <- ret.coef

coefMean <- colMeans(ret.coef)
coefVar <- (1/(bootstraps-1))*colSums((ret.coef-coefMean)^2)
boot$se <- sqrt(coefVar)
boot$bcoefficients <- coefMean

boot$errorList <- as.vector(do.call(rbind, lapply(1:bootstraps, function(x) ret[[x]]$error)))

class(boot) <- "mmKDEboot"
saveRDS(boot, paste(bootName, ".rds", sep=""))
}

finalize()
}

#' @describeIn mmKDEboot print method for mmKDEboot.
#'
#' @param x An mmKDEboot object.
#'
#' @export
print.mmKDEboot <- function(x, ...)
{
cat("Call:\n")
print(x$call)
cat("\nCoefficients:\n")
print(x$coefficients, digits=5)
}

#' @describeIn mmKDEboot summary method for mmKDEboot.
#'
#' @param object An mmKDEboot object.
#'
#' @return summary.mmKDEboot: A list of class summary.mmKDEboot with the following components:
#' \itemize{
#' \item call: Original call to \code{mmKDEboot} function.
#' \item coefficients: A matrix with estimates, estimated errors, and 95\% parameter confidence intervals (based on the inverse empirical distribution function).
#' \item moments: A matrix of the first \eqn{n} moments of the dependent and independent variables that were matched.  The final row corresponds to the estimated bandwidth parameters for each, i.e. \code{h_y} and \code{h_X}, respectively.
#' \item r.squared: The \eqn{r^{2}} coefficient.
#' \item adj.r.squared: The adjusted \eqn{r^{2}} coefficient.
#' \item sigma: The residual standard error.
#' \item df: Degrees of freedom for the model.
#' \item error: Value of the objective function.
#' \item time: Min, mean and max time incurred by the computation, as obtained from \code{\link{comm.timer}}.
#' \item residSum: Summary statistics for the distribution of the residuals.
#' \item errorSum: Summary statistics for the distribution of the value of the objective function.
#' }
#'
#' @export
summary.mmKDEboot <- function(object, ...)
{
ci <- do.call(rbind, lapply(1:ncol(object$coefDist), function(j) quantile(object$coefDist[,j], probs=c(0.025, 0.975), names=FALSE, type=1)))

pval <- do.call(rbind, lapply(1:length(object$coefficients), function(i) {
fn <- ecdf(object$coefDist[,i]-object$bcoefficients[i])
1-fn(object$coefficients[i])
}))

TAB <- cbind(Estimate = coef(object),
StdErr=object$se,
ci=ci,
b.value = object$bcoefficients,
p.value = pval)

rownames(TAB) <- names(object$coefficients)
    colnames(TAB) <- c("Estimate", "StdErr", "LCI", "UCI", "b.value", "p.value")

momTAB <- cbind(LHS = c(object$MOMy, object$h_y), RHS = c(object$MOMX, object$h_X))

rownames(momTAB) <- c(1:length(object$MOMy), "BW")
    colnames(momTAB) <- c("LHS", "RHS")

f <- object$fitted.values
r <- object$residuals
mss <- if(object$intercept) sum((f-mean(f))^2) else sum(f^2)
rss <- sum(r^2)

r.squared <- mss/(mss+rss)
    df.int <- if (object$intercept) 1L else 0L
n <- length(f)
rdf <- object$df
adj.r.squared <- 1-(1-r.squared)*((n-df.int)/rdf)

mom <- list(call=object$call,
coefficients=TAB,
moments = momTAB,
r.squared=r.squared,
adj.r.squared=adj.r.squared,
sigma=sqrt(sum((object$residuals)^2)/rdf),
df=object$df,
error=object$error,
time=object$time,
residSum=summary(object$residuals, digits=5)[-4],
errorSum=summary(object$errorList, digits=5)[-4])

class(mom) <- "summary.mmKDEboot"
mom
}

#' @describeIn mmKDEboot print method for summary.mmKDEboot.
#'
#' @return print.summary.mmKDEboot: The object passed to the function is returned invisibly.
#'
#' @export
print.summary.mmKDEboot <- function(x, ...)
{
cat("Call:\n")
print(x$call)
cat("\nResiduals:\n")
print(x$residSum)

cat("\nKDE Moments:\n")
print(x$moments)
cat("\nValue of objective function:\n")
print(x$errorSum)
cat("\n")

printCoefmat(x$coefficients, P.values=TRUE, has.Pvalue=TRUE)
digits <- max(3, getOption("digits")-3)

cat("\nResidual standard error: ", formatC(x$sigma, digits=digits), " on ",
formatC(x$df, digits=0, format="f"), " degrees of freedom\n", sep="")
cat("Multiple R-squared: ", formatC(x$r.squared, digits=digits),
",\tAdjusted R-squared: ",formatC(x$adj.r.squared, digits=digits),
"\nValue of objective function: ",formatC(x$error,digits=digits, format="f"), "\n", sep="")

cat("\nElapsed time:\n", sep="")
print(x$time)
invisible(x)
}

#' @describeIn mmKDEboot formula method for mmKDEboot.
#' @export
mmKDEboot.formula <- function(formula, data=list(), lower, upper, itermax, type, bootstraps, bootName, ...)
{
mf <- model.frame(formula=formula, data=data)
x <- model.matrix(attr(mf, "terms"), data=mf)
y <- model.response(mf)

mom <- mmKDEboot.default(formula, data=data, lower=lower, upper=upper, itermax=itermax, type=type, bootstraps, bootName, ...)
mom$call <- match.call()
mom$formula <- formula
mom$intercept <- attr(attr(mf, "terms"), "intercept")
mom
}

#' @describeIn mmKDEboot predict method for mmKDEboot.
#'
#' @param newdata The data on which the estimated model is to be fitted.
#'
#' @return predict.mmKDEboot: A vector of predicted values resulting from the estimated model.
#'
#' @export
predict.mmKDEboot <- function(object, newdata=NULL, ...)
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

y <- as.vector(x%*%coef(object))
}
y
}