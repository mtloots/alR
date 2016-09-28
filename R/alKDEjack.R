#' Arc length matching for kernel density estimators.
#'
#' Bias corrected jackknife estimates, along with standard  errors and confidence intervals, of a linear model, resulting from arc length matching of kernel density estimates.
#'
#' @param formula An LHS ~ RHS formula, specifying the linear model to be estimated.
#' @param data A data.frame which contains the variables in \code{formula}.
#' @param lower,upper Numeric vectors, of length equal to the number of independent variables, for the lower and upper bounds for the parameters to be estimated.
#' @param q1,q2 Numeric vectors, of length equal to the number of independent variables, for the lower and upper bounds of the intervals over which arc lengths are to be computed.
#' @param itermax Number of iterations for the Differential Evolution algorithm.
#' @param type An integer specifying the bandwidth selection method used, see \code{\link{bw}}.
#' @param jackName The name of the .rds file to store the alKDEjack object.  May include a path.
#' @param ... Arguments to be passed on to \code{DEoptim.control()} of the Differential Evolution algorithm.
#'
#' @return A generic S3 object with class alKDEjack.
#'
#' @import pbdMPI
#' @importFrom DEoptim DEoptim DEoptim.control
#' @importFrom stats coef ecdf fitted model.frame model.matrix model.response printCoefmat quantile
#'
#' @export
alKDEjack <- function(formula, data=list(), lower, upper, q1, q2, itermax, type, jackName, ...) UseMethod("alKDEjack")

#' @describeIn alKDEjack default method for alKDEjack.
#'
#' @return alKDEjack.default: A list object (saved using \code{saveRDS} in the specified location) with the following components:
#' \itemize{
#' \item intercept: Did the model contain an intercept TRUE/FALSE?
#' \item coefficients: A vector of estimated coefficients.
#' \item coefDist The jackknife parameter distribution.
#' \item jcoefficients: A vector of jackknife coefficients, resulting from jackknife estimation.
#' \item df: Degrees of freedom of the model.
#' \item se: The standard errors for the estimates resulting from jackknife estimation.
#' \item error: The value of the objective function.
#' \item errorList: A vector of values of the objective function at jackknife points.
#' \item fitted.values: A vector of estimated values.
#' \item residuals: The residuals resulting from the fitted model.
#' \item call: The call to the function.
#' \item h_y: The KDE bandwidth estimator for the dependent variable.
#' \item h_X: The KDE bandwidth estimator for the independent variables, i.e. \eqn{\mathbf{X}\underline{\hat{\beta}}}.
#' \item ALy: \eqn{n} arc length segments of the KDE cast over the dependent variable, where $\eqn{n} is the number of columns in the design matrix.
#' \item ALX: \eqn{n} arc length segments of the KDE cast over the independent variables \eqn{\mathbf{X}\underline{\hat{\beta}}}, where $\eqn{n} is the number of columns in the design matrix.
#' \item time: Min, mean and max time incurred by the computation, as obtained from \code{\link{comm.timer}}.
#' p1: The vector of quantiles in the domain of \eqn{y} corresponding to \code{q1}.
#' p2: The vector of quantiles in the domain of \eqn{y} corresponding to \code{q2}.
#' }
#' 
#' @export
alKDEjack.default <- function(formula, data=list(), lower, upper, q1, q2, itermax, type, jackName, ...)
{
comm.set.seed(123, diff=TRUE)
N <- nrow(data)

ret.time <- comm.timer({
ret <- task.pull(1:nrow(data), function(jid) alKDEshort(formula, data[-jid], lower, upper, q1, q2, itermax, type, ...))
})

if(comm.rank() == 0)
{
mf <- model.frame(formula=formula, data=data)
X <- model.matrix(attr(mf, "terms"), data=mf)
y <- model.response(mf)

jack <- alKDE(formula, data, lower, upper, q1, q2, itermax, type, ...)

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

class(jack) <- "alKDEjack"
saveRDS(jack, paste(jackName, ".rds", sep=""))
}

finalize()
}

#' @describeIn alKDEjack print method for alKDEjack.
#'
#' @param x An alKDEjack object.
#'
#' @export
print.alKDEjack <- function(x, ...)
{
cat("Call:\n")
print(x$call)
cat("\nCoefficients:\n")
print(x$coefficients, digits=5)
}

#' @describeIn alKDEjack summary method for alKDEjack.
#'
#' @param object An alKDEjack object.
#'
#' @return summary.alKDEjack: A list of class summary.alKDEjack with the following components:
#' \itemize{
#' \item call: Original call to the \code{alKDEjack} function.
#' \item coefficients: A matrix with estimates, estimated errors, and 95\% parameter confidence intervals (based on the inverse empirical distribution function).
#' \item arclengths: A matrix of the arc length segments that were matched, for the dependent and independent variables.  The final row corresponds to the estimated bandwidth parameters for each, i.e. \code{h_y} and \code{h_X}, respectively.
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
summary.alKDEjack <- function(object, ...)
{
ci <- do.call(rbind, lapply(1:ncol(object$coefDist), function(j) quantile(object$coefDist[,j], probs=c(0.025, 0.975), names=FALSE, type=1)))

pval <- do.call(rbind, lapply(1:length(object$coefficients), function(i) {
fn <- ecdf(object$coefDist[,i]-object$jcoefficients[i])
1-fn(object$coefficients[i])
}))

TAB <- cbind(Estimate = coef(object),
StdErr=object$se,
ci=ci,
j.value = object$jcoefficients,
p.value = pval)

rownames(TAB) <- names(object$coefficients)
    colnames(TAB) <- c("Estimate", "StdErr", "LCI", "UCI", "j.value", "p.value")

alTAB <- cbind(LHS = c(object$ALy, object$h_y), RHS = c(object$ALX, object$h_X))

rownames(alTAB) <- c(paste("[", round(object$p1, 5), ", ", round(object$p2, 5), "]", sep=""), "BW")
    colnames(alTAB) <- c("LHS", "RHS")

f <- object$fitted.values
r <- object$residuals
mss <- if(object$intercept) sum((f-mean(f))^2) else sum(f^2)
rss <- sum(r^2)

r.squared <- mss/(mss+rss)
    df.int <- if (object$intercept) 1L else 0L
n <- length(f)
rdf <- object$df
adj.r.squared <- 1-(1-r.squared)*((n-df.int)/rdf)

al <- list(call=object$call,
coefficients=TAB,
arclengths = alTAB,
r.squared=r.squared,
adj.r.squared=adj.r.squared,
sigma=sqrt(sum((object$residuals)^2)/rdf),
df=object$df,
error=object$error,
time=object$time,
residSum=summary(object$residuals, digits=5)[-4],
errorSum=summary(object$errorList, digits=5)[-4])

class(al) <- "summary.alKDEjack"
al
}

#' @describeIn alKDEjack print method for summary.alKDEjack.
#'
#' @return print.summary.alKDEjack: The object passed to the function is returned invisibly.
#'
#' @export
print.summary.alKDEjack <- function(x, ...)
{
cat("Call:\n")
print(x$call)
cat("\nResiduals:\n")
print(x$residSum)

cat("\nKDE Arc Lengths:\n")
print(x$arclengths)
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

#' @describeIn alKDEjack formula method for alKDEjack.
#' @export
alKDEjack.formula <- function(formula, data=list(), lower, upper, q1, q2, itermax, type, jackName, ...)
{
mf <- model.frame(formula=formula, data=data)
x <- model.matrix(attr(mf, "terms"), data=mf)
y <- model.response(mf)

al <- alKDEjack.default(formula, data=data, lower=lower, upper=upper, q1=q1, q2=q2, itermax=itermax, type=type, jackName, ...)
al$call <- match.call()
al$formula <- formula
al$intercept <- attr(attr(mf, "terms"), "intercept")
al
}

#' @describeIn alKDEjack predict method for alKDEjack.
#'
#' @param newdata The data on which the estimated model is to be fitted.
#'
#' @return predict.alKDEjack: A vector of predicted values resulting from the estimated model.
#'
#' @export
predict.alKDEjack <- function(object, newdata=NULL, ...)
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