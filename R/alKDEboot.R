#' Arc length matching for kernel density estimators.
#'
#' Bootstrap estimates, along with standard  errors and confidence intervals, of a linear model, resulting from arc length matching of kernel density estimates.
#'
#' @param formula An LHS ~ RHS formula, specifying the linear model to be estimated.
#' @param data A data.frame which contains the variables in \code{formula}.
#' @param xin Numeric vector of length equal to the number of independent variables, of initial values, for the parameters to be estimated.
#' @param q1,q2 Numeric vectors, for the lower and upper bounds of the intervals over which arc lengths are to be computed.
#' @param type An integer specifying the bandwidth selection method used, see \code{\link{bw}}.
#' @param bootstraps An integer giving the number of bootstrap samples.
#' @param bootName The name of the .rds file to store the alKDEboot object.  May include a path.
#' @param ... Arguments to be passed on to the control argument of the \code{\link{optim}} function.
#'
#' @return A generic S3 object with class alKDEboot.
#'
#' @import pbdMPI
#' @importFrom stats coef ecdf fitted model.frame model.matrix model.response optim printCoefmat
#'
#' @export
alKDEboot <- function(formula, data=list(), xin, q1, q2, type, bootstraps, bootName, ...) UseMethod("alKDEboot")

#' @describeIn alKDEboot default method for alKDEboot.
#'
#' @return alKDEboot.default: A list object (saved using \code{saveRDS} in the specified location) with the following components:
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
#' \item ALy: Arc length segments of the KDE cast over the dependent variable.
#' \item ALX: Arc length segments of the KDE cast over the independent variables \eqn{\mathbf{X}\underline{\hat{\beta}}}.
#' p1: The vector of quantiles in the domain of \eqn{y} corresponding to \code{q1}.
#' p2: The vector of quantiles in the domain of \eqn{y} corresponding to \code{q2}.
#' \item time: Min, mean and max time incurred by the computation, as obtained from \code{\link{comm.timer}}.
#' }
#' 
#' @export
alKDEboot.default <- function(formula, data=list(), xin, q1, q2, type, bootstraps, bootName, ...)
{
comm.set.seed(123, diff=TRUE)
N <- nrow(data)

ret.time <- comm.timer({
ret <- task.pull(1:bootstraps, function(jid)
{
id <- sample(1:N, N, replace=TRUE)
alKDEshort(formula, data[id,], xin, q1, q2, type, ...)
})
})

if(comm.rank() == 0)
{
mf <- model.frame(formula=formula, data=data)
X <- model.matrix(attr(mf, "terms"), data=mf)
y <- model.response(mf)

boot <- alKDE(formula, data, xin, q1, q2, type, ...)

boot$call <- match.call()

boot$time <- ret.time
ret.coef <- do.call(rbind, lapply(1:bootstraps, function(x) ret[[x]]$coefficients))

boot$coefDist <- ret.coef

coefMean <- colMeans(ret.coef)
coefVar <- (1/(bootstraps-1))*colSums((ret.coef-coefMean)^2)
boot$se <- sqrt(coefVar)
boot$bcoefficients <- coefMean

boot$errorList <- as.vector(do.call(rbind, lapply(1:bootstraps, function(x) ret[[x]]$error)))

class(boot) <- "alKDEboot"
saveRDS(boot, paste(bootName, ".rds", sep=""))
}

finalize()
}

#' @describeIn alKDEboot print method for alKDEboot.
#'
#' @param x An alKDEboot object.
#'
#' @export
print.alKDEboot <- function(x, ...)
{
cat("Call:\n")
print(x$call)
cat("\nCoefficients:\n")
print(x$coefficients, digits=5)
}

#' @describeIn alKDEboot summary method for alKDEboot.
#'
#' @param object An alKDEboot object.
#'
#' @return summary.alKDEboot: A list of class summary.alKDEboot with the following components:
#' \itemize{
#' \item call: Original call to the \code{alKDEboot} function.
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
summary.alKDEboot <- function(object, ...)
{
ci <- do.call(rbind, lapply(1:ncol(object$coefDist), function(j) c(qsamp(object$coefDist[,j], 0.025), qsamp(object$coefDist[,j], 0.975))))

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

class(al) <- "summary.alKDEboot"
al
}

#' @describeIn alKDEboot print method for summary.alKDEboot.
#'
#' @return print.summary.alKDEboot: The object passed to the function is returned invisibly.
#'
#' @export
print.summary.alKDEboot <- function(x, ...)
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

#' @describeIn alKDEboot formula method for alKDEboot.
#' @export
alKDEboot.formula <- function(formula, data=list(), xin, q1, q2, type, bootstraps, bootName, ...)
{
mf <- model.frame(formula=formula, data=data)
x <- model.matrix(attr(mf, "terms"), data=mf)
y <- model.response(mf)

al <- alKDEboot.default(formula, data=data, xin=xin, q1=q1, q2=q2, type=type, bootstraps, bootName, ...)
al$call <- match.call()
al$formula <- formula
al$intercept <- attr(attr(mf, "terms"), "intercept")
al
}

#' @describeIn alKDEboot predict method for alKDEboot.
#'
#' @param newdata The data on which the estimated model is to be fitted.
#'
#' @return predict.alKDEboot: A vector of predicted values resulting from the estimated model.
#'
#' @export
predict.alKDEboot <- function(object, newdata=NULL, ...)
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