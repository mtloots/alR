#' Arc length matching for kernel density estimators.
#'
#' Estimate parameters of a linear model by matching the arc lengths of kernel density estimators.
#'
#' @param formula An LHS ~ RHS formula, specifying the linear model to be estimated.
#' @param data A data.frame which contains the variables in \code{formula}.
#' @param lower,upper Numeric vectors, of length equal to the number of independent variables, for the lower and upper bounds for the parameters to be estimated.
#' @param q1,q2 Numeric vectors, of length equal to the number of independent variables, for the lower and upper bounds of the intervals over which arc lengths are to be computed.
#' @param itermax Number of iterations for the Differential Evolution algorithm.
#' @param type An integer specifying the bandwidth selection method, see \code{\link{bw}}.
#' @param ... Arguments to be passed on to \code{DEoptim.control()} of the Differential Evolution algorithm.
#'
#' @return A generic S3 object with class alKDE.
#' @importFrom DEoptim DEoptim DEoptim.control
#' @importFrom stats coef fitted model.frame model.matrix model.response printCoefmat
#'
#' @export
alKDE <- function(formula, data=list(), lower, upper, q1, q2, itermax, type, ...) UseMethod("alKDE")

#' @describeIn alKDE default method for alKDE.
#'
#' @return alKDE.default: A list with all components from \code{\link[DEoptim]{DEoptim}}, as well as:
#' \itemize{
#' \item intercept: Did the model contain an intercept TRUE/FALSE?
#' \item coefficients: A vector of estimated coefficients.
#' \item df: Degrees of freedom of the model.
#' \item error: The value of the objective function.
#' \item fitted.values: A vector of estimated values.
#' \item residuals: The residuals resulting from the fitted model.
#' \item call: The call to the function.
#' \item h_y: The KDE bandwidth estimator for the dependent variable.
#' \item h_X: The KDE bandwidth estimator for the independent variables, i.e. \eqn{\mathbf{X}\underline{\hat{\beta}}}.
#' \item ALy: \eqn{n} arc length segments of the KDE cast over the dependent variable, where $\eqn{n} is the number of columns in the design matrix.
#' \item ALX: \eqn{n} arc length segments of the KDE cast over the independent variables \eqn{\mathbf{X}\underline{\hat{\beta}}}, where $\eqn{n} is the number of columns in the design matrix.
#' p1: The vector of quantiles in the domain of \eqn{y} corresponding to \code{q1}.
#' p2: The vector of quantiles in the domain of \eqn{y} corresponding to \code{q2}.
#' }
#' 
#' @examples
#' x <- 1:10
#' y <- x+rnorm(10)
#' alKDE(y~x, lower=c(-2,2), upper=c(2,2), q1=c(0.025,0.5), q2=c(0.5,0.975), itermax=50, type=1)
#'
#' @export
alKDE.default <- function(formula, data=list(), lower, upper, q1, q2, itermax, type, ...)
{
mf <- model.frame(formula=formula, data=data)
X <- model.matrix(attr(mf, "terms"), data=mf)
y <- model.response(mf)

h_y <- bw(y, type)

p1 <- sapply(1:length(q1), function(i) qkdeGauss(q1[i], y, h_y)$result)
p2 <- sapply(1:length(q2), function(i) qkdeGauss(q2[i], y, h_y)$result)

ALy <- kdeGaussInt2(y, h_y, p1, p2, FALSE)

al <- DEoptim(alrKDE, lower=lower, upper=upper, control=DEoptim.control(trace=FALSE, itermax=itermax, strategy=6, ...), gamma=X, aly=ALy, q1=p1, q2=p2, type=type)

al$intercept <- if(attr(attr(mf, "terms"), "intercept") == 1) TRUE else FALSE

al$coefficients <- al$optim$bestmem

labels <- if(al$intercept) c("Intercept", attr(attr(mf, "terms"), "term.labels")) else attr(attr(mf, "terms"), "term.labels")

names(al$coefficients) <- labels

al$df <- nrow(X)-ncol(X)

al$error <- al$optim$bestval

al$fitted.values <- as.vector(X%*%al$coefficients)
al$residuals <- y-al$fitted.values
al$call <- match.call()

al$h_y <- h_y
al$ALy <- ALy
al$h_X <- bw(al$fitted.values, type)
al$ALX <- kdeGaussInt2(al$fitted.values, al$h_X, p1, p2, FALSE)
al$p1 <- p1
al$p2 <- p2

class(al) <- "alKDE"
al
}

#' @describeIn alKDE print method for alKDE.
#'
#' @param x An alKDE object.
#'
#' @export
print.alKDE <- function(x, ...)
{
cat("\nCall:\n")
print(x$call)
cat("\nCoefficients:\n")
print(x$coefficients, digits=5)
}

#' @describeIn alKDE summary method for alKDE.
#'
#' @param object An alKDE object.
#'
#' @return summary.alKDE: A list of class summary.alKDE with the following components:
#' \itemize{
#' \item call: Original call to the \code{alKDE} function.
#' \item coefficients: A vector with parameter estimates.
#' \item arclengths: A matrix of the \eqn{n} arc length segments of the dependent and independent variables that were matched.  The final row corresponds to the estimated bandwidth parameters for each, i.e. \code{h_y} and \code{h_X}, respectively.
#' \item r.squared: The \eqn{r^{2}} coefficient.
#' \item adj.r.squared: The adjusted \eqn{r^{2}} coefficient.
#' \item sigma: The residual standard error.
#' \item df: Degrees of freedom for the model.
#' \item error: Value of the objective function.
#' \item residSum: Summary statistics for the distribution of the residuals.
#' }
#'
#' @export
summary.alKDE <- function(object, ...)
{
TAB <- cbind(Estimate = coef(object))

rownames(TAB) <- names(object$coefficients)
    colnames(TAB) <- c("Estimate")

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
residSum=summary(object$residuals, digits=5)[-4])

class(al) <- "summary.alKDE"
al
}

#' @describeIn alKDE print method for summary.alKDE.
#'
#' @return print.summary.alKDE: The object passed to the function is returned invisibly.
#'
#' @export
print.summary.alKDE <- function(x, ...)
{
cat("\nCall:\n")
print(x$call)
cat("\nResiduals:\n")
print(x$residSum)

cat("\nKDE Arc Lengths:\n")
print(x$arclengths)
cat("\n")

printCoefmat(x$coefficients, P.values=FALSE, has.Pvalue=FALSE)
digits <- max(3, getOption("digits") - 3)
cat("\nResidual standard error: ", formatC(x$sigma, digits=digits), " on ",
formatC(x$df, digits=0, format="f"), " degrees of freedom\n", sep="")
cat("Multiple R-squared: ", formatC(x$r.squared, digits=digits),
",\tAdjusted R-squared: ",formatC(x$adj.r.squared, digits=digits),
"\n", sep="")
cat("\tValue of objective function: ",formatC(x$error, digits=digits, format="f"), "\n", sep="")
invisible(x)
}

#' @describeIn alKDE formula method for alKDE.
#' @export
alKDE.formula <- function(formula, data=list(), lower, upper, q1, q2, itermax, type, ...)
{
mf <- model.frame(formula=formula, data=data)
x <- model.matrix(attr(mf, "terms"), data=mf)
y <- model.response(mf)

al <- alKDE.default(formula, data=data, lower=lower, upper=upper, q1=q1, q2=q2, itermax=itermax, type=type, ...)
al$call <- match.call()
al$formula <- formula
al$intercept <- attr(attr(mf, "terms"), "intercept")
al
}

#' @describeIn alKDE predict method for alKDE.
#'
#' @param newdata The data on which the estimated model is to be fitted.
#'
#' @return predict.alKDE: A vector of predicted values resulting from the estimated model.
#'
#' @examples
#' u <- 11:20
#' v <- u+rnorm(10)
#' al <- alKDE(y~x, lower=c(-2,2), upper=c(2,2), q1=c(0.025,0.5), q2=c(0.5,0.975), itermax=50, type=1)
#' predict(al, newdata=data.frame(y=v, x=u))
#'
#' @export
predict.alKDE <- function(object, newdata=NULL, ...)
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