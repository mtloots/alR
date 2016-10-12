#' Moment matching for kernel density estimators.
#'
#' Estimate parameters of a linear model by matching the moments of kernel density estimators.
#'
#' @param formula An LHS ~ RHS formula, specifying the linear model to be estimated.
#' @param data A data.frame which contains the variables in \code{formula}.
#' @param xin Numeric vector of length equal to the number of independent variables, of initial values, for the parameters to be estimated.
#' @param type An integer specifying the bandwidth selection method used, see \code{\link{bw}}.
#' @param ... Arguments to be passed on to the control argument of the \code{\link{optim}} function.
#'
#' @return A generic S3 object with class mmKDE.
#' @importFrom stats coef fitted model.frame model.matrix model.response optim printCoefmat
#'
#' @export
mmKDE <- function(formula, data=list(), xin, type, ...) UseMethod("mmKDE")

#' @describeIn mmKDE default method for mmKDE.
#'
#' @return mmKDE.default: A list with all components from \code{\link{optim}}, as well as:
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
#' \item MOMy: The first \eqn{n} non central moments of the dependent variable, where $\eqn{n} is the number of columns in the design matrix.
#' \item MOMX: The first \eqn{n} non central moments of the independent variables \eqn{\mathbf{X}\underline{\hat{\beta}}}, where $\eqn{n} is the number of columns in the design matrix.
#' }
#' 
#' @examples
#' x <- 1:10
#' y <- x+rnorm(10)
#' XIn <- lm(y~x)
#' mmKDE.default(y~x, xin=coef(XIn), type=-1)
#'
#' @export
mmKDE.default <- function(formula, data=list(), xin, type, ...)
{
mf <- model.frame(formula=formula, data=data)
X <- model.matrix(attr(mf, "terms"), data=mf)
y <- model.response(mf)

h_y <- bw(y, type)
MOMy <- kdeGaussMom(ncol(X), y, h_y)

mom <- optim(xin, momKDE, gamma=X, momy=MOMy, kdeGaussMom=kdeGaussMom, type=type, control=list(...))

mom$intercept <- if(attr(attr(mf, "terms"), "intercept") == 1) TRUE else FALSE

mom$coefficients <- mom$par

labels <- if(mom$intercept) c("Intercept", attr(attr(mf, "terms"), "term.labels")) else attr(attr(mf, "terms"), "term.labels")

names(mom$coefficients) <- labels

mom$df <- nrow(X)-ncol(X)

mom$error <- mom$value

mom$fitted.values <- as.vector(X%*%mom$coefficients)
mom$residuals <- y-mom$fitted.values
mom$call <- match.call()

mom$h_y <- h_y
mom$MOMy <- MOMy
mom$h_X <- bw(mom$fitted.values, type)
mom$MOMX <- kdeGaussMom(ncol(X), mom$fitted.values, mom$h_X)

class(mom) <- "mmKDE"
mom
}

#' @describeIn mmKDE print method for mmKDE.
#'
#' @param x An mmKDE object.
#'
#' @export
print.mmKDE <- function(x, ...)
{
cat("\nCall:\n")
print(x$call)
cat("\nCoefficients:\n")
print(x$coefficients, digits=5)
}

#' @describeIn mmKDE summary method for mmKDE.
#'
#' @param object An mmKDE object.
#'
#' @return summary.mmKDE: A list of class summary.mmKDE with the following components:
#' \itemize{
#' \item call: Original call to \code{mmKDE} function.
#' \item coefficients: A vector with parameter estimates.
#' \item moments: A matrix of the first \eqn{n} moments of the dependent and independent variables that were matched.  The final row corresponds to the estimated bandwidth parameters for each, i.e. \code{h_y} and \code{h_X}, respectively.
#' \item r.squared: The \eqn{r^{2}} coefficient.
#' \item adj.r.squared: The adjusted \eqn{r^{2}} coefficient.
#' \item sigma: The residual standard error.
#' \item df: Degrees of freedom for the model.
#' \item error: Value of the objective function.
#' \item residSum: Summary statistics for the distribution of the residuals.
#' }
#'
#' @export
summary.mmKDE <- function(object, ...)
{
TAB <- cbind(Estimate = coef(object))

rownames(TAB) <- names(object$coefficients)
    colnames(TAB) <- c("Estimate")

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
residSum=summary(object$residuals, digits=5)[-4])

class(mom) <- "summary.mmKDE"
mom
}

#' @describeIn mmKDE print method for summary.mmKDE.
#'
#' @return print.summary.mmKDE: The object passed to the function is returned invisibly.
#'
#' @export
print.summary.mmKDE <- function(x, ...)
{
cat("\nCall:\n")
print(x$call)
cat("\nResiduals:\n")
print(x$residSum)

cat("\nKDE Moments:\n")
print(x$moments)
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

#' @describeIn mmKDE formula method for mmKDE.
#' @export
mmKDE.formula <- function(formula, data=list(), xin, type, ...)
{
mf <- model.frame(formula=formula, data=data)
x <- model.matrix(attr(mf, "terms"), data=mf)
y <- model.response(mf)

mom <- mmKDE.default(formula, data=data, xin=xin, type=type, ...)
mom$call <- match.call()
mom$formula <- formula
mom$intercept <- attr(attr(mf, "terms"), "intercept")
mom
}

#' @describeIn mmKDE predict method for mmKDE.
#'
#' @param newdata The data on which the estimated model is to be fitted.
#'
#' @return predict.mmKDE: A vector of predicted values resulting from the estimated model.
#'
#' @examples
#' u <- 11:20
#' v <- u+rnorm(10)
#' XIn <- lm(y~x)
#' mom <- mmKDE(y~x, xin=coef(XIn), type=-1)
#' predict(mom, newdata=data.frame(y=v, x=u))
#'
#' @export
predict.mmKDE <- function(object, newdata=NULL, ...)
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