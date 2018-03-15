#' Arc length estimation.
#'
#' A framework for arc length estimation.
#'
#' \itemize{
#' \item Estimate distributional parameters using the method of arc lengths.
#' \item Simulate bootstrap distributions for parameter estimates, resulting from sample arc length statistics.
#' }
#'
#' This method is currently only implemented for the normal and generalised Pareto distributions.  The underlying C code for the Nelder-Mead method of the optim function is used for optimising the objective function.  The tolerance level is set at 1e-15, and a maximum number of 1000 iterations is allowed.  The maximum likelihood estimates are used as initial values for the Nelder-Mead algorithm.
#'
#' @param X A vector of sample values.
#' @param q1,q2 Vectors specifying the quantiles over which arc length segments are to be computed.
#' @param dc TRUE/FALSE:  Should the discrete or continuous sample statistic be used.
#' @param type The type of bandwidth estimator for the underlying KDE; see \code{\link{bw}}.
#' @param bootstraps An integer specifying the size of the parametric bootstrap.
#' @param distribution The distribution to be fitted, 1=normal (default), 2=generalised Pareto.
#' @param ... Additional arguments passed to \code{alEfit} (not currently used).
#'
#' @return alEfit: A generic S3 object with class alEfit.
#'
#' @export
alEfit <- function(X, q1, q2, dc, type, bootstraps, distribution=1, ...) UseMethod("alEfit")

#' @describeIn alEfit default method for alEfit.
#'
#' @return alEfit.default: A list with all components from \code{\link{alE}}, as well as :
#' \itemize{
#' \item dc: TRUE/FALSE Was the discrete or continuous sample arc length statistic used?
#' \item q1, q2: The segments over which the arc length(s) were calculated.
#' \item bw: The bandwidth for the kernel density estimator.
#' \item distribution: The distribution fitted to the data.
#' \item dist: A numeric matrix whose columns represent a bootstrap distribution for the corresponding parameter estimate.
#' \item se: A numeric vector with standard errors, obtained by a parametric bootstrap.
#' \item bootstraps: Number of bootstrap samples.
#' }
#' 
#' @examples
#' alEfit(x, q1=0.025, q2=0.975, dc=TRUE, type=-1, bootstraps=50)
#' alEfit(x, q1=0.025, q2=0.975, dc=FALSE, type=-1, bootstraps=50)
#'
#' @export
alEfit.default <- function(X, q1, q2, dc, type, bootstraps, distribution=1, ...)
{
al <- alE(X, q1, q2, dc, type, distribution)

al$distribution <- switch(distribution, "normal distribution", "generalised Pareto distribution")

al$dist <- alEfitdist(X, q1, q2, dc, type, bootstraps, distribution)

coefMean <- colMeans(al$dist)
coefVar <- (1/(bootstraps-1))*colSums((al$dist-coefMean)^2)
al$se <- sqrt(coefVar)

al$bw <- bw(X, type)
al$dc <- dc
al$q1 <- q1
al$q2 <- q2
al$bootstraps <- bootstraps

class(al) <- "alEfit"
al
}

#' @describeIn alEfit print method for alEfit.
#'
#' @param x An alEfit object.
#'
#' @export
print.alEfit <- function(x, ...)
{
cat(paste("Fitting of the ", x$distribution, " by the method of arc lengths\n", sep=""))
    cat("Parameters:\n")
out <- cbind.data.frame("estimate" = x$par, "Std. Error" = x$se)
if (x$distribution=="normal distribution") rownames(out) <- c("mu","sigma")
if (x$distribution=="generalised Pareto distribution") rownames(out) <- c("mu","sigma","alpha")
print(out)
}