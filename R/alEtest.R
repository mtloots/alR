#' Arc length estimation.
#'
#' Goodness-of-fit using arc lengths.
#'
#' First the distributional parameters of a sample is estimated using the continuous arc length sample statistic (see \code{\link{alE}}).  The calculated sample arc length statistic is then compared to the distribution of that particular sample statistic, obtained by a parametric bootstrap, using the estimated parameters (see \code{\link{alEdist}}).  This finally leads to the calculation of a p-value for a goodness-of-fit test, based on the simulated distribution.
#'
#' This method is currently only implemented for the normal distribution, and for a single arc length segment.
#'
#' @param X A vector of sample values.
#' @param mu A real value specifying the mean of the normal distribution.
#' @param sigma A positive real number specifying the scale parameter of the normal distribution.
#' @param q1,q2 Vectors specifying the quantiles (or points if quantile=FALSE) over which arc length segments are to be computed.
#' @param quantile TRUE/FALSE whether q1 and q2 are quantiles, or elements of the domain of \code{x}.
#' @param dc TRUE/FALSE:  Should the discrete or continuous sample statistic be used.
#' @param type The type of bandwidth estimator for the underlying KDE; see \code{\link{bw}}.
#' @param n An integer specifying the sample size.
#' @param bootstraps An integer specifying the size of the parametric bootstrap.
#' @param ... Additional arguments passed to \code{alEtest} (not currently used).
#'
#' @return alEtest: A generic S3 object with class alEtest.
#'
#' @export
alEtest <- function(X, mu, sigma, q1, q2, type, bootstraps, ...) UseMethod("alEtest")

#' @describeIn alEtest default method for alEtest.
#'
#' @return alEtest.default: A list with the following components:
#' \itemize{
#' \item q1, q2: The segment over which the arc length was calculated.
#' \item mu: A real value specifying the mean of the normal distribution.
#' \item sigma: A positive real number specifying the scale parameter of the normal distribution.
#' \item bw: The bandwidth for the kernel density estimator.
#' \item dist: A numeric matrix whose columns represent a bootstrap distribution for the corresponding sample arc length statistic.
#' \item statistic: The value of the observed sample statistic.
#' \item pvalue: The p-value for the test based on a parametric bootstrap sample.
#' \item bootstraps: Number of bootstrap samples.
#' }
#'
#' @examples
#' \dontrun{
#' x <- rnorm(1000)
#' s1 <- alE(x, 0.025, 0.975, TRUE, -1)
#' alEtest(x, mu=s1$par[1], sigma=s1$par[2], q1=0.025, q2=0.975,
#' type=-1, bootstraps=50)
#' s2 <- alE(x, 0.025, 0.975, FALSE, -1)
#' alEtest(x, mu=s2$par[1], sigma=s2$par[2], q1=0.025, q2=0.975,
#' type=-1, bootstraps=50)
#' }
#'
#' @export
alEtest.default <- function(X, mu, sigma, q1, q2, type, bootstraps, ...)
{
if (length(q1) > 1 | length(q2) > 1)
{
stop("This test is only implemented for single arc length segments.")
}
else
{
al <- list()
p1 <- qsamp(X, q1)
p2 <- qsamp(X, q2)
al$dist <- alEdist(length(X), bootstraps, mu, sigma, p1, p2, FALSE, FALSE, type)
al$mu <- mu
al$sigma <- sigma
al$bw <- bw(X, type)
al$q1 <- q1
al$q2 <- q2
al$bootstraps <- bootstraps
al$statistic <- kdeGaussInt(X, al$bw, p1, p2, FALSE)$value

lpv <- mean(al$dist <= al$statistic)
al$pvalue <- 2*min(lpv, 1-lpv)

class(al) <- "alEtest"
al
}
}

#' @describeIn alEtest print method for alEtest.
#'
#' @param x An alEtest object.
#'
#' @export
print.alEtest <- function(x, ...)
{
cat("Arc length test of goodness-of-fit\n")
cat("Null hypothesis: N(", x$mu, ", ", x$sigma, ") distribution\n")
cat("\n s2 = ", x$statistic, ", p-value = ", x$pvalue, "\n")
}