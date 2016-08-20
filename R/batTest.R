#' Bhattacharryya test for comparing  two samples.
#'
#' Use the multinomial distribution to test the hypothesis that two samples come from the same distribution.
#' 
#' It is assumed that the two samples come from the same kernel density distribution.  The support of the KDE of the first sample is divided into \eqn{k} equally spaced quantiles, and then compared to the corresponding proportions of the second.
#'
#' @param y,x Samples to be compared.
#' @param k Number of proportions.
#'
#' @return bat.test: A list with the following components:
#' \itemize{
#' \item k: number of proportions used.
#' \item chi: The test statistic for the Bhattacharrayya test.
#' \item p.value: The p-value of the test.
#' }
#' 
#' @examples
#' y <- rnorm(1000)
#' x <- rnorm(1000)
#' bhatt.test(y,x,10)
#'
#' @export
bhat.test <- function(y, x, k)
{
h_y <- Silverman(y)
pts <- seq(0,1,length.out=k+1)[2:k]
CDF <- sapply(pts, function(i) qkdeGauss(i, y, h_y))
yecdf <- ecdf(y)
xecdf <- ecdf(x)

yprop <- numeric(k)
xprop <- numeric(k)

for (i in 1:length(CDF))
{
yprop[i] <- yecdf(CDF[i])
xprop[i] <- xecdf(CDF[i])
}

bhattout <- list(k,
yprop,
xprop)
}