#' Bhattacharryya test for comparing  two samples.
#'
#' Use the multinomial distribution to test the hypothesis that two samples come from the same distribution.
#' 
#' It is assumed that the two samples come from the same kernel density distribution.  The support of the KDE of the first sample is divided into \eqn{k} equally spaced quantiles, and then compared to the corresponding proportions of the second.
#'
#' @param y,x Samples to be compared.
#' @param k Number of proportions.
#'
#' @importFrom stats ecdf pchisq
#'
#' @return bat.test: A list with the following components:
#' \itemize{
#' \item df=2*k: where k is the number of proportions used.
#' \item y.prop, x.prop: Vectors of proportions.
#' \item D2: Measure of divergence between samples.
#' \item test.stat: The test statistic for the Bhattacharrayya test.
#' \item p.value: The p-value of the test.
#' }
#' 
#' @examples
#' y <- rnorm(1000)
#' x <- rnorm(1000)
#' bhatt.test(y,x,10)
#'
#' @export
bhatt.test <- function(y, x, k)
{
h_y <- Silverman(y)
pts <- seq(0,1,length.out=k+1)[2:k]
CDF <- sapply(pts, function(i) qkdeGauss(i, y, h_y)$result)
yecdf <- ecdf(y)
xecdf <- ecdf(x)

yprop <- numeric(k+1)
xprop <- numeric(k+1)

for (i in 1:length(CDF))
{
yprop[i+1] <- yecdf(CDF[i])
xprop[i+1] <- xecdf(CDF[i])
}

yprop[k+1] <- 1
xprop[k+1] <-1

Yprop <- diff(yprop)
Xprop <- diff(xprop)

numerator <- (Yprop-Xprop)^2
denominator <- Yprop+Xprop
chi <- length(y)*sum(numerator/denominator)
p.value <- pchisq(chi, 2*k, lower.tail=FALSE)

bhattout <- list(df=2*k,
y.prop=Yprop,
x.prop=Xprop,
D2=chi/(2*length(y)),
test.stat=chi,
p.value=p.value)
bhattout
}