library(alR)

context("General functionality of kdeGaussInt")

test_that("numerical result is valid", {
set.seed(1)
mu <- rnorm(1000)
h <- bw(mu, type=1)

integrand <- function(y, mu, h)
{
res <- 0
for (i in 1:length(mu))
{
res = res+(y-mu[i])*exp(-0.5*((y-mu[i])/h)^2)
}
sqrt(1+(1/(2*pi*((length(mu))^2)*(h^6)))*(res)^2)
}

int_r <- numeric(6)
int_r[1] <- integrate(integrand, qkdeGauss(0.025, mu, h)$result, qkdeGauss(0.975, mu, h)$result, mu=mu, h=h)$value
int_r[2] <- integrate(integrand, qkdeGauss(0.025, mu, h)$result, qkdeGauss(0.5, mu, h)$result, mu=mu, h=h)$value
int_r[3] <- integrate(integrand, qkdeGauss(0.5, mu, h)$result, qkdeGauss(0.975, mu, h)$result, mu=mu, h=h)$value
int_r[4] <- integrate(integrand, -1.96, 1.96, mu=mu, h=h)$value
int_r[5] <- integrate(integrand, -1.96, 0, mu=mu, h=h)$value
int_r[6] <- integrate(integrand, 0, 1.96, mu=mu, h=h)$value

int_c <- numeric(6)
int_c[1] <- kdeGaussInt(mu, h, 0.025, 0.975, TRUE)$value
int_c[2] <- kdeGaussInt(mu, h, 0.025, 0.5, TRUE)$value
int_c[3] <- kdeGaussInt(mu, h, 0.5, 0.975, TRUE)$value
int_c[4] <- kdeGaussInt(mu, h, -1.96, 1.96, FALSE)$value
int_c[5] <- kdeGaussInt(mu, h, -1.96, 0, FALSE)$value
int_c[6] <- kdeGaussInt(mu, h, 0, 1.96, FALSE)$value

expect_equal(round(int_r[1],0), round(int_c[1],0))
expect_equal(round(int_r[2],0), round(int_c[2],0))
expect_equal(round(int_r[3],0), round(int_c[3],0))
expect_equal(round(int_r[4],0), round(int_c[4],0))
expect_equal(round(int_r[5],1), round(int_c[5],1))
expect_equal(round(int_r[6],1), round(int_c[6],1))
})

test_that("quantile T/F should yield similar results", {
set.seed(1)
mu <- rnorm(1000)
h <- bw(mu, type=1)

X1 <- qkdeGauss(0.05, mu, h)$result
X2 <- qkdeGauss(0.2, mu, h)$result

expect_equal(kdeGaussInt(mu, h, q1=0.05, q2=0.2, quantile=TRUE)$value, kdeGaussInt(mu, h, q1=X1[[1]], q2=X2[[1]], quantile=FALSE)$value)
})