library(alR)

context("General functionality of GaussInt")

test_that("numerical result is valid", {
set.seed(1)
mu <- 2
sigma <- 3.5

integrand <- function(y, mu, sigma)
{
sqrt(1+(1/(2*pi*sigma^6))*((y-mu)^2)*exp(-((y-mu)/sigma)^2))
}

int_r <- numeric(6)
int_r[1] <- integrate(integrand, qnorm(0.025, mu, sigma), qnorm(0.975, mu, sigma), mu=mu, sigma=sigma)$value
int_r[2] <- integrate(integrand, qnorm(0.025, mu, sigma), qnorm(0.5, mu, sigma), mu=mu, sigma=sigma)$value
int_r[3] <- integrate(integrand, qnorm(0.5, mu, sigma), qnorm(0.975, mu, sigma), mu=mu, sigma=sigma)$value
int_r[4] <- integrate(integrand, -1.96, 1.96, mu=mu, sigma=sigma)$value
int_r[5] <- integrate(integrand, -1.96, 0, mu=mu, sigma=sigma)$value
int_r[6] <- integrate(integrand, 0, 1.96, mu=mu, sigma=sigma)$value

int_c <- numeric(6)
int_c[1] <- GaussInt(mu, sigma, 0.025, 0.975, TRUE)$value
int_c[2] <- GaussInt(mu, sigma, 0.025, 0.5, TRUE)$value
int_c[3] <- GaussInt(mu, sigma, 0.5, 0.975, TRUE)$value
int_c[4] <- GaussInt(mu, sigma, -1.96, 1.96, FALSE)$value
int_c[5] <- GaussInt(mu, sigma, -1.96, 0, FALSE)$value
int_c[6] <- GaussInt(mu, sigma, 0, 1.96, FALSE)$value

expect_equal(int_r[1], int_c[1])
expect_equal(int_r[2], int_c[2])
expect_equal(int_r[3], int_c[3])
expect_equal(int_r[4], int_c[4])
expect_equal(int_r[5], int_c[5])
expect_equal(int_r[6], int_c[6])
})

test_that("quantile T/F should yield similar results", {
mu <- 2
sigma <- 3.5

X1 <- qnorm(0.05, mu, sigma)
X2 <- qnorm(0.2, mu, sigma)

expect_equal(GaussInt(mu, sigma, q1=0.05, q2=0.2, quantile=TRUE)$value, GaussInt(mu, sigma, q1=X1, q2=X2, quantile=FALSE)$value)
})