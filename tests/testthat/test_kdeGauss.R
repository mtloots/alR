library(alR)

context("General functionality of kdeGauss")

test_that("numerical result is valid", {
set.seed(1)
x <- rnorm(100)
h_x <- bw(x, type=1)
expect_equal(dkdeGauss(0, x, h_x), (1/(length(x)*h_x*sqrt(2*pi)))*sum(exp(-(1/2)*((0-x)/h_x)^2)))
expect_equal(dkdeGauss(-1, x, h_x), (1/(length(x)*h_x*sqrt(2*pi)))*sum(exp(-(1/2)*((-1-x)/h_x)^2)))
expect_equal(dkdeGauss(2.5, x, h_x), (1/(length(x)*h_x*sqrt(2*pi)))*sum(exp(-(1/2)*((2.5-x)/h_x)^2)))
})