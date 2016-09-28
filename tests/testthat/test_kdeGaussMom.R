library(alR)

context("General functionality of kdeGaussMom")

test_that("numerical result is valid", {
set.seed(1)
x <- rnorm(100)
h_x <- bw(x, type=1)
mom.theo <- numeric(3)
mom.theo[1] <- mean(x)
mom.theo[2] <- (1/100)*(sum((x)^2+(h_x)^2))
mom.theo[3] <- (1/100)*(sum((x)^3+3*x*(h_x)^2))
mom.calc <- kdeGaussMom(3, x, h_x)
expect_equal(mom.calc[1], mom.theo[1])
expect_equal(mom.calc[2], mom.theo[2])
expect_equal(mom.calc[3], mom.theo[3])
})