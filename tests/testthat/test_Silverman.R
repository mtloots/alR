library(alR)

context("General functionality of Silverman")

test_that("numerical result is valid", {
x1 <- c(1,2,3)
expect_equal(Silverman(x1), sd(x1)*(4/(3*3))^(1/5))
x2 <- c(10,20,30)
expect_equal(Silverman(x2), sd(x2)*(4/(3*3))^(1/5))
})