library(alR)

context("General functionality of qsamp")

test_that("numerical result is valid", {
x <- rnorm(100)
expect_equal(qsamp(x, 0.1), quantile(x, 0.1, type=1)[[1]])
expect_equal(qsamp(x, 0.5), quantile(x, 0.5, type=1)[[1]])
expect_equal(qsamp(x, 0.9), quantile(x, 0.9, type=1)[[1]])
})