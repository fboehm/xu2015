library(xu2015)
context("calc_allocation_prob helper function")

test_that("output length of calc_allocation_prob is same length as y vector input", {
  expect_equal(length(calc_allocation_prob(rnorm(5), 1/2, 0, 1)), 5)
  expect_equal(length(calc_allocation_prob(rnorm(10), 1/10, 0, 1)), 10)
})
