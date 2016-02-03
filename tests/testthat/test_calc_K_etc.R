library(xu2015)
context("calc_allocation_prob helper function")

test_that("output length of calc_allocation_prob is same length as y vector input", {
  expect_equal(length(calc_allocation_prob(rnorm(5), 1/2, 0, 1)), 5)
  expect_equal(length(calc_allocation_prob(rnorm(10), 1/10, 0, 1)), 10)
})

context("define_extra_parameters helper function for update_K")

test_that("output is length 3", {
  expect_equal(length(define_extra_parameters()), 3)
})

test_that("all three entries are between 0 and 1 (since they're draws from beta distributions)", {
  expect_equal(define_extra_parameters() < 1, rep(TRUE, 3))
  expect_equal(define_extra_parameters() > 0, rep(TRUE, 3))
})

