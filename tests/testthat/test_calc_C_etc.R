library(xu2015)
context("calc_little_c")

test_that("output length of calc_little_c is 1", {
  expect_equal(length(calc_little_c(rnorm(1), rnorm(1), theta = 1 )), 1)
})

test_that("output of calc_little_c is non-negative", {
  expect_gte(calc_little_c(rnorm(1), rnorm(1), theta = 1), 0)
})

context("calc_little_q")

test_that("output of calc_little_q is non-negative", {
  expect_gte(calc_little_q(rnorm(1), tau = 1), 0)
})

context("calc_C")

test_that("output is a square matrix", {
  expect_equal(nrow(calc_C(rep(0, 4), 1, 1)), ncol(calc_C(rep(0,4), 1, 1)))
  expect_equal(nrow(calc_C(0, 1, 1)), ncol(calc_C(0, 1, 1)))
})
