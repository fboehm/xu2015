library(xu2015)
context("calc_allocation_prob helper function")

test_that("output length of calc_allocation_prob is same length as y vector input", {
  expect_equal(length(calc_allocation_prob(rnorm(5), 1/2, 0, 1)), 5)
  expect_equal(length(calc_allocation_prob(rnorm(10), 1/10, 0, 1)), 10)
})

##################################################
context("define_extra_parameters helper function for update_K")

test_that("output is length 3", {
  expect_equal(length(define_extra_parameters()), 3)
})

test_that("all three entries are between 0 and 1 (since they're draws from beta distributions)", {
  expect_equal(define_extra_parameters() < 1, rep(TRUE, 3))
  expect_equal(define_extra_parameters() > 0, rep(TRUE, 3))
})
##################################################
context("define_big_w helper function for update_K")

test_that("length of output equals 1 plus length of input w vector", {
  expect_equal(length(define_big_w(w = 1:3 / 6, alpha = runif(n = 1, 0, 1), ind1 = 1, ind2 = 4)), 4)
  expect_equal(length(define_big_w(w = 1:4 / 10, alpha = runif(n = 1, 0, 1), ind1 = 1, ind2 = 5)), 5)
  expect_equal(length(define_big_w(w = 1:10 / sum(1:10), alpha = runif(n = 1, 0, 1), ind1 = 1, ind2 = 11)), 11)
})

test_that("output vector sums to 1", {
  expect_equal(sum(define_big_w(w = 1:3 / 6, alpha = runif(n = 1, 0, 1), ind1 = 1, ind2 = 4)), 1)
  expect_equal(sum(define_big_w(w = 1:4 / 10, alpha = runif(n = 1, 0, 1), ind1 = 1, ind2 = 5)), 1)
  expect_equal(sum(define_big_w(w = 1:10 / sum(1:10), alpha = runif(n = 1, 0, 1), ind1 = 1, ind2 = 11)), 1)
})


test_that("all entries of output are positive", {
  expect_equal(define_big_w(w = 1:3 / 6, alpha = runif(n = 1, 0, 1), ind1 = 1, ind2 = 4) > 0, rep(TRUE, 4))
  expect_equal(define_big_w(w = 1:4 / 10, alpha = runif(n = 1, 0, 1), ind1 = 1, ind2 = 5) > 0 , rep(TRUE, 5))
  expect_equal(define_big_w(w = 1:10 / sum(1:10), alpha = runif(n = 1, 0, 1), ind1 = 1, ind2 = 11) > 0, rep(TRUE, 11))
})
