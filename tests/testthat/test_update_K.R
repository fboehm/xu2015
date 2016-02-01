library(xu2015)
context("Lengths of components outputted by update_K function")

test_that("outputted mu is same length as outputted w", {
  expect_equal(length(update_K(rnorm(15), mu = 1:3, sigma = c(1.1, 1.3, 1.5), w = c(1/2, 1/3, 1/6), s = 1 + rbinom(n=15, size = 2, prob = 1/2), tau = 1, theta = 1, delta = 1)$mu), length(update_K(rnorm(15), mu = 1:3, sigma = c(1.1, 1.3, 1.5), w = c(1/2, 1/3, 1/6), s = 1 + rbinom(n=15, size = 2, prob = 1/2), tau = 1, theta = 1, delta = 1)$w))
  expect_equal(length(update_K(rnorm(25), mu = 1:3, sigma = c(1.1, 1.3, 1.5), w = c(1/2, 1/3, 1/6), s = 1 + rbinom(n=25, size = 2, prob = 1/2), tau = 1, theta = 1, delta = 1)$mu), length(update_K(rnorm(25), mu = 1:3, sigma = c(1.1, 1.3, 1.5), w = c(1/2, 1/3, 1/6), s = 1 + rbinom(n=25, size = 2, prob = 1/2), tau = 1, theta = 1, delta = 1)$w))
})
