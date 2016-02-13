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
###########################################################################
context("define_big_mu helper for update_K")

test_that("output has length one greater than inputted mu", {
  expect_equal(length(define_big_mu(mu = c(1,2), w = c(0.5, 0.5), ind1 = 1, ind2 = 3, sigma = c(0.5, 0.4), r = 1/2)), 3)
  expect_equal(length(define_big_mu(mu = c(1,2,3), w = c(0.3, 0.3, 0.4), ind1 = 1, ind2 = 4, sigma = c(0.5, 0.4, 0.3), r = 1/2)), 4)
  expect_equal(length(define_big_mu(mu = c(1,2, 3, 4), w = c(0.25, 0.25, 0.25, 0.25), ind1 = 1, ind2 = 5, sigma = c(0.5, 0.4, 0.3, 0.3), r = 1/2)), 5)
})

test_that("all but two components of big mu are in little mu", {
  expect_identical(define_big_mu(mu = c(1,2), w = c(0.5, 0.5), ind1 = 1, ind2 = 3, sigma = c(0.5, 0.4), r = 1/2)[-c(1,3)], 2)
})
###########################################################################
context("define_big_s helper for update_K")

test_that("length of output equals length of inputted s",{
  sample_size <- 300
  mu_true <- c(-10, 10)
  w_true <- c(1/2, 1/2)
  means <- sample(mu_true, size = sample_size, replace = TRUE, prob = w_true)
  sd_true <- 0.01
  sigma_true <- rep(sd_true, 2)
  y <- rnorm(n = sample_size, mean = means, sd = sd_true)
  s <- sample(1:3, size = sample_size, replace = TRUE)
  mu <- c(-9, 9, 11)
  sigma <- rep(1, 3)
  s_tab <- table(s)
  s_tab[1:3 ] -> s_tab2
  names(s_tab2) <- NULL
  w <- s_tab2 / sample_size
  w_new <- define_big_w(w, 1/2, 2)
  s_big <- define_big_s(s, 2, 4, y, w_big = rep(0.25, 4), mu_big = c(mu, 5), kappa_big = rep(1,4))
  ###
  expect_equal(length(s_big), sample_size)
  })

test_that("exactly 1 group is split into two groups",{
  sample_size <- 300
  mu_true <- c(-10, 10)
  w_true <- c(1/2, 1/2)
  means <- sample(mu_true, size = sample_size, replace = TRUE, prob = w_true)
  sd_true <- 0.01
  sigma_true <- rep(sd_true, 2)
  y <- rnorm(n = sample_size, mean = means, sd = sd_true)
  s <- sample(1:1, size = sample_size, replace = TRUE)
  mu <- c(-9)
  sigma <- rep(1, 1)
  s_tab <- table(s)
  s_tab[1:1 ] -> s_tab2
  names(s_tab2) <- NULL
  w <- s_tab2 / sample_size
  w_big <- define_big_w(w, 1/2, 2)
  s_big <- define_big_s(s, 2, 4, y, w_big = w_big, mu_big = c(mu, 5), kappa_big = rep(1,2))
  big_tab<- table(s_big)
  small_tab <- table(s)
  ###
  expect_equal(sum(big_tab %in% small_tab), length(mu) - 1)
})
