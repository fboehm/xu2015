#' Calculate rho, the acceptance probability for the update K step
#'
#' @param omega_small list containing the current state of parameter vector
#' @param omega_big list containing the proposed state of parameter vector
#' @param y data vector
#' @param ind1 indicator denoting which cluster is affected
#' @param ind2 indicator denoting other affected cluster
#' @param a hyperparameter scalar
#' @param b a scalar hyperparameter
#' @param alpha scalar parameter for dimension-matching
#' @param beta scalar parameter for dimension-matching
#' @param r scalar parameter for dimension-matching
#' @param delta scalar hyperparameter
#' @export
calc_rho_old <- function(y, omega_small, omega_big, ind1, ind2, a, b, alpha, beta, r, delta, theta, tau){
  # work in one-dimension
  ##########
  # unpack omega
  K_big <- omega_big$K
  kappa_big <- omega_big$kappa
  s_big <- omega_big$s # use full vector, s
  w_big <- omega_big$w
  mu_big <- omega_big$mu
  #unpack omega_small
  kappa_small <- omega_small$kappa
  s_small <- omega_small$s
  w_small <- omega_small$w
  K_small <- omega_small$K
  mu_small <- omega_small$mu
  # calc kappa ratio
  #kappa_ratio <- (1 / gamma(a / 2)) * (kappa_big[ind1]) ^ (1 - a / 2) * kappa_big[ind2] * (b / (2 * kappa_big[ind2])) ^ (a / 2) * exp(- 0.5 * b * (1 / kappa_big[ind1] + 1 / kappa_big[ind2] - 1 / kappa_small[ind1])) / kappa_small[ind1] ^ (1 - a / 2)
  log_kappa_ratio <- - lgamma(a / 2) + (1 - a / 2) * log(kappa_big[ind1]) + log(kappa_big[ind2]) + (a / 2) * log(b / (2 * kappa_big[ind2])) - 0.5 * b * (1 / kappa_big[ind1] + 1 / kappa_big[ind2] - 1 / kappa_small[ind1]) - (1 - a / 2) * log(kappa_small[ind1])
  # calc w ratio
  n_big <- numeric(length = 2)
  n_big[1] <- sum(s_big == ind1)
  n_big[2] <- sum(s_big == ind2)
  #w_ratio <- w_big[ind1] ^ (delta - 1 + n_big[1]) * w_big[ind2] ^ (delta - 1 + n_big[2]) / (w_small[ind1] ^ (delta - 1 + n_big[1] + n_big[2]) * beta(delta, K_small * delta))
  log_w_ratio <- (delta - 1 + n_big[1]) * log(w_big[ind1]) + (delta - 1 + n_big[2]) * log(w_big[ind2]) - (delta - 1 + n_big[1] + n_big[2]) * log(w_small[ind1]) - lbeta(delta, K_small * delta) # note use of lbeta() function to return log of beta
  # calc mu ratio
  mu_ratio <- det(calc_C(omega_big$mu, theta, tau)) / det(calc_C(omega_small$mu, theta, tau))
  # calc lik_ratio
  sd_big <- sqrt(1/omega_big$kappa)
  sd_small <- sqrt(1/omega_small$kappa)
  foo_diff <- numeric(length = length(y))
  for (n in 1:length(y)){
    # we're only interested in those entries that belong to the components/clusters of interest
    if (s_big[n] %in% c(ind1, ind2)){
      foo_diff[n] <- dnorm(y[n], mean = mu_big[s_big[n]], sd = sd_big[s_big[n]], log = TRUE) - dnorm(y[n], mean = mu_small[s_small[n]], sd = sd_small[s_small[n]], log = TRUE)
    }
  }
  log_lik_ratio <- sum(foo_diff, na.rm = TRUE)
  ############
  log_posterior_ratio <- log_lik_ratio + log_kappa_ratio + log_w_ratio + log(mu_ratio)
  #### define constants
  q_K_big_d <- 0.5
  q_K_big_c <- 1 / (K_big * (K_big-1))
  qKu <- (K_small == 1) + (K_small > 1) / 2
  qKs <- 1 / K_small
  qu <- dbeta(alpha, 1, 1) * dbeta(beta, 1, 1) * dbeta(r, 2, 2)
  detJ <- (w_small[ind1] ^ (3 + 1) / (w_big[ind1] * w_big[ind2]) ^ (3 / 2)) * kappa_small[ind1] ^ 1.5 * (1 - r ^ 2)
  ###
  ratios <- c(log_lik_ratio, log_kappa_ratio, log_w_ratio, log(mu_ratio), log_posterior_ratio)
  posterior_ratio <- exp(log_posterior_ratio)
  acc_ratio <- posterior_ratio * q_K_big_d * q_K_big_c * detJ / ( (K_big) * qKu * qKs * qu)
  return(list(acc_ratio = acc_ratio, ratios = ratios))
}

