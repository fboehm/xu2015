#' Calculate the allocation probability
#'
#' @param y scalar, ie, one element of data vector
#' @param w a scalar weight
#' @param mu a scalar class mean
#' @param kappa a scalar inverse variance
#' @export
calc_allocation_prob <- function(y, w, # w a scalar
                                 mu, # mu a scalar
                                 kappa # scalar
){
  w * sqrt(kappa) * dnorm(y, mean = mu, sd = 1 / sqrt(kappa))
}


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
calc_rho <- function(y, omega_small, omega_big, ind1, ind2, a, b, alpha, beta, r, delta, theta, tau){
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



#' Update $K$ in Gibbs sampling
#'
#' @param y data vector
#' @param mu vector of class means
#' @param w vector of class weights
#' @param sigma vector of inverse variances, one entry per class
#' @param tau hyperparameter
#' @param theta hyperparameter
#' @param delta hyperparameter
#' @export
update_K <- function(y, mu, w, sigma, s, tau, theta, delta){
  kappa <- 1/sigma^2
  a <- 1 / (4 * tau ^ 2)
  b <- 1 / (2 * theta ^ 2)
  q_down <- 0
  if (length(mu) > 1){
    q_down <- 0.5
  }
  ### define extra parameters
  alpha <- rbeta(n = 1, 1, 1)
  beta <- rbeta(n = 1, 1, 1)
  r <- rbeta(n = 1, 2, 2)
  # define s_new
  s_new <- s
  ### decide to split or combine
  split <- as.logical(rbinom(n = 1, size = 1, prob = 1 - q_down))
  ##########
  if (split){
    omega_small <- list(K = length(mu), mu = mu, kappa = kappa, w = w, s = s)
    #sampling_vec <- as.integer(names(table(s)))
    sampling_vec <- 1:length(mu)
    # we introduce sampling_vec because there's a chance that one or more clusters has no observations assigned to it.
    ind1 <- sample(sampling_vec, size = 1, replace = FALSE)
    ind2 <- length(mu) + 1
    ## edit w; make K+1 the 'new' component
    w_new <- w
    w_new[ind1] <- alpha * w[ind1]
    w_new[ind2] <- (1 - alpha) * w[ind1]
    ## edit mu
    mu_new <- mu
    mu_new[ind1] <- mu[ind1] - sqrt(w_new[ind2] / w_new[ind1]) * r / sigma[ind1]
    mu_new[ind2] <- mu[ind1] + sqrt(w_new[ind1] / w_new[ind2]) * r / sigma[ind1]
    ## edit kappa (sigma)
    kappa_new <- kappa
    kappa_new[ind1] <- beta * (1 - r) ^ 2 * (w[ind1] / w_new[ind1]) * kappa[ind1]
    kappa_new[ind2] <- (1 - beta) * (1 - r) ^ 2 * (w[ind1] / w_new[ind2]) * kappa[ind1]
    ## Yanxun told me to see Richardson & Green 1997
    # to get the method for
    # allocation, ie, to assign values to n_new
    for (i in 1:length(s)){
      if (s[i] == ind1){
        foo_p1 <- calc_allocation_prob(y = y[i], w = w_new[ind1], mu = mu_new[ind1], kappa = kappa_new[ind1])
        foo_p2 <- calc_allocation_prob(y = y[i], w = w_new[ind2], mu = mu_new[ind2], kappa = kappa_new[ind2])
        foo_bin <- rbinom(n = 1, size = 1, prob = foo_p1 / (foo_p1 + foo_p2))
        s_new[i] <- foo_bin * ind1 + (1 - foo_bin) * ind2
      }
    }

    ############
    omega_big <- list(K = length(mu) + 1, mu = mu_new, kappa = kappa_new, w = w_new, s = s_new)
    foo <- calc_rho(y, omega_small, omega_big, ind1, ind2, a, b, alpha, beta, r, delta = 1, theta = theta, tau = tau)
    #print(foo)
    u <- runif(n = 1, min = 0, max = 1)
    # compare u to acceptance ratio & decide to accept or reject
    if (u <- foo$acc_ratio) {out <- list(w = w_new, mu = mu_new, kappa = kappa_new, s = s_new, ar = foo, u = u, split = split)} else {out <- list(w = w, mu = mu, kappa = kappa, s = s, ar = foo, u = u, split = split)}
  }else { ## combine
    sampling_vec <- as.integer(names(table(s)))
    # we introduce sampling_vec because there's a chance that none of the y's are assigned to some of our clusters.
    indices <- sample(sampling_vec, size=2, replace=FALSE)
    ind1 <- min(indices)
    ind2 <- max(indices)
    # edit w
    w_new <- w
    w_new[ind1]<- w[ind1] + w[ind2]
    w_new <- w_new[-ind2]
    # edit mu
    mu_new <- mu
    mu_new[ind1] <- (w[ind1] * mu[ind1] + w[ind2] * mu[ind2])/w_new[ind1]
    mu_new <- mu_new[-ind2]
    # edit kappa
    kappa_new <- kappa
    kappa_new[ind1] <- (w[ind1] / w_new[ind1]) * kappa[ind1] + (w[ind2]/w_new[ind1])*kappa[ind2] + (w[ind1]*w[ind2]/(w_new[ind1])^2)*(mu[ind1] - mu[ind2])^2
    kappa_new <- kappa_new[-ind2]
    # edit s
    s_new[s_new == ind2]<- ind1
    for (k in 1:(length(mu)-2)){
      s_new[s_new == (k+ind2)] <- k + ind2 - 1
    }
    #################
    # calculate acceptance ratio
    omega_small <- list(K = length(mu) - 1, mu = mu_new, kappa = kappa_new, w = w_new, s = s_new)
    omega_big <- list(K = length(mu), mu = mu, kappa = kappa, w = w, s = s)
    bar <- calc_rho(y, omega_small = omega_small, omega_big = omega_big, ind1, ind2, a, b, alpha, beta, r, delta = 1, theta = theta, tau = tau)
    #print(bar)
    acc_ratio <- bar$acc_ratio
    u <- runif(n = 1, min = 0, max = 1)
    # compare u to acceptance ratio & decide to accept or reject
    if (u < acc_ratio) {out <- list(w = w_new, mu = mu_new, kappa = kappa_new, s=s_new, ar = bar, u = u, split = split)} else {out <- list(w = w, mu = mu, kappa = kappa, s = s, ar = bar, u = u, split = split)}
  }
  return(out)
}
